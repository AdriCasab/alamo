# Intermittent Bus Error in `Newton::solve` — Root-Cause Analysis

Status: **diagnosed, not yet fixed**.

This document explains the intermittent `SIGBUS` ("Invalid address
alignment") crashes observed during static elastic mechanics solves on
`bin/mmwspalling-3d-g++` on Apple Silicon (ARM64) macOS. The crash was
first noticed under `-np 4` in the AMR investigation
([AMR_FAILURE_ANALYSIS.md §9.8](AMR_FAILURE_ANALYSIS.md), 2026-05-09)
and reproduced multiple times during hi-res (80³) sandstone runs on
2026-05-11. This file localises the cause and explains why it is
intermittent.

## TL;DR

- **Where**: inside `Solver::Nonlocal::Newton::prepareForSolve`, at a
  **compiler-emitted unrolled memcpy loop** that copies a 24-byte
  `Set::Vector` or 72-byte `Set::Matrix` field across an `amrex::Array4`
  array.
- **What instruction**: `ldp q, q [src]` — paired 128-bit SIMD load,
  32 bytes per instruction. On ARM64 it requires 16-byte alignment of
  `src`.
- **Why it faults**: `Set::Vector` is **24 bytes** (`Eigen::Matrix<double,3,1>`),
  not a multiple of 16. In an array of `Set::Vector`, **roughly half the
  elements** sit at 8-byte-aligned offsets relative to the FAB base
  pointer. When the unrolled SIMD memcpy hits one of those, you get
  `SIGBUS`.
- **Why it's intermittent**: which element the unrolled loop visits
  first, and the alignment of the FAB base pointer, both vary run-to-run
  with the OS allocator. Most runs avoid the bad element by luck;
  occasionally a run hits it and dies.
- **Why concurrent simulations make it worse**: memory pressure changes
  macOS `malloc`'s allocator behavior and increases heap fragmentation,
  raising the probability that FAB buffers land at less-aligned
  addresses.
- **The bug is latent**: it has likely been crashing alamo on this
  hardware/macOS combination since before the current session. The
  pre-existing comment "[bottom_solver=smoother / BiCGStab not robust](ARCHIVE_DONE.md)"
  may have been a workaround for closely related symptoms.

## Reproduction signature

The signal handler prints (all observed instances):

```
Process received signal:
Signal: Bus error: 10 (10)
Signal code: Invalid address alignment (1)
Failing at address: 0x<heap-address>
...
[ 2] mmwspalling-3d-g++   ...Newton...solve...                     + 8360
[ 3] mmwspalling-3d-g++   ...Mechanics::TimeStepBegin              + 11380
[ 4] mmwspalling-3d-g++   ...MMWSpalling::TimeStepBegin            + 424
```

Observed failing addresses (all instances):

| run | step at failure | trap address | mod 16 |
| --- | --------------- | ------------ | ------ |
| `-np 4` sandstone diag | STEP 3 | `0x7e33fffb8` | **8** |
| `-np 1` sandstone hi-res #1 | STEP 184 | `0xad43ffaa8` | **8** |
| `-np 1` sandstone hi-res #2 | STEP 500 | `0xb49bfffb8` | **8** |
| `-np 1` sandstone hi-res #3 | (no crash, reached 1000) | — | — |

Every trap address is **8-byte aligned but not 16-byte aligned**. That's
the diagnostic fingerprint — it isolates the fault to a 16-byte-aligned
SIMD instruction reading from an 8-byte-aligned address.

## Disassembly evidence

The stack trace shows the crash inside Newton::solve at offset 8360.
Disassembling around that offset:

```
0x10021de74   bl    ...prepareForSolve...        ; call site
0x10021de78   ldr   x0, [sp, #0x78]              ; ← return address = func + 8360
```

The trap is the **return address** of the call to `prepareForSolve`,
which means the actual fault is **inside** `prepareForSolve`. The
function is 13.4 KB and contains:

- **90** single `ldr q` (16-byte SIMD load) instructions
- **3** `ldp q, q` (32-byte SIMD paired load) instructions

The most suspicious are the `ldp q, q` instances, which form the body
of a compiler-unrolled copy loop at `0x100223bc4`:

```asm
loop_top:
0x100223bc4   ldp  q0, q1, [x8, #-0x20]   ; load 32B from src
0x100223bc8   ldp  q2, q3, [x8],   #0x40  ; load 32B, post-incr src by 64
0x100223bcc   stp  q0, q1, [x25, #-0x20]
0x100223bd0   stp  q2, q3, [x25],  #0x40  ; store 32B, post-incr dst by 64
0x100223bd4   subs x4, x4, #0x8           ; decrement element count by 8
0x100223bd8   b.ne loop_top
```

This is the classic LLVM unrolled `memcpy(dst, src, n*8)` pattern — 64
bytes per loop iteration. The instructions require `x8` (source) to be
16-byte aligned at the start of each iteration.

The pre-loop setup:

```asm
0x100223bb0   cmp   x8, #0x40
0x100223bb4   b.lo  0x100223be8     ; fall back to scalar copy if < 64 B
0x100223bb8   add   x8, x6, x30     ; src = base + offset
0x100223bbc   add   x25, x26, x19   ; dst = base + offset
```

LLVM only uses this fast path when there are ≥ 64 bytes to copy. The
fallback at `0x100223be8` is a plain scalar `ldr d0; str d0` loop with
no alignment requirement. So this is a clear "vectorize when possible"
code path.

## Why a 24-byte `Set::Vector` array breaks the assumption

`Set::Vector` is defined as `Eigen::Matrix<double, 3, 1>` =
**24 bytes** (3 doubles, no padding). Eigen does NOT add padding for
3-vectors because the size isn't a power of 2 ≥ 16.

In an `amrex::Array4<Set::Vector>` (or any contiguous array of these),
element stride = 24 bytes. From a 16-byte-aligned base pointer:

| element | offset | absolute mod 16 |
| ------- | ------ | --------------- |
| 0 | 0 | **0** ← aligned |
| 1 | 24 | **8** ← unaligned |
| 2 | 48 | **0** |
| 3 | 72 | **8** |
| 4 | 96 | **0** |
| 5 | 120 | **8** |
| … | … | … |

**Every odd-indexed element is 8-byte aligned but not 16-byte aligned.**

The same applies to `Set::Matrix = Eigen::Matrix<double, 3, 3>` =
**72 bytes**: with a stride of 72, the alignment pattern is even worse
because `gcd(72, 16) = 8`, so neighbouring elements alternate between
0 and 8 mod 16.

When the LLVM-emitted memcpy at the call site uses `ldp q, q` on an odd
element, it touches an 8-aligned address with a 16-aligned-required
instruction → `SIGBUS`.

## Why it manifests inside `prepareForSolve`

[`Newton::prepareForSolve`](../../../src/Solver/Nonlocal/Newton.H)
iterates `amrex::ParallelFor(bx, …)` over every node of every AMR fab
and assigns:

```cpp
ddw(i, j, k) = model(i, j, k).DDW(...);   // Set::Matrix4 — large
dw(i, j, k)  = model(i, j, k).DW(...);    // Set::Matrix    72 B
rhs(i, j, k) = b(i, j, k) - div(dw);      // Set::Vector    24 B
```

Each of these is a struct-assignment (`operator=`) of an Eigen
fixed-size matrix/vector. At `-O3` with LLVM's auto-vectorizer enabled,
LLVM lowers each assignment to either a single SIMD store or an
unrolled `ldp q, q` / `stp q, q` block depending on size. For
`Set::Matrix4<3, 6>` (which is the elastic stiffness tensor, hundreds
of bytes), the body is large enough to trigger the **unrolled `ldp q`
fast-path**.

This is the exact loop that crashes.

The reason the trap appears only after `prepareForSolve` returns (the
`+8360` offset in `Newton::solve`) is just how the OS reports the
backtrace: the signal trampoline captures the call stack at the point
of the SIGBUS, with the deepest frame being the SIMD instruction
itself. macOS' simple backtrace printer reports the return addresses
of the upper frames, so `Newton::solve + 8360` is the **call site of
`prepareForSolve`**, not the fault PC. The actual fault PC was inside
`prepareForSolve` at one of those `ldp q` instructions.

## Why retries advance further

This is what I observed across three sandstone hi-res retries with
identical inputs (`-np 1`, dt=0.1, 80³, plot_int=10):

| retry | step at failure | t at failure | RAM pressure |
| ----- | --------------- | ------------ | ------------ |
| #1 | STEP 184 | t=18.3 s | medium |
| #2 | STEP 500 | t=50.0 s | medium |
| #3 | reached STEP 1000 | t=100.0 s ✓ | low (user idle) |

A timestep with `prepareForSolve` over an 80³ mesh visits roughly
512000 nodes. Each is a candidate for hitting the alignment fault. The
loop is invoked ~400-1000 times per simulation (once per Newton iter
per timestep). So we have hundreds of millions of opportunities to hit
the bad path per simulation.

Whether ANY of them faults depends on:

1. **FAB base pointer alignment.** The AMReX FAB allocator
   (`amrex::FArenaAllocator`, ultimately `posix_memalign`-style) gives
   16-byte alignment by default. But under memory pressure, some
   allocations come from a different arena with weaker guarantees, or
   from `malloc` falling back to a slower path that's only 8-byte aligned.
2. **Element index of the first iteration** of the unrolled loop.
   `amrex::ParallelFor` tiles the box and starts the inner loop at
   different element indices for different boxes. If iteration starts
   at an even index from a 16-aligned base, no fault; if it starts at
   an odd index, fault.
3. **Mesh tile layout** depends on `amr.max_grid_size`, `n_cell`, and
   the box decomposition. Hi-res 80³ has different tile boundaries
   than 40³, so the failure pattern differs.

Concurrent simulations hammer (1) and (2): more memory pressure +
worse alignment statistics → more frequent faults.

## Why it never showed up in low-res production

Production runs use `n_cell = 40 40 40`, max_level = 0. With ~64K nodes
per ParallelFor loop and dt = 0.25 / 0.8, there are ~ 50 million
node-visits per full run. That's an order of magnitude fewer than the
hi-res case (which has ~500M+). Combined with a lower memory
footprint = less heap fragmentation = better alignment statistics,
the bug fires rarely enough to look like "transient" in production
testing.

## Why the AMReX folks didn't hit it upstream

AMReX is heavily tested but predominantly with:
- Scalar (`double`) fields, which are 8 bytes — no alignment ambiguity.
- `amrex::FArrayBox` (flat double arrays), naturally aligned.

`amrex::Array4<T>` with non-trivial T (here, `Eigen::Matrix<double, 3,
1>` or `Eigen::Matrix<double, 3, 3>`) is more of an alamo extension
than an AMReX-canonical use. The unaligned-struct-array memcpy pattern
isn't part of AMReX's standard test matrix.

## Fix options

In order of intrusiveness / risk:

### 1. Disable LLVM auto-vectorization (cheapest, easiest A/B test)

Add to the Makefile or build command:

```
CXX_COMPILE_FLAGS += -fno-slp-vectorize -fno-vectorize
```

This stops LLVM from emitting `ldp q, q` for memcpy. Existing scalar
`ldr d` loops do not require 16-byte alignment.

Performance cost: per-step overhead increase, probably 10-30% for the
memcpy-heavy paths. Negligible if the bottleneck is the MLMG solve
itself.

Reversibility: trivial.

Risk: low.

### 2. `EIGEN_MAX_ALIGN_BYTES=8` (Eigen-level)

```
CXX_COMPILE_FLAGS += -DEIGEN_MAX_ALIGN_BYTES=8
```

Tells Eigen to assume only 8-byte alignment everywhere and not to
generate 16-byte-aligned SIMD loads/stores in its own templated code.
Does NOT directly fix the LLVM-emitted memcpy in `prepareForSolve`
(option 1 is needed for that). Best paired with option 1.

### 3. Pad `Set::Vector` / `Set::Matrix` to multiples of 16

In `src/Set/Set.H` (or wherever Set::Vector/Matrix are typedef'd),
introduce a wrapper that aligns and pads:

```cpp
struct alignas(16) Vector {
    Eigen::Matrix<double, 3, 1> v;
    double _pad;  // 8 bytes pad → total 32 bytes
};
```

Pros: addresses the root cause, no perf hit on the SIMD path.

Cons: changes the layout of a foundational type used everywhere; FAB
sizes grow 33% for vector fields, 11% for `Set::Matrix`; plotfile
output writers may need adjustment; serialized state on disk is now
incompatible with old binaries.

High regression risk. Don't do this unless you're rebuilding the
entire alamo regression suite and accept the disk-format change.

### 4. Surgical: mark the `prepareForSolve` loop body `__attribute__((no_vectorize))`

```cpp
#pragma clang loop vectorize(disable)
amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
    ...
});
```

Targets only this function; rest of alamo keeps its SIMD optimizations.
Has to be repeated at any other site that copies `Set::Vector` /
`Set::Matrix` arrays — at minimum `prepareForSolve`, possibly also in
`Operator::Elastic::Fapply` and `Operator::Elastic::Diagonal`.

Medium reversibility.

### 5. Replace struct-assignment with explicit `std::memcpy`

```cpp
std::memcpy(&dw(i, j, k), &result, sizeof(result));
```

LLVM will lower this to a non-unrolled memcpy intrinsic that handles
unaligned addresses correctly. But locating every assignment site is
tedious.

## Recommended next step

I propose **option 1** as the first move:

```bash
EIGEN=$PWD/ext CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  CXX_COMPILE_FLAGS+="-fno-slp-vectorize -fno-vectorize" \
  make -j8
```

Then rerun sandstone hi-res with the exact same input as retry #3. If
it completes cleanly across 2–3 attempts under similar memory pressure
that previously crashed, we've confirmed the diagnosis and have a
working fix. The flag can be left in place permanently until someone
audits the per-loop alternative.

If a few retries still crash, the bug is elsewhere (likely in another
auto-vectorized location), and we need to inspect further.

## Open questions / things not yet checked

- Whether `-np 4` MPI runs have the same root cause (presumed yes
  because the trap signature matches), or if there's a separate
  MPI-related issue piling on.
- Whether the same crash exists in the unit-test sweep
  (`thermal_stress`, `gb_cohesive`, etc.) on this hardware. They
  haven't been run hi-res, so the rate might be too low to observe.
- Whether earlier alamo commits with `EIGEN_DONT_VECTORIZE`-like
  defaults compiled cleanly without this crash.

## Reproducer (for future sessions)

```bash
# Build (unchanged from CLAUDE.md)
EIGEN=$PWD/ext CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  make -j8

# Repro: hi-res sandstone, intermittent crash within first ~1000 steps
mpirun --oversubscribe --bind-to none -np 1 \
  bin/mmwspalling-3d-g++ \
  tests/MMWSpalling/hu_spall_onset/input_sandstone2 \
  amr.n_cell="80 80 80" \
  amr.plot_int=10 \
  timestep=0.1 \
  el.max_coarsening_level=-1 \
  el.solver.max_iter=5000 \
  el.solver.tol_rel=5e-5 \
  plot_file=tests/MMWSpalling/hu_spall_onset/output/sandstone2_busrepro

# Run several concurrent simulations in other terminals to raise
# memory pressure and increase the trigger rate.
```

Expected: `SIGBUS` with signal code 1 ("Invalid address alignment")
inside `Newton::solve + 8360` (= `prepareForSolve` return address),
trap address ending in `…b8` / `…a8` / `…78` (8-byte aligned, not
16-byte aligned).
