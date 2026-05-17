# AMR Failure Analysis — Hu Spall-Onset (Step 11b) Multi-Level Mechanics

Investigation of why the post-regrid composite (multi-AMR-level) MLMG mechanics
solve fails to converge in `tests/MMWSpalling/hu_spall_onset/input_{granite2,sandstone2}`
when run with `amr.max_level = 1`.

**Status:** root cause narrowed to the multi-AMR-level *operator algebra* in
[`src/Operator/Operator.cpp`](../../../src/Operator/Operator.cpp) and
[`src/Operator/Elastic.cpp`](../../../src/Operator/Elastic.cpp). Black-box
probing has hit its limit; the next step requires opening the operator.

## 1. Problem fingerprint

After all the operator-level fixes from the prior session (verified correct
in this session's review):

- **Single-level (`max_level=0`) runs converge** for both materials in 30–60
  MLMG iterations to ratio ~1e-6.
- **First multi-level solve (STEP 4)** — immediately after `regrid` creates
  `Level_1`, before any thermal advance has happened on it — **converges**
  identically to single-level.
- **Second multi-level solve (STEP 5)** — after `Level_1` has been advanced
  once and now has fine-resolution thermal/eigenstrain content the coarse
  level doesn't see — **fails**:
  - Residual reaches a fixed-magnitude *absolute floor* of ~7e9.
  - V-cycle has *zero contraction factor* beyond that floor (verified via
    per-iteration MLMG verbose log: ratio moves 0.08944 → 0.08943 over
    iterations 25–50).
  - AMReX repeatedly prints `MLMG: Bottom solve failed.` 3× per V-cycle.

The floor is *bit-identical* across initial guesses: `~7e9` whether starting
from `u=0`, from a warm-started displacement, with `fixed_iter`, or with any
solver tuning. That signature is **irreducible-residual** (a near-null-space
mode the V-cycle cannot eliminate), not slow convergence.

### Spatial signature

From the partial-residual plotfile
(`output/diag_sandstone2_fixedit/00005`), `res_z` dominates and concentrates
at the top of `Level_1`:

| z (Level_1 cell-center)     | L2(`res_z`) over xy slab |
| --------------------------- | ------------------------ |
| 0.090625 (CF interface)     | 9.7e8                    |
| 0.094375                    | 1.7e10                   |
| 0.096875                    | 2.4e11                   |
| 0.098125                    | 4.8e11                   |
| **0.099375 (top fine cell)**| **1.7e12**               |

The residual grows ~1700× over the 8 fine cells of `Level_1`. It is *not* at
the coarse-fine ring (z = 0.09), and *not* at the supports (zlo). It tracks
where the eigenstrain `F0_zz` is largest (peak 1.18e-3 on `Level_1` at zhi).

## 2. What has been ruled out

Each row is a falsified hypothesis with the probe and result.

| Hypothesis                                         | Probe                                            | Result                                                |
| -------------------------------------------------- | ------------------------------------------------ | ----------------------------------------------------- |
| Material-specific (sandstone ν=0.34 hard)          | Run Granite with same setup                      | Granite fails identically                             |
| Reflux is broken                                   | Bypass `Operator<Grid::Node>::reflux`            | Made it WORSE (residual blew up to 8e30); reflux helps|
| AMR-level coefficient inconsistency                | `el.solver.average_down_coeffs=1`                | 0.02% change in stall ratio — no effect               |
| Bottom solver iteration-limited                    | `el.solver.bottom_max_iter=5000`                 | Same plateau, 24× longer bottom timer                 |
| Multigrid depth                                    | `el.max_coarsening_level = 0,1`                  | Both worse than default `2` even pre-regrid           |
| Coarse-fine strategy                               | `el.solver.cf_strategy=none`                     | Catastrophic; resid0 explodes 14 orders to 6e24       |
| u=0 initial guess insufficient                     | `el.zero_out_displacement=0` warm-start          | Same ~7e9 absolute residual floor                     |
| Time evolution sharpens F0                         | `timestep=0.005` (50× smaller dt)                | Failed *worse* (ratio 0.42 vs 0.088); dt isn't trigger|
| `temp_mf` not averaged-down across levels          | yt comparison Level_0 covered vs Level_1 sub-avg | Match to 0.000K — `temp_mf` is correctly averaged     |
| Per-level Newton rhs disagrees at AMR overlap      | Throwaway `average_down(rhs_mf)` in `Newton.H`   | Confirmed real (25% max change) but doesn't help      |
| BC at zlo (rank-deficient roller pin)              | Replace `zlo_roller_321` with full clamp         | Same residual floor                                   |
| BC stencil at zhi top face (traction-free w/ F0)   | Clamp ALL 6 faces (zero traction-free DOFs)      | Same residual floor                                   |

Material, BC type, solver tuning, initial guess, time-step magnitude, reflux,
coefficient averaging, rhs averaging, coarsening depth, CF strategy: all
ruled out.

## 3. What's left

The bug is in the **multi-AMR-level operator algebra at the refinement
boundary** in
[`src/Operator/Operator.cpp`](../../../src/Operator/Operator.cpp) /
[`src/Operator/Elastic.cpp`](../../../src/Operator/Elastic.cpp), specifically
how the eigenstrain forcing `F0` is treated when Level_1 has fine-resolution
content the coarse level doesn't see.

The trigger condition: **`Level_1` has been advanced through one thermal
step.** After that one step, Level_1's local cell-temperature solver
develops a sharper near-surface profile than Level_0 ever had. The
nodal model rebuilt from that sharper temp gives a per-level `F0` whose
discrete divergence (`div(C·F0)`) on Level_1 differs from Level_0's
divergence-of-the-averaged-down-`F0`. The composite solver cannot reconcile
the two levels' equations at the AMR overlap region.

We confirmed (via diagnostic prints) that the per-level rhs values *do*
disagree by ~25%, but enforcing rhs-consistency by hand was insufficient.
That suggests the inconsistency lives deeper — in either:

- (A) `Operator::Elastic::Fapply` — the per-level `Apply` operation. With
  per-level `m_ddw_mf` (stiffness only; `F0` is *not* in `m_ddw_mf` —
  it's in the user's `model_mf` which becomes `dw` via Newton), Apply
  computes `div(C·grad u)` per-level. The consistency between levels
  depends on AMReX's CFStrategy::ghostnodes filling Level_1's coarse-fine
  ghost displacement values from Level_0 — and on the discrete operator's
  symmetry/SPD properties at the AMR boundary.
- (B) The restriction or interpolation of state between MG levels *inside*
  Level_1's V-cycle. AMReX MLNodeLaplacian has special handling for nodal
  scalar Laplacian-like operators; ALAMO's `Operator<Grid::Node>` is a
  custom subclass and may not have ported all of that machinery for the
  vector elastic operator.
- (C) The **null space of the multi-level operator under Dirichlet
  ghosts**. With `CFStrategy::ghostnodes`, Level_1's local problem is
  Dirichlet at the CF boundary (from Level_0) and free at the physical
  zhi face. When `F0_zz` has a sharp surface-layer pattern that grows
  monotonically toward zhi, there's no smooth Level_1 displacement that
  *both* (a) satisfies the local fine-grid PDE *and* (b) matches the
  coarse-level displacement that Level_0 would predict from the averaged
  forcing. The composite system is over-constrained at exactly the cell
  layer where the residual is observed.

Hypothesis (C) is consistent with *all* the evidence: it explains why
fresh-Level_1 (FillPatched smooth) works (per-level F0 effectively the
same on both), why advanced-Level_1 fails (F0 patterns diverge between
levels), why no BC choice fixes it (the conflict is between Levels'
forcing, not the boundary), why no solver tuning fixes it (the operator's
range doesn't span the full rhs), and why the residual concentrates at zhi
top (where F0 is sharpest).

## 4. Recommended path forward

### 4.1 First — confirm the operator-level diagnosis cheaply

Before changing operator code, do one more black-box test that takes ~30
minutes:

**Disable the per-level F0 evolution on Level_1.** In
[`src/Integrator/MMWSpalling.H`](../../../src/Integrator/MMWSpalling.H)
`UpdateModel`, force `F0` on `Level_1` (only) to be the FillPatch-from-Level_0
interpolated value rather than a fresh `CellToNodeAverageInDomain` of
Level_1's evolved temp. If the AMR smoke now passes, hypothesis (C) is
confirmed and the fix is localised to the per-level forcing build, not the
operator. If it still fails, the operator itself needs surgery.

This is a ~10-line throwaway edit guarded by an env var (e.g.
`ALAMO_LEV1_F0_FROM_COARSE`), parallel to the reflux-bypass and
`average_down_rhs` probes already exercised.

### 4.2 If 4.1 confirms — fix the per-level forcing build

If forcing-build is the issue, the cleanest fix is to enforce that
Level_1's `F0` is the consistent fine-resolution restriction of Level_0's
view. Two implementations to consider, in order of intrusiveness:

1. After `UpdateModel` populates per-level `model_mf`, restrict (average-up)
   from Level_0 to Level_1 at coarse-fine ghosts, and then `FillBoundary`
   so the operator sees consistent F0 across the AMR boundary. This is the
   complement of the standard "average-down for state."
2. Compute `F0` *only* on the finest level, then `average_down` to coarser
   levels. Per-level `Apply` then operates on a coherent multi-level
   forcing field by construction.

### 4.3 If 4.1 fails — open the operator

The investigation now moves into
[`src/Operator/Operator.cpp`](../../../src/Operator/Operator.cpp) and
[`src/Operator/Elastic.cpp`](../../../src/Operator/Elastic.cpp). Concrete
reading list, in order:

1. **`Operator<Grid::Node>::reflux`** ([Operator.cpp:549](../../../src/Operator/Operator.cpp#L549)).
   Read end-to-end. Note that the restriction stencils are written
   explicitly per face/edge/corner type — compare against AMReX's
   `MLNodeLaplacian::reflux` for shape and weight differences.
2. **`Operator<Grid::Node>::Diagonal`**, **`Fsmooth`**, **`applyBC`**,
   **`solutionResidual`/`correctionResidual`**, **`fixUpResidualMask`** —
   each should be checked against the corresponding AMReX
   `MLNodeLaplacian` routine. The previous session fixed several real bugs
   here (3D coloring, ghost-init, nodalSync); there may be more.
3. **`Elastic<SYM>::averageDownCoeffsDifferentAmrLevels`**
   ([Elastic.cpp:611](../../../src/Operator/Elastic.cpp#L611)) — currently
   gated behind `m_average_down_coeffs` (default false; we tested it on
   and saw no effect). Read the restriction it does for `m_ddw_mf` and
   verify it actually does what it claims, especially across the AMR
   boundary.
4. **`Elastic<SYM>::Fapply`**
   ([Elastic.cpp:147](../../../src/Operator/Elastic.cpp#L147)) — the line
   `Set::Matrix sig = (DDW(i, j, k) * gradu) * psi_avg;` confirms `F0` is
   *not* in the operator's stress computation. Confirm where `F0` enters
   the multi-level system: `Newton::prepareForSolve` builds `dw =
   model.DW(kinvar(grad u))` per level (which contains `−C·F0` at u=0)
   and feeds `b − div(dw)` as the linearised rhs to MLMG. Verify that
   this linearisation is consistent across AMR levels — that's the heart
   of the suspected bug.

### 4.4 Diagnostic experiments to run *while* reading the operator

Each takes minutes; pair them with the relevant code section.

- **Project the residual onto rigid-body modes.** The fine-level Level_1
  patch (free-floating except for CF Dirichlet ghosts) has a 6-mode null
  space if those Dirichlet values are not sufficiently strong. Compute
  the inner product of the failing residual with each of 6 unit
  rigid-body vectors. If the residual lives mostly in the null space,
  the problem is null-space deflation; if it lives in the range, the
  problem is operator inconsistency.
- **Compare `Apply(0)` per level.** Run a one-shot Apply on `u = 0` and
  inspect what each level produces. If `Apply(0)` on Level_0 and on
  Level_1 disagree at coarse-fine corresponding nodes, the operator
  itself is inconsistent (independent of forcing).
- **Synthetic constant-F0 problem.** Build a throwaway input where
  `F0_xx = F0_yy = F0_zz = 1e-3` everywhere (uniform thermal expansion,
  no spatial structure), force AMR via thermal IC, and run the multi-
  level solve. A uniform F0 has zero divergence so the operator should
  satisfy `Apply(u) = 0` with `u` matching the boundary conditions. If
  this fails, the operator has a bug independent of any fine/coarse F0
  mismatch.

### 4.5 Maintainer / community check

Before deep operator surgery, the AMReX MLMG nodal vector-elastic path may
have known limitations. Worth searching the AMReX issue tracker and asking
the ALAMO maintainers whether multi-AMR-level static elasticity with
spatially-varying eigenstrain is expected to work, and if so, whether
there's a reference test case that exercises it.

## 5. Diagnostic infrastructure already in place — reuse, don't reinvent

| Mechanism                                  | Where                                                      | Purpose                                                       |
| ------------------------------------------ | ---------------------------------------------------------- | ------------------------------------------------------------- |
| `el.print_residual = 1`                    | [Base/Mechanics.H:128](../../../src/Integrator/Base/Mechanics.H#L128) | Writes `res_mf` as plotfile vector (`res_x/y/z`)        |
| `el.solver.verbose = 2` or `3`             | [Solver/Nonlocal/Linear.H:273](../../../src/Solver/Nonlocal/Linear.H#L273) | MLMG per-V-cycle (2) or per-iteration (3) norms        |
| `el.solver.fixed_iter = N`                 | [Solver/Nonlocal/Linear.H](../../../src/Solver/Nonlocal/Linear.H) | Force exactly N V-cycles regardless of convergence; lets `res_mf` plotfile capture partial state |
| `el.solver.cf_strategy = ghostnodes\|none` | [Solver/Nonlocal/Linear.H:334](../../../src/Solver/Nonlocal/Linear.H#L334) | Toggle AMReX MLMG CF strategy                          |
| `el.solver.average_down_coeffs = 1`        | [Solver/Nonlocal/Linear.H:365](../../../src/Solver/Nonlocal/Linear.H#L365) | Enable cross-AMR coefficient averaging in `Operator::Elastic` |
| `mmw.debug_amr_mechanics = 1`              | [Integrator/MMWSpalling.H:2440](../../../src/Integrator/MMWSpalling.H) | Finite-field checks with valid / boundary-ghost / CF-ghost classification |
| `Util::RealFillBoundary`                   | various                                                    | Standard ghost-fill helper                                    |

## 6. Throwaway diagnostic inputs

All under [`tests/MMWSpalling/hu_spall_onset/diag/`](./diag/) — never touch
the real Hu inputs. Reuse for new probes; delete when the investigation
ends.

| Input file                                | Purpose                                                             |
| ----------------------------------------- | ------------------------------------------------------------------- |
| `input_granite2_base`                     | Granite AMR baseline (`max_level=1`, `regrid_int=3`, verbose, debug)|
| `input_sandstone2_base`                   | Sandstone AMR baseline                                              |
| `input_sandstone2_clamp`                  | Sandstone with zlo full Dirichlet (no roller pins)                  |
| `input_sandstone2_clampall`               | Sandstone with all 6 faces Dirichlet u=0                            |

Recommended override pattern when adding a new probe: pass overrides on the
command line rather than editing inputs in place, e.g.
`mpirun ... ./bin/mmwspalling-3d-g++ diag/input_sandstone2_base
el.solver.fixed_iter=10
plot_file=tests/MMWSpalling/hu_spall_onset/output/diag_<probe_name>`.

## 7. Empirical numbers to reproduce

Reference values for Sandstone with `regrid_int=3`, `dt=0.25`,
`stop_time=5.0`, all current operator fixes in place:

```
STEP 1 (single-level)        : rhs 4.20e10, converged in 46 iter, ratio 1.9e-6
STEP 2                       : rhs 4.21e10, converged in 48 iter, ratio 1.8e-6
STEP 3                       : rhs 4.22e10, converged in 37 iter, ratio 1.7e-6
STEP 4 (first multi-level)   : rhs 4.21e10, converged in 39 iter, ratio 2.0e-6
STEP 5 (second multi-level)  : rhs 7.89e10, FAILED at 1000 iter, ratio 0.0882
                                final residual 6.96e9 (the irreducible floor)
```

For Granite: STEP 5 fails at ratio 0.123, residual 6.4e9 — same floor.

## 8. Practical recommendation in the meantime

If `max_level=1` is not strictly required for the Hu spall-onset
validation, run single-level (`max_level=0`) to unblock the validation
campaign and document the AMR limitation in `ROADMAP.md`. The fix is real
work and may take days; single-level produces correct physics, just
without adaptive resolution.

## 9. Session 2 (2026-05-09 / 2026-05-10) findings

This session executed §4.1 (probe), §4.2 (forcing fix), part of §4.3
(operator surgery), and §4.4 (diagnostic instrumentation). Net result:
the bug is localised but not yet fixed. The key data points below are
intended to make the next investigation cheap to resume.

### 9.1 Probe (`ALAMO_LEV1_F0_FROM_COARSE=1`) — partial fix

A new helper `ApplyAMRF0Consistency()` in
[`src/Integrator/MMWSpalling.H`](../../../src/Integrator/MMWSpalling.H)
runs after `UpdateModel` and either (env-var probe) overwrites Level_1's
F0 with `InterpFromCoarseLevel(Level_0)`, or (default) does an
`average_down_nodal` fine→coarse on F0.

| Mode | STEP 5 rhs | Final residual | Ratio | Verdict |
| ---- | ---------- | -------------- | ----- | ------- |
| Baseline | 7.892e10 | 6.917e9 | 0.0876 | stalls |
| Probe (Level_1 F0 = interp coarse) | 7.751e10 | **4.758e9** | **0.0614** | -30% residual, still stalls |
| Production fix (`average_down_nodal` fine→coarse) | 7.892e10 | 6.917e9 | 0.0876 | **bit-identical no-op** |

The production-mode `average_down_nodal` is a no-op because `temp_mf` is
already correctly averaged-down between levels (the doc's earlier
finding: "Match to 0.000K"). With temp consistent, the cell→node
averaging produces nodal F0 values that are also consistent at AMR
overlap nodes — there is no fine→coarse F0 mismatch left to fix.

The probe direction (smoothing Level_1) helps because it removes the
*sharper* fine-level F0 gradient near the zhi free surface, not because
it improves AMR consistency.

### 9.2 §4.3 candidate: nodalSync in solutionResidual/correctionResidual

Pattern audit of [`src/Operator/Operator.cpp`](../../../src/Operator/Operator.cpp)
showed that every other use of `realFillBoundary` is paired with
`nodalSync` (lines 74–75, 96–97, 154–155, 385–386, 455–456, 487–488,
657–658). The two exceptions are `solutionResidual` ([line 672](../../../src/Operator/Operator.cpp#L672))
and `correctionResidual` ([line 684](../../../src/Operator/Operator.cpp#L684)).

Adding `nodalSync(amrlev, mglev, resid)` to both made convergence
slightly **worse** (ratio 0.0876 → 0.0952). Reverted. The asymmetry is
intentional: residual is computed deterministically at every node from
already-synced `apply` inputs; an extra sync at the end *averages* away
correct values at shared nodes.

Cross-check: AMReX `MLNodeLinOp::solutionResidual`
(`ext/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeLinOp.cpp:122`) doesn't
sync residual either.

### 9.3 §4.4 diagnostic: `DiagnoseAMRResidual()`

Added a method in
[`src/Integrator/MMWSpalling.H`](../../../src/Integrator/MMWSpalling.H)
that, when `ALAMO_AMR_DIAG=1` and `el.print_residual=1`, prints after
each post-mechanics-solve `TimeStepBegin`:

- per-level L2(res_mf) and max|res_mf| with location;
- per-z-slab L2 on the finest level;
- 6-mode rigid-body projection of the finest-level residual (TX, TY, TZ,
  RX, RY, RZ contributions to ‖res‖²).

Run with `el.solver.fixed_iter=N` to suppress MLMG's abort-on-failure so
the diagnostic can fire on the stalled state.

### 9.4 Diagnostic results — sandstone STEP 5, `fixed_iter=50`

**Baseline:**

| z (Level_1) | L2(res) at slab |
| ----------- | --------------- |
| 0.085 (CF interface) | 2.5e8 |
| 0.0925 | 4.6e11 (mid-L1 ring peak) |
| 0.09875 (top fine cell) | **1.38e12 (dominant peak)** |
| 0.1 (physical zhi) | 5.0e9 (drop) |

L0 ‖res‖₂ = 3.6e11; L1 ‖res‖₂ = 1.8e12.

RBM projection on Level_1: TZ 8.8%, RX 1.8%, RY 1.9%, others <0.001%.
**Total RBM = 12.5% of ‖res‖²**, so 87.5% lives in the range space.

**Probe (`ALAMO_LEV1_F0_FROM_COARSE=1`):**

| z (Level_1) | L2(res) | vs baseline |
| ----------- | ------- | ----------- |
| 0.0925 | 1.13e11 | **−75%** |
| 0.09875 | 1.39e12 | **unchanged** |
| 0.1 | 5.0e9 | unchanged |

Total RBM ≈ 13.6%.

### 9.5 Localised hypothesis

The mid-Level_1 ring at z=0.0925 was an F0-sharpness artefact (probe
removes it). The dominant peak at z=0.09875 is **invariant to F0
smoothing** — therefore not a forcing-level bug.

Looking at [`src/Operator/Elastic.cpp:200-206`](../../../src/Operator/Elastic.cpp#L200-L206):
at any node on the physical domain boundary of the *current AMR level*,
`Fapply` returns `f = (*m_bc)(u, gradu, sig, ...)` directly, replacing
the divergence-of-stress with the BC's evaluation. For traction-free
zhi, this is `sig · ẑ = sigma_zz`. With `Newton::prepareForSolve`
encoding `rhs = b - sigma_zz_initial = (λ+2μ)·F0_zz` at zhi, the V-cycle
must drive `apply(u) = (λ+2μ)·F0_zz` at every fine zhi node.

The remaining stall is not in the forcing or in `nodalSync`; it is in
how the multi-level operator drives this fine-level zhi-face equation
to zero. A working hypothesis: the BC value at fine zhi nodes
(`sigma_fine · ẑ`) is computed from a fine-stencil `gradu` whose ghost
cells are filled via AMReX `realFillBoundary` (line 473–476, alamo's
which currently ignores the geometry argument and only does
`MultiFab::FillBoundary` — i.e. internal box-to-box ghosts only, no
physical-BC ghost fill). At the *physical* zhi face on Level_1, the
ghost row `k = hi_z + 1` is therefore left at whatever `MultiFab::define`
initialised it to (zeros, typically), poisoning the fine `gradu` stencil
in a way the V-cycle cannot correct via reflux alone.

### 9.6 Hypothesis test — ghost-fill probe (partially refuted)

The `DiagnoseAMRResidual` was extended to dump `disp_mf` at Level_1's
k=80 (valid zhi top, n=1521 nodes) and k=81 (ghost above zhi, n=1521).
Result on stalled STEP 5:

| Location | L2(disp) | max\|disp_z\| |
| -------- | -------- | ------------- |
| k=80 (valid zhi top) | 0.00333 | 1.37e-4 |
| k=81 (ghost above zhi) | 0.00153 | 6.2e-5 |

The ghost row is **not** zero — alamo's `RealFillBoundary` (or upstream
operations like `ParallelCopy` between MG levels) leave some data
there. It's also not a proper traction-free mirror (which would need
`u(k=81) = u(k=79) + 2·dz·F0_zz`). The simple "ghost is unset" version
of 9.5 is wrong; whatever fills the ghost is doing so non-deterministically
relative to the BC.

A second observation that complicates the diagnosis: revisiting the
z-slab residual table, the peak at **z=0.09875** is at **nodal k=79**
(one node *below* the physical zhi at k=80), not at the boundary node
itself. The boundary node k=80 is handled by the BC short-circuit path
in [`Fapply` (Elastic.cpp:200-206)](../../../src/Operator/Elastic.cpp#L200-L206)
which evaluates `sigma·n` directly and doesn't use ghost values. The
bug is in the **interior stencil at k=79**, which reads `u(k=78,79,80)`
— all valid data, no ghost involvement. So the persistent residual
isn't a boundary-ghost issue at all.

### 9.7 Fapply per-level dump — root cause localized

`Elastic::Fapply` was instrumented (env-var-guarded
`ALAMO_FAPPLY_DUMP=1`) to log `u`, `gradu`, and `f` at probe nodes
(20,20,39|40) on Level_0 and (40,40,78|79|80) on Level_1 — physical
column at the surface_patch center, near the zhi boundary. STEP 5 with
`fixed_iter=2`. Key data at one representative iteration:

```
amrlev=0 (20,20,39)  u_z=1.99e-5  gradu_z[2,2]=0.00330  f_z=4.07e10
amrlev=0 (20,20,40)  u_z=3.29e-5  gradu_z[2,2]=0.00521  f_z=1.47e8 (boundary BC)
amrlev=1 (40,40,78)  u_z=1.99e-5  gradu_z[2,2]=0.00237  f_z=2.92e10
amrlev=1 (40,40,79)  u_z=2.37e-5  gradu_z[2,2]=0.00521  f_z=8.98e10  ← residual peak
amrlev=1 (40,40,80)  u_z=3.29e-5  gradu_z[2,2]=0.00737  f_z=2.03e8 (boundary BC)
```

Interpretation:

1. The V-cycle correctly synchronises `u` at coincident coarse-fine
   nodes: `u_lev0(20,20,39) == u_lev1(40,40,78) == 1.99e-5` and
   `u_lev0(20,20,40) == u_lev1(40,40,80) == 3.29e-5`. So the
   coarse↔fine interpolation/restriction of `u` is fine.
2. **At the same physical point z=0.0975**, the operator computes
   `f_z = 4.07e10` on Level_0 and `f_z = 2.92e10` on Level_1.
   These **disagree by ~30%** because the discrete stencil for
   `gradu_z` is different on the two levels (coarse reads
   `u(k=38), u(k=40)` over `2·dx_coarse=0.005`, fine reads
   `u(k=77), u(k=79)` over `2·dx_fine=0.0025`). The linear-interpolated
   intermediate value `u_lev1(40,40,79)=2.37e-5` is not the average of
   the two coincident coarse values (which would be 2.64e-5) — the
   V-cycle has shifted it, and the fine stencil's `gradu_z` reflects
   that adjusted value.
3. The peak residual lives at the **intermediate fine node**
   `(40,40,79)` (z=0.09875), which has **no coincident coarse node**.
   This node is filled by the V-cycle's linear interpolation between
   the two coincident coarse values, then only Level_1's Jacobi
   smoother can adjust it. Coarse-level smoothing has no direct way to
   reduce `f_z` at this specific node.

### 9.8 Diagnosis

The remaining floor is **a fundamental multigrid smoother limitation,
not a code bug**:

- The discrete operator must give different values at coincident
  coarse-fine nodes (different stencil width) — that's expected for
  any nodal AMR multigrid.
- AMReX MLMG handles this by *reflux*: the fine residual is restricted
  to the coarse residual at coarse-fine boundary nodes, so the coarse
  level "sees" the fine residual after reflux. Alamo implements this
  in `Operator<Grid::Node>::reflux` ([Operator.cpp:549](../../../src/Operator/Operator.cpp#L549)),
  and the reflux logic is correct.
- The problem is that **intermediate fine nodes** (e.g. `(40,40,79)`)
  have no coarse counterpart and are filled by interpolation, then
  must be relaxed by fine-level smoothing alone. With sharp F0 at the
  free surface, the local rhs at intermediate fine nodes is large,
  and Jacobi (alamo's default smoother — see
  [`Operator<Grid::Node>::Fsmooth`](../../../src/Operator/Operator.cpp#L100))
  is too weak to drive these to zero in 1000 V-cycles.

This matches the rigid-body-mode projection observation in §9.4: 12.5%
of `‖res‖²` is in RBM (translations + rotations), 87.5% in the range
space. The range-space residual is the discrete "what fine Jacobi
can't smooth" component.

### 9.9 Suggested fix paths (multi-day work, not attempted)

Pick one in increasing order of intrusiveness:

1. **Stronger fine smoother**. Replace alamo's pure Jacobi
   ([`Operator<Grid::Node>::Fsmooth`](../../../src/Operator/Operator.cpp#L100))
   with a red-black Gauss-Seidel or line-Jacobi smoother in z. This is
   the most common fix for stalled multigrid at high-frequency
   boundary modes. AMReX has SOR-style smoothers in MLNodeLaplacian
   that can be ported. ~2-3 days.
2. **Increase smoothing iterations near the boundary**. Add a `presmooth`
   / `postsmooth` count override in `Solver::Nonlocal::Linear` and
   parametrize it in the input. The cost is per-V-cycle time, but if
   the issue is just slow convergence, more smoothing helps. ~half a
   day.
3. **Change boundary stencil at zhi free surface to use a one-sided
   second-order derivative**. The current centered stencil at k=79
   reads `u(k=78,79,80)` where `u(k=80)` is the boundary node value
   (large because of thermal expansion). A higher-order one-sided
   stencil might be better-conditioned. ~1 day, risk of regressions on
   other tests.
4. **Avoid the issue at the system level**: don't let AMR refinement
   reach the zhi physical boundary. Configure refinement to leave a
   buffer of coarse cells between Level_1 and the physical free
   surface. This eliminates the intermediate-fine-near-boundary mode
   entirely. ~1 day; physically reasonable for the Hu spall-onset
   problem.

**Practical recommendation:** option 4 (refinement buffer at zhi) is
the most pragmatic for unblocking AMR in this validation. The bug isn't
in alamo's elastic operator per se; it's in the interaction of nodal
multigrid + traction-free physical boundary + sharp surface forcing,
which Jacobi smoothing can't resolve. Avoiding refinement at the free
surface sidesteps this entirely.

Also worth running but not yet attempted (cheaper diagnostics for
later sessions):

- `apply(0)` per-level dump (single Apply on `u=0` everywhere).
- Synthetic uniform-F0 experiment (§4.4 of this doc).

### 9.7 Files changed in session 2

| File | Change | Status |
| ---- | ------ | ------ |
| [`src/Integrator/MMWSpalling.H`](../../../src/Integrator/MMWSpalling.H) | Added `ApplyAMRF0Consistency()` (probe + average-down) and `DiagnoseAMRResidual()`. Hooked into `UpdateModel` and `TimeStepBegin`. | Kept (no-op default; probe + diag gated on env vars) |
| [`src/Operator/Operator.cpp`](../../../src/Operator/Operator.cpp) | Tested `nodalSync` in solutionResidual/correctionResidual. Hurt convergence. | **Reverted** |
| [`tests/MMWSpalling/hu_spall_onset/AMR_FAILURE_ANALYSIS.md`](AMR_FAILURE_ANALYSIS.md) | This section 9. | Kept |

### 9.8 Reproducer commands

```bash
# Build (-np 1 only; -np 4 has a separate bus error not yet investigated)
EIGEN=$PWD/ext \
  CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  make -j8

# Baseline failure
mpirun --oversubscribe --bind-to none -np 1 \
  bin/mmwspalling-3d-g++ \
  tests/MMWSpalling/hu_spall_onset/diag/input_sandstone2_base \
  el.solver.verbose=2

# Probe (smoothed Level_1 F0): partial fix, ratio 0.061
ALAMO_LEV1_F0_FROM_COARSE=1 mpirun --oversubscribe --bind-to none -np 1 \
  bin/mmwspalling-3d-g++ \
  tests/MMWSpalling/hu_spall_onset/diag/input_sandstone2_base \
  el.solver.verbose=2

# Diagnostic on the stalled residual
ALAMO_AMR_DIAG=1 mpirun --oversubscribe --bind-to none -np 1 \
  bin/mmwspalling-3d-g++ \
  tests/MMWSpalling/hu_spall_onset/diag/input_sandstone2_base \
  el.solver.verbose=2 el.solver.fixed_iter=50
```

## 10. Session 3 (2026-05-10) — five hypotheses refuted, doc retraction

This session followed the Session 2 plan (probe-then-RBM-deflation),
plus a new "smoother bandwidth via input override" Phase 0. **Net
result: every actionable hypothesis from the prior plan was refuted.**

### 10.1 Phase 0 — smoother bandwidth (no source change)

| `el.solver.pre_smooth = post_smooth` | STEP 5 ratio | Wall |
| ------------------------------------ | ------------ | ---- |
| 2 (default)                          | 0.08764      | 23.3s |
| 4                                    | 0.08770      | 37.3s |
| 8                                    | 0.08773      | 58.6s |
| 16                                   | 0.08774      | 92.3s |

The ratio is constant within rounding across an 8× increase in
fine-level Jacobi sweeps. This **refutes the §9.4 working hypothesis
that pure Jacobi was too weak** — the offending mode has effective
smoothing factor μ_J ≈ 1.0 (null-space-like behaviour), not just slow
contraction.

### 10.2 Probes — bottom solver and CF interface placement

| Probe | STEP 5 ratio | Wall |
| ----- | ------------ | ---- |
| `el.solver.bottom_solver=bicgstab` | 0.08765 | 105s |
| `el.solver.bottom_solver=smoother` | 0.08764 | 10.7s |
| `amr.n_error_buf=4` (vs 1)         | 0.08608 | 88s |

Bottom solver type is irrelevant (CG and BiCGStab both fail at
iteration limit on the coarsest matrix; `smoother` removes the
warning but doesn't move the ratio). Bigger refinement buffer gives
1.8% — discretization noise, not a fix. Combined with the persistent
`Bottom solve failed` print, this confirms the **coarsest matrix is
genuinely near-singular in the AMR composite**.

### 10.3 RBM deflation experiment — catastrophic divergence (refutes §9.5 leading hypothesis)

Implemented Gram-matrix RBM deflation in
[`Solver::Nonlocal::Newton::solve`](../../../src/Solver/Nonlocal/Newton.H)
(env-var `ALAMO_RBM_DEFLATE=1`), with `el.solver.nriters=5` so the
deflation feeds back into subsequent MLMG calls within a single
timestep. STEP 5 trace:

| Newton iter | rhs        | norm(ddisp) |
| ----------- | ---------- | ----------- |
| 1           | 7.9e10     | 1.1e-4      |
| 2           | 1.7e11     | 1.3e-3      |
| 3           | 9.5e12     | 0.43        |
| 4           | 3.2e15     | 142         |
| 5           | 1.0e18     | 46037       |

The deflation **diverges** in the multi-level case (single-level STEP 1
did not diverge). The deflation code was reverted after this experiment
to keep the tree clean per CLAUDE.md guidance.

**Why the hypothesis was wrong:** The `zlo_roller_321` BC with corner
pins at xloylozlo and xhiyhizlo **already constrains all 6 RBM** in
the underlying operator. The composite operator has **no actual null
space**. The 12.5% "RBM contribution" measured in §9.4 was the
coincidental shape of the unconverged residual, not a true null mode.

When the experiment removes RBM from `u`, it destroys legitimate
displacement that the multi-level system needs. In particular, the
Level_1 patch's CF Dirichlet ghosts (interpolated from Level_0)
contain spatial content with a non-trivial RBM projection, because
the matched-CF condition `u_fine(CF interface) = u_coarse(CF
interface)` IS RBM-content-bearing. Stripping RBM from Level_1's `u`
violates this matching, and the next iteration's MLMG sees a
catastrophically wrong rhs.

The *single-level* test (STEP 1) didn't diverge because there is no
CF coupling there. This is consistent with the explanation.

### 10.4 Combined diagnosis after Session 3

| Hypothesis | Tested | Result |
| ---------- | ------ | ------ |
| Smoother bandwidth | Yes (10.1) | refuted |
| Smoother type (Jacobi vs RB-GS) | Theoretical (10.1) | argued against |
| Bottom solver type | Yes (10.2) | refuted |
| CF interface placement | Yes (10.2) | refuted (1.8% noise) |
| Forcing-level F0 inconsistency | Yes (§9.1) | partial — gives 30%, not full |
| RBM / null-space deflation | Yes (10.3) | refuted (catastrophic) |

What we know:
- V-cycle contraction factor is ≈ 1.0 for the offending mode.
- Coarsest grid matrix is near-singular in the AMR composite.
- The mode is **not** in the operator's null space (RBM deflation
  destabilises rather than helps).
- That combination is consistent with a **CF-coupling-induced
  near-null mode in the multi-level composite operator** that does
  not exist in either level's standalone problem.

What we have not investigated (next-session candidates from §4.3):

- **(A) Per-level operator coefficient consistency** at the AMR
  boundary. The Fapply-per-level dump in §9.7 showed expected
  discretization-induced disagreement at coincident nodes. Reflux is
  correct. The remaining question is whether the resulting
  multi-level coarse-grid operator is a sufficiently good
  approximation of the fine-level operator that the V-cycle's
  contraction factor stays bounded away from 1.
- **(B) Vector-elastic-specific multigrid machinery in
  alamo's `Operator<Grid::Node>`** vs AMReX `MLNodeLaplacian`.
  Specifically Galerkin-vs-rediscretised coarse-grid operator
  construction for vector problems with CF-Dirichlet coupling.

Both are multi-day surgery with regression risk.

### 10.5 Session 3 disposition

- Diagnostic infrastructure (§9 helpers) **kept**, env-var-guarded.
- RBM-deflation code **reverted** after experiment.
- [`AMR_SMOOTHER_LIMITATION.md`](AMR_SMOOTHER_LIMITATION.md) **rewritten**
  to retract the §9.5 smoother-weakness diagnosis and document the
  five negative results.
- Production inputs remain at `amr.max_level=0`. Hu spall-onset
  validation continues at single resolution.
- Per user direction: investigation halted; return control to Step 12
  (Mie-Gruneisen EOS).

### 10.6 Reproducer commands (session 3 specific)

```bash
# Phase 0 smoother sweep (already-merged plumbing, no source change)
mpirun --oversubscribe --bind-to none -np 1 \
  bin/mmwspalling-3d-g++ \
  tests/MMWSpalling/hu_spall_onset/diag/input_sandstone2_base \
  el.solver.verbose=3 el.solver.fixed_iter=200 \
  el.solver.pre_smooth=N el.solver.post_smooth=N
# (try N ∈ {2, 4, 8, 16})

# Probe — bottom solver swap
mpirun --oversubscribe --bind-to none -np 1 \
  bin/mmwspalling-3d-g++ \
  tests/MMWSpalling/hu_spall_onset/diag/input_sandstone2_base \
  el.solver.verbose=3 el.solver.fixed_iter=200 \
  el.solver.bottom_solver=bicgstab

# Probe — bigger refinement buffer
mpirun --oversubscribe --bind-to none -np 1 \
  bin/mmwspalling-3d-g++ \
  tests/MMWSpalling/hu_spall_onset/diag/input_sandstone2_base \
  el.solver.verbose=3 el.solver.fixed_iter=200 \
  amr.n_error_buf=4

# RBM deflation experiment (code reverted; reinstate from session 3 git history if reproducing).
```
