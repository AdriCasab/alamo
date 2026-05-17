# Multi-Level Mechanics Solve Stalls — Post-Experiment Diagnosis

This is a companion to
[AMR_FAILURE_ANALYSIS.md](AMR_FAILURE_ANALYSIS.md). That document is a
chronological investigation log; this one is the **executive summary**
of three sessions of work, current as of 2026-05-10.

> **Earlier draft of this document (2026-05-09) confidently named pure
> Jacobi smoothing as the root cause. That diagnosis turned out to be
> wrong, and is retracted below.** The real diagnosis is more
> uncomfortable: we have systematically refuted five of the most
> obvious hypotheses, and the surviving candidates require deep
> operator surgery.

## TL;DR

- **Implementation correctness checks are all green.** Alamo's elastic
  operator, Newton/MLMG forcing construction, reflux logic, BC
  treatment, and input-file parameters all check out against AMReX
  conventions and the documented physics.
- **The composite multi-AMR-level static-elastic V-cycle has a
  contraction factor ≈ 1.0 for some specific mode.** Residual ratio
  plateaus at ~0.087 across 200+ V-cycles regardless of how we tune
  the standard knobs.
- **Five hypotheses ruled out by experiment.** Smoother bandwidth,
  smoother type (theoretical), bottom-solver type, CF interface
  placement, and rigid-body-mode deflation. Each was the most natural
  next thing to test; none moved the needle (deflation actively
  diverged).
- **The single-level (`amr.max_level=0`) case converges normally** and
  is the production configuration today. Hu spall-onset validation is
  unblocked at single resolution.
- **Remaining candidates require multi-day operator surgery** with
  meaningful regression risk: per-level operator coefficient
  consistency, or vector-elastic-specific multigrid restriction /
  interpolation machinery that alamo's `Operator<Grid::Node>` may not
  have fully ported from AMReX `MLNodeLaplacian`.

## Setup

The Hu spall-onset test heats a circular patch at the top of a
0.1 m × 0.1 m × 0.1 m cube and drives a thermoelastic response with a
`static` elastic solve every step. With `amr.max_level = 1` and
refinement triggered by the near-surface temperature gradient, AMReX
creates a Level_1 patch covering the upper portion of the domain
(z ∈ [~0.085, 0.1]).

```
            z=0.1  (physical zhi, traction-free)
   +======================+    ← Level_1 zhi face
   |                      |
   |     Level_1 patch    |    fine cells, dx = 0.00125
   +----------------------+    ← coarse-fine interface, z ≈ 0.085
   |                      |
   |     Level_0 only     |    coarse cells, dx = 0.0025
   |                      |
   +======================+
            z=0.0  (zlo, roller_321 BC)
```

The composite MLMG solve must satisfy elastic equilibrium
$\nabla \cdot \boldsymbol\sigma = \mathbf{0}$ on both levels with
consistent CF coupling, where
$\boldsymbol\sigma = \mathbb{C}:(\nabla\mathbf{u} - \mathbf{F}_0)$ and
$\mathbf{F}_0$ is the prescribed thermal eigenstrain. The first
multi-level solve (immediately after regrid, before Level_1 has been
advanced) converges normally. The **second** multi-level solve fails
to converge.

## Five experiments, five negative results

All numbers below from `tests/MMWSpalling/hu_spall_onset/diag/input_sandstone2_base`,
STEP 5 final residual ratio with `el.solver.fixed_iter=200`.

### 1. Smoother bandwidth (Phase 0 of the fix plan)

Raised the outer pre-smooth / post-smooth count via input key:

| `el.solver.pre_smooth = post_smooth` | Final ratio | Wall (s) |
| ------------------------------------ | ----------- | -------- |
| 2 (default)                          | **0.08764** | 23.3     |
| 4                                    | 0.08770     | 37.3     |
| 8                                    | 0.08773     | 58.6     |
| 16                                   | 0.08774     | 92.3     |

Doubling, quadrupling, octupling Jacobi sweeps does **literally
nothing** — ratio is constant within the 4th digit. This is consistent
with the operator's residual living in a mode whose Jacobi smoothing
factor is exactly 1.0 — *not* the "Jacobi is too weak, GS would be
better" regime where more sweeps would help even slightly.

### 2. Smoother type (theoretical, based on #1)

Pure Jacobi having effective μ_J = 1 for the offending mode means the
mode is in (or very near) the iteration-matrix's null space. Red-black
Gauss–Seidel has the same problem for true null-space modes; GS_RB
beats Jacobi only for high-frequency *non-null* modes. So based on the
#1 result, we predicted GS_RB would not help either, and skipped
implementing it.

This is the part of the earlier draft of this doc that was **wrong**
in the original framing. The literature on Jacobi-vs-GS-vs-Polynomial
smoothers (Yang & Henson 2004; Adams et al. LLNL-JRNL-473191; Kushida
2024; AMReX MAESTROeX) genuinely does say GS_RB usually beats Jacobi
*for high-frequency boundary modes that aren't null-space-like*. Our
mode is null-space-like, so the literature doesn't apply.

### 3. Bottom solver type

`el.solver.bottom_solver=bicgstab` (vs default `cg`):

| Bottom solver | Final ratio | "Bottom solve failed" |
| ------------- | ----------- | --------------------- |
| `cg` (default) | 0.08764    | 3× per V-cycle |
| `bicgstab`     | 0.08765    | 3× per V-cycle |
| `smoother`     | 0.08764    | none (smoother can't fail) |

CG and BiCGStab both fail at the iteration limit on the coarsest grid.
That message means the coarsest matrix is genuinely near-singular —
both Krylov methods hit a (near-)null direction they can't deflate.
Switching to a smoother-only bottom solver removes the warnings but
does not move the ratio. The bottom solver is **not** the limiter.

### 4. CF interface placement

`amr.n_error_buf=4` (vs default 1) — pushes Level_1 patches 3 extra
cells away from the steepest gradient region:

| `n_error_buf` | Final ratio | Wall (s) |
| ------------- | ----------- | -------- |
| 1 (default)   | 0.08764    | 23.3 |
| 4             | 0.08608    | 88.1 |

Marginal 1.8% improvement — not enough to attribute to anything beyond
discretization noise. The CF interface placement is **not** the
limiter.

### 5. Rigid-body-mode deflation

The diagnostic in
[`DiagnoseAMRResidual()`](../../../src/Integrator/MMWSpalling.H)
showed 12.5% of `‖res‖²` projected onto the 6 rigid-body modes
(mostly TZ at 8.8%, RX/RY at ~1.8% each). The hypothesis was that this
RBM content was a near-null mode that no Krylov / smoother / coarse
correction could deflate.

We implemented Gram-matrix RBM deflation in `Newton::solve`
(env-var-guarded `ALAMO_RBM_DEFLATE=1`) and tested with
`el.solver.nriters=5`. **Result: catastrophic divergence** in the
multi-level case:

| Newton iter | rhs        | norm(ddisp) |
| ----------- | ---------- | ----------- |
| 1           | 7.9e10     | 1.1e-4      |
| 2           | 1.7e11     | 1.3e-3      |
| 3           | 9.5e12     | 0.43        |
| 4           | 3.2e15     | 142         |
| 5           | 1.0e18     | 46037       |

The deflation code was reverted after this experiment.

**Why this refuted the hypothesis:** the BCs (`zlo_roller_321` with
corner pins at xloylozlo and xhiyhizlo) **already constrain all 6
rigid-body modes** in the underlying operator. The composite operator
has **no actual null space**. The 12.5% "RBM contribution" in the
diagnostic was the coincidental shape of the unconverged residual —
not a true null mode of the discrete system.

When we deflate `u` against RBM, we destroy legitimate displacement
content that the multi-level system needs (especially the CF Dirichlet
coupling — Level_1's `u` has to *match* Level_0's at the CF interface,
and that matching has RBM-like spatial content). The next iteration's
rhs is then huge because `u` is far from any equilibrium, and the
system blows up. Single-level (STEP 1) didn't blow up — there's no CF
coupling there — which is consistent with this explanation.

### What the data actually says

Putting items #1, #3, and #5 together:

- The V-cycle has effectively **zero contraction** for the offending
  residual mode (#1).
- The coarsest grid matrix is **genuinely near-singular** in the AMR
  composite (#3).
- The mode is **not** in the operator's null space (#5 — deflating
  destabilises).

This combination is consistent with the **multi-level composite
operator** having a near-null mode that does *not* exist in either
level's standalone operator. That mode is a CF-coupling artifact, not
a property of the discrete elastic equations on either level.

## What's left to investigate

Two candidates from
[AMR_FAILURE_ANALYSIS.md §4.3](AMR_FAILURE_ANALYSIS.md), neither yet
attempted:

- **(A) Per-level operator coefficient consistency** at the AMR
  boundary. The Fapply-per-level dump in §9.7 showed the operator
  gives different `f` at coincident coarse-fine nodes (4.07e10 vs
  2.92e10 at z=0.0975). This is *expected* discretization, not a bug
  per se — different stencil widths give different finite-difference
  values. AMReX MLMG handles this via reflux, which alamo implements
  correctly. But the *strength* of the discretization-induced
  inconsistency relative to the multi-level operator's smallest
  eigenvalue may be what produces the near-null composite mode. Fixing
  this requires either (i) modifying the stencils to produce more
  compatible per-level operators, or (ii) adding explicit deflation of
  the specific multi-level near-null mode.

- **(B) Vector-elastic-specific multigrid machinery** that alamo's
  `Operator<Grid::Node>` may not have ported from AMReX
  `MLNodeLaplacian`. AMReX has special handling of the dense nodal
  stencil (the "nodal scalar Laplacian" lineage); alamo's custom
  vector-elastic subclass may be missing some of that. Specifically
  worth checking: how `MLNodeLaplacian` constructs its coarse-grid
  operators (Galerkin vs rediscretised) for vector problems with
  CF-Dirichlet coupling at non-trivial physical boundaries.

Either is genuinely multi-day work with regression risk on every
elastic test. The single-level configuration produces correct physics
without it.

## Practical recommendation

Run production at `amr.max_level=0`. The Hu spall-onset validation
campaign is unblocked at single resolution. AMR for surface-driven
thermoelastic spallation problems is an open issue tracked in this
file pair and `AMR_FAILURE_ANALYSIS.md`.

The diagnostic infrastructure added during these sessions is kept in
the tree (env-var-guarded, default off) for the next investigator:

- [`DiagnoseAMRResidual()` in src/Integrator/MMWSpalling.H](../../../src/Integrator/MMWSpalling.H)
  — env var `ALAMO_AMR_DIAG=1` + `el.print_residual=1`. Reports
  per-level L2(res), per-z-slab L2 on the finest level, 6-mode RBM
  projection, and a boundary-ghost probe at zhi.
- [`ApplyAMRF0Consistency()` in same file](../../../src/Integrator/MMWSpalling.H)
  — env var `ALAMO_LEV1_F0_FROM_COARSE=1`. Forces Level_1's F0 to be
  the coarse-interpolated version (~30% residual reduction; useful for
  isolating F0-sharpness effects from operator effects).
- See [AMR_FAILURE_ANALYSIS.md §9](AMR_FAILURE_ANALYSIS.md) for the
  full chronological log including reproducer commands.

## In one sentence

Three sessions of investigation have ruled out five plausible causes
without finding the root; the implementation is correct everywhere we
audited; the surviving hypothesis is a CF-coupling-induced near-null
mode in the composite multi-AMR-level operator that requires deep
operator surgery to address, and the campaign continues at
`max_level=0` in the meantime.
