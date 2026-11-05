# MMWSpalling code census

Written 2026-10-08 from a complete read of `src/Integrator/MMWSpalling.H` (11,033 lines at
the start of the read) and `src/Integrator/MMWSpalling/Removal.H` (3,654 lines). The
`/implement` session for P2-2b was editing both files while this was written, so they
ended the day about 80 lines longer. Every line number below is therefore approximate.
Find things by function name; the numbers only say roughly where to look.

This document changes no code. It is the map the restructure will work from.

---

## How to read this document

The document has four parts.

**Part A** explains what the integrator does and how one timestep flows through it. Read
it first; everything else refers back to it.

**Part B** is the list of the biggest problems, ranked. Each one says what the problem is,
where it lives, why it matters, what fix is proposed, and how we would know the fix is
safe. If you read only one part, read this one.

**Part C** is the tour of the code: every function in MMWSpalling.H grouped by subsystem,
then Removal.H broken into its stages, then the data members, then every input key with a
verdict on whether it is live.

**Part D** is the proposed target layout and the order in which to get there.

A few terms used throughout:

- A **plotfile** is the AMReX output directory a run writes every N steps. It contains every
  registered field with the `plot = true` flag, in the order the fields were registered.
- **thermo.dat** is the one-row-per-step text file of scalar diagnostics. Its columns are
  the "integrated variables" the code registers, in registration order.
- **The gate** (or **the identity sweep**) is the test under
  `tests/MMWSpalling/studies/tree_sweep_0922/`. It runs a fixed set of inputs and hashes
  every one of the 2,639 plotfiles and thermo.dat files they produce. If all 2,639 hashes
  match the saved list, the code's behaviour has not changed. It is the oracle for the
  whole restructure.
- A **column** is one (i, j) stack of cells. Almost all the surface machinery works per
  column, with a host-side vector of length `ncols` holding one number per column.
- **Reduce** means an MPI reduction: combining a per-rank value into one value every rank
  holds. The code does this by hand for every per-column vector and every diagnostic scalar.
- The **criterion** is the spallation test `Sp = K_I / K_Ic >= 1`: does the stress intensity
  at a surface flaw exceed the rock's fracture toughness at that temperature.

---

# Part A. What the integrator is

## A.1 The class

`Integrator::MMWSpalling` inherits virtually from `HeatConduction` (which supplies the
temperature field `temp_mf` with its boundary conditions and ghost cells) and from
`Mechanics<Affine::Isotropic>` (which supplies the elastic solve). On top of these it adds:

- an **enthalpy** field `H_mf` as the conserved thermal state, with temperature recovered
  from enthalpy every step through a table (this is how latent heat is handled);
- a **beam source** (Gaussian or flat-top, Beer-Lambert attenuated) and surface losses;
- a **mineral microstructure**: per-cell phase id, grain id, and per-cell copies of the
  thermal and elastic properties;
- a **surface patch**: the flame-jet boundary condition with its stand-off, feet, wall-jet
  enthalpy march and side-face heating;
- the **spall criterion** and the **removal pipeline** that turns a firing into voided cells;
- a large set of **diagnostics**: an energy ledger, a beam closure, and three CSV logs.

The removal pipeline is in Removal.H, which is not a header in the usual sense. It is
`#include`d inside the class body at the position of its one function, so it is a class-
body fragment. That is the result of an earlier "R1 refactor" that cut the file in two
without changing anything else.

## A.2 The four ways heat is conducted

`Advance` picks one of four conduction kernels, and a run never mixes them:

| Kernel | Chosen when | Handles voids | Handles the surface patch | Handles side faces | Used by |
|---|---|---|---|---|---|
| `AdvanceConstAlpha` | no `material` block and no microstructure | no | no | no | the Step 2-4 regressions (stefan, beam, equilibrium) |
| `AdvanceMaterial` | `material.type` given, `solver = explicit` | yes | prescribed-T or constant-h Robin | no | Zhang and Hu thermal validation |
| `AdvanceMaterialImplicit` | `material.type` given, `solver = implicit` | no | no | no | Zhang implicit |
| `AdvanceMicrostructure` | `microstructure.nphases` given | yes | every form | yes | every Meier, Kant and study run |

Only the fourth matters for the physics programme. The first three are kept for validation
tests and share about sixty lines of copy-pasted setup with it (see issue B6).

## A.3 One timestep, in order

This is `Advance` with the calls it makes unrolled. The order is not arbitrary; several
bugs in the record came from a consumer running before its producer.

1. **`BuildSurfaceTColumns`** reconstructs one surface temperature per column from the top
   two cell averages, assuming an exponential profile. Both criterion sites and the
   Robin boundary condition read this one number, which was the point of packet P0b.
2. **`UpdateWallLawState`** computes the per-step state of the derived wall-jet law
   (bore radius, gas state, duct coefficient). It is a no-op unless
   `surface_patch.side_face_h_law = derived`.
3. **`UpdatePropertyFields`** rewrites the per-cell thermal and elastic properties from the
   Grimsel temperature laws. It is a no-op unless a phase sets `law = grimsel`.
4. **`Mechanics::Advance`** runs the elastic solve. Every Meier and Kant deck sets
   `el.type = disable`, so on the production path this does nothing and the stress is
   reconstructed algebraically from temperature instead.
5. **`UpdateDamageAfterMechanics`** and **`UpdateCohesiveAfterMechanics`** are the Step 8
   and Step 9 damage-law machinery. Both are off on the production path.
6. **`UpdateSpAfterMechanics`** is **criterion site A**. It writes the field `Sp_field` at
   the top cell of every column. This is what the Kant onset validation scores.
7. **`UpdateSpClusters`** labels connected firing regions on the top plane. Only used when
   `weibull.detach_mode = connected_cluster`, which two inputs set.
8. **`UpdateRemovalAfterCohesive`** is all of Removal.H and contains **criterion site B**:
   the depth scan that decides how thick a flake is, the winner between spall, vapour and
   bit, the void flip, and the enthalpy remap that takes the flake's heat out of the
   column.
9. **`PerCellRemoval`** flips mid-column cells with `Sp >= 1`. The key is set by no input.
10. The enthalpy and temperature fields are **swapped** into their `_old` copies.
11. **One conduction kernel** runs. Inside `AdvanceMicrostructure` this also runs, in order:
    `BuildFlameColumns` (per-column flame coefficient and gas temperature), which itself
    calls `FeetDescent` (where the burner sits), `JetEnthalpyMarch` (how the gas cools as
    it spreads over the floor) and `WallStepDiagnostic`; then `BuildPinEffective` (the
    pinned surface closure); then the per-cell kernel; then `BeamEnergyBalance`,
    `EnergyLedger` and `PatchDiagnostics`.
12. **`InvertHtoT`** recovers temperature from enthalpy everywhere.

One consequence worth stating plainly: the criterion (step 8) runs before the flame is
rebuilt (inside step 11). Anything the criterion reads from `flame_h`, `flame_Tg` or
`side_T_eff` is therefore one step stale. The comments rediscover this at four separate
places and a unit test (`lateral_flux_lock`) pins the sites that must agree.

## A.4 The three kinds of state

Understanding the state explains half the bugs in the record.

**AMReX fields.** About 45 of them, declared as `Set::Field`, registered in `Parse`,
checkpointed, and plotted when registered with `true`. Thirty-four are plotted and so
appear in the hashed plotfiles.

**Per-column host vectors.** About 35 of them, each a `std::vector` or `amrex::Vector`
of length `ncols`. Some are rebuilt every step (`flame_h`, `flame_Tg`, `flame_s`,
`pin_T_eff`, the seven `surf_*_col` arrays); some are carried across steps (`pin_T`,
`pin_t`, `surf_x_hold`, the D2z clock `col_t_last`/`col_dt`/`col_nvoid`, the lateral
accumulator `lat_accum_col`, `foot_clear_k`, `jet_fs_z0`, `wall_corner`). **None of these
is checkpointed.** A restarted run silently starts them from zero. The comments say so at
each declaration; nothing enforces it and no test restarts a run.

**Thermo scalars.** About 130 of them, registered as integrated variables with
`extensive = false`, which means the Integrator does not reduce them: every one is
reduced by hand before being assigned. The record contains at least four incidents of
getting this wrong in one direction or the other (a per-column sum reduced when it was
already rank-identical, giving 4x on 4 ranks; or a per-box sum not reduced). The rule
"reduce what a per-box loop accumulated, never what a per-column loop did" is written in
comments six times and enforced nowhere.

## A.5 The two criterion sites and the four side-face sites

Two patterns recur so often that they deserve naming here.

**Two criterion sites.** The same fracture test is implemented twice. Site A
(`UpdateSpAfterMechanics`) evaluates `Sp` at the surface of every column and writes a
field; it is what the Kant onset anchors score. Site B (the depth scan inside Removal.H)
walks downward from the surface asking at what depth the crack would stop, which sets the
flake thickness; it is what decides Meier firings. Both compute the thermoelastic stress
from the same reconstructed exponential profile, integrate `K_I` with the same
`SpCriterion::K_I`, evaluate `K_Ic` at the crack tip with the same toughness law, and apply
the same relief factor, confining offset and melt skeleton. They are two copies, held
together by comments saying "BOTH firing sites must change together" at six places. P2-2b,
in flight right now, is editing both.

**Four side-face sites.** Once a hole exists, the vertical faces of the staircase wall are
exposed to the gas. The heat they receive is computed in four places: the conduction
kernel (which actually heats the rock), the jet enthalpy march (which charges that heat to
the jet's budget so the gas cools accordingly), the energy ledger (which re-computes it to
check conservation), and the lateral criterion (which asks whether a wall face spalls).
Each of the four resolves the heat-transfer coefficient per face (explicit pair, radius
gate, derived law, cap) and the first three each contain the full three-way area-law
chain (D2z step-height with the P1e staircase factor, or D2x exact riser, or the D2k
uniform factor). Every bug in the P1f to P2-1 line was one copy drifting from another.

---

# Part B. The biggest issues, ranked

Each issue has the same five headings so they can be compared. "Where" gives function
names and approximate lines. "Safe because" says what proves a fix did not change
behaviour; in almost every case that is the gate.

## B1. The order of field registration is a hidden contract

**What.** The plotfile layout and the thermo.dat column order are decided by the order in
which `Parse` calls `RegisterNewFab` and `RegisterIntegratedVariable`. The gate hashes
both. So reordering two registration lines, with no change to any number, breaks 2,639
hashes. The code knows this: fifteen comments say "registered LAST" or "APPENDED so
nothing above moves", and the P1g-3 counters carry a fifteen-line comment about the four
places they were tried before one worked.

**Where.** `Parse`, roughly lines 480 to 1180, where about 45 field registrations and
about 130 thermo registrations are nested inside feature-key `if` blocks and interleaved
with about twenty cross-key abort checks.

**Why it matters.** It makes `Parse` the one function in the file that cannot be
reorganised by ordinary means, and it has shaped the code in ways that have nothing to do
with physics. Counters for D2aa, D2ab and P2-1 are registered inside
`if (side_face_step_height)` not because they depend on that key but because every study
leg at the time happened to set it. The guards express history, not dependency.

**Proposed fix.** Two parts.

1. During the restructure, never reorder a registration. Split `Parse` by moving blocks
   whole into per-feature functions called in the same sequence, so the registration
   order is preserved by construction. The gate confirms it on every keys-off input.
2. Afterwards, make the contract explicit. Put every registration in one function,
   `RegisterFields()`, in a flat ordered list with a comment saying this order is the
   plotfile layout. Then anyone who adds a field appends to that list and sees why. If
   the team ever wants to clean up the order itself, that is a deliberate one-time
   re-pin of the gate with a header-only diff, done as its own packet.

**Safe because.** The gate hashes exactly the thing this contract protects. If a move
preserves order, 2,639 of 2,639 match. If it does not, the failure is immediate and the
diff is a header line.

## B2. Side-face heating is implemented four times

**What.** The heat delivered to an exposed vertical face is computed independently in
the kernel, the jet march, the ledger and the lateral criterion. Three of the four carry
the whole area-law chain. The ledger exists partly to catch the other copies drifting,
which it does by `ledger_err` blowing up.

**Where.**

- Kernel: `AdvanceMicrostructure`, the `col_in_patch` block, about lines 2070 to 2220.
- Jet march: `JetEnthalpyMarch`, the `SideHit` gather (about 3160 to 3260) and the
  `q_col` lambda (about 3400 to 3470).
- Ledger: `EnergyLedger`, about lines 4280 to 4560, plus the opposite-law diagnostic.
- Lateral criterion: Removal.H stage 5, about lines 560 to 640.

**Why it matters.** The P1f unit bug (pascals compared with megapascals), the P1g-0 length
bug (the Newton's length taken as the normal distance instead of the cell width), the
P1g-0b cap bug (the criterion capped the coefficient unconditionally while the kernel
capped it only under one key) and the P2-0 review's warning about "the fourth instance of
the family" were all the same defect: one copy read something the others did not. Every
packet that touches side faces has to edit four places and prove they agree.

**Proposed fix.** One function that answers "what does this face receive" and is the only
place that knows:

```
struct SideFaceResult { Set::Scalar h_face, area, q_in, q_out, q_robin, q_loss, T_face; };
SideFaceResult SideFaceFlux(int col, int k, const FaceOpenness& faces,
                            const SideFaceContext& ctx) const;
```

`SideFaceContext` is built once per step in `Advance` (it holds the explicit pair, the
radius gates, the cap, the wall-law state, the D2z phases, which area law is on). The
kernel calls `SideFaceFlux` and adds `q_in * area / V` to its divergence. The march calls
it and adds `q_robin * area` to the bin's draw on the jet. The ledger calls it and sums
`q_robin * area` and `q_loss * area`. The criterion calls it for `T_face` and `h_face`.
The ledger then stops being a mirror and becomes a sum, which is what a ledger should be.

The kernel's call is inside a device-capable `ParallelFor` lambda. The context struct
must therefore be trivially copyable (plain scalars and raw pointers to the per-column
arrays), which the current hoisted locals already are.

**Safe because.** With the four copies replaced by one call, the keys-off inputs take
exactly one branch of one function. The gate proves identity; `ledger_err` on the
side-face unit tests proves the three live roles still agree.

## B3. The column-top surface is computed four ways and about eight times per step

**What.** "Which cell is the top solid cell of each column" is the most basic question
the surface code asks. It is answered by `ColumnTopSolid` (sweep `removed_mf`, reduce by
max, sentinel `INT_MIN` for an empty column), by `GatherColumnSurface` (same sweep, same
rule, different argument list, falls back to the domain top when removal is off), by
Removal.H PASS 1 (same sweep inline), and by `RepairRemovalStateAfterRegrid` (same sweep
inline). Within one step `ColumnTopSolid` is then called from `AdvanceMicrostructure`,
`UpdateWallLawState`, `PerCellRemoval`, `CheckSurfaceCovered`, `TagCellsForRefinement`,
Removal.H's mechanical clearance, and `BuildSurfaceTColumns` goes through
`GatherColumnSurface`. Each call is a full sweep of the level plus an MPI reduce.

**Where.** `ColumnTopSolid` about 5539; `GatherColumnSurface` about 6871; Removal.H
226 to 258; `RepairRemovalStateAfterRegrid` about 4830.

**Why it matters.** Four copies of the rule means the rule can diverge (it did once: the
old form also required `phi >= 0`, and D2q-1 had to remove that from each copy). Eight
sweeps per step is pure cost. And every function that needs the surface has to be handed
it or recompute it, which is why `k_top_col` is threaded through so many argument lists.

**Proposed fix.** Compute it once at the top of `Advance`, after removal has settled from
the previous step, into a small per-step structure:

```
struct ColumnState {
    amrex::Vector<int> k_top;          // top solid cell, INT_MIN if none
    amrex::Vector<int> solid_run;      // contiguous solid cells below the top
    amrex::Vector<int> solid_above;    // only when the beam is on
};
```

Every consumer takes a `const ColumnState&`. Removal recomputes `k_top` after it has
voided cells (it is the one producer of a change), and the next step's `Advance` rebuilds
the structure. `GatherColumnSurface` and the two inline copies are deleted.

**Safe because.** The rule is identical in all four copies today, so the result is
bit-identical. The gate confirms it.

## B4. The firing criterion is written twice

**What.** Described in A.5. Site A evaluates `Sp` at the surface; site B scans downward.
The shared physics is the stress profile (three branches: a prescribed expression, the
local thermoelastic reconstruction, or the mechanics stress field), with the relief
factor, the confining offset, the melt skeleton, the reference temperature from
`T_init`, the flaw size with its volume factor, the `K_I` integral and `K_Ic` at the tip.
Each site has its own lambda implementing all of this.

**Where.** Site A: `UpdateSpAfterMechanics`, the `sigma_at_depth` and `T_recon_at` lambdas,
about lines 7120 to 7260. Site B: Removal.H stage 8, the lambdas `T_at_didx`,
`T0_at_didx`, `LL_at_didx`, `T_at_zv`, `kic_pa_at_cell`, `kic_pa_at_z`, `k_i_at_zcrack`,
about lines 1500 to 1850.

**Why it matters.** Any change to the criterion is two edits that must be identical, and
the only check is reading them side by side. P0 "closed the gap deliberately"; P2-2a
re-applied change A at both; P2-2b is doing it again now. The difference between the two
sites' surface handling (site A lerps inside the top half-cell, site B does too but with
a different fallback) is exactly the kind of thing that drifts.

**Proposed fix.** One object that holds a column's profile and answers the two questions:

```
struct ColumnCriterion {
    // built once per column from the gathered T, T0 and Lambda_L profiles
    Set::Scalar sigma(Set::Scalar z_vertical) const;   // the one stress formula
    Set::Scalar K_I(Set::Scalar z_crack, Set::Scalar a_f) const;
    Set::Scalar K_Ic_at(Set::Scalar z_vertical) const;
};
```

Site A builds it from the cell's own column (it already has the gathered `k_top`) and asks
`K_I(0, a) / K_Ic_at(a)`. Site B builds it from the gathered profiles it already has and
runs its integer scan and bisection on `K_I(z) >= K_Ic_at(z + a)`. The three stress
branches live inside `sigma` once.

**Safe because.** The arithmetic is moved, not changed; the comments at site B note where
an expression must stay "verbatim" for bit identity, and those expressions move intact.
The gate covers the Kant onset anchors (site A) and the Meier legs (site B).

## B5. Removal.H is a single 3,654-line function

**What.** `UpdateRemovalAfterCohesive` is one function with about 25 nested lambdas, 20
stage banners in comments, and roughly 1,500 lines of comment. It cannot be tested below
the whole-run level, a reader has to hold a dozen per-column buffers in their head at
once, and the stage boundaries are not enforced by anything.

**Where.** All of Removal.H. The stage map is in C.5.

**Why it matters.** The function does at least nine separable things: calibrate, clear the
rock under the burner feet, find the surface, run the criterion, pick a winner per column,
apply the void, remap the enthalpy, assert the invariants, write the logs. A bug in the
log writer is as hard to find as a bug in the criterion because they are the same
function. And about a third of the function (stages 5, 14's banking half, 17 and 18, about
1,250 lines) is the P1f/P1g lateral line, which no input file switches on.

**Proposed fix.** Cut it at the stage boundaries into member functions that share one
per-step state object:

```
struct RemovalStep {
    int lev; Set::Scalar time, dt;
    const ColumnState& cols;             // from B3
    amrex::Vector<int> spall_fires, vapor_fires, regime;
    amrex::Vector<Set::Scalar> h_s_cand, h_v_cand, h_chosen, ...;
    // the rev_*, rc_*, lat_* buffers as named members
};
void RemovalCalibrate(RemovalStep&);
void RemovalMechanicalClearance(RemovalStep&);
void RemovalCriterionScan(RemovalStep&);        // stage 8, using B4
void RemovalChooseWinner(RemovalStep&);         // stage 13
void RemovalApply(RemovalStep&);                // stage 15
void RemovalAssertInvariants(const RemovalStep&);
void RemovalWriteEvents(const RemovalStep&);
```

The lateral stages become their own file (`Lateral.H`) or are deleted under B7. The
lambdas that only exist because everything is in one scope (`col_index`, `T_at_didx` and
friends) become small private helpers or move into `ColumnCriterion`.

**Safe because.** Nothing in the arithmetic moves; only the braces do. The gate covers
every removal path the studies use, and `unit/spall_event`, `unit/removal_increment`,
`unit/seed_void` cover the damage-law and debug paths it does not.

## B6. The four conduction kernels share copy-pasted setup

**What.** Each `Advance*` kernel begins with the same ten-line block that copies the beam
parameters into locals, the same six-line block for the loss constants, and in the kernel
body the same radiation-plus-convection loss expression, the same harmonic-mean face
conductivity lambda, and the same three lines computing the cell centre. The loss
expression also appears in `BuildSurfaceTColumns` and the lateral criterion; the
`kc = kappa * exp(-a_D * D)` expression appears seven times (kernel, ledger, march twice,
`WallStepDiagnostic`, `PatchDiagnostics`, `RefreshKappaEff`, the lateral criterion).

**Where.** `AdvanceConstAlpha` 1534, `AdvanceMaterial` 1575, `AdvanceMaterialImplicit`
1704, `AdvanceMicrostructure` 1850 and the diagnostics listed.

**Why it matters.** Low risk, but it is about 250 lines of pure repetition, and the
"perf: skip exp when a_D == 0" optimisation lives in only one of the seven `kc` copies.

**Proposed fix.** Three small device-safe helpers:

```
struct BeamParams { /* P, w0, lam, x0, y0, z0, R, profile, radius, alpha_f */
                    static BeamParams At(const Numeric::MMWBeam&, Set::Scalar t); };
AMREX_GPU_HOST_DEVICE Set::Scalar SurfaceLoss(Set::Scalar T, const LossParams&);
AMREX_GPU_HOST_DEVICE Set::Scalar DamagedKappa(Set::Scalar kappa, Set::Scalar D, Set::Scalar a_D);
AMREX_GPU_HOST_DEVICE Set::Scalar HarmonicFace(Set::Scalar a, Set::Scalar b);
```

Separately, the team should decide whether `AdvanceConstAlpha` and
`AdvanceMaterialImplicit` are still wanted. Both exist for early validation inputs; the
implicit path is single-level, hard-codes Neumann boundaries and ignores voids, and
CLAUDE.md already lists it as something to revisit.

**Safe because.** The helpers are the same expressions in the same order. The Step 2-4
and Zhang regressions are unit tests with their own reference checks (they are not in
the 2,639-file sweep), so those paths are covered by `unit/stefan`, `unit/beam`,
`unit/equilibrium` and the Zhang validation rather than by the gate.

## B7. About fifteen study switches are still in the code after their studies closed

**What.** Every study in the D2k to P2-1 line added a default-off key. Most keys that
produced a negative or inert result were left in: the Parse block, the members, the
thermo registration under a nested guard, the kernel branch, the ledger mirror, and a
comment explaining why the key exists. Nothing is wrong with any of them; they are simply
code for experiments whose answer is in the archive.

**Where.** Listed in C.7.3. The largest bodies are the P1f/P1g lateral line (Removal.H
stages 5, 14, 17, 18, about 1,250 lines, plus `PerCellRemoval` and about 40 members), the
P1e staircase factor at four sites, the D2x exact-riser area law at three sites, the
D2c/D2e jet options (`rec_fixed`), the D2v above-plane pair, the D2aa pad shield, and the
P0b `side_face_closure = cell` control.

**Why it matters.** Each dead key makes every live site longer by one branch and every
new packet has to reason about one more interaction. The `if / else if / else` area-law
chain in B2 is three laws because two studies each added one.

**Proposed fix.** Delete whole branches, one key per commit, in this order: Parse block,
member declarations, thermo registrations, kernel branch, march branch, ledger branch,
criterion branch, the unit test that only exists to exercise the key. Record the deletion
in `ARCHIVE_DONE.md` with the study that closed it, so the result is not lost. The keys
and their evidence are in C.7.3; the first refactor packet should confirm each verdict
against the current `CTRL7`/`CTRL8` argv set before deleting, because the production
configuration lives in argv (see B12).

**Safe because.** A default-off key is off in every one of the 2,639 gate inputs, so
deleting its branch cannot change a hash. The unit tests that set the key
(`unit/lateral_drive`, `unit/lateral_flux_lock`, `unit/side_staircase`,
`unit/side_face_area`) are deleted with it or re-pointed.

## B8. Diagnostics are interleaved with the physics

**What.** `EnergyLedger` is 699 lines, `WallStepDiagnostic` 180, `PatchDiagnostics` 120,
`BeamEnergyBalance` 170, `AssertStepSplit` 50, plus three CSV writers (jet profile, side
power, removal events, about 250 lines together) and `UpdateSpClusters` at 213. That is
about 1,700 lines whose only output is a thermo.dat column, a CSV row or an abort. They
are called from inside the kernel's epilogue and from the middle of Removal.H, and the
ledger duplicates the kernel's side-face arithmetic (B2) in order to check it.

**Where.** Listed; see C.4 for lines.

**Why it matters.** The production path cannot be read without reading the diagnostics
around it, and the ledger's mirror is the single largest source of "change it in four
places" (B2). The 130 thermo scalars are declared wherever their packet put them, so the
member list is scattered (B10).

**Proposed fix.** Move everything listed into `MMWSpalling/Diagnostics.H`. Once B2 is
done, `EnergyLedger` shrinks to the enthalpy sums plus a loop calling `SideFaceFlux`
and summing; it no longer needs to know what an area law is. Make `Advance` call a
single `RunDiagnostics(lev, time, dt, cols)` at the end of the step so the kernel's
epilogue holds no diagnostics at all.

**Safe because.** Diagnostics write thermo.dat and CSVs, which the gate hashes. Moving a
call does not change its output as long as it runs after the same state; the gate
confirms it.

## B9. Comments carry packet history instead of invariants

**What.** A third of MMWSpalling.H and two fifths of Removal.H are comment lines, 317 of
them flagged with a warning glyph. Most explain why a guard sits where it does, what the
previous revision got wrong, which packet rejected what, and the per-rank reduction rule.
`ParseSpallSettings` is 508 lines of which about 300 are comment. The same pad-ring
shielding comment is pasted twice in `BuildFlameColumns`; the same "READ 1 of 2" comment
is pasted twice in `UpdateSpAfterMechanics`; the "reduce per box, not per column" rule
appears six times.

**Why it matters.** The history is already in `ARCHIVE_DONE.md`, the daily logs and the
review findings. In the source it buries the one line that states the invariant. A
reader who wants to know what `side_face_h_max` does has to read why P1g-0b found it
capped unconditionally, three times, before reaching the value.

**Proposed fix.** For each comment block over ten lines: keep the invariant ("the kernel,
march, ledger and criterion must apply the same cap; see `unit/lateral_flux_lock`"), and
replace the history with a one-line pointer to the archive entry. Where a rule is
restated more than twice, turn it into code: the reduction rule becomes two helpers,
`ReducePerBox(x)` and `AssertRankIdentical(x)`, and the per-column loops call the second.

**Safe because.** Comments do not change hashes. The helpers for the reduction rule are
wrappers around the existing `ParallelDescriptor` calls.

## B10. Members are declared wherever their packet left them

**What.** The data members nominally start at about line 9010, but from 8766 onward they
are interleaved with function definitions: the D2z clock is declared before `ColPhi`,
`side_cols_above` and the `wall_step_*` group right after `StepConductance`, the `wall_*`
state after `side_faces`, the `lat_*` counters inside a run of P1 comments at about 10300,
and the `SurfacePatch` struct (110 members and its own 400-line `Parse`) sits in the
middle. Nested types (`SurfaceRecon`, `MineralPhase`, seven enums) are declared where
they were first needed, including `SurfaceRecon` mid-class at 1468 because Removal.H
needed it declared before its use.

**Why it matters.** Finding a member means grepping. Related state (the pin vectors, the
reconstruction arrays, the wall-law state) is not visibly related.

**Proposed fix.** One ordered member section with sub-headings by kind: fields, per-step
column vectors, carried column vectors (with the "not checkpointed" note once, at the
top), thermo scalars grouped by the function that writes them, parsed settings grouped by
prefix, nested types first. Pure reordering of declarations.

**Safe because.** Member order inside a class has no effect on behaviour or on output.

## B11. Carried host state is not checkpointed

**What.** The 35 per-column vectors and a dozen scalars (`nozzle_z`, `feet_init`,
`jet_T_rec_next`, `side_T_eff`, `ledger_started` and its baselines, `jet_profile_last`,
the CSV header latches) persist across steps but are not written to checkpoints. A
restart re-lands the burner, forgets every pin, re-baselines the ledger, resets the
D2z clock and the lateral accumulator, and re-writes CSV headers.

**Why it matters.** Today the long runs are not checkpointed at all (the memory note says
host sleep loses them), so this costs nothing yet. The moment a 490-second Meier run is
checkpointed to survive a crash, every one of these resets silently changes the
trajectory after the restart.

**Proposed fix.** Not a refactor task, but the refactor should make it tractable: once the
carried state is in one place (B10), either register it (as fields where it is per
column, or through the Integrator's checkpoint hook for scalars) or add a single
`Util::Abort` on restart with removal enabled, so the limitation is explicit. Add one
restart test.

**Safe because.** A restart abort changes no non-restart behaviour.

## B12. The production configuration is not an input file

**What.** The Meier feet deck (`validation/meier/sp_meier_pilot/input_feet_d2c`) sets the
D2c jet options, the feet, `robin_form = pinned` and the Weibull flaws. Every closure
adopted since then is added by the study harness on the command line: `robin_form =
surface`, `spall.surface_T = reconstructed`, `spall.removal_energy = flake`,
`spall.flake_h_max`, the whole `side_face_*` line, `side_face_h_law = derived`, and now
the P2-2 keys. The `CTRL7`/`CTRL8` definitions in `ACTIVE_STEP.md` are the only place
the current model's configuration is written down.

**Why it matters.** A reader of `tests/` cannot tell what the model runs with. Several
keys that look unused by any input (`flake_h_max`, `side_face_factor_above`,
`side_face_h_law`) are in fact on every current leg. The grep in C.7 undercounts for
exactly this reason.

**Proposed fix.** Commit one input file per production configuration (for example
`input_feet_p2` with the full CTRL8 set) and have the harness reference it with argv
reserved for the one key each leg varies. Then B7's verdicts can be read off a file.

**Safe because.** Adding an input file changes nothing that runs.

## B13. Multi-phase machinery for a single-phase criterion

**What.** The microstructure layer supports Voronoi grains, per-grain minerals, grain-
boundary flags, per-facet Weibull flaws, grain-boundary contact conductance and
per-phase enthalpy tables (about 700 lines). The criterion reads `phases[0]` directly in
`BuildSurfaceTColumns`, `UpdateSpAfterMechanics`, Removal.H stages 4 and 5, `RockLaw`
and the confining offset, and aborts if there is no phase 0. Every production deck is
`microstructure.mode = expression` with one phase, so `is_grain_boundary` is zero
everywhere and the grain-boundary flaw branch never runs (a memory note records the GB
narrative as fiction for exactly this reason).

**Why it matters.** This is a decision, not a bug. The grand research direction names an
adjoint-supervised microstructure surrogate, which needs the multi-phase layer. If that
remains the plan, the single-phase reads in the criterion are the debt. If it does not,
the Voronoi, grain-boundary and per-facet code is 700 lines that can go.

**Proposed fix.** Decide. If multi-phase stays: make the criterion read the per-cell
property fields (`kappa_mf`, `E_phase_mf`, `beta_mf`, `T_ref_phase_mf`) instead of
`phases[0]`, which is what P2-2b's "convert every site to one style" is already doing for
rho and c_p. If it goes: delete the Voronoi mode, the grain id, the GB flags and the
per-facet flaw branch, and `InitializeMicrostructure` becomes a property-field fill.

**Safe because.** Reading the field instead of the scalar is bit-identical under
`law = constant` (the fields hold the scalars); the gate confirms. Deleting the Voronoi
mode affects only `unit/voronoi`, `unit/grain_topology` and `unit/amr_microstructure_regrid`.

## B14. Feature interlocks are enforced in three places

**What.** Cross-key rules ("`side_face_flux` requires `follow_mask`", "`side_face_h` and
`side_face_T` both or neither", "`lateral_drive` requires `side_face_staircase`") live in
`SurfacePatch::Parse` (requires-chains, both-or-neither, and a `CheckPadRing` that had to
be called on both exits because of an early return), in `MMWSpalling::Parse` (an abort
matrix at about 960 to 1180), and at runtime (`AdvanceMicrostructure`'s cubic-cell abort,
Removal.H's P1g-0b assertions, and the P1f parse-time self-test that calls
`SurfaceCellFlux` twice with sentinel arguments to check that the right branch writes
the face temperature).

**Why it matters.** When a rule is in three places, a new key has to be checked against
all three, and the record has several "guard placed where it never fires" incidents.

**Proposed fix.** After the parsing is split by feature (B1), one `ValidateSettings()`
function called once at the end of `Parse`, with all cross-key rules in it, each a
one-line condition and an abort. The runtime asserts that only check configuration
(cubic cells, the P1g-0b form equality) move there too. The parse-time self-test of
`SurfaceCellFlux` becomes a unit test.

**Safe because.** An abort that fires on the same input in a different function is the
same behaviour; an abort that never fired on any gate input still does not.

## B15. Several things are recomputed every step that never change

**What.** `T0_col_prof`, the per-column profile of the stress-free reference
temperature, is gathered and reduced every step although `T_init_mf` is written once in
`Initialize`. `ColumnSolidAboveCount` reduces `ncols x nz` integers per step when the
beam is on. `PatchAzimuthFraction` evaluates 2,048 trig samples per step under the wall
law. `AssertStepSplit` runs every step under `side_face_step_height`. Under AMR with
mechanics, `InitializeMicrostructure` (including the O(cells x grains) Voronoi search)
runs twice per step because both `TimeStepBegin` and `TimeStepComplete` call
`RefreshMicrostructureDerivedFields`.

**Why it matters.** The profile memory note put 45 percent of a Meier step in the
conduction kernel and 14 to 18 percent in the enthalpy inversion; these gathers are
smaller but they are also avoidable. The double microstructure refresh is not on the
production path (mechanics is disabled there) but it is why AMR mechanics runs are slow.

**Proposed fix.** Cache `T0_col_prof` in `ColumnState` (B3) and invalidate on regrid, which
`Regrid` already hooks. Compute `PatchAzimuthFraction` once at init (the domain does not
move). Make `AssertStepSplit` debug-only or every N steps. Remove one of the two refresh
calls after checking which one the AMR mechanics tests need.

**Safe because.** A cached value that equals the recomputed value is bit-identical; the
gate and `unit/amr_microstructure_regrid` confirm.

---

# Part C. The tour

## C.1 Status codes used in the tables

- **core**: on the Meier/Kant production path.
- **val**: exists for a validation case or unit test only.
- **diag**: writes a diagnostic; no effect on the trajectory.
- **study**: behind a default-off key whose study has closed; see C.7.3 for the verdict.
- **dead**: unreachable, or computed and then discarded.

## C.2 MMWSpalling.H, region by region

### C.2.1 Lifecycle and parsing

`Parse` (line 86, 1,100 lines) is the entry point for every input key. It delegates to
`HeatConduction`, `Mechanics`, `MMWBeam` and `SurfacePatch::Parse` through `queryclass`,
calls six `Parse*Settings` helpers, and parses inline the `phase.*`, `losses.*`,
`material.type`, `solver`, `alpha`, the whole `microstructure.*` block (nphases, mode,
seed, per-phase properties with the P2-2 `law = grimsel` option, phase change, alpha-beta
transition), `conductivity.*`, `damage.ic`, `surface.follow_mask`,
`energy_ledger.enabled`, and the AMR band keys. The second half is the registration
section described in B1. Note that `material.type` can select only `Zhang` or `Constant`;
`Grimsel` and `Table` are included at the top of the file but are not selectable here
(Grimsel is reached through the per-phase law instead).

| Function | Line | Size | Status | Notes |
|---|---|---|---|---|
| `Parse` | 86 | 1,100 | core | see B1 |
| `Initialize` | 1188 | 35 | core | builds the enthalpy table, calls the base inits, stores the initial temperature, initialises phi, seeds a void if asked, builds the mask and microstructure, converts T to H |
| `StoreInitialTemperature` | 1225 | 10 | core | the only write of `T_init_mf` |
| `TimeStepBegin`, `TimeStepComplete` | 1235, 1247 | 10 each | core | both refresh the microstructure under AMR with mechanics (B15) |
| `Regrid` | 1258 | 20 | val | repairs phi/removed, mask, microstructure, inverts H to T, repairs mechanics |
| `UpdateModel` | 1278 | 206 | core only with mechanics on | builds the nodal elastic model from cell properties and free strain; two near-identical loops (heterogeneous and homogeneous); carries a prescribed-T override for the Hu case |
| `Advance` | 1484 | 50 | core | A.3 |
| `TagCellsForRefinement` | 5851 | 30 | val | base tagging plus the surface band |
| `ParseDamageSettings` | 5882 | 100 | val | |
| `ParseCohesiveSettings` | 6142 | 50 | val | |
| `ParseSpallSettings` | 6269 | 508 | core | about 300 lines of comment; includes the P1f self-test executed at parse time (B14) |
| `ParseVaporSettings` | 6777 | 40 | val | |
| `ParseBitSettings` | 6818 | 30 | val | |
| `ParseSpallationModel` | 10835 | 200 | core | model, confining, `sp.*`, toughness law, prescribed overrides, `weibull.*` |
| `SurfacePatch::Parse` | about 9500 | about 400 | core | the side-face block is pasted one indentation level too far left inside the `side_face_flux` branch; it compiles because braces balance, but it reads as if it were outside the branch |
| `SurfacePatch::CheckPadRing` | about 9470 | 40 | core | three aborts that must run on both parse exits |

### C.2.2 Thermal core

`AdvanceMicrostructure` (line 1850, 583 lines) is the production kernel. Its first 190
lines hoist about 45 locals (beam, losses, every side-face setting, the flame arrays, the
pin arrays, the Robin form flags) and call `ColumnTopSolid`, the D2z `UpdateColumnPhase`
and `AssertStepSplit`, `ColumnSolidAboveCount` (beam only), `BuildFlameColumns` and
`BuildPinEffective`. The per-cell lambda then builds six face conductivities with
`k_face` (harmonic mean of the damage-reduced conductivity, zero across a void face,
grain-boundary contact resistance in series, the D2z `StepConductance` correction at a
staircase step on lateral faces), forms the divergence, adds the top-face patch flux
through `SurfaceCellFlux`, adds the side-face flux (B2), adds the beam source with
`SolidFacePath`, and updates H. Four locals are computed and discarded: the phase ids
`pc` and `p_n`, and the `q_rob`/`q_los` outputs of the top-face call.

| Function | Line | Size | Status | Notes |
|---|---|---|---|---|
| `AdvanceConstAlpha` | 1534 | 40 | val | Step 4 path |
| `AdvanceMaterial` | 1575 | 130 | val | explicit kappa(T), harmonic faces, void-aware, patch |
| `AdvanceMaterialImplicit` | 1704 | 146 | val | backward Euler through MLABecLaplacian; single level; Neumann hard-coded; ignores voids |
| `AdvanceMicrostructure` | 1850 | 583 | core | above |
| `UpdatePropertyFields` | 5381 | 236 | core (P2-2) | per-step Grimsel refresh of kappa, E, mu, rho, Cp and the secant beta |
| `ColumnTopSolid` | 5539 | 25 | core | B3; its header comment is stranded 500 lines above, just before the wall-law banner |
| `ColumnSolidAboveCount` | 5567 | 45 | core, beam only | |
| `ColumnSolidRunBelowTop` | 5611 | 45 | core | |
| `BeamEnergyBalance` | 5617 | 170 | diag | E1 closure; carries a retired invariant as a comment |
| `BeamSource`, `BeamIncidentIntensity` | 5787, 5824 | 40, 30 | core | clean and pure |
| `BuildEnthalpyTable`, `SetPhaseFractionsConstAlpha` | 7572, 7581 | 5, 10 | val | |
| `ConvertTtoH`, `InvertHtoT` | 7590, 7644 | 55, 65 | core | each has the same three-way dispatch (phase table, global table, closed form) |
| `RefreshKappaEff` | 5981 | 15 | diag | fills a plotted field the kernel does not read |

### C.2.3 Surface boundary and reconstruction

`SurfaceCellFlux` (line 3968, 130 lines) is the one routine every surface flux goes
through. It has four forms selected by its arguments: `surface` (the reconstructed
surface temperature is supplied, so no solve is needed), `face` (a ten-iteration
safeguarded Newton on the half-cell balance to find the face temperature), `pinned`
(the face form, then the surface is held at the column's last firing temperature if the
Newton would exceed it), and `cell` (flux on the cell average). It has 20 positional
parameters, four of them outputs, plus two optional output pointers. Eight call sites.

`BuildSurfaceTColumns` (line 3610, 150 lines) is the P0b reconstruction: gather the top
two cell temperatures per column, take `x = log((T_top - T0)/(T_1 - T0))`, hold the last
valid `x` on a flat column, optionally compute the P1g-3 flux-based length, and store
`T_surf`, `ell_eff`, `h_f`, `T_bar` and a flag per column. It reads `phases[0].T_ref`
and `phases[0].kappa` directly (B13). Its `dt` argument is unused since rev5.

| Function | Line | Size | Status | Notes |
|---|---|---|---|---|
| `BuildSurfaceTColumns` | 3610 | 150 | core | |
| `SharedSurfaceRecon` | 3762 | 20 | core | rebuilds the struct from the arrays |
| `SurfaceRobinT` | 3785 | 5 | core | returns the column surface T only under `robin_form = surface`; a five-line orphaned comment sits above it |
| `BuildPinEffective`, `PinRule` | 3792, 3814 | 20, 20 | core (pinned) | the idle rule |
| `PatchDiagnostics` | 3838 | 120 | diag | |
| `RiserWeight` | 3962 | 5 | study (D2x) | |
| `SurfaceCellFlux` | 3968 | 130 | core | above; the `AMREX_GPU_HOST_DEVICE` at 3952 attaches to nothing |
| `RawSurfaceX`, `ReconFromX`, `ReconstructSurfaceT` | 10688, 10708, 10799 | 15, 90, 10 | core | the reconstruction algebra, pure static |
| `SurfaceGasTempCap` | 10815 | 10 | core | |
| `InitSurfaceMask`, `UpdateSurfaceMaskFromPhi` | 7708, 7739 | 15 each | core | the phi-window mask; D2q-1 "site 7", the one site that could not be made removed-mask authoritative without changing 24 plotfiles |
| `SurfaceOwner` | 4815 | 1 | core | |

### C.2.4 Flame, jet and feet

`BuildFlameColumns` (line 2433, 170 lines) produces, once per step, the per-column flame
coefficient `flame_h`, gas temperature `flame_Tg` and stand-off `flame_s`. On the
constant path it assigns `h_conv` and `T_flame`. On the expression path it gets the
nozzle height from `FeetDescent` (or the prescribed feed), computes the stand-off,
aborts on a collision, and then takes one of two branches: the jet-enthalpy branch
(coefficient from `h_expr`, gas temperature left for the march) or the plain expression
branch (both from expressions). The twenty-line pad-ring shielding block appears verbatim
in both branches. It then calls `JetEnthalpyMarch` and `WallStepDiagnostic`.

`JetEnthalpyMarch` (line 3041, 569 lines) is the D2a2 wall-jet closure. It is serial and
replicated on every rank. It gathers the top-cell temperature and conductivity of every
column, gathers the exposed side faces sparsely (the `SideHit` structure, B2 copy 2),
bins the columns by radius, computes the stagnation state (mixing temperature from fixed,
exhaust or a fixed recirculation value; decay diameter from the nozzle or a momentum
equivalent; the core-length decay; entrained mass), then marches outward bin by bin: each
bin draws its Robin flux from the gas, the gas cools by that amount divided by `m cp`,
optionally dilutes at the free surface, and writes a profile row. It ends by writing
thirteen `jet_*` scalars and aborting if the enthalpy budget does not close to 1e-9.

`FeetDescent` (line 2873, 170 lines) is the D2b burner-on-feet model: collect the
columns under the foot annulus, split them into pads by azimuth, take a nearest-rank
quantile of each pad's floor height, average, add the stand-off, and let the nozzle
descend to that height at a bounded rate without ever rising. Under D2c it also arms
`foot_clear_k` so Removal.H can cut rock standing proud of a pad.

| Function | Line | Size | Status |
|---|---|---|---|
| `BuildFlameColumns` | 2433 | 170 | core |
| `FeetDescent` | 2873 | 170 | core (feet) |
| `JetEnthalpyMarch` | 3041 | 569 | core (enthalpy closure) |
| `UpdateColumnPhase` | 2602 | 40 | core (D2z) |
| `AssertStepSplit` | 2641 | 50 | diag (D2z) |
| `WallStepDiagnostic` | 2693 | 180 | diag |

### C.2.5 Side faces and the derived wall law

The P2-1 wall law is about 370 lines: nine `static constexpr WJ_*` constants,
`WallMu` (Sutherland), `WallZFill`, `WallJetBranch`, `WallJetH`, the per-face resolver
`SideFaceH` (returns the typed coefficient unless the law is on and ready),
`ColumnWallCorner` (the height datum by a ray march to the axis), `PatchAzimuthFraction`
and `UpdateWallLawState` (bore radius from a histogram of column tops, gas state,
momentum product, Gnielinski duct branch). It is live and review-accepted.

The D2z helpers at 8766 to 8888 (`ColPhi`, `PairDelta`, `StepCellWeight`,
`StaircaseFactor`, `AreaOfDelta`, `StepConductance`) are about 90 lines of pure static
code under 150 lines of comment. `StaircaseFactor` and the `wetted` branch of
`AreaOfDelta` are closed study legs (P1e, D2z's second law).

### C.2.6 Column geometry and AMR repair

| Function | Line | Size | Status |
|---|---|---|---|
| `RepairRemovalStateAfterRegrid` | 4824 | 90 | val (AMR) |
| `CheckSurfaceCovered` | 4917 | 50 | val (AMR) |
| `GatherColumnSurface` | 6871 | 50 | core; the second copy of B3 |
| `ColumnSurfaceCos` | 6922 | 35 | study (Step 20) |

### C.2.7 Microstructure

`InitializeMicrostructure` (line 8066, 401 lines) builds the per-phase materials and
enthalpy tables (Grimsel-aware, copying the latent-heat scalars across first), generates
the Voronoi grain centres lazily with a seeded generator or evaluates the phase and grain
expressions, assigns each cell its phase and grain by nearest centre (an O(cells x grains)
search), fills the per-cell property fields from the per-phase scalars (the one
legitimate scalar-to-field writer), fills ghost cells, computes the phase- and
grain-boundary flags, applies the damage initial condition, refreshes the effective
conductivity, resets the cohesive fields, computes a bulk conductivity diagnostic, and
samples the flaws. `InitializeFlaws` (line 8467, 299 lines) draws a Weibull flaw length
per cell from a hash of the cell index and the seed, per facet on grain-boundary cells and
per cell elsewhere, so the field is identical for any rank count.

| Function | Line | Size | Status |
|---|---|---|---|
| `InitializeMicrostructure` | 8066 | 401 | core |
| `InitializeFlaws` | 8467 | 299 | core (weibull) |
| `RefreshMicrostructureDerivedFields` | 7937 | 10 | val (AMR) |
| `UpdateMicrostructureFreeStrain`, `ThermalFreeStrain`, `AlphaBetaFraction`, `CellToNodeAverageInDomain` | 7980 to 8066 | 70 | core with mechanics on |
| `WeibullSizeFactor` | 9009 | 5 | core (P2-2) |
| `RockLaw`, `LawSigma` | 7610, 7597 | 5 each | core (P2-2b) |
| `SplitMix64` | 10827 | 5 | core |

### C.2.8 The damage-law lineage

Steps 8 to 11 built a Drucker-Prager damage evolution, a grain-boundary cohesive law and a
spall rule based on a damage threshold and a principal-stress ratio. `sp_weibull`
replaced all of it for every production run. It survives for `unit/dp_yield`,
`unit/gb_cohesive`, `unit/spall_event` and the Hu validation. The cohesive feedback into
the elastic operator was deferred in Step 9 and never implemented.

| Function | Line | Size |
|---|---|---|
| `DamageGrowthRHS`, `IntegrateDamageRK4` | 5999, 6013 | 15 each |
| `UpdateDamageAfterMechanics` | 6031 | 110 |
| `UpdateCohesiveAfterMechanics`, `ResetCohesiveDiagnostics` | 6195, 6260 | 65, 10 |
| `PrincipalStressRatio` | 6848 | 20 |
| `RepairMechanicsAfterRegrid`, `FiniteMatrix` | 7946, 7970 | 25, 5 |

### C.2.9 Criterion and removal

`UpdateSpAfterMechanics` (line 6956, 403 lines) is criterion site A (B4). It chooses
the stress source, computes the confining offset and the Step 20 slope cosines, reads the
P0b reconstruction with seven witness counters, and then for every column's top cell
(or every exposed cell under `per_cell_removal`) evaluates the stress at depth, the flaw
with its size factor, `K_I`, `K_Ic` at the tip, and writes `Sp`.

| Function | Line | Size | Status |
|---|---|---|---|
| `UpdateSpAfterMechanics` | 6956 | 403 | core |
| `UpdateSpClusters` | 7359 | 213 | val (two inputs) |
| `UpdateRemovalAfterCohesive` | Removal.H | 3,654 | core; C.5 |
| `PerCellRemoval` | 7784 | 85 | study (D2q-2b) |
| `ApplySeedVoid` | 7869 | 30 | debug |
| `InitPhi` | 7902 | 35 | core |
| `MeltSkeleton` | 10211 | 5 | study (P1 melt) |

### C.2.10 The energy ledger

`EnergyLedger` (line 4104, 699 lines, level 0 only) sums `H * V` over all cells and over
solid cells with compensated summation, re-evaluates every top-face and side-face flux
through `SurfaceCellFlux` to get the Robin power and the losses, and closes
`ledger_err = (dE_H - E_src) / E_src`. On the way it computes about forty thermo
scalars: the D2aa bore/pads/wall split, the P2-1 coefficient histogram quantiles, the
D2x realised riser weight, the D2z step-height statistics, the P1e staircase statistics,
the P0b side-power shares, and the per-column side-power CSV. See B2 and B8.

## C.3 Removal.H, stage by stage

One function. Host loops only. Line ranges from the 3,654-line snapshot.

**Stage 0, lines 20 to 115. Prologue.** The damage-law flake thickness
`h_spall = C_h sqrt(alpha t)`, the vapour constants, the per-step reset of the event and
regime fields, and the `col_index` lambda.

**Stage 1, lines 158 to 225. Mechanical clearance under the feet (D2c).** For each
annulus column standing proud of its pad, flip every cell above the pad index to void with
regime 4 and shift the column's phi. Live on the feet deck.

**Stage 2, lines 226 to 258. PASS 1.** Find each column's top solid cell (B3, copy three)
and, under Step 20, the slope cosines.

**Stage 3, lines 260 to 340. PASS 1.5.** Per-column maximum damage over the `spall.sample`
window. Damage-law path only.

**Stage 4, lines 342 to 412. Captures.** Allocate the per-column candidate buffers; a
twenty-line comment about the P1b support rule whose code was deleted; capture phase 0's
E, beta, nu, T_ref, kappa and the prescribed-temperature callable (B13).

**Stage 5, lines 413 to 1024. The lateral criterion (P1f).** Six hundred lines. Two runtime
assertions that the criterion's hard-coded `SurfaceCellFlux` arguments match the kernel's.
Then for every solid cell on the hole wall with an open lateral face: riser counts on each
axis, the coefficient through `SideFaceH` with the cap (B2 copy four), `SurfaceCellFlux`
for the face temperature, the two-cell slope `ell` along the face normal, the P1g-3 flux
length, the flake fraction, the flaw with its size factor, a `sigma_lat` lambda, `K_I`,
`K_Ic` at the tip (in megapascals, converted), the score, and if it fires, thirteen logged
quantities and the cos-squared split weights. Fifteen `ReduceRealMax` calls to pick the
per-column winner. No input file enables this; see B7.

**Stage 6, lines 1024 to 1100. Gates.** The Step 19 stress-source selector and its abort,
the confining offset, the `HColEvent` struct, and the booleans `pin_on`, `rev_on`,
`inc_on`, `inc_active`, `recon_on`, `lat_csv` that gate everything below.

**Stage 7, lines 1100 to 1250. Profile gathers (PASS 1b).** For every column, a window
of `Dmax + 1` cells below the top: temperature, the stress-free reference (re-gathered
every step though it never changes, B15), enthalpy and phi (for the remap), liquid
fraction (for melt). Each is an `ncols x colw` reduction.

**Stage 8, lines 1359 to 2120. PASS 2, the per-column criterion (site B).** The owner rank
of each column's top cell records the event-log values and the per-cell rho and c_p,
then takes one of three branches. Under `sp_weibull`: the detach gate (cluster id or
`Sp_top >= 1`), the lambdas `T_at_didx`, `T0_at_didx`, `LL_at_didx`, `T_at_zv` (the
surface lerp), `kic_pa_at_cell`, `kic_pa_at_z`, `k_i_at_zcrack` (with `sigma_at_depth`
inside it: melt skeleton, then either the mechanics stress or the thermoelastic law or
the verbatim pre-P2-2b expression), an integer scan downward with the crack-event edge
trigger, an eight-step bisection, the flake thickness `h_col`, the running-max field, and
the first-firing CSV event. Under damage law: the threshold test and the principal
stress ratio. Then the vapour branch with its three rate forms. A second P1b comment with
no code.

**Stage 9, lines 2120 to 2185. Reductions** of every candidate buffer and the
step-to-step change of `ell_eff`.

**Stage 10, lines 2185 to 2225. Flake coherence cap.** A firing column may recede at most
`coh_cells` below its shallowest neighbour. Bypassed under Step 20. The Meier drilling
input sets 4 mm; the feet deck sets 0.

**Stage 11, lines 2225 to 2277. The `h_col_events` CSV**, written rank by rank with
barriers.

**Stage 12, lines 2277 to 2315. The bit target** per column (Step 19).

**Stage 13, lines 2315 to 2474. The winner loop.** Per column: regime 1 (spall) or 2
(vapour) with spall priority under `sp_weibull`, 3 (bit controls), 6 (lateral, with zero
phi movement); then the D2z-1 whole-cell override; then the P0b flake thickness from the
reconstruction when `removal_energy = flake`.

**Stage 14, lines 2474 to 2600. Remap build and banking.** For each firing column, the
new enthalpy profile shifted up by the flake thickness using the exponential blend, the
energy the flake carries, and the predicted void count; the P1g-1 banking decision; the C1
pin update.

**Stage 15, lines 2600 to 2870. PASS 3, apply.** Every rank shifts phi in its own cells,
flips cells whose phi crossed zero to void (writing regime, RoP, event and thickness
diagnostics, and deliberately not resetting the void cell's temperature), voids or banks
the regime-6 cell, and writes the remapped enthalpy with the window-boundary energy term.

**Stage 16, lines 2870 to 2990. Invariants.** On the IO rank: every spall firing voided
exactly the predicted number of cells, with a message that distinguishes "ran out of rock
at the domain floor" from "the phi tie at a cell boundary"; and the Step 20
`sn_recession_superseded` witness.

**Stage 17, lines 2990 to 3400. The lateral thermal shift (P1g-2).** For each banking
column, gather two rows of enthalpy along the face normal (32 cells, refuse under 6),
build the shifted rows with the same exponential rule and the cos-squared split, compute
an independent energy estimate, record every touched cell in a `std::map` so two rows
writing one cell compose their deltas (P1g-2c), write back, count orphans. Twelve
diagnostics. No input enables it; see B7.

**Stage 18, lines 3400 to 3470. Lateral and banking diagnostics**, with four comment
blocks about which counters must and must not be reduced.

**Stage 19, lines 3470 to 3600. The `removal_events` CSV.** A header with three
conditional column groups and one row per fired column, plus mechanical-clearance rows.

**Stage 20, lines 3600 to 3654.** `CheckSurfaceCovered` under AMR; `InvertHtoT` and a
ghost refill because the remap wrote enthalpy only.

## C.4 Diagnostics, collected

For B8, the diagnostic code and where it is called from:

| Diagnostic | Lines | Called from | Output |
|---|---|---|---|
| `EnergyLedger` | 699 | kernel epilogue | about 40 thermo scalars, `side_power` CSV |
| `BeamEnergyBalance` | 170 | kernel epilogue | 8 thermo scalars |
| `PatchDiagnostics` | 120 | kernel epilogue | 7 thermo scalars |
| `WallStepDiagnostic` | 180 | `BuildFlameColumns` | 10 thermo scalars |
| `AssertStepSplit` | 50 | kernel prologue | abort only |
| `UpdateSpClusters` | 213 | `Advance` | a field, `clusters` CSV |
| jet profile CSV | 40 | `JetEnthalpyMarch` | file |
| `h_col_events` CSV | 50 | Removal.H stage 11 | file |
| `removal_events` CSV | 130 | Removal.H stage 19 | file |
| P0 witness counters | 30 | `UpdateSpAfterMechanics` | 8 thermo scalars |
| lateral diagnostics | 70 | Removal.H stage 18 | 30 thermo scalars |

## C.5 Data members, by kind

| Kind | Count | Where declared | Notes |
|---|---|---|---|
| `Set::Field` | 45 | 9010 onward, scattered | 34 plotted and therefore hashed |
| per-column host vectors | 35 | scattered | none checkpointed (B11) |
| thermo scalars | 130 | beside their function or packet | hand-reduced (A.4) |
| parsed settings (scalars, bools, enums, parsers) | 230 | scattered; 110 inside `SurfacePatch` | |
| nested types | `SurfacePatch`, `MineralPhase`, `SurfaceRecon`, 7 enums | wherever first needed | |
| `static constexpr` | 14 | with the wall law and the flake fraction | |

## C.6 The gather-and-reduce pattern

Because a column spans ranks in z, every per-column quantity is produced by the same
idiom: allocate a vector of `ncols` zeros, let the rank that owns the relevant cell write
its value, then `ReduceRealSum` (or `ReduceIntMax` for indices) so every rank holds the
same vector. It appears about thirty times across the two files. It is correct and
decomposition-independent, which is why it was adopted after the 0-ghost-cell seam bug.
Its cost is one reduction per vector per step and the discipline, stated in comments
six times, that a value every rank already computed identically must not be reduced
again. A `ColumnGather<T>` helper taking the owner-write lambda would make the pattern
one line and could carry the "already rank-identical" flag as a type.

## C.7 Input key census

The evidence columns are: **deck** means the key is set in the Meier feet deck
(`validation/meier/sp_meier_pilot/input_feet_d2c`) or the Kant onset deck; **inputs** is
how many of the 143 non-output `input*` files under `tests/` set the key to a non-default
value; **refs** is how many harness scripts, inputs and notes under `tests/MMWSpalling`
name the key at all. Study legs set most keys on the command line, so "inputs" badly
undercounts and "refs" is the honest reach (B12). The verdicts come from `ROADMAP.md`,
`ARCHIVE_DONE.md` and the memory notes, not from a re-run.

### C.7.1 Live on the production path

| Key | deck | inputs | refs | Note |
|---|---|---|---|---|
| `microstructure.*` (expression mode, one phase) | yes | many | | single phase on every production deck |
| `spall.enabled`, `spallation.model = sp_weibull`, `spall.stress_source = local_thermoelastic` | yes | 53 | | |
| `weibull.enabled`, `a0_gb`, `a0_ig`, `m`, `seed`, `V0` | yes | 31 | | grain-boundary branch inert on a single grain |
| `confining.p` | yes | 28 | | |
| `spall.sample = top_cell`, `spall.h_spall_max` | yes | 46 | | the other two sample modes are damage-law only |
| `surface.follow_mask`, `energy_ledger.enabled` | yes | 24, 21 | | |
| `surface_patch` convective flame, `h_expr`, `jet_closure = enthalpy`, `jet_entrained_mass`, `jet_T_ent_mode = exhaust`, `jet_decay_diameter = momentum`, `jet_bin_update = exponential`, `jet_negative_flux = count`, `jet_free_surface` | yes | 1 to 6 | 19 to 20 | the "D2c options" are all on in the feet deck |
| `surface_patch.nozzle_descent = feet`, `foot_*`, `foot_body_clearance` | yes | 4, 1 | 9 | |
| `surface_patch.robin_form` | `pinned` in the deck; `surface` on every study since P0b | 17 | | `cell` and `face` are controls |
| `spall.surface_T = reconstructed`, `spall.removal_energy = flake` | argv since P0b | 6, 6 | | the current rate model |
| `spall.flake_h_max` | argv since P0b rev7 | 0 | 21 | its own comment calls it a stopgap |
| `surface_patch.side_face_flux`, `side_face_h` and `side_face_T` (live), `side_face_h_r_min`, `side_face_factor_above`, `side_face_step_height` | argv (CTRL legs) | 0 to 2 | 22 to 28 | the D2k to D2z side-face line; step height review-accepted at D2z-2 |
| `surface_patch.side_face_h_law = derived` | argv (P2-1) | 0 | 4 | review-accepted |
| `microstructure.phase0.law = grimsel`, `spallation.sp.stress_relief`, `weibull.size_volume = flux` | argv (P2-2) | 0 | 0 to 1 | being implemented now |
| `spall.removal_events_csv`, `jet_profile_interval`, `side_power_interval`, `step_diag_radius` | yes / yes / argv / argv | | 17 to 19 | diagnostics the studies score |

### C.7.2 Validation and unit-test only

`material.type` (Zhang or Constant, 11 inputs), `solver = implicit`, the constant-alpha
path keys, `damage.*` (11), `cohesive.*` (2), `spall.prescribed_ratio`, `spall.alpha`,
`spall.t_spall`, `spall.xi_s`, `spall.C_h`, `vapor.*` (15), `bit.*` (2 inputs, 7 unit
references), `spallation.sp.prescribed_stress`, `prescribed_T_surface`,
`weibull.detach_mode = connected_cluster` with `A_crit` (2), `spall.clusters_csv`,
`spall.h_col_events_csv`, `surface.refine_depth` and `amr.max_level > 0` (3 references,
no input), `spall.seed_void` (debug, 3 unit references), `surface_patch.mode =
prescribed_T` (Hu), `alpha_beta.*`.

### C.7.3 Study legs the record closed as negative or inert

These are the candidates for B7.

| Key | inputs | refs | What the record says |
|---|---|---|---|
| `spall.lateral_drive`, `spall.lateral_removal` | 0 | 14 | P1f: Leg A inert by scope; Leg B accepted after the unit bug; P1g-0 and P1g-0b fixes accepted; never entered a production deck |
| `spall.lateral_subcell` | 0 | 13 | P1g-1 accepted; P1g-2b found "banking is 0.00 mm, the drift is the binary" |
| `spall.lateral_thermal` | 0 | 9 | P1g-2 accepted; P1g-2c closed the leak, worth 0.16 mm of bore |
| `spall.flake_length = flux`, `spall.flake_fraction = weibull` | 0 | 7, 8 | P1g-3 rejected on provenance |
| `spall.per_cell_removal` | 0 | 5 | D2q-2b; the P1b, P1c, P1e support-rule line all negative |
| `surface_patch.side_face_staircase` | 1 | | P1e accepted negative: cures the stall, fails the physics bar |
| `surface_patch.side_face_area_exact` | 0 | 2 | D2x: no side-face treatment closes the wall band |
| `surface_patch.side_face_area_law = wetted` | 0 | 21 | D2z: both laws held, `step` is the one used |
| `surface_patch.pad_shield_open` | 0 | 2 | D2aa bounded; P1 confirmed, P2 refuted |
| `surface_patch.jet_above_h`, `jet_T_above` | 0 | 3 | D2v: the funnel was the map, not this |
| `surface_patch.side_face_closure = cell` | 0 | 1 | P0b rev7 control |
| `spall.surface_normal` | 4 | | Step 20; the recession half is superseded by whole-cell removal and its own parse message says so; the criterion half remains |
| `spall.melt_aware` | 9 | | P1 melt skeleton; pins the surface near 850 K; not on the feet deck |
| `spall.flake_coherence_length` | 1 | 49 | superseded for the feet line (set to 0 there) |
| `surface_patch.jet_T_ent_mode = rec_fixed` | 0 | | D2e diagnostic |
| `surface_patch.robin_form = cell`, `face` | controls | | |
| damage-law spall (`spall.sample = column_max`, `shallow_max`) | Hu only | | |

### C.7.4 Dead code and stranded text

- `AdvanceMicrostructure`: `pc`, `p_n`, `q_rob`, `q_los` computed then discarded.
- `BuildSurfaceTColumns`: the `dt` parameter.
- Removal.H PASS 3: the temperature array bound and discarded.
- `SurfacePatch::Parse`: the legacy `type` key, consumed and checked for equality only.
- `BeamEnergyBalance`: invariant (b) exists only as a comment.
- Removal.H stages 4 and 8: two P1b support-rule comment blocks, about 45 lines, with no
  code.
- `SurfaceRobinT`: an orphaned "THE REFERENCE TEMPERATURE" comment above it.
- `ColumnTopSolid`: its comment 500 lines above the body, under the wall-law banner.
- `SurfaceCellFlux`: the `AMREX_GPU_HOST_DEVICE` at 3952 attaches to nothing.
- `Parse`: "P1b: the acceptance test for the support rule is that", unfinished.
- `BuildFlameColumns`: the pad-ring block pasted twice.
- `UpdateSpAfterMechanics`: the "READ 1 of 2" comment pasted twice.

---

# Part D. Where to go and in what order

## D.1 Target layout

Each partial header stays a class-body fragment, exactly as Removal.H is today, so the
inheritance graph and the virtual plumbing do not change (CLAUDE.md lesson 1).

| File | Content | Comes from |
|---|---|---|
| `MMWSpalling.H` | the class, members by kind (B10), `Parse` as a sequence of calls, `RegisterFields` (B1), lifecycle hooks, `Advance` | lines 1 to 1534 |
| `MMWSpalling/Settings.H` | `SurfacePatch`, `MineralPhase`, the enums, the `Parse*Settings` functions, `ValidateSettings` (B14) | 5882 to 6848, 9400 to 9900, 10835 to end |
| `MMWSpalling/Thermal.H` | the four kernels with the shared helpers (B6), `ConvertTtoH`, `InvertHtoT`, `BeamSource`, `UpdatePropertyFields` | 1534 to 2433, 5381 to 5850, 7572 to 7710 |
| `MMWSpalling/Surface.H` | `SurfaceCellFlux`, the reconstruction, the pin rule, `SideFaceH` and the wall law, `SideFaceFlux` (B2), `ColumnState` and its builders (B3) | 3610 to 4104, 5018 to 5380, 5539 to 5620, 10688 to 10835 |
| `MMWSpalling/Jet.H` | `BuildFlameColumns`, `FeetDescent`, `JetEnthalpyMarch` | 2433 to 3610 |
| `MMWSpalling/Criterion.H` | `ColumnCriterion` (B4), `UpdateSpAfterMechanics`, the depth scan from Removal.H stage 8 | 6956 to 7360, Removal.H 1359 to 2120 |
| `MMWSpalling/Removal.H` | the stage functions of B5 | Removal.H |
| `MMWSpalling/Lateral.H` | stages 5, 14 (banking), 17, 18 and `PerCellRemoval`, or deleted under B7 | Removal.H, 7784 |
| `MMWSpalling/Microstructure.H` | `InitializeMicrostructure`, `InitializeFlaws`, free strain, `UpdateModel` | 1278 to 1484, 7937 to 8766 |
| `MMWSpalling/DamageLaw.H` | damage, cohesive, `PrincipalStressRatio`, stage 3 | 5981 to 6270 |
| `MMWSpalling/Diagnostics.H` | everything in C.4 | |
| `MMWSpalling/AMR.H` | regrid repair, `CheckSurfaceCovered`, tagging | 4824 to 5020, 5851 |

## D.2 Order of work

Each item is one or more packets. Each packet's acceptance criterion is "2,639 of 2,639
hashes match and every unit test passes", with no physics number. One semantic move per
commit; the gate between commits.

1. **Make the gate one command** and pin the current binary's hashes as the baseline.
   Nothing else starts until this runs green on the unchanged tree.
2. **Confirm the C.7.3 verdicts** against the current `CTRL7`/`CTRL8` argv set, and commit
   the production configuration as an input file (B12).
3. **Delete closed study keys** (B7), one key per commit, Parse block to unit test.
4. **Remove dead code and stranded comments** (C.7.4) and trim the history comments (B9).
5. **Build the single-copy pieces**: `ColumnState` (B3), the thermal helpers (B6),
   `SideFaceFlux` (B2), `ColumnCriterion` (B4). Each replaces its copies in one commit.
6. **Cut Removal.H at its stages** (B5) and move the diagnostics out (B8).
7. **Reorder the members** (B10) and move the interlocks into `ValidateSettings` (B14).
8. **Split the files** along D.1. `Parse` last, by moving registration blocks whole (B1).
9. Then, as separate decisions outside the refactor: the implicit path's future (B6),
   checkpointing the carried state (B11), single-phase versus multi-phase (B13), and the
   caches in B15.

Physics changes stay out of every one of these packets. The gate cannot tell a cleanup
from a behaviour change, so a mixed commit that fails the gate cannot be diagnosed.

## D.3 Caveats

- Line numbers are from the 11,033 and 3,654-line snapshot and have already drifted by
  about 80 lines. Re-anchor by function name.
- Key reach counts are greps, not an execution trace. A harness that names a key in a
  comment counts as a reference; a key passed by argv does not appear in any input file.
- The status verdicts reflect the review record as written in `ROADMAP.md`, the archive
  and the memory notes, not a re-run. Step 2 of D.2 is where they are confirmed.
