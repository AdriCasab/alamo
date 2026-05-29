# Active Step: R1 - Extract UpdateRemovalAfterCohesive into a partial header

Status: completed

## Context

This is a **behavior-preserving refactor**, not a physics feature. It begins
decomposing `src/Integrator/MMWSpalling.H` (~4377 lines), which has grown into a
monolith owning thermal, microstructure, two damage models (continuous-DP
`damage_law` and LEFM `sp_weibull`), cohesive, and removal physics on one class.

The decomposition rule (set with the user): **keep one integrator class and do
not touch its virtual multiple-inheritance graph** (`MMWSpalling : virtual public
HeatConduction, virtual public Mechanics<...>`). CLAUDE.md "Lessons Learned" #1:
earlier attempts to restructure the base classes caused vtable/thunk crashes
during AMR setup. So the split must live *below* the integrator, never in its
inheritance.

R1 extracts the single biggest, most self-contained piece: the
`UpdateRemovalAfterCohesive` method (~957 lines, lines ~2100–3053, ≈22% of the
file). It is shared by BOTH damage models — `damage_law` (D-threshold spall) and
`sp_weibull` (LEFM cluster spall) — plus vaporisation, so it is correctly a
shared removal stage, not a model-specific one.

**Mechanism (lowest-risk possible): physical code-motion only.** Move the method
*body verbatim* into a new partial header that is `#include`d back **inside the
class body** at the method's original location. The method stays an inline class
member; nothing about its semantics, call sites, or member access changes. The
git diff must show the method deleted from `MMWSpalling.H` and added to the new
header with **zero content change**.

Deliberately deferred to a later refactor step (R2+): relocating the removal
*data members* (`phi_mf`, `removed_mf`, `spall_event/thickness`, `RoP_spall`,
`vapor_*`, `regime_mf`, `RoP_mf`, `regime_field_mf`, `h_spall_field_mf`,
`crack_event_count_mf`) and their config scalars, and moving removal-only helpers
(`PrincipalStressRatio`). Keeping R1 to a single method body is what makes it
provably safe and trivially reviewable. Do NOT move data members in R1.

## Sources

All facts below were verified against the working tree during planning.

- **Method to move**: `void UpdateRemovalAfterCohesive(int lev, Set::Scalar
  time, Set::Scalar dt)` at
  [src/Integrator/MMWSpalling.H:2100](src/Integrator/MMWSpalling.H#L2100). Its
  leading doc comment starts at line ~2092 (`// pass populates each rank's local
  maximum k_top ...`). The method body ends at the `}` on line ~3053 (the lines
  just above it are `UpdateSurfaceMaskFromPhi(lev); RefreshKappaEff(lev);
  (void)dz; }`). The very next member is `void BuildEnthalpyTable()` at
  [line ~3057](src/Integrator/MMWSpalling.H#L3057). The member just above the
  doc comment is the end of `UpdateSpClusters` (its closing `}` at line ~2090).
- **Sole call site**: `Advance()` at
  [src/Integrator/MMWSpalling.H:737](src/Integrator/MMWSpalling.H#L737)
  (`UpdateRemovalAfterCohesive(lev, time, dt);`). This does NOT move and does NOT
  change.
- **Includes already present** (lines 31–58): the method needs nothing beyond
  what `MMWSpalling.H` already includes before the class
  (`Numeric/Stencil.H`, `Numeric/SpCriterion.H`, AMReX headers, `<algorithm>`,
  `<cmath>`, `<limits>`, etc.). The partial header therefore needs NO `#include`
  of its own.
- **Class span**: `class MMWSpalling :` opens at line 62; the class closes at the
  `};` shortly before `} // namespace Integrator` / `#endif` at the end of the
  file. The `#include` of the partial header must sit **between** those, at the
  method's original spot.
- **Include-path convention**: the build compiles with `-I./src/`
  (Makefile rule lines 161/169), and existing includes use src-relative paths
  (e.g. `#include "Integrator/HeatConduction.H"`). So the new partial header is
  included as `#include "Integrator/MMWSpalling/Removal.H"`. A subdirectory under
  `src/Integrator/` is already an established pattern (`src/Integrator/Base/`).
- **Build safety (verified in Makefile lines 53–63)**: objects are built only
  from `*.cpp` at `find src/ -mindepth 2` and top-level `*.cc`
  (`src/mmwspalling.cc`). `*.H` files are **never** compiled as standalone
  translation units — they are include-only and tracked for incremental rebuilds
  via the generated `.cc.d`/`.cpp.d` dependency files (`-MM`). A partial header
  that is only meaningful when textually included inside the class is therefore
  safe: it will not be standalone-compiled, and editing it will still correctly
  trigger a rebuild of `mmwspalling.cc` because it appears in the `.d` deps.
- **C++ complete-class context** (why this works with zero forward-declaration
  fuss): inside an inline member-function body the entire class is treated as
  complete, so the moved method may freely call members declared *later* in the
  class (`RefreshKappaEff`, `UpdateSurfaceMaskFromPhi`, `PrincipalStressRatio`,
  `removal_enabled()`) and reference fields declared at lines 3773+ — exactly as
  it does today. Relocating the body to the same point via `#include` preserves
  this.

## Goal

Numbered deliverables:

1. **Create** `src/Integrator/MMWSpalling/Removal.H` containing **only** the
   `UpdateRemovalAfterCohesive` member-function definition — its leading doc
   comment (from line ~2092) through its closing `}` (line ~3053) — copied
   **verbatim** (preserve indentation and blank lines exactly). The file must
   contain **no** `#include`, **no** `namespace`, and **no** `class` wrapper: its
   contents are pasted directly inside `class MMWSpalling { ... }`. An optional
   include guard (`#ifndef INTEGRATOR_MMWSPALLING_REMOVAL_H` / `#define` /
   `#endif`) wrapping the method is allowed (preprocessor only) but not required.
2. **Edit** `src/Integrator/MMWSpalling.H`: delete the method body that now lives
   in the partial header (lines ~2092–3053) and replace it, **at the same
   location inside the class body**, with:
   `#include "Integrator/MMWSpalling/Removal.H"`. Keep one short comment line
   above the include naming what it pulls in (e.g.
   `// UpdateRemovalAfterCohesive — extracted to keep this file navigable.`).
   The surrounding members (`UpdateSpClusters` above, `BuildEnthalpyTable` below)
   must be untouched.
3. **Do not** move or edit any data members, config scalars, other methods, or
   the `Advance()` call site. R1 is a single-method code-motion.
4. **Build**: `make -j8` (env per `## Commands`) completes and produces
   `bin/mmwspalling-3d-g++`.
5. **Regression-verify byte-identity** on the removal-exercising tests (see
   `## Commands`). Record results in the completion notes.

## Guardrails

- **Inheritance graph is untouchable.** No edits to the base-class list, virtual
  inheritance, or any `HeatConduction`/`Mechanics` wiring (CLAUDE.md Lessons #1).
- **Verbatim move.** Not one character of the method body may change — same
  whitespace, same comments. The diff is the proof of behavior preservation:
  `git diff` should show a pure deletion from `MMWSpalling.H` and a pure addition
  to `Removal.H` with identical content. If the implementer is tempted to "clean
  up" anything inside the method, STOP — that belongs to a later step.
- **Include placement.** The `#include "Integrator/MMWSpalling/Removal.H"` goes
  **inside the class body**, at the method's original position — NOT at the top
  of the file with the other `#include`s, and NOT outside the class. Placing it
  anywhere else will fail to compile or change semantics.
- **Partial header content.** No `#include`, no `namespace`, no `class` in
  `Removal.H` — it is a class-body fragment. All headers the method needs are
  already included by `MMWSpalling.H` before the class.
- **No data-member / helper relocation in R1.** `phi_mf`, `removed_mf`,
  `spall_*`, `vapor_*`, `regime*`, `RoP*`, `h_spall_field_mf`,
  `crack_event_count_mf`, their config scalars, and `PrincipalStressRatio` stay
  exactly where they are. (Deferred to R2.)
- **Do Not Touch** (CLAUDE.md): `ext/`, `bin/`, `obj/`, `build/`,
  `compile_commands.json`, `configure`, `LICENSE`, and unrelated integrators.
  Do not touch `src/Numeric/*` or `src/BC/Operator/Elastic/ZloRoller321.H`.
- **Acceptance bar = byte-identical regressions.** Because this is pure code
  motion, every listed test must reproduce its existing baseline:
  - `sp_rossi_damage_profile` — **canonical witness** for the `sp_weibull`
    prescribed_T removal path (P1–P4: 528 cells, 96.0% GB, max h_spall 1.5564
    mm; pre-existing R1 Rossi-band FAIL is expected/deferred, unchanged).
  - `hu_end_to_end` — `damage_law` removal (granite onset 40.0 s, sandstone
    90.5 s; both fire spall_event).
  - `spall_event` and `regime_low_high_power` — synthetic spall/vapor removal
    paths (closed-form `h_spall` / phi-shift checks).
  - `sp_v_n_regime` — `sp_weibull` vapor (local-Q) path through removal
    (524 firing-cell observations, Pearson r ≈ +0.988).
  - `sp_kant_onset` verification case — confining/onset witness (461.28 °C).
  A build that links `bin/mmwspalling-3d-g++` plus these tests at baseline is the
  full acceptance condition. No physics numbers should move at all.

## Commands

```bash
# Build from /Users/tzetze20/amr_tools/alamo.
EIGEN=$PWD/ext \
  CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  make -j8

# Confirm the binary exists and is freshly built.
ls -l bin/mmwspalling-3d-g++

# Removal-path regressions (4 MPI ranks by default on this machine).
# Each test runs its own mpirun internally; use the venv python for yt/numpy.
VENV=/Users/tzetze20/Desktop/code/.venv/bin/python
$VENV tests/MMWSpalling/sp_rossi_damage_profile/test
$VENV tests/MMWSpalling/hu_end_to_end/test
$VENV tests/MMWSpalling/spall_event/test
$VENV tests/MMWSpalling/regime_low_high_power/test
$VENV tests/MMWSpalling/sp_v_n_regime/test
$VENV tests/MMWSpalling/sp_kant_onset/test   # verification case is the witness

# If you prefer to drive a sim manually:
mpirun --oversubscribe --bind-to none -np 4 \
  bin/mmwspalling-3d-g++ tests/MMWSpalling/spall_event/input
```

## Claude completion notes

### Files changed

- **New** [src/Integrator/MMWSpalling/Removal.H](src/Integrator/MMWSpalling/Removal.H)
  (975 lines): the `UpdateRemovalAfterCohesive` member-function definition
  (doc comment + body, original lines 2092–3053) moved **verbatim**, wrapped in
  an `#ifndef INTEGRATOR_MMWSPALLING_REMOVAL_H` guard plus an 8-line header
  comment explaining it is a class-body fragment (no `#include` / `namespace` /
  `class` of its own). The method body is byte-identical to HEAD.
- **Edited** [src/Integrator/MMWSpalling.H](src/Integrator/MMWSpalling.H):
  deleted the method body at lines 2092–3053 and replaced it, at the same
  position **inside the class body** (now line 2095), with a 3-line comment +
  `#include "Integrator/MMWSpalling/Removal.H"`. Nothing else in this file was
  touched by R1. The sole call site (`Advance()`, line 737) is unchanged.

### Byte-exactness proof (before build)

Verified by explicit range-diffs of the rebuilt working tree against
`git show HEAD`:
- head lines 1–2091 identical to HEAD (`HEAD_OK`);
- tail (new 2096–3419) identical to HEAD 3054–4377 (`TAIL_OK`);
- moved block (`Removal.H` lines 12–973) identical to HEAD 2092–3053
  (`BLOCK_OK`);
- main-header brace balance 260/260; include present at line 2095 inside the
  class. So the change is provably *remove block + insert 4-line include*, zero
  content change.

### Build

`make -j8` (project env) links `bin/mmwspalling-3d-g++`. Warnings are all
pre-existing: three `[=, this]` C++20-capture warnings (`MMWSpalling.H:1010,
2139, 2188`, unchanged) and one `unused variable 'D'` now reported at
`Removal.H:334` — the same warning from the verbatim-moved body, just at its
new file location, NOT a new issue.

### Tests run (against the clean-R1 binary, before the concurrent edit below)

| Test | Result | Notes |
|---|---|---|
| `spall_event` | **PASS** | h_spall closed-form exact (0.25), phi errors 0.0, no-event state matches IC. |
| `regime_low_high_power` | **PASS** | phi shift errors 0.0 (both regimes); RoP_vap/h_vap mean 0.25 exact. |
| `sp_v_n_regime` | **PASS** | regime_field==regime on all firing cells; radial spread 82.84%, Pearson r +0.988 (both byte-identical to Step 19 baseline). The aggregate "firing-cell observations" count is MPI-reduction-order nondeterministic — same binary gave 513 then 524 on back-to-back runs; the 524 reproduces the Step 19 number. Not a regression signal. |
| `sp_rossi_damage_profile` | **byte-identical on all witnesses** | P1 528/1024 columns, P2 528 cells removed, P3 96.0% GB, P4 max h_spall 1.5564 mm, peak bin centre 705 µm, GB/IG 507:21 — all exact vs Step 19 baseline. The script's overall "FAIL" verdict is solely the pre-existing, deferred Rossi-R1 band (705 µm ∉ [100,200] µm; ROADMAP Known Note + rossi-validation-diagnostic-design.md §11), unchanged. Canonical witness for the `sp_weibull` removal path. |
| `hu_end_to_end` | **PASS** | granite onset 40.0 s / LRST 787.82 K / 236 cells; sandstone 90.5 s / 907.27 K / 300 cells; ordering correct — exact Step 19 baseline. Canonical `damage_law` removal witness. |
| `sp_kant_onset` (verification) | **PASS (461.28 °C)** | Ran the single `input_verification` sim through the new binary (480 steps to t=48 s, no abort, 97 plotfiles) and scored onset with the test's own `onset_homogeneous` definition: t=43.21 s, dT=461.28 °C — exact Step 19 baseline. The full 25-case sweep (killed for time in Step 19) was not run; this case also leaves `spall.enabled` unset, so `removal_enabled()` is false and the moved method early-returns — it does not exercise the R1-touched code, and R1 left the Sp/onset path untouched. |

### Concurrent out-of-scope edit by the user (flag for the reviewer)

While R1 was building/testing, the user manually added a separate **Weibull
mesh-objectivity** feature to `MMWSpalling.H` (new `weibull.V0` parser key,
`weibull_V0` / `weibull_A_crit_warned` members, a `vol_factor =
(V_cell/V0)^(1/m)` weakest-link size-effect scaling of the per-cell flaw length
in `InitializeFlaws`, and a `connected_cluster` A_crit resolvability warning).
All of it lands in the **tail** of the file (`InitializeFlaws`, member
declarations, `Parse`) — outside the moved block — so it does not conflict with
the R1 include or `Removal.H`, and R1's integrity is intact. **This feature is
NOT part of R1**: the reviewer should attribute those diff hunks to the user,
not the extraction. At the user's request I rebuilt the **combined** tree (R1 +
their edit) and it compiles/links cleanly (same 4 pre-existing warnings). I did
NOT re-run the regressions against the combined binary and make no byte-identity
claim for it: their feature intentionally rescales flaw length at meshes ≠ 1 mm
(`V0` default = (1 mm)³, so it is exactly unity only at 1 mm), so the
`sp_weibull` flaw-dependent tests (e.g. `sp_rossi` at 1.5625 mm) are *expected*
to move once their feature is exercised — that is their feature working, not an
R1 regression. The R1 byte-identity results above stand against the clean-R1
binary that predates their edit. Verifying/own the Weibull feature is out of R1
scope.

## Implementation takeaways

1. **Partial-header-inside-the-class-body is a clean, zero-risk split for a
   header-only ALAMO integrator.** The build (`Makefile` lines 53–63) only
   compiles objects from `*.cpp` (mindepth ≥2) and top-level `*.cc`; every
   `*.H` is include-only and tracked for incremental rebuilds via the generated
   `.cc.d`/`.cpp.d` deps (`-MM`). So a fragment that is only meaningful when
   `#include`d inside `class MMWSpalling { ... }` is never standalone-compiled,
   yet editing it still correctly triggers a `mmwspalling.cc` rebuild. The
   fragment carries an include guard but no `#include`/`namespace`/`class`.
2. **Scripted line-range extraction beats hand-retyping for a ~960-line move.**
   I built `Removal.H` and the new main header with `sed -n` range copies, then
   *proved* byte-exactness with three range-diffs against `git show HEAD`
   (head, tail, moved block) plus a brace-balance check — far safer than a
   manual `Edit` of a 960-line `old_string`, and the diffs are the acceptance
   evidence.
3. **`git diff --stat` is misleading for a large block move.** The diff showed
   "57 insertions / 963 deletions" even though R1 added only 4 lines — git's
   alignment heuristic re-pairs short common lines (`}`, blanks) around a big
   deletion. Trust explicit range-diffs over the stat. (Here the apparent extra
   insertions turned out to be the user's concurrent edit; the lesson stands
   regardless.)
4. **Complete-class context makes the move trivially safe.** Inside an inline
   member body the class is treated as complete, so the moved method's calls to
   members declared later (`RefreshKappaEff`, `UpdateSurfaceMaskFromPhi`,
   `PrincipalStressRatio`, `removal_enabled()`) and references to fields at
   3773+ resolve exactly as before — relocating to the same point via `#include`
   preserves it.
5. **Some MMWSpalling aggregate diagnostics are MPI-reduction-order
   nondeterministic.** `sp_v_n_regime`'s "firing-cell observations" count varied
   513↔524 across two runs of the *same* binary. When checking byte-identity,
   anchor on the discrete/exact gates (cell counts, GB %, interpolated onset
   °C, closed-form errors), not on threshold-sensitive aggregate sums.
6. **R2 (next refactor step) is now unblocked**: relocate the removal *data
   members* + config scalars (and `PrincipalStressRatio`) alongside the method,
   or proceed to the next field-ownership cluster. Note the new
   `weibull_V0`/`weibull_A_crit_warned` members the user added belong to the
   LefmSp cluster, not Removal.
