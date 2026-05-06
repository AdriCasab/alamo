# Active Step: 10 - Surface Advancement And Spall Detachment

Status: ready-for-implementation.

Step 9 is complete and archived. It added the default-off bilinear
grain-boundary cohesive-zone law, history, and diagnostics, but deliberately
deferred traction feedback into the elastic operator. Step 10 should now add
the first material-removal geometry update: a level-set-style surface field,
spall-detachment detection, field reset in detached material, and a focused
unit regression.

## Goal

Add a default-off spallation removal capability and `spall_event` regression:

- maintain a moving surface/depth field derived from a level-set convention
- detect surface-connected fully damaged zones that also satisfy the
  Hoek-Brown spalling-ratio cutoff
- compute `h_spall = C_h * sqrt(alpha * t_spall)`
- advance the local surface by that thickness
- reset detached cells to ambient enthalpy/temperature and `D = 0`
- update the existing `surface` mask from the moving surface
- record enough diagnostics for a deterministic yt regression

This is a geometry/removal verification step, not the Hu onset validation yet.
Keep the feature default-off so existing Hu/Zhang/topology/mechanics
regressions behave exactly as before when `spall.enabled = 0`.

## Sources To Use

- Plan section: `in-main-tex-you-will-quizzical-treasure.md`, Step 10.
- `main.tex`, Hoek-Brown spalling limit:
  `sigma_1 / sigma_3 >= xi_s ~= 8`, applied only after `D -> 1` in a
  connected near-surface zone.
- `main.tex`, spall scaling:
  `h_spall = C_h * sqrt(alpha * t_spall)` and
  `RoP_spall = h_spall / t_spall`.
- Current code paths in `src/Integrator/MMWSpalling.H`:
  - `surface_mf` is currently static, set by `InitSurfaceMask()` to the top
    z-row.
  - `Advance()` order is mechanics -> damage -> CZM -> enthalpy advance ->
    `InvertHtoT()`.
  - `D_mf`, `damage_F_mf`, and `damage_Psi_mf` already live on the
    microstructure path.
  - Step 8's `damage.prescribed_stress.*` hook exists for deterministic
    damage tests, but it does not by itself provide a spall ratio diagnostic.

## Current Code State

- `surface_mf` is plotted as `surface` and consumed by surface losses. It is
  always the top domain row today.
- There is no `phi_mf`, removed-material mask, spall event diagnostic, or
  spall RoP diagnostic.
- The heat state is conserved enthalpy `H_mf`; `Temp` is recovered through
  `InvertHtoT()`. Reuse existing H/T conversion logic when resetting detached
  cells to ambient state.
- Damage is active only with `microstructure.nphases`; passive `damage.ic.*`
  still works when `damage.enabled = 0`.
- Step 9 CZM fields are diagnostics only. Do not make Step 10 depend on CZM
  traction feedback.

## Implementation Direction

1. Add default-off parser/settings.
   - Suggested keys:
     - `spall.enabled = 0|1`
     - `spall.damage_threshold` default near `0.99`
     - `spall.xi_s` default `8.0`
     - `spall.C_h` default `1.0`
     - optional `spall.alpha` override for the unit test; otherwise use the
       available thermal diffusivity (`therm_alpha` or per-cell
       `kappa_eff/(rho*Cp)` on the microstructure path)
     - optional `spall.t_spall` override for deterministic tests; otherwise
       use the event time convention you document, likely `time + dt`
     - optional deterministic ratio hook such as
       `spall.prescribed_ratio.enabled` and `spall.prescribed_ratio.expr`
       for `tests/MMWSpalling/spall_event/`
   - Reject invalid non-finite/negative parameters when enabled.

2. Add fields and a clear level-set convention.
   - Prefer a plotted scalar `phi` with convention
     `phi = z_surface(x,y) - z_cell_center`.
   - Initial surface is the domain `zhi` plane, so top cells have
     `phi ~= 0.5*dz`, bulk material has `phi > 0`, and detached/void cells have
     `phi < 0`.
   - Existing `surface` mask should be derived from `phi`: surface cells are
     material cells nearest the current surface, e.g. `0 <= phi < dz`.
   - When `spall.enabled = 0`, keep current static top-row behavior.
   - Useful diagnostics for the test: `spall_event`, `spall_thickness`,
     `RoP_spall`, and removed enthalpy/energy tally if practical.

3. Detect spall events after mechanics/damage.
   - Candidate cells must be material cells, have `D >= damage_threshold`, and
     pass the stress-ratio check.
   - For production mode, compute principal stresses from the symmetric stress
     tensor using the existing Eigen pattern (`SelfAdjointEigenSolver` on
     `Set::Matrix`). Document the sign convention because ALAMO stores tensile
     stress positive and compression negative.
   - For the unit test, prefer a prescribed ratio hook so the geometry/removal
     logic is deterministic and does not depend on elastic solver details.
   - Require a surface-connected component. A CPU flood fill / 6-neighbour
     component scan per level is acceptable for Step 10; the regression can be
     single-level. If AMR component connectivity is deferred, state that in
     completion notes.

4. Advance the surface and reset detached cells.
   - On event, compute `h_spall`, subtract it from `phi` in the affected
     surface-connected columns/region, and recompute `surface`.
   - Detached cells are those crossed by the advancing surface, e.g. material
     cells whose old `phi >= 0` and new `phi < 0`.
   - Reset detached cells: `D = 0`, `H = H(T_amb)`, `Temp = T_amb`, and refresh
     derived fields such as `kappa_eff` when present.
   - Be careful with the current `Advance()` swaps. If the reset happens before
     `std::swap(*H_mf, *H_old_mf)` / `std::swap(*temp_mf, *temp_old_mf)`, make
     sure the reset state becomes the old state consumed by the thermal update.

5. Add `tests/MMWSpalling/spall_event/`.
   - Use a small deterministic 3D single-level input with microstructure
     enabled so `D_mf` exists.
   - Prescribe a top-connected slab/patch with `D = 1` using `damage.ic.*`.
   - Use the prescribed ratio hook above threshold for the event case.
   - Include a no-event companion case: ratio below threshold or a damaged zone
     disconnected from the surface.
   - Choose `dt`, `alpha`, and `C_h` so `h_spall` is easy to check against the
     mesh, ideally one or two cell widths.
   - The Python test should verify one event, expected `phi` shift, expected
     `surface` mask movement, reset `D/H/Temp` in detached cells, no reset
     outside the detached zone, and no event in the companion case.
   - Write `output/comparison.png`, for example before/after `phi` or surface
     masks with the expected detached thickness annotated.

## Pass Criteria

- With `spall.enabled = 0`, existing behavior is unchanged and `surface` stays
  the static top z-row.
- Event case triggers exactly one spall event in the prescribed connected zone.
- Computed thickness matches `C_h * sqrt(alpha * t_spall)` within a tight
  deterministic tolerance.
- `phi` advances by the expected thickness only in the detected region.
- `surface` mask moves to the first remaining material cells.
- Detached cells have `D = 0`, ambient `Temp`, ambient-compatible `H`, and
  refreshed derived conductivity.
- No-event companion leaves `phi`, `surface`, `D`, `H`, and `Temp` unchanged
  except for normal background evolution expected by the input.
- Existing damage, CZM, topology, and thermal-stress regressions still pass.

## Verification Target

Build:

```bash
EIGEN=$PWD/ext \
  CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  make -j8
```

Run the new simulation(s) with 4 MPI ranks unless the input metadata says
otherwise:

```bash
mpirun --oversubscribe --bind-to none -np 4 \
  bin/mmwspalling-3d-g++ tests/MMWSpalling/spall_event/input
/Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/spall_event/test
```

Nearby regressions to run:

- `tests/MMWSpalling/dp_yield/test`
- `tests/MMWSpalling/gb_cohesive/test`
- `tests/MMWSpalling/grain_topology/test`
- `tests/MMWSpalling/amr_microstructure_regrid/test` if derived-field repair
  or AMR initialization changes
- `tests/MMWSpalling/thermal_stress/test`
- `tests/MMWSpalling/alpha_beta_transition/test`

If any shared surface-loss, `surface_patch`, `Advance()` ordering, or default
thermal behavior changes outside `spall.enabled`, also rerun `hu_conduction`,
`hu_breakage_index`, and the Hu thermoelastic regressions.

## Guardrails

- Keep spall removal default-off and additive.
- Do not change Step 8b's Hu indicator validation thresholds or top-node
  `surface_patch` eigenstrain fix.
- Do not reinterpret `is_gb`; CZM/damage/spall microstructure logic should use
  true `is_grain_boundary` where grain-boundary specificity is needed.
- Preserve passive `damage.ic.*` behavior when `damage.enabled = 0`.
- Keep the Step 10 regression deterministic; do not rely on a fragile elastic
  solve just to test geometry bookkeeping.
- Do not implement vaporisation or unified RoP here; that is Step 11.
- If cross-level AMR spall connectivity is not complete in Step 10, document
  the limitation clearly and keep the verification single-level.

## End-of-Task Handoff For Claude

When Step 10 is complete, update this file before ending the session:

- Mark status as completed, blocked, or completed-with-caveat.
- Add `Claude completion notes`: files changed, tests run, pass/fail, and tests
  not run.
- Add `Implementation takeaways`: exact parser keys, field names, level-set
  sign convention, event-time convention, stress-ratio convention, AMR
  limitations, and any tests or validations that should be rerun by Codex.
