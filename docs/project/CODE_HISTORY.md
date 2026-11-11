# Code history: comments moved out of the MMWSpalling sources

Moved **verbatim** by P2-30 (2026-11-10), per `tests/MMWSpalling/studies/p2_30_structure/COMMENT_SORT.txt`. Each block is listed under the slice it left and the member it stood at; the source no longer carries it. The packet that wrote a block is in its first words; the code it describes may be long gone.

## MMWSpalling.H

### TagCellsForRefinement

```cpp
        // P2-28: the AMR surface band (surface.refine_depth / refine_height) is retired with the column path
```

## Parse.H

### Parse

```cpp
        // P2-29: `solver` (explicit) and `alpha` (0 0) are retired at their defaults (the table above): one heat
        // path, forward Euler; the thermal expansion is each phase's beta.
```

```cpp
            // P2-27 stage 3 (census C1b): D, kappa_eff and k_eff -- the damage
            // field, kappa_phase under another name, and a uniform bulk average --
            // are no longer registered: every later plotfile component moves up
            // three places (a declared re-pin of both gates, component-audited).
```

```cpp
            // P2-28: surface.follow_mask is retired at 1 (the table above): the mask column top is the
            // face list's `top` class, and removal requires the list, so the switch IS the list's.
```

```cpp
                // P2-28: surf_cols_robin / surf_q_mean (the column reconstruction's witnesses) are deleted
```

```cpp
            // D2z-1 Goal 4 (revised): NO firings_per_void column. Under
            // whole-cell it is 1.000 by construction, so it carries no
            // information, and it appended after the eight jet_* columns and
            // broke jet_enthalpy criterion (a). The invariant lives in the
            // n_voided == 1 assertion in the column path (deleted, P2-28) instead.
```

## Model.H

### UpdateModel

```cpp
        // P2-29: the homogeneous branch (the global alpha[] x eta eigenstrain) is deleted: every deck states
        // its rock as phases (microstructure.nphases is required), so the branch above is the only one.
```

## Heat.H

### Advance

```cpp
    // (P2-29: the const-alpha, material and implicit kernels are deleted; every deck states its rock as phases).
```

### AdvanceMicrostructure

```cpp
                        // P2-27: h_gb0 (1 - D_face) with D = 0 on every deck since
                        // the damage lineage was retired
```

## Flame.H

### FeetDescent

```cpp
    //   max / mean: the whole-annulus reduction (one "pad").
```

```cpp
        // P2-29: pads is the feet rule (max / mean, closed with no live deck, are deleted)
```

## BeamProps.H

### BeamEnergyBalance

```cpp
    //   (b) RETIRED in D2q-1. It asserted phi - n*dz in [0, dz), i.e. that phi
    //       and the `removed` mask agree -- that phi is a valid void oracle,
    //       which is exactly the assumption D2q-1 removed.
```

## ParseRemoval.H

### ParseSpallSettings

```cpp
        // ------------------------------------------------------------------
        // P2-28: spall.surface_T is retired at `reconstructed` (the table in Parse): the face list's
        // R1 IS the reconstruction; site A's cell-average and diagnostic forms are deleted.
```

```cpp
        // P2-28: spall.removal_energy is retired at `flake` (the table in Parse): the face list's R3
        // IS the flake remap.
```

```cpp
        // P2-27 (census C1a): the damage-law removal's keys (spall.damage_threshold,
        // xi_s, C_h, alpha, t_spall, prescribed_ratio.*, sample, shallow_depth) are
        // retired with it -- see the abort list and the retired-at-default table.
```

## Shared.H

### (file head)

```cpp
//
// WHAT LIVES HERE AND WHY. The pieces BOTH removal paths call: the column path
// (the column path (deleted, P2-28), `surface_patch.face_list = 0`, retiring) and the line path
// (FaceRemoval.H, `face_list = 1`, production). Until P2-18 they were declared
// at the top of the column path (deleted, P2-28), so deleting the column file would have taken the
// production path's own helpers with it. Each helper body is P2-9's, moved
// verbatim (clang contracts a*b + c only within ONE expression under
// -ffp-contract=on, so a textual move keeps every frozen tree byte-identical),
// plus the two prologue pieces the face list silently relied on the column path (deleted, P2-28) for:
// the step guards and the per-step reset of the regime / rate fields.
    // The column path (the column removal step (deleted, P2-28), key off) and the line path
    // (FaceRemoval.H, surface_patch.face_list = 1) call these and nothing else,
    // so on a +z line (mu = 1) the two run the SAME reconstruction, scan and
    // remap on the same inputs (P2-9 a2 (i)). Each body is the column path's
    // former lambda moved expression for expression: clang contracts a*b + c
    // only within ONE expression (-ffp-contract=on), so moving whole
    // expressions keeps the key-off tree byte-identical (a1).
```

## Fields.H

### ConvertTtoH

```cpp
        // P2-29: one path -- the per-cell phase table (this->phase_tables[pid]); the global material's table and
        // the const-alpha closed form are deleted
```

## SurfacePatch.H

### SurfacePatch

```cpp
        // P0b deliverable 1: `surface` is a FOURTH form, not a branch that
        // overrides the other three. Revision 1 selected it whenever
        // spall.surface_T was applied, which would have made the `face` and
        // `cell` control legs identical to the scored leg by construction and
        // destroyed the comparison the gate depends on.
        // P2-28: `pinned` (the column pin rule) is retired -- it aborts by value.
```

```cpp
                    // P2-29: surface_patch.foot_rule is retired at `pads` (the table in Parse)
```

## State.H

### MineralPhase

```cpp
    // (`damage_law` was retired at P2-27 and aborts by value.)
```

```cpp
    // ⚠ THE DEFECT THIS REPLACES. Under `spall.stress_source =
    // local_thermoelastic` (the Meier and Kant configuration) BOTH criterion
    // sites read the TOP CELL AVERAGE as the surface temperature:
    //   * site A (deleted, P2-28) -- sigma_xx is flat at T_a(i,j,k) over the
    //     whole cell (a0 = 4-20 um << dz, so every K_I sample lands in it) and
    //     K_Ic is evaluated at the same value;
    //   * the column path (deleted, P2-28)'s depth scan (T_at_zv) -- its surface lerp has
    //     T_face = T_top, i.e. the top HALF-CELL is isothermal at the average.
    // ⚠ THREE modes, not two, and the third is what makes the second
    // FALSIFIABLE. A threshold criterion fires when ITS OWN temperature
    // variable reaches the threshold, so "the reconstruction fires at 462 K on
    // every mesh" is nearly a tautology on its own -- the cell-average
    // criterion fires at 462 K of CELL AVERAGE on every mesh just as flatly.
    // What actually differs is the PHYSICAL STATE at firing: the cell-average
    // criterion fires the rock at a true surface temperature that moves with
    // dz. `diagnostic` computes and logs the reconstruction while leaving the
    // criterion on the cell average, so that number can be MEASURED on the
    // control leg instead of inferred.
    // P2-28: the three P0 flags (sp_surface_T_recon / _apply, sp_removal_energy_flake) are deleted:
    // the face list always reconstructs (R1) and always removes a flake (R3).
```

### ParseSpallationModel

```cpp
            // P2-28: weibull.detach_mode is retired at per_face, weibull.A_crit aborts (Parse)
```


# Comment review, 2026-11-10 (after P2-30 was accepted)

Every comment in the class, the beam file and the rock-property files was read against the code beside it and
corrected where it was no longer true. The texts below are the comment lines that were **replaced or removed**,
verbatim, per file (the replacement text is in the source). No code changed: the comment-stripped text of all 20
files and the preprocessed token stream (2,096,272 tokens) are identical to the pre-image
`tests/MMWSpalling/studies/p2_30_structure/src_pre_5/`, and the six source-reading tests print the same 258 lines.

### src/Integrator/MMWSpalling.H

```
// Step 4b: + material-driven T-dependent ρ(T), κ(T), c_p(T) and time-varying
//           beam power. When a `material` block is supplied (Zhang granite or
//           any future subclass of Numeric::Material::Material), Advance uses
//           a conservative flux-divergence form `∂H/∂t = ∇·(κ(T)∇T) + Q − q_loss`
//           and inverts H↔T through a one-time tabulated map that bakes in
//           the smoothed latent-heat plateaux. Without a `material` block the
//           integrator falls back to step 4's constant-α path.
// plateaux — see CLAUDE.md). H is the conserved variable; T is recovered by
// inverting H each step. T is held at T_m or T_v across phase plateaux, so
// adjacent cells in a plateau have ∇T = 0 and latent heat is correctly
// absorbed without conduction "leaking" through the front.
// face list, the ONE removal path since P2-28 deleted the column path).
```

### src/Integrator/MMWSpalling/BeamProps.H

```
// the BeamSource call in AdvanceMicrostructure). With the removal
// invariant n_above == floor(phi/dz), but the packet formula
// (floor(phi/dz) + 1/2)*dz is NOT evaluated from phi: at a round-off tie
// the column shift leaves e.g. phi = -5e-18 (voided) above phi = dz
// exactly (top solid), floor gives 1, and the whole top slice (98% of the
// column power at alpha*dz = 4) is lost -- seen on the Meier dev harness,
// where it stalled one column for 22 s. Counting solid cells from the
// `removed` mask is exact by construction. Do NOT use
// ⚠⚠ THE PROPERTIES ARE STATIC PER-CELL FIELDS FILLED AT INIT.
// `kappa_mf` and `beta_mf` are each written at exactly ONE site
// (`InitializeMicrostructure`), so temperature dependence is NOT a
// parameter swap -- without a refresh site the criterion fires on the rock
// the fields were initialised with, whatever the laws say. (a5) verified
// that premise by counting the writers by DECLARED TYPE.
//
// Called FIRST in Advance (where P2-1's wall-law state also ran), so the
// criterion and the kernel/march/ledger all see ONE set by
// construction -- the ordering P2-1 validated.
// NOT alpha(T). That is what makes the seven stress sites compute
// E(T) * INTEGRAL(alpha) / (1-nu) while still reading one field: every site
// already forms E * beta * (T - T_ref), so writing the secant here makes
// the product the accumulated strain identically. ⚠⚠ Writing alpha(T) here
// instead would give every site alpha(T)*dT -- the planner's wrong reading
// -- and R would silently absorb the factor. That is gate (b6).
// P2-27 stage 4: also called at the end of InitializeMicrostructure with
// and independent of the z box decomposition. D2q-1: every column-top
// scan in the code now uses this rule -- `removed` alone -- so the column path (deleted, P2-28)
// PASS 1 no longer additionally requires phi >= 0.
// solid (the material paths, and the microstructure path without removal);
// Cost: one (ncols*nz + nz)-double allreduce per step, ONLY on a deck whose
// P2-24: the ROCK PHASE's solidus on the microstructure path (phase.T_m left the beam deck)
// cells in (k, k_top] (D2q-1: a solid count, no longer k_top - k):
//       RESTATED in D2q-2b. It used to read "cell above is void or k is the
//       domain top <=> n == 0", a CONTIGUITY test, which an interior void
//       breaks: the cell under the void satisfies the left side with solid
//       still above it. The current form is exact at any topology because
//       k_top is itself solid, so the count vanishes only at the top cell.
// reason: the old contiguity form counted incident power twice per column
// once an interior void existed (measured 14.29 % closure error).
// Single level: lev == 0 only (multi-level unsupported; composite sums
// would double-count covered cells).
// follow (studies/p2_21_loss_table/step0). Unbuilt: the closed form, as before.
// (P2-20: alpha no longer depends on T, so no temperature is read
// here -- the call is valid on every advance path.)
// tally below. This replaces the old `is_top_solid` test --
// "the cell directly above is void, or this is the domain top"
// -- which is a contiguity test and counts MORE THAN ONE cell
// per column the moment an interior void exists: the true top,
// and again the cell sitting just under the void. That
// double-counted the incident beam power (measured: a 14.29 %
// closure error against a predicted 140/120 over-count on the
// seeded probe). Identical while void is contiguous.
// (a) RESTATED (D2q-2b item 0.3): `n_above == 0` iff
// `k == k_top`.
// The old form compared against `is_top_solid` -- "the cell
// directly above is void" -- which is a CONTIGUITY test: a
// while still having solid above, so it would be counted as
// a violation the moment per-cell removal creates one.
// The restated form is exact at any topology, because
// n_above counts SOLID cells in (k, k_top] and k_top is
// itself solid: below the top there is always at least one,
// so the count vanishes exactly at the top cell. With void
// contiguous the two forms agree, which keeps identity.
// D2q-1: invariant (b) -- phi - n*dz in [0, dz) -- is
// RETIRED. It asserted that phi and the `removed` mask
// agree, i.e. that phi is a valid void oracle, which is
// exactly the assumption this step removes. It reads 0
// violations under both forms today, and n_viol feeds the
// `beam_invariant_viol` thermo.dat column, so dropping it
// is byte-identical.
// with n_above = solid cells above in the column (== floor(phi/dz), since
// removal keeps the top solid cell at phi in [0, dz)).
// still evolves with the geometric distance z_0 - z regardless, since the
// beam propagates geometrically through the void cavity.
// uniform disk. The Gaussian waist evolves with the geometric z0 - z.
```

### src/Integrator/MMWSpalling/FaceList.H

```
// the conductivity exactly as the kernel forms it (kappa_phase;
// the damage factor exp(-a_D D) is retired, P2-27)
// BEFORE the march, which bills every face at this h.
// (byte-identical to the pre-S1 inline expressions).
```

### src/Integrator/MMWSpalling/FaceRemoval.H

```
// the class body of src/Integrator/MMWSpalling.H, right after the column path (deleted, P2-28). A
// R2  ONE criterion: the depth scan (the column path (deleted, P2-28)'s DepthScanH and its
//     pieces, the SAME functions the column path calls) along the line,
// way; it is not heat-driven spallation) -- the column path's block in
// the column path (deleted, P2-28) (which now skips the face list) carried onto the face list's
// exactly the vapour rule's booking. foot_mech_vol / foot_mech_cols as on
// the column path. One regime-4 row per cleared column in the event log:
// path (the column path is deleted). The prologue both paths shared
// (Shared.H), the column tops, then the lines. Called by Advance once per
// step, BEFORE the thermal swap.
// ReduceRealSum-collapsed (PASS 1b pattern): each cell has exactly
// fallback's depth and the logged diagnostics read it; nothing else)
// profile -- no longer ell_q = k_c (T_s - T0) / q_net, which is
// the window has no depth: no recession (fr_nodepth). q_net and the
// old tangent depth ell_q are kept for the log only. Nothing defaults
// to one cell -- ReconFromX's flat h_f = dz is NOT read. (o.fb is set
// on every fallback evaluation, as before, so fr_fb_* is unchanged.)
// the old depth (log only)
// nothing to recede (P2-19: was "no heated skin", q_net <= 0)
// phi follows a +z line's column, as the column path moves it
// the whole column's phi, the column path's PASS 3 shift
// the tail every removal step runs (the column path (deleted, P2-28)'s own): T back from H,
// ghosts refilled, the derived masks rebuilt
// solver's own BC (RefillTempGhosts, the column path's call too).
```

### src/Integrator/MMWSpalling/Fields.H

```
// Build the global-material H<->T table once. Microstructure tables
// are built per phase in InitializeMicrostructure(); see phase_tables.
// Surface cells are material cells (phi >= 0) that are the highest
// material cell in their column, equivalently 0 <= phi < dz. Implement
// the simple per-column scan; AMR is single-level for Step 10.
// D2q-1: this is the ONE site that could NOT be made rem-authoritative
// without changing behaviour, so it is left on phi and ESCALATED rather
// than approximated. Replacing the phi window with "solid here, void
// directly above" changed 24 plotfiles in `dev2d` -- and changed the
// trajectory, not merely the written `surface` field: the differing files
// migrate between the 40000cell / 48000cell / 56000cell output
// directories. `input_2d_dev` sets no `surface.follow_mask`, so this mask
// is a live decision path feeding SurfaceCellFlux, not a diagnostic.
// See studies/d2q1_percellvoid/RESULTS.md section "Site 7".
// D2q-2b item 0: flip the seeded index-space box to void, once, right
// inverse, with phase fractions inferred from H.
// only `T_ref`, so `T_m`/`fus_width`/`L_m`/`T_vap_lo`/`T_vap_hi`/
// `L_v` sit at Material's base defaults and
// `effective_rhocp`'s boxcars would contribute NOTHING -- the
// latent heat vanishes silently. ⚠ Inert on both production
// decks (neither sets a phase-change scalar), so the ledger
// cannot see it; `unit/latent_heat` is the only witness.
// ⚠⚠ P2-2b SITE 5 -- THE ONE LEGITIMATE SCALAR -> FIELD WRITER, and
// the ONLY place in the rock path that may read a per-phase scalar.
// Every consumer reads the FIELD; this is where per-phase values
// become per-cell ones.
// `law = grimsel`, and skips non-Grimsel phases. Under the
// shipped default the fields would then be NEVER WRITTEN, and the
// idle clock divides by `pin_rho*pin_Cp`. Together the two cover
// every phase under every law: this at init, that per step.
// ⚠ LAW-AWARE so the fields hold the law's values at EVERY instant,
// not only from the first Advance. Without this a grimsel deck
// carries scalar values between Initialize and the first
// UpdatePropertyFields -- a window nothing reads today, which is
// exactly the kind of latent inconsistency that later becomes a bug.
// above (no-op unless some phase sets law = grimsel; no counter moves).
// sweep. is_gb is preserved as a phase-boundary alias so existing
// tests (Step 5/6/6b) keep their behaviour; later CZM/damage steps
// should consume is_grain_boundary instead.
// a0_ig otherwise. The RNG is seeded from a splitmix64 hash of the
// GLOBAL cell index (i, j, k) and weibull.seed, so flaw_a_mf is
// cumulative. Microstructure path, level 0 only.
// S1: surface.follow_mask applies the convective_flame patch and the
// surface losses at the `removed`-mask column top (ColumnTopSolid)
// instead of the fixed domain top / the phi-window surface_mf.
// jet_T_exhaust; jet_T_ent on the first step). Host state, NOT checkpointed.
// sides of the ratio. See the column path (deleted, P2-28)'s pad-clearance block for why this is a
// stress at ALL SEVEN sites.
```

### src/Integrator/MMWSpalling/Flame.H

```
// D2c (foot_body_clearance = 1, pads only): each annulus column whose
// removal pass (the column path (deleted, P2-28)). Clipping values above a nearest-rank
// pad height at the next removal pass (the column path (deleted, P2-28), before PASS 1).
// q_robin is SurfaceCellFlux's h*(T_gas - T_s) with the same inputs as the
// H update (top-cell T_old and k_c, robin_form, and, for pinned, the pin
// from the pin rule (deleted, P2-28) at this h and T_gas), so the surface temperature and branch
// match the kernel: the pin rule (deleted, P2-28), which runs next, reproduces the pin
// exactly from the final flame vectors. Columns hotter than the gas
// jet_r_reach = max column r over bins with P_bin > 0. Every rank gathers the same top-cell data and runs the same
// P2-28: the march bills the face list's faces (its column flux, and the top-cell gather it
// read, are deleted with the column path; Parse requires the list under jet_closure = enthalpy)
// T_stag = T_ent + (T_nozzle - T_ent)*phi; march mass m = mdot, or
// mdot/phi under jet_entrained_mass = 1.
// update, negative flux, invariants); free surface and the pin are
// refused under path (Parse).
// the gas leaving the floor, which is the gas leaving the domain (P2-10:
// the wall march that carried it on up the bore is deleted).
//     T_gas = flame_Tg[c].
//     bin) and the line rule read F1-F4 and nothing else.
// (di, dj, dk, index into FrPathLen); |o|^2 <= 16
// weights), grad by central differences of F_s. P2-8's +-1 indicator
// gradient read every step edge as a 45 deg corner: +5 to +10 % area at
```

### src/Integrator/MMWSpalling/Heat.H

```
// temperature BEFORE the criterion and the conduction dispatch read
// Step 6b: every path evolves H_mf as the conserved state, so
// dispatch a single InvertHtoT here. The microstructure variant
// selects the per-phase Material::Table by cell phase ID.
//   ∂H/∂t = ∇·(κ_eff(x,D) ∇T) + Q − q_loss/dz
// where κ_eff(x,D) = κ_phase(x) · exp(−α_D · D). Face conductivities are
// built from the harmonic mean of the cell-centred κ_eff values, with an
// optional grain-boundary contact-conductance term `h_gb(D) = h_gb0·(1−D)`
// phase IDs.
// Numeric::Material::Table (built from the phase's Constant material).
// For phases without phase change this reduces to T_new = T_old +
// dt·(div+Q-loss)/(rho*Cp), matching the Step 6 behaviour exactly.
// Perf (2026-09-18): phase, grain_id and kappa_phase are written only
// every step. (D, refilled every step until P2-27, is deleted with the damage lineage.)
// Step 10: if spall is active, removed cells are void and must not
// through ParmParse-style ghost reads.
// `removed` reads 0, the list and its column tops are built as with removal)
// or (S1) with surface.follow_mask / the energy ledger on.
// patch, the side faces and the losses (coarser levels: conduction
// Step 11b: the Hu onset validation uses a homogeneous one-phase
// microstructure (until P2-27 so its DP damage fields existed). Apply the same
// surface patch here as the deleted material kernel did (P2-29): prescribed_T injects the
// kernel, the pin rule, the ledger and the diagnostics.
// S1: with surface.follow_mask the patch and the losses sit on
// the mask column top; otherwise the fixed domain top and the
// phi-window surface_mf (unchanged).
// divergence term is 2 kc (T_patch - T_c) / dz^2. Until P2-28 this read
// order, so refinement hid it). Under the face list the term is the
// P2-8 F4: the column-top +z face keeps the surface form, x w
// (gain AND losses: both scale with the face's area).
// mask. Removal keeps solid <=> phi >= 0 with phi linear down
// the column, so the top solid cell has phi in [0, dz) and a
// cell m cells below it has n_above = floor(phi/dz) == m; its
// path interval [m*dz, (m+1)*dz] tiles the column with no gap
// and no overlap for ANY phi_top. Passing raw phi instead left
// n_above is taken from the mask (k_top - k), which equals
// floor(phi/dz) except at round-off ties phi_top ~ 0 / ~ dz,
// where the floor picks the wrong cell (see SolidFacePath).
// H stays frozen at its value when voided (removal never
// writes H).
//   constants: flame_h = h_conv, flame_Tg = T_flame (bit-identical path).
// patch radius) has s <= 0 (nozzle in the rock), or if an in-patch column
// gets a non-finite / non-positive h or T_gas. patch_min_standoff is the
// final there). P2-26 (census B3): the D2v above-plane pair and
// the D2aa pad shield that reached this site are deleted.
```

### src/Integrator/MMWSpalling/Ledger.H

```
// give (38 fields, none of them a side flux). Sized only when the key
// P2-8 F4: the top face x w, exactly as the kernel weights it.
// P0b witnesses for robin_form = surface -- a FLAME patch's witnesses
// (registered only with the patch on; they read flame_h, which is empty
// without it: P2-23)
// rev7 d2: the per-column side-power log. Deliverable 2 is explicit
// that the wall/riser split "cannot be a re-analysis" -- the
// plotfiles carry 38 fields and not one of them is a side flux, and
```

### src/Integrator/MMWSpalling/Model.H

```
// elastic eigenstrain below reads cell-centred T through
// CellToNodeAverage. With Neumann zhi on temp_mf, the ghost row equals
// the topmost real cell, so the surface node sees a half-cell-lagged
// temperature instead of the prescribed Dirichlet sp_T(t). Override
// the top-node temperature inside the patch with sp_T(a_time) so that
// the eigenstrain at z=z_top reflects the BC Hu prescribes.
// reads solved temp_mf through CellToNodeAverage with no override.
// P0b: declared here, ahead of its first use as a parameter type in
// the column reconstruction (deleted, P2-28). An in-class member function BODY is a
// the face list's old tangent depth k_c (T_s - T0) / q_net, logged
// only (FaceRemoval.H); 0 everywhere else
// flake thickness min(ell_eff, dz) [m]
// no usable gradient; T_surf == T_top
// per path when the loss table varies and the beam is on. Indexed
// Declared here, ahead of its first use as a parameter type (as SurfaceRecon).
// P2-21: a cell's (tau_top, tau_bot) for the kernels, captured by value in a
// ParallelFor. Unbuilt => -1, -1 => BeamSource's closed form.
```

### src/Integrator/MMWSpalling/Parse.H

```
// description per rock (beam.eps_r, beam.tan_delta) now gives both, alpha =
// P2-29 (census D): a deck without microstructure.nphases ran a separate kernel; there is none now
// Step 6d: optional prescribed-temperature surface patch (Hu LRST
// conduction validation). When `surface_patch.enabled = 1`, a
// circular patch on a single domain face injects a Dirichlet-like
// flux  kappa * (T_patch(t) - T_cell) / dz^2  on in-patch top-row
// cells. Default disabled keeps every prior MMWSpalling test
// unchanged. Currently only the explicit material path consumes
// this term; the implicit path hardcodes Neumann domain BCs.
// Step 10: optional surface advancement and spall detachment.
// Disabled by default so existing static-surface behavior is
// preserved. When enabled, registers a level-set `phi`, a moving
// `removed` mask, and per-step event diagnostics, and reuses the
// microstructure path for resetting detached cells.
// phi / removed with the spall path. When enabled, registers
// RoP_vap, vapor_event, vapor_thickness, regime, and unified RoP
// diagnostics.
// Surface mask: 1.0 on the topmost z-row of each level, 0.0 elsewhere.
// Step 11 will repurpose this as a level-set-driven moving surface.
// Step 5/6: mineral microstructure. Activated when
// `microstructure.nphases` is present, so existing tests (Stefan, beam,
// equilibrium, Zhang) keep their homogeneous behaviour. Each
// target volume fraction and the per-phase thermophysical / mechanical
// properties used by later steps (heterogeneous κ, ρ, Cp, E for steps 6
// and 7; β, A_D for step 8; α_attenuation for the Beer-Lambert source).
// or a user-supplied amrex::Parser expression of (x,y,z) (added
// in Step 6 so the heterogeneous-kappa test can prescribe a
// deterministic laminate without random seeding).
// object. When this was inside the `grimsel` branch, a leg
// got a DEFAULT-CONSTRUCTED law, so
// `<pfx>.grimsel.alpha_cap_theta` was silently ignored and
// (c7)'s two alpha fits would have come back IDENTICAL,
// Defaults disable phase change so existing Step 6 input
// files (e.g. heterogeneous_kappa) remain pure linear-H
// materials. Latent heats are physical specific heats
// a Util::RealFillBoundary call; properties are pure cell-local so
// they are written with 0 ghosts.
// Util::RealFillBoundary. `is_phase_boundary` is the renamed
// semantics of the old `is_gb` (kept as an alias);
// `is_grain_boundary` is the new true-grain neighbour flag.
// never updated. It replaces `phases[0].T_ref` -- a DECLARED input
// (`unit/gradient_blind` types 293.15 against an 820 K IC, so the
// model booked 527 K of superheat that no heating produced).
//   * `T_ref_phase_mf` is read by the EIGENSTRAIN (:1298), so
// P2-28: and under the face list's evaluate-only mode (list_state)
// phi is a level-set state: phi = z_surface(x,y) - z_cell.
// Material cells have phi >= 0; void cells (already removed)
// have phi < 0. Standard AMReX FillPatch preserves it across
// Step 17: cell-resolved regime field under sp_weibull. The
// scalar `regime_mf` remains the per-cell single-regime
// indicator (registered before regime_field; kept); `regime_field_mf` is
// the canonical per-cell diagnostic for the sp_weibull
// per-cell v_n logic — distinct so plotfile consumers can
// tell the two diagnostics apart.
// surface cell). P2-28: site A's fields (Sp_field, Sp_cluster_id, h_spall_field,
// crack_event_count, T_surf_field, ell_eff_field) and the recon_* columns are deleted with it.
// P2-27 (Oct census B14): every cross-key ABORT of this block now lives
// S1: surface-resolution study switches (default-off; microstructure
// path, single level).
// P0b deliverable 1: robin_form = surface. Same prerequisites as
// pinned -- it is the same closure with the surface temperature
// computed rather than remembered -- PLUS the reconstruction
// itself, without which the form has no surface to sit on.
// P0b: appended last, and ALSO gated on the reconstruction --
// a tree input with the ledger on but no P0b key must not gain
// header alone whatever the numbers do. P2-28: the reconstruction's one home is the face list.
// D2a: the spatial-flame diagnostic, after every pre-D2a variable (P2-28: the pinned
// branch counters and patch_P_robin -- the column patch diagnostics -- are deleted).
// (its requirements -- the face list, the reconstruction -- are checked in
// CheckInterlocks).
// already been parsed. Each condition and message is VERBATIM from where it
// used to sit, and they run in the order they used to fire: the sub-parse
// checks (ParseSpallSettings, ParseBitSettings, in Parse's call order; the
// two ParseSpallationModel ones died with the damage lineage, P2-27 stage 2),
// then Parse's own end-of-parse block, then the face
// list's cubic-cell requirement (a per-call runtime check in BuildFaceList
// until P2-27; geom is known at Parse). ONLY ABORTS live here: every key
// ⚠ A deck that breaks TWO rules can now meet a different one first: a
// sub-parse that used to abort early now runs on to its later key queries
// (bit.enabled with spall and vapour off still reads bit.feed_rate, which
// is required). One broken rule fires exactly as before.
// P2-27 (census C1a): with the damage-law removal deleted, the only
// P2-28 (census H): ONE REMOVAL PATH. The column path is deleted, so removal runs on the face
// P2-28: the jet march bills the face list's faces; its column flux is deleted
// P0b deliverable 1: robin_form = surface, the same prerequisites as
// pinned PLUS the reconstruction.
// P2-28: the form's surface temperature is the face list's (R1): the column
// reconstruction it used without the list is deleted
// P2-15: the buckling condition lives on the face list only (the column path never
// reads the key -- retire, don't stack) and reads the reconstructed skin length, so
// either missing would leave the key SILENTLY INERT: both abort, naming the key.
// P2-24: regime 2 IS carried now -- the face list's vapour rule (a LISTED cell whose enthalpy has
// P2-27 stage 5 (census G): regime 4 IS carried now -- the pad clearance
// Until P2-27 this was checked in BuildFaceList on every call; the
```

### src/Integrator/MMWSpalling/ParseRemoval.H

```
// spall.enabled early-return because it also governs the Sp_field
// ONSET diagnostic (site A (deleted, P2-28)), which runs under
// spallation.model = sp_weibull even when material removal is off
// (e.g. the Kant onset verification runs). Default `mechanics`
// (stress_mf from the elastic MLMG solve) preserves existing
// behaviour. `local_thermoelastic` reconstructs sigma_xx from the
// temperature field (Kant 1-D confinement), letting el.type = disable.
// A1: every-removal event log (validated at the end of Parse, once
// spallation.model is known). P2-28: the clusters and h_col CSVs are retired with the column path.
// clamp. ⚠ The min(., dz) clamp P0 used is the measured culprit for the
// 20.6 % residual in its mesh ladder: at 0.5 mm, ell_eff (0.992 mm) is
// already larger than dz, so the "flake" was the whole cell again and
// T_bar collapsed back onto the cell average.
// ⚠ c_f IS THE ONE EMPIRICAL CONSTANT THIS PACKET INTRODUCES. Under
// `flake` the thickness is IMPOSED (the column path (deleted, P2-28) bypasses the depth scan's
// h_col), so it is not derived from fracture here. 0.3-0.5 is what
// P2-26 (census B1): parsed HERE, before the flake_h_max check below and
// before the spall-off return, so the check reads the DECK's value -- it
// used to read the member's initialiser, because this parse ran after it.
// the column reconstruction (deleted, P2-28) reads it on a spall-off onset deck as well.
// spall.h_spall_max, is parsed above the flake_h_max check (P2-26).
// space, `ilo jlo klo ihi jhi khi`, flipped to void once after
// Initialize. Absent = off, and nothing in the code path runs.
// depth-scan clamp, the three k_top scans, ColumnSolidAboveCount and
// before any mechanism creates such a state. D2q-1 shipped that
// handling with nothing able to produce an interior void, so it was
// proved inert but never proved correct.
// P1b: SEVERAL boxes, 6 ints each. Six ints is exactly what it
// was, so every existing input is unchanged. ⚠ The extension is
// not cosmetic: a standing column with void on two sides CANNOT be
// built from one box -- if (i+1, j) and (i, j+1) are both inside a
// box then so is (i, j) -- so a single box can express a pit, a
// ridge or a step, but never the tower the support rule exists to
// remove. The unit test needs that geometry prescribed, not grown.
// ⚠ `static` and taking the law by POINTER on purpose. Every caller --
// the `Sp_field` diagnostic site (an `amrex::LoopOnCpu` with a `[&]`
// capture today; it was a `ParallelFor` under `[=, this]` when this was
// written) and the depth scans -- hoists one pointer before its loop; a
// `SpCriterion::K_I`'s convention and the double-negation the pre-P2-2b
// expressions perform. The column scan's stress_mf branch (the column path (deleted, P2-28))
// returns `-sg(0,0)` for the same
// reason. A reversed bracket here would discard every hot flake and stop
// all drilling -- it is what rev 4 of this packet would have shipped.
// deck's `E/(2 mu) - 1`. On both production decks those agree at 0.30; if
// a deck ever disagrees, the law wins on the law path by construction.
// The rock law for the criterion, or nullptr under the shipped `constant`
// default. ⚠ ONE place decides this, so "is the law on?" cannot be
// answered differently at two of the three stress sites -- which is the
// failure mode `MMWSpalling.H`'s two-firing-site invariant exists to stop.
// ⚠ Phase 0 by construction: the T-derived sigma_xx path is already
// phase-0-only (it aborts on an empty `phases`), and widening that is not
// in this packet.
// Gating this on `phases[0].grimsel_law` -- which is what the first
// implementation did -- ties the criterion to conduction, so the Kant
```

### src/Integrator/MMWSpalling/Shared.H

```
// immediately before the column path (deleted, P2-28). It is a class-body fragment, NOT a standalone
// nu is outside (-1, 1). `ok` = phase 0 exists with E > 0 and mu > 0 (nu
// is 0 otherwise); each caller keeps its own abort. The four former copies
// (two in the column path (deleted, P2-28), one in MMWSpalling.H, one in FaceRemoval.H) computed
// these same expressions.
// inert. Without the law, the scalar pair as before.
// `ok` no longer asks for them: every consumer takes LawSigma then.
// kernel read the PRE-remap temperature across every domain side. This is
// scan's two former blocks (at a cell, at a depth) were textually
// identical; they are this ONE function now, so they cannot disagree.
// K_I/K_Ic at the surface (s = 0) -- the line path's score. ⚠ K_Ic is
// evaluated BEFORE K_I at every sample, as it always was: the stress
// lambda carries counters.
// along the line (see the column path's note at its call):
// P2-18 A3: THE STEP PROLOGUE both paths share (formerly the first lines
// of the column removal step (deleted, P2-28)).
// spallation_model == SpWeibull && (spall_on || vapor_on).
```

### src/Integrator/MMWSpalling/State.H

```
//   constant (DEFAULT) -- the scalar kappa/beta/E/mu/rho/Cp above, the
//                         shipped behaviour; every frozen study is
//                         byte-identical.
// ⚠ The scalars stay parsed and are what `constant` uses, so turning
// the law on is the only thing that changes behaviour.
// index-aligned with `phases`. Built once in Initialize().
// heterogeneous thermoelastic update. Cell-centred, 1-ghost; sampled by
// CellToNodeAverage in UpdateModel.
// `phi = z_surface(x,y) - z_cell` is a level-set state; cells with
// phi >= 0 are material, phi < 0 are void. `removed` is the binary
// void mask. `spall_event`, `spall_thickness`, and `RoP_spall` are
// per-step diagnostics reset at the start of each spall update.
// D2z-1: whole-cell removal is UNCONDITIONAL -- there is no key and no
// half-cell mode. Under computed-T the depth scan's surface region is flat
// across the whole top half-cell (the T_face = T_top fallback), so at a
// marginal firing the K_I crossing landed just past dz/2 and a column fired
// TWICE per cell. Now a firing voids the top cell (h = dz), so an event in
// removal_events.csv IS a removal.
//   * The depth scan is RETAINED under prescribed_T. That is CONSERVATISM,
//     not a justified split: D2z-0 takeaway 7 found the same dz/2 landing
//     in a prescribed_T control with a real gradient, which points at the
//     GATE LOCATION (Sp_top >= 1 tested at the cell CENTRE while depth is
//     measured from the FACE) rather than the lerp. Keeping it untouched
//     keeps Hu and Kant unchanged; fixing it is a Step 16d question.
//   * Applies to the SPALL regime only. Vapor sets its own h_v and the bit
//     is drawworks-set; neither had the half-cell artefact.
// in index space. Empty = off (the default), and ApplySeedVoid returns
// `spall.h_spall_max` caps the downward column scan in metres (safety
// stop; in practice the K_I < K_Ic crossing always fires first). The
// when set to 0 (default), so the depth scan uses the same Simpson
// resolution as `site A (deleted, P2-28)`'s surface Sp.
// Stress source for the Sp = K_I/K_Ic criterion (both the Sp_field surface
// diagnostic/gate AND the h_spall depth scan). Default `mechanics` samples
// the elastic MLMG solve's stress_mf. `local_thermoelastic` reconstructs
// sigma_xx algebraically from the temperature field via the 1-D-confined
// relation sigma_xx = -E*beta*(T - T_ref)/(1 - nu) (phases[0]), so the
// solve, no convergence cliff, no max_level = 0 resolution ceiling. Valid
// with the spall path. RoP_vap follows main.tex:
//   RoP_vap = (1 - R) * P0(t) / (rho * L_tot * A_beam)
//   L_tot   = Cp * (T_vap_lo - T0) + L_m + L_v
// The trigger is H >= H_vap = phase_table.H_of_T(T_vap_lo). If both
// spall and vapor want the same column in the same step, the larger
// candidate thickness wins and `regime` records the winning mode.
// sp_weibull. {0 none, 1 spall, 2 vap}. Registered only when
// listed line; nothing is removed. Until P2-28 the list refused removal off (Parse), so list_state() equals
// removal_enabled() on every deck that ran before it.
// surface flaw length a0 (m); Weibull randomization arrives in Step 15.
// compressive). Applied as +p·nu/(1-nu) offset in both sigma_xx lambdas
// (site A (deleted, P2-28) and the k_i_at_zcrack depth scan), NOT as an
// under robin_form = surface). P2-28: surf_cols_robin / surf_q_mean, the
// column reconstruction's witnesses, are deleted.
// (spall.removal_events_csv; sp_weibull + spall.enabled, level 0). One
// row per (step, column) whose surface moved: time, col_i, col_j, regime,
// k_top, T_top, a_f, Sp_top, h_scan, h_applied, n_voided. Unlike the
// column path's retired edge-triggered CSVs, it is a complete removal record.
// THE FIX (spall.surface_T = reconstructed). For an exponential profile
// ⚠ T0 IS phases[0].T_ref. That is the temperature at which the
// thermoelastic driving stress is zero, and in every input this programme
// uses it equals the undisturbed rock temperature (293.15 K in both
// validation/kant/* and validation/meier/*). It is NOT a new constant and
// it is NOT tunable here. If a future input sets T_ref away from the
// `nsolid_fb` / `flat` witnesses in the event log are what would show it.
// The ONE place the reconstruction is written. Both criterion sites call
// this, so they cannot drift apart. Pure function of its arguments --
// no state, no fields, no MPI -- which is also what makes it unit-testable
//   h_f     = min(ell_eff, dz) = dz / max(x, 1)
//                                                    = min(x, 1)
// no usable decaying gradient ("flat"), which the CALLER resolves -- under
// P0b it holds the column's last valid smoothed x (deliverable 1.3) rather
// than dropping to the cell average, because dropping gives each such
// column a x1.3 flux kick and 12.8 % of Meier firings are flat.
// Everything from a given x. Split out of the reconstruction so the
// per-column builder can substitute a smoothed or held x (deliverables
// 1.3, 5.1, 5.2) without duplicating a line of the algebra.
// is a MESH-INDEPENDENT multiple of the surface rise (P6). With P0's
// min(ell_eff, dz) the multiple was mesh-dependent and collapsed to the
// cell average wherever ell_eff exceeded dz.
// block) off a two-cell ratio of 1.0014. The assertion caught it.
// It binds on nearly-flat columns only -- h_f is 0.4 mm against a
// 15 mm cap on every 1-D column leg -- and its bind fraction is
// reported (flag bit 8).
// in this model says a preheated slab should spall less. The real fix
// is a gradient-aware criterion (buckling / energy release with a
// flake thickness), which is a P1 DESIGN ITEM. This key buys the
// Meier comparison; it does not earn it.
// BuildFlameColumns -- which runs AFTER the criterion (see the call-order
// note in the column path (deleted, P2-28)). Returns 0 when no gas temperature is defined, which
// P2-27: the default `none` is what the 48 decks that never named the key
// always had -- no criterion loaded (no Sp_field, sp.a0 not required).
// runs; the offset is only applied in the sp_weibull sigma_xx
// lambdas (site A (deleted, P2-28) + Step 16d's k_i_at_zcrack).
// only on one branch would be silently inert on the other. Its cross-key requirements
// are checked in Parse after the face list and spall.surface_T are both parsed.
// Checked against the span in Parse, beside the span's own face-list check.
// runs at :206 and `value.phases` is not resized until :293, so a
// legal leg. The check sits with the other cross-key aborts after the
// microstructure block instead.
// reason. Inside it, a run with `weibull.enabled = 0` never executed the
// parse and silently took the MEMBER default 15.0 (:9608) -- which is
// exactly what P1f's lateral-drive fixture did, while Meier's
// `input_drilling:124` sets 20.0. The flake fraction below reads `m`, so
// the fixture and the validation leg would have derived thicknesses
// 33 % apart for no stated reason.
// equals the member default at :9608, so any run that does not set the
```

### src/Integrator/MMWSpalling/SurfacePatch.H

```
//     2*kappa*(T_patch - T_cell)/dz^2 on the in-patch top row (P2-28
//     corrected the factor 2; P2-29 deleted the material kernel that still
//     had the half term). Under the face list the top faces carry it.
//     q = h_fl*(T_flame - T_surface), injected as the divergence
// keeps every Step 6d/7b/8b/11b/14 input byte-for-byte.
// P2-8 group A: `side_face_closure` (SideTForm) DELETED -- `cell` was
// set only by an abort test; the lateral closure is always `follow`.
// list (BuildFaceList / FaceListHeat / FaceWeight). 0 => every site
// runs its pre-P2-8 code.
// D2a: spatial flame expressions (x, y, r, s, t); use_expr = both given.
// patch ON they stay required exactly as before; with it OFF the mask
// S1: Robin discretisation for convective_flame (microstructure
// path only; see SurfaceCellFlux). `cell` (default) applies
// P2-28: pinned_idle_cycles governed the pin rule, retired with robin_form = pinned.
// T_ext / h_conv / T_flame. Under convective_flame, h_conv is the
// constants, never a mix. (surface.follow_mask = 1 is checked in
// MMWSpalling::Parse.)
// simulated time between rows. 0 = off, so no existing leg gains a
// file. P2-26 (census B1): parsed under EVERY jet_closure -- it sat
// inside the enthalpy branch, so a jet_closure = none deck (CAL)
// could not set it although its consumer is closure-independent.
// which must stay finite and > 0.
```

### src/Numeric/MMWBeam.H

```
// `beam.R` are retired (MMWSpalling::Parse aborts naming them). Surface
// position z_0 is constant in step 3; surface advancement is step 11.
// Beam power [W]; default zero so the integrator reduces to step 2 if
// no beam block is supplied.
// of the two derived members below is parsed before they are formed).
// aborts when they differ from the no-op and the loss table is flat (no sum to carry them).
// skipped -- a long measured schedule (oglesby's 159 intervals) no longer
// sits on two 159-token lines. The file and the arrays are mutually
```

### src/Numeric/Material/Constant.H

```
// Constant-property material: rho, c_p, kappa do not depend on T. Used by
// the microstructure path to give every mineral phase the same H<->T
// machinery as the global Zhang material (Step 6b). Phase-change scalars
```

### src/Numeric/Material/Grimsel.H

```
//   exposes lambda and lets the integrator form the diffusivity, exactly as
//   the reference demands. Its three fits agree to 5.0 % against a 6.5 %
// ⚠ H(theta) here is ANALYTIC and is the exact integral of c_p, so the H<->T
//   inversion is consistent by construction rather than by agreement.
// `unit/lateral_flux_lock` asserts that all seven stress sites take it from
// here and from nowhere else. If anyone writes `alpha(T)*dT` at one site, the
// uniform 1 K grid from 20 to 1000 degC, linearly interpolated between
// nodes. ⚠ 1 K spacing makes the table agree with a dense offline
// ⚠ Called from the constructor AND from the end of `Parse`: both
// `alpha_cap_theta` (which reshapes the strain table) and `E0` are parsed
// ignore them. `(c7)` is exactly the `alpha_cap_theta` sensitivity, so a
```

### src/Numeric/Material/Material.H

```
// Abstract base for temperature-dependent thermophysical materials. The
// ALAMO MMW spalling integrator dispatches κ(T), ρ(T), c_p(T) through this
// interface; phase-change scalars (T_m, T_v, latent heats, smoothing widths)
// are stored on the base so they are uniform across material implementations.
// Step 4b only needs Zhang/Oglesby granite. New rock types or per-mineral
// models slot in by deriving a new subclass and adding it to MMWSpalling's
// `select_default<...>("material", ...)` factory call.
// [W/(m K)]
// ρ · c_p with the smoothed latent-heat boxcars baked in. This is what
// implicit time-steppers want as their thermal-mass coefficient, since
// it is the slope dH/dT used by the H↔T table. Equal to ρ·c_p outside
```

### src/Numeric/Material/Table.H

```
// base (fusion / vaporisation widths and specific latent heats). The table
// is built once per material instance and used by:
//
//   - the global `material` path (single Zhang/whatever instance), and
//   - the microstructure path (one Table per phase, indexed by phase ID).
// fusion plateau (fus_width = 0) where T jumps discontinuously across the
// latent interval still resolves the partial-melt fraction correctly.
// by T_of_H above its top. Until P2-29 only the cached band top carried it, so a liquid cell read L/c_p too
// hot (100 K on the Stefan fixture) and ConvertTtoH dropped rho*L. No live deck had a sharp band.
// sharp-vap reference T (unused if width>0)
// For sharp fusion (fus_width == 0) the integrator does not put a
// m_H[lo+1] that the old bisection returned, so the interpolated
// T is bit-identical. Replaces ~14 data-dependent bisection steps
// per cell per step.
```

### src/Numeric/Material/Zhang.H

```
// rho(T)   = a_rho * theta + b_rho                       (T <= T_rho_cap)
//          = rho_cap                                     (T >  T_rho_cap)
// Phase-change parameters (stored on the Material base).
```
