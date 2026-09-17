# Active Step: D2b - Burner-on-feet descent and the scored Meier ROP / hole-diameter test

Status: completed

Revised 2026-09-16 (before implementation), after the user's review of the
first D2b packet:
- tripod pad rule replaces the column max;
- nominal foot annulus corrected to [0.028, 0.040] m;
- exhaust-recirculation entrainment is the blind-hole baseline;
- hole diameter is scored depth-mean, with a pre-registered depth profile;
- Stage A has a pre-decided tie-break;
- the mesh check is mandatory and runs before the sensitivities;
- Stage C is cut to three runs;
- the meaning of a stall is stated.

## Context

**This is the step that turns ROP into a prediction.** D2a and D2a2 moved
the nozzle on a **prescribed** path. Under a prescribed feed the centre ROP is
an identity: the face settles at whatever stand-off supplies
`rho*Cp*dT_fire*v_feed`, so it says nothing about the jet. Meier's burner
instead **rests on feet** on the rock. It can drop only as the rock under the
feet is removed, so **ROP = the recession rate of the face at the foot
ring**, where the jet flux is lowest.

**What D2a2 left (accepted 2026-09-16):**
- The energy-conserving wall-jet closure (`jet_closure = enthalpy`) with the
  free-jet stagnation decay; D2a2 recommended `jet_stagnation = decay` +
  `jet_entrained_mass = 1` before D2b. This step keeps both.
- In 3-D the centre **leads** the nozzle (`s_c` 60-90 mm, `T_stag`
  1000-1270 K at 1900 K with 293 K entrainment). Ring ROP at r = 40 mm was
  **0.5-1.3 m/h, below Meier's band**, from an early transient window; no run
  reached a steady bowl.
- All 1900 K holes were Ø 84-103 mm at 10 mm depth; at 1436 K, Ø 71-79 mm.
- Open flag: `input_drilling`'s `spall.flake_coherence_length = 0.004` caps
  walls at 63-65 deg, and capped cells overheat under the pinned flux.
- 6-23 never-fired rim columns in the bracket decay runs.

**Four physics points fixed before any run (user review, 2026-09-16):**

1. **The feet are a rigid tripod.** It rests on the plane through its three
   contact points: the highest rock under each pad, averaged across the pads.
   A max over every 2 mm column in the annulus would let one never-fired
   column freeze the burner for the whole run - a discretisation artefact
   recorded as physics. **Baseline: per-pad high quantile, mean across pads.**

2. **Hole diameter under the feet rule is mostly geometry.**
   - Any column receding slower than the ring falls behind the burner. Once
     it is `foot_standoff` (50 mm) behind, it is at the nozzle plane, gets no
     flux (the D2a2 `s <= 0` rule) and freezes.
   - A column at radius `r` therefore freezes at depth
     `d(r) ~= foot_standoff * v / (v - v(r))`.
   - With a flux falling roughly as `1/r`, `v(r)` drops a few % per mm at the
     foot radius (user estimate ~2.5 %/mm). That puts the Ø 93 wall at a few
     hundred mm depth and the Ø 85 wall below the 0.5 m block.
   - So **the hole narrows with depth toward the 80 mm burner**, and Meier's
     85-93 mm follows from the feet, the 50 mm stand-off and the profile shape
     almost regardless of `h`.
   - **A hole-diameter pass is therefore a check of the profile shape and the
     above-nozzle rule, not support for the anchor.** Score the depth-mean
     diameter and pre-register the depth profile.

3. **In a blind hole the jet entrains its own exhaust, not 293 K air.**
   - Below the first few centimetres the only gas available is the exhaust
     flowing back up. The drawing makes this stronger: the jet expands for its
     whole 50 mm inside the burner's feet collar (flange bore Ø 45, tube ID
     Ø 56), and the exhaust leaves through three slots.
   - Feeding the march's own exhaust temperature back as the entrained
     temperature is **mass and energy conservation, not a fit**. D2a2's slab
     measured ~30 % ring-ROP sensitivity to this temperature.
   - **This is the scored baseline** (new `jet_T_ent_mode = exhaust`).
   - **Stated risk:** this choice raises ROP, the direction of D2a2's
     shortfall. It is decided here, **before** any D2b run, on the confinement
     geometry alone. The 293 K case is the first Stage C run, and it is
     reported in the same verdict table.

4. **A stall is where Meier's operator pushed.** The crane drove the burner
   with "intermittent and rough axial displacement" (thesis pp. 219-220). A
   model stall marks where a real operator would have pushed weakly attached
   rock away - the old Step 19 Sp-gated bit was exactly that push. It stays
   out of the baseline. **A stall must be reported as "the jet alone does not
   clear the ring here", not as "the physics says it cannot drill".**

**No tuning.** The scored anchor is **J-M**: Martin `h` at the measured mass
flow, `T_nozzle` = 1900 K (adiabatic), exhaust recirculation. J-5 and J-10 (the
reference bracket) are reported as the bracket's implication. **If J-M misses
Meier's band, report the miss.** Do not change the anchor, `h_ref`,
`T_nozzle`, the entrainment mode, the foot geometry or rule, or the wall
treatment (beyond the pre-stated Stage A rule) to pass.

**Caveats carried forward:**
- (C1) Pin state is not checkpointed; `pinned_idle_cycles` stays **2**.
- (C1/D2a) Score by time windows and fits (nozzle height,
  `removal_events.csv` `h_applied`), never by end depth or plotfile
  staircases.
- (D2a/D2a2) **2 mm is not shown mesh-converged.** Most of the march is
  face-form uptake by unfired rock (pinned share 0.01-0.12), which S1/C1
  showed is mesh-dependent at 2 mm. If J-M misses low, that is the leading
  suspect, and only Stage D tests it.
- (D2a) Flux is applied on horizontal cell area; no `1/cos(theta)`, no
  side-wall Robin.
- (D2a) **AMReX 25.12 parser bug:** numeric arguments first in every
  `max`/`min`; `unit/jet_enthalpy` (g) guards it.
- (S1b) `V0 = V_cell` switches the Weibull size effect off.
- (D2a2) `jet_mdot` is the modelled sector's share: quarter domain = mdot/4.
- (D2a2) **The planner's hand tables were wrong twice** (a non-computable
  criterion, then decay rows at D = 7.1 mm with `s_c = SOD`). Hand predictions
  here are computed by the implementer with `walljet.py` at D = 7.5 mm,
  **before** any run (deliverable 5). This packet gives no planner numbers.

## Sources

Line numbers are approximate. The tree carries uncommitted
P1 + R1 + E1 + S1 + A1 + C1 + D2a + D2a2 edits.

### Burner geometry (drawings; full record in `ROADMAP.md` Known Notes "Burner drawings")

- **Sheet 15/17, "Abstand Halter F VT5", section N-N (re-read 2026-09-16 at
  higher resolution):**
  - an 8 mm flange (OD Ø 80, bore Ø 45 +0.2/+0.05, Ø 6.2 bolt holes with Ø 10
    counterbores on a Ø 68 circle);
  - a **tube OD Ø 80 / ID Ø 56** running 58 - 8 = **50 mm** below it;
  - **three slots** cut into the tube at the bolt positions.
  - The tube's lower rim is the foot.
- **Foot annulus: r in [0.028, 0.040] m** (Ø 56-80). This corrects the first
  packet's [0.034, 0.040] "nominal", which misread the Ø 68 bolt circle as a
  wall. The slots remove part of the rim; their angular width is not
  dimensioned.
- **`foot_standoff = 0.050 m`** (sheet 1/17: nozzle exit to foot tips);
  nozzle bore **D = 7.5e-3 m**, `SOD/D = 6.67`; burner Ø at the drilling end
  **80 mm**.
- **Confinement:** the jet leaves the Ø 7.5 bore and travels its 50 mm inside
  the Ø 45 / Ø 56 collar; the exhaust leaves through the three slots and, once
  rock under the rim is removed, under it.
- **Model limitations to state, not fix:**
  - in reality rock under a steel pad is covered and gets no gas; the model
    heats every ring column with the full wall-jet flux;
  - a quarter domain with two symmetry planes is 4-fold symmetric, so three
    "pads" in it are a proxy for a tripod, not a tripod.

### Meier targets (`validation/meier/meier-2017-pilot-660kg-validation-reference.md`, `Claude_markdowns/2026-09-15b.md` §1, §5)

| target | criterion |
|---|---|
| ROP | 1.3-1.6 m/h +/-20 %, i.e. **[1.04, 1.92]** |
| hole diameter (depth-mean over the 0.5 m block) | **85-93 mm**, and >= 80 mm everywhere |
| volume rate | 2.47 cm3/s +/-20 %, i.e. **[1.98, 2.96]** |
| `dT_fire` | 500-560 K |
| ROP flatness | steady-window sub-fits within +/-10 % of the mean |
| ledger | closed to round-off |
| mesh | ring ROP and `T_fire` within 5 % between 2 and 1 mm, on a steady window |

- 3.42 L / 0.50 m is exactly Ø 93 mm: **the volume is the depth-mean
  diameter restated.** Since the simulation drills only part of the block, the
  0-0.5 m depth-mean comes from the hand depth profile (deliverable 5)
  **checked against** the simulated profile over the depth actually drilled.
- The consistent data triples are (2.47 cm3/s, Ø 93, 1.3 m/h) and
  (2.47, Ø 85, 1.57). Rock-side removal power from the data: **2.83 kW**.
- Wall temperature is **not** a target.

### Code this step touches (`src/Integrator/MMWSpalling.H`)

- `BuildFlameColumns(lev, time, k_top_col)` ~1499-1527: per column `r`,
  `z_face = PLO_z + (k_top_col+1)*dz`, and `s = z_n - z_face` with
  **`z_n = nozzle_z0 - nozzle_feed*time`** (the line the feet rule
  generalises). Then the D2a expression path or, under `enthalpy`,
  `JetEnthalpyMarch`.
- **`JetEnthalpyMarch`** ~1564-1721:
  - gathers top-cell T_old / k_c;
  - bins by r and takes `s_c` from the innermost non-empty bin;
  - `phi = min(1, core*D/s_c)`,
    `T_stag = T_ent + (T_nozzle - T_ent)*phi`, `m = mdot` or `mdot/phi`;
  - marches `T_gas <- max(T_ent, T_gas - P_bin/(m*cp))`;
  - checks the invariant (abort at 1e-9);
  - sets the 8 `jet_*` thermo values.
  - **Deliverable 2 changes this function, and only under the new mode.**
- `PinRule` ~1752 (shared with `BuildPinEffective` ~1731); do not change.
- `SurfacePatch::Parse` ~4561-4608 (D2a2 `jet_*` keys), plus the D2a keys
  `nozzle_z0` (required with expressions), `nozzle_feed` (default 0) and
  `nozzle_collision_radius`.
- Thermo registration ~602 (`jet_*` appended after every older column).
- **Pattern to mirror:** the Step 19 bit feed, `bit.*`, parse ~2725-2751,
  members ~4606-4619. **Do not change `bit.*`.**
- **Wall treatment** (selected, not edited):
  - `spall.flake_coherence_length` (parse ~2914; applied at `Removal.H:1062`
    only when `surface_normal` is off);
  - `spall.surface_normal` (parse ~2788; per-column `cos(theta)` at
    `Removal.H` ~169, floor 0.15). Under it a firing recedes the column by
    `h/cos(theta)` **with no matching increase** in the flux per horizontal
    area, so its energy ratio exceeds 1 on slopes by construction.
- `k_top_col` is the removed-mask top (`ColumnTopSolid`); per-column
  reductions must be domain-wide MPI reductions.
- `removal_events.csv`:
  `time,col_i,col_j,regime,k_top,T_top,a_f,Sp_top,h_scan,h_applied,n_voided`;
  `k_top` is the **pre-removal** top.

### Study helpers

- `tests/MMWSpalling/studies/d2a2_walljet/walljet.py` (D = 7.5e-3):
  `martin`, `h_anchor`, `h_expr`, `h_py`, `phi_of`, `t_stag`, `loss`,
  `rop_closed`, `march(h_ref, T_nozzle, cp, s, ..., mdot)`, `centre_rop`,
  `s_equilibrium(h_ref, T_nozzle, feed_mh, stag)`.
- `run.py` there builds the CLI overrides (`flame_keys`, quarter domain,
  `jet_mdot = MDOT/4`, dt from the stagnation flux).
- `d2a_jet_face/jet.py` keeps D = 7.1e-3 and `walljet.py` overrides it. Do not
  edit `jet.py`.

## Goal

### Code (default-off; microstructure path; single level)

1. **Feet descent rule.** New key
   `surface_patch.nozzle_descent = prescribed | feet`, default `prescribed`
   (D2a/D2a2 behaviour, bit-identical).

   Under `feet`, required:
   - `surface_patch.foot_r_inner`, `surface_patch.foot_r_outer` [m], the
     bearing annulus about `(x0, y0)`;
   - `surface_patch.foot_standoff` [m].

   Optional:
   - `surface_patch.foot_rule = pads | max | mean`, default `pads`;
   - `surface_patch.foot_npads` (default 3);
   - `surface_patch.foot_pad_quantile` (default 0.9);
   - `surface_patch.nozzle_feed_max` [m/s], default `+inf` (no operator cap).

   **`pads` rule**, evaluated at lev 0, once per step, in `BuildFlameColumns`
   **before** `s` is computed:
   - live columns with `r in [foot_r_inner, foot_r_outer]` are split into
     `foot_npads` equal azimuth sectors of the annulus's **observed** azimuth
     range (`atan2` about `(x0, y0)`, so a quarter domain splits its 90 deg);
   - each pad's height is the `foot_pad_quantile` quantile of `z_face` over its
     columns, using a deterministic rank rule (state it);
   - `z_foot` is the mean of the pad heights.
   - `max` and `mean` are the whole-annulus reductions, kept for unit tests
     and comparison.
   - All reductions are domain-wide and rank-independent (gather the annulus
     `z_face` values, then sort).

   Descent, for every rule:
   - `z_target = z_foot + foot_standoff`;
   - `z_n = min(z_n_prev, max(z_target, z_n_prev - nozzle_feed_max*dt))`; the
     burner never rises;
   - first step: `z_n = z_target`;
   - `z_n` is host state and **not checkpointed** (same class as the C1 pin
     state); say so in the takeaways.

   Aborts:
   - parse, under `feet`:
     - `nozzle_z0` or `nozzle_feed` given;
     - bad radii (`foot_r_inner >= foot_r_outer`, non-finite or negative);
     - `foot_standoff <= 0`;
     - `foot_npads < 1`;
     - quantile outside `(0, 1]`;
     - flame expressions not enabled;
   - runtime:
     - the annulus has no live column;
     - a pad has no live column.

   A **stall** is a result: log it, do not abort, do not add a fallback
   descent.

2. **Exhaust-recirculation entrainment.** New key
   `surface_patch.jet_T_ent_mode = fixed | exhaust`, default `fixed`
   (D2a2 behaviour, bit-identical).
   - Under `exhaust`, the temperature of the entrained gas is
     `T_rec = jet_T_exhaust` from the **previous** step (first step:
     `jet_T_ent`). The stagnation mixture is
     `T_stag = T_rec + (T_nozzle - T_rec)*phi` with `m = mdot/phi`.
   - **Requires `jet_entrained_mass = 1`** (abort otherwise). Only then is the
     mixture energy-consistent: the stagnation enthalpy is the nozzle flow
     plus `(m - mdot)` of recirculated gas.
   - **The march floor stays at the fixed `jet_T_ent`, never `T_rec`.**
     Flooring at the lagged exhaust temperature would make `T_exhaust`
     non-decreasing step to step (a ratchet), and the gas could never cool.
     Unit test (j) checks this.
   - Invariants, all referenced to the fixed `jet_T_ent`:
     - `P_face + P_exhaust = P_cap = m*cp*(T_stag - T_ent)`, as today;
     - `P_face + P_exhaust = mdot*cp*(T_nozzle - T_ent) + P_recirc`, with
       `P_recirc = (m - mdot)*cp*(T_rec - T_ent)`.
     - At steady state this implies
       `P_face = mdot*cp*(T_nozzle - T_exhaust)`: the net enthalpy leaving the
       hole is the nozzle flow at `T_exhaust`.
   - `T_rec` is host state, not checkpointed.
   - **Limitation to state:** near the surface (the first few cm) the real
     entrained gas is ambient air; the mode applies from t = 0.

3. **Collision radius under `feet`.** The default `nozzle_collision_radius`
   becomes `foot_r_inner`. Keep it overridable and document it.

4. **Thermo**, registered **after** every pre-D2b column:
   - under `feet`:
     - `nozzle_z` [m], `foot_z` [m];
     - `foot_cols` (live annulus columns);
     - `foot_carry_cols` (columns with `z_face >=` their pad's height, i.e.
       the ones carrying the burner);
     - `foot_stalled_steps` (cumulative steps with `z_target >= z_n_prev`);
   - under `jet_T_ent_mode = exhaust`:
     - `jet_T_rec` [K], `jet_P_recirc` [W].

### Hand prediction (before any run)

5. **Pre-registered feet prediction, computed with `walljet.py`** and written
   into `RESULTS.md` §0 before any D2b run, together with the code that
   produced it. Setup: exhaust recirculation (iterate `T_rec` to its fixed
   point), entrained mass, `cp = 1250`, D = 7.5 mm, full jet.
   - **Feet equilibrium:**
     - the ring (the annulus, reduced like the pad rule) sits at
       `s = foot_standoff`, and the centre at its own `s_c`;
     - find `(v, s_c)` such that the centre ROP at `s_c` equals `v` **and**
       the ring ROP from the march started at `T_stag(s_c, T_rec)` equals `v`;
     - state how `s` is interpolated between centre and ring.
   - **Depth profile:**
     - freeze depth `d(r) = foot_standoff * v / (v - v(r))` for `r` beyond
       the ring;
     - wall profile `r_wall(z)`;
     - **depth-mean hole Ø over 0-0.5 m**, and over the depth Stage B will
       actually drill.
   - Report `v`, `s_c`, `T_stag`, `T_rec`, `r_wall(z)` and the depth-mean Ø
     for J-M, J-5 and J-10 at 1900 K, plus J-M with fixed 293 K entrainment
     and J-M at 1436 K.
   - Feed the model D2a2's measured `s_c` and say whether it reproduces D2a2's
     early-window ring ROPs.

### Tests

6. **New `tests/MMWSpalling/unit/robin_feet/`** (`input_*` + `test`, the
   `robin_jet` / `jet_enthalpy` pattern, 1 rank unless noted). **Re-run the
   test after any edit to it and keep the passing log** in
   `output/test_pass.log`.
   - **(a) Defaults are bit-identical.** `prescribed` + `fixed` reproduce a
     `jet_enthalpy` input's Cell_D, `removal_events.csv` and every thermo
     column exactly.
   - **(b) Kinematics.** On every row, `nozzle_z - foot_z == foot_standoff`
     (print precision) whenever the feet are down, and `nozzle_z` is
     non-increasing.
   - **(c) Pad rule.** Held column heights give `foot_z` equal to the hand
     value of the pad quantiles' mean, and `foot_carry_cols` equals the hand
     count. Check `max` and `mean` the same way.
   - **(d) Operator cap.** With `nozzle_feed_max` below the ring recession
     rate, `-d(nozzle_z)/dt` equals the cap and the feet lift off.
   - **(e) Stall.** One ring column that never fires:
     - under `max` it freezes `nozzle_z` and increments
       `foot_stalled_steps`, with no abort;
     - under `pads` with quantile 0.9 the burner keeps descending.
   - **(f) Aborts.** Each parse abort (including exhaust mode without
     entrained mass) and each runtime abort returns `rc != 0` with its
     message.
   - **(g) Rank identity.** A `feet` + `exhaust` run on 3 ranks gives
     identical `removal_events.csv`, identical new thermo text, and identical
     Level_0 FABs compared **through `Cell_H` FabOnDisk offsets**.
   - **(h) Closed form.** A single-ring-column feet run reaches a steady ROP
     matching `rop_closed` at `s = foot_standoff` and the march's `T_gas`,
     within 2 %.
   - **(i) Exhaust identities.** On every row, `jet_T_rec` equals the previous
     row's `jet_T_exhaust` (at thermo interval 1) and the `P_recirc` identity
     holds to round-off. At steady state,
     `P_face = mdot*cp*(T_nozzle - T_exhaust)` within 1 %.
   - **(j) No ratchet.** A run in which the face cools the gas more over time
     (e.g. a growing unfired area) shows `jet_T_exhaust` **decreasing** on some
     rows.

### Study and scored test (3-D only)

7. **Harness `tests/MMWSpalling/studies/d2b_feet_rop/`**: `run.py`,
   `analyze.py`, `RESULTS.md`, `README.md`. **No 2-D slab part.**

   **Common settings** (`input_drilling` + CLI overrides):
   - **Domain:** quarter domain with the axis at the xlo/ylo corner
     (Neumann-0 = symmetry); 2 mm; 0.12 x 0.12 x **0.30 m** (60x60x150);
     initial top z = 0.30.
   - **Off:** `beam.P0 = 0`, `bit.enabled = 0`.
   - **Removal and ledger:** `weibull.V0 = V_cell`; `robin_form = pinned`,
     `pinned_idle_cycles = 2`; `surface.follow_mask = 1`; energy ledger and
     `removal_events_csv` on.
   - **Jet:** `jet_closure = enthalpy`, `jet_stagnation = decay`,
     `jet_entrained_mass = 1`, **`jet_T_ent_mode = exhaust`**,
     `jet_T_ent = 293.15`, `jet_mdot = mdot/4`, `jet_cp = 1250`,
     `jet_D = 7.5e-3`, `jet_T_nozzle = 1900`, `jet_profile_interval` on.
   - **Feet:** `nozzle_descent = feet`, annulus **[0.028, 0.040]**,
     `foot_standoff = 0.050`, `foot_rule = pads`, `foot_npads = 3`,
     `foot_pad_quantile = 0.9`.
   - **Time and output:** dt by the overshoot rule (4 ms cap);
     **`amr.plot_int` sparse** (every 50-100 s).
   - **Stall watchdog (in `run.py`):** stop a run cleanly and mark it
     `stalled` if `nozzle_z` has not descended for 60 s of simulated time.
     Keep the data up to the stop.

   **Stage A: wall-treatment decision** (J-M, ~150 s each), comparing:
   - **A1** `flake_coherence_length = 0.004` (as in `input_drilling`);
   - **A2** `flake_coherence_length = 0`;
   - **A3** `spall.surface_normal = 1`.

   **Decision rule, written into `RESULTS.md` before Stage A runs:**
   1. **Energy ratio over the whole hole:**
      `R = rho*Cp*dT_fire*(removed volume rate) / (absorbed face power)`
      must be <= 1.05. Also report `R` per radial bin; a rim bin can
      legitimately exceed 1 through lateral conduction from the hotter bowl.
   2. **Overheating:** removal `T_top` p99 no more than 50 K above the pinned
      `T_fire` band.
   3. **Artefacts:** no needles (single-column spikes), no never-fired columns.
   4. **If a treatment passes 1-3,** choose among passers in the order
      A2 > A3 > A1.
   5. **Tie-break, if none passes 1-3:** ring recession is the scored
      quantity. Accept the treatments whose violations of 1-3 do **not**
      touch the foot-annulus bins (ring-bin `R <= 1.05`, no overheated or
      never-fired ring columns, no needles in the ring), choose among them in
      the order A2 > A3 > A1, and carry the violation as a caveat into the
      verdict.
   6. **Stop and report only if every treatment violates inside the ring.**
   - **Do not look at ROP proximity to Meier.**

   **Then run the stages in the order B, D, C.**

   **Stage B: scored runs** (chosen wall treatment; stop at 600 s, at the
   watchdog, or when the centre is 40 mm above the domain bottom):
   - **J-M** (scored), **J-5**, **J-10**.
   - Run the cost probe first. D2a2 measured 20-38 min per 300 s at depth
     0.20 m. If a run would exceed ~2 h on 4 ranks, shorten it but keep a
     steady window of >= 150 s, and say so.
   - **Steady** means `nozzle_z` descent-rate sub-fits within +/-10 % and
     `jet_s_c` drift under 0.02 mm/s.

   **Stage D: mesh check (required)**, J-M at 1 mm.
   - **Trimmed domain:** the quarter domain cut to **r <= 0.06 m**, with depth
     sized to the steady window. The march is outward, so rock beyond the
     60 mm edge cannot change the gas that reaches the 40 mm ring; the trim
     changes only the outer hole shape and the lateral boundary.
   - Report the trimmed 2 mm counterpart as well, so the mesh comparison is
     like for like.
   - Stop once >= 60 s of steady window exists.
   - **Budget:** probe first. If the probe exceeds ~4 h on 4 ranks, shorten to
     the minimum steady window. If it still exceeds 4 h, report the probe and
     mark the mesh target "not evaluated", with the numbers. **Do not
     substitute a 2-D check.**
   - **Score:** ring ROP and `T_fire`; also report the pinned share and
     face-form uptake between 2 and 1 mm (the D2a2 review's open question).

   **Stage C: sensitivities** (J-M, three runs):
   - `jet_T_ent_mode = fixed` (293 K ambient entrainment; D2a2's treatment);
   - `jet_entrained_mass = 0` with `jet_T_ent_mode = fixed` (exhaust mode
     requires entrained mass);
   - `jet_T_nozzle = 1436 K`.

8. **`RESULTS.md`.**
   - **§0, written before the runs it governs:** the hand prediction
     (deliverable 5) and the Stage A decision rule.
   - **Per run:**
     - **ROP:** least-squares fit of `-d(nozzle_z)/dt` over the steady
       window, cross-checked against the ring columns' `h_applied` fits.
     - **Flatness:** fits over 50 s sub-windows.
     - **Stall:** `foot_carry_cols` over time, `foot_stalled_steps`, and the
       watchdog status.
     - **Hole wall profile** `r_wall(z) = max{r : that column has receded
       below z}` over the drilled depth, compared with the hand `r_wall(z)`.
       Also the depth-mean Ø over the drilled depth, the hand-extrapolated
       depth-mean Ø over 0-0.5 m, and the minimum `2*r_wall` vs the 80 mm
       burner.
     - **Volume rate:** `d/dt` of `integral(pi*r_wall(z)^2 dz)` over the
       steady window (quarter-domain `r_wall` is already a radius), and the
       0-0.5 m volume from the extrapolated profile.
     - **Temperatures and jet:** `T_fire` mean/p10/p90 and `dT_fire`;
       `jet_s_c`, `jet_T_stag`, `jet_T_rec`, `nozzle_z - foot_z`.
     - **Power:** face power `jet_P_face` x4 against `jet_P_cap` x4, against
       `mdot*cp*(T_nozzle - T_exhaust)` and against Meier's 2.83 kW (in kW,
       not % of 38 kW).
     - **Checks:** `R` over the hole and per radial bin, branch-reason
       fractions, rim check, ledger, maximum face slope, and in-hole area
       outside Martin's validity range.
   - **Scored table for J-M**, pass or fail as it comes out, next to the hand
     prediction. Then the same table for J-5 and J-10 as the bracket's
     implication.
   - **Discussion:**
     - (i) The verdict, stated plainly. **A hole-diameter or volume pass is
       read as a check of the profile shape and the above-nozzle rule, not as
       support for the anchor. The ROP line is the test of the jet.**
     - (ii) Hand prediction vs simulation, and what the difference is made
       of.
     - (iii) Sensitivities (%/unit) and which could move the verdict. The
       fixed-293 K run is shown **in the verdict table**, not only here.
     - (iv) Stalls, interpreted as in Context point 4.
     - (v) What D3 must say:
       - the slope area factor;
       - side-wall and exhaust heating;
       - pad-covered rock heated in the model;
       - the quarter-domain "tripod";
       - face-form mesh sensitivity and the Stage D result;
       - entrainment near the surface;
       - the out-of-range correlation area;
       - the size effect switched off.
     - (vi) The sentence **"no parameter was adjusted toward Meier's ROP,
       hole diameter or volume rate"**, plus every number taken from outside
       the repo with its source.

9. **Retarget the stale Meier scored test.**
   - Add `validation/meier/sp_meier_pilot/input_feet` (the Stage B J-M
     configuration, shortened if needed) and a new scored `test_feet`.
   - `test_feet` checks the J-M table's criteria **as they came out**. A
     failing criterion is recorded as an expected-fail with its numbers, not
     loosened.
   - Leave `input_drilling` and `test` untouched, and add a note at the top of
     `validation/meier/README.md` that the uniform-beam test is superseded.
   - If the configuration is too expensive for a regression, make
     `test_feet` a check-only script over the Stage B output, and say so.

## Guardrails

- **No tuning.** J-M is scored whatever it gives. The following are all set
  as stated **before** any run:
  - the anchor grid, `h_ref`, `T_nozzle`, `cp`;
  - the closure treatment (decay + entrained mass + exhaust recirculation);
  - the foot geometry, pad rule and quantile;
  - `pinned_idle_cycles`.

  The wall treatment is chosen **only** by the Stage A rule.
- **Bit identity**, with `nozzle_descent = prescribed` and
  `jet_T_ent_mode = fixed` (the defaults). Use targeted witnesses:
  - unit: `robin_face`, `robin_pinned`, `robin_jet`, `jet_enthalpy`,
    `beam_void_closure`, `scalar_flaw`, `spall_event`;
  - Meier `input_2d_dev`;
  - S1 `run_sweep.py --smoke --force`.

  Hash Level_0 Cell_D and CSVs pre/post for these. New modes only add thermo
  columns. **The full 26k-file sweep is deferred to the commit** that
  includes D2b.
- **Do not change:**
  - `SurfaceCellFlux`, `PinRule`, `BuildPinEffective`, `EnergyLedger`,
    `PatchDiagnostics`;
  - `JetEnthalpyMarch` other than the `exhaust` branch of deliverable 2;
  - the removal cadence, phi shift and void flip, K_I scan;
  - the coherence cap and surface-normal code (they are selected, not
    edited);
  - E1, S1, A1, C1, D2a, D2a2 keys and defaults;
  - `bit.*`;
  - the `weibull.V0` default;
  - any existing input or assertion (deliverable 9 adds files only).
- **Top-cell logic** uses the `removed` mask, never `floor(phi/dz)`.
- **Column-top patch only.** Flux on `dx*dy`; no `1/cos(theta)`; no flux for
  `s <= 0`.
- **Out of scope:**
  - a slope area factor;
  - side-wall or exhaust heating of the wall;
  - pad shielding of the rock under the feet;
  - a depth-dependent blend of ambient and exhaust entrainment;
  - an operator push in the baseline;
  - AMR, `alpha(T)`, MMW;
  - D3;
  - an overdraw counter (deferred).
- **Do Not Touch** (CLAUDE.md): `ext/`, `bin/`, `obj/`, `build/`,
  `compile_commands.json`, `configure`, `LICENSE`, unrelated integrators,
  the inheritance graph.
- **Git:** build on the uncommitted tree (do not revert or commit). If the
  user has committed, use that commit as the baseline.
- **Re-read `ACTIVE_STEP.md` (check its mtime) before the long Stage B/D runs
  and before writing notes** (D2a2 takeaway 1).

## Commands

```bash
# Build from /Users/tzetze20/amr_tools/alamo. Save a pre-change binary first.
cp bin/mmwspalling-3d-g++ <scratchpad>/mmwspalling-preD2b
EIGEN=$PWD/ext \
  CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current \
  make -j8

VENV=/Users/tzetze20/Desktop/code/.venv/bin/python
RUN="mpirun --oversubscribe --bind-to none -np 4"

# New unit test
$VENV tests/MMWSpalling/unit/robin_feet/test

# Targeted bit-identity witnesses (hash Level_0 Cell_D_* and CSVs pre/post)
for t in unit/robin_face unit/robin_pinned unit/robin_jet unit/jet_enthalpy \
         unit/beam_void_closure unit/scalar_flaw unit/spall_event; do
  $VENV tests/MMWSpalling/$t/test
done
$RUN bin/mmwspalling-3d-g++ tests/MMWSpalling/validation/meier/sp_meier_pilot/input_2d_dev
$VENV tests/MMWSpalling/studies/s1_surface_resolution/run_sweep.py --smoke --force

# Study (order: predict, A, probe, B, D, C)
$VENV tests/MMWSpalling/studies/d2b_feet_rop/run.py --predict     # deliverable 5, RESULTS §0
$VENV tests/MMWSpalling/studies/d2b_feet_rop/run.py --stage A
$VENV tests/MMWSpalling/studies/d2b_feet_rop/run.py --probe
$VENV tests/MMWSpalling/studies/d2b_feet_rop/run.py --stage B
$VENV tests/MMWSpalling/studies/d2b_feet_rop/run.py --probe-1mm
$VENV tests/MMWSpalling/studies/d2b_feet_rop/run.py --stage D
$VENV tests/MMWSpalling/studies/d2b_feet_rop/run.py --stage C
$VENV tests/MMWSpalling/studies/d2b_feet_rop/analyze.py           # writes RESULTS.md

# Scored Meier test (deliverable 9)
$VENV tests/MMWSpalling/validation/meier/sp_meier_pilot/test_feet
```

## Claude completion notes

Implementer, 2026-09-16. The packet was read at 16:49 and unchanged at
close-out (mtime 16:49:55). Stage A was decided at 18:20; Stages B, D and C
ran 18:23–21:27; the final regressions ran 21:28–21:38.

### Files changed

- `src/Integrator/MMWSpalling.H`: the only source file.
  - `nozzle_descent = feet` (`FeetDescent`) and its keys.
  - `jet_T_ent_mode = fixed | exhaust` (in `JetEnthalpyMarch`).
  - The collision-radius default under `feet`.
  - New thermo columns: `nozzle_z`, `foot_z`, `foot_cols`,
    `foot_carry_cols`, `foot_stall_time`, `jet_T_rec`, `jet_P_recirc`.
  - `BuildFlameColumns` now takes `dt`.
- `tests/MMWSpalling/unit/robin_feet/` (new): `input_feet`, `input_column`,
  `test`; the passing log is `output/test_pass.log`.
- `tests/MMWSpalling/studies/d2b_feet_rop/` (new): `feetmodel.py`, `run.py`,
  `analyze.py`, `README.md`, `RESULTS.md`, `stageA.png`, `stageB.png`,
  `stageC.png`, `stageD.png`, and `output/` (6.4 GB; runs A_*, B_*, C_*,
  D_*, probe, probe_1mm).
- `tests/MMWSpalling/validation/meier/sp_meier_pilot/input_feet` (new):
  generated by `run.py --write-input-feet --wall A2 --stop-feet 340`.
- `tests/MMWSpalling/validation/meier/sp_meier_pilot/test_feet` (new):
  check-only scoring of the study's B_JM_A2 output.
- `tests/MMWSpalling/validation/meier/README.md`: a "superseded
  uniform-beam test" note at the top. `input_drilling` and `test` are
  untouched.
- `Claude_markdowns/2026-09-16_imp.md`: the status note the user asked for.

### Tests run

All on 4 ranks, against the final binary (built 21:28:19).

| test | result |
|---|---|
| `unit/robin_feet/test` | **PASS 20/20**, checks (a)–(j) plus the aborts |
| `unit/robin_face`, `robin_pinned`, `robin_jet`, `jet_enthalpy`, `beam_void_closure`, `scalar_flaw`, `spall_event` | **PASS** (all rc 0) |
| Meier `input_2d_dev` run and the S1 smoke (`run_sweep.py --smoke`) | rc 0 |
| **Hash witness** over the targeted outputs (Level_0 Cell_D, events, clusters, jet profile, thermo) | **2501 / 2501 files identical** to the pre-D2b baseline |
| `validation/meier/sp_meier_pilot/test_feet` | **PASS**: every outcome as recorded, 4 of them expected-fail (below) |
| `input_feet` smoke (0.2 s) | runs; its t = 0.1 s thermo row equals B_JM_A2's |
| Study | Stage A (3 runs), cost probes (2 mm, 1 mm), Stage B (3), Stage D (2), Stage C (3); all completed, statuses in the `.done` files |

**Scored J-M result (B_JM_A2, stopped at the bottom watchdog at 341.6 s;
window 171–342 s):**

| criterion | result | verdict |
|---|---|---|
| ROP | 1.36 m/h | PASS |
| hole Ø, depth-mean over 0–0.5 m | 104.3 mm | FAIL |
| hole Ø, minimum | 84.9 mm | PASS |
| volume rate inside Ø 93 | 3.27 cm³/s | FAIL |
| ΔT_fire | 528 K | PASS |
| flatness | sub-fits 2.5 %, but `jet_s_c` drifts 0.39 mm/s | FAIL |
| ledger | 1.8e-14 | PASS |
| mesh | ring −4.8 %, T_fire 0.00 %, but no steady window | FAIL |

J-5 and J-10 are above the band (2.43 and 2.82 m/h). J-M at 1436 K is below
it (0.67 m/h). Both fixed-293 K runs stalled. The full tables and the
discussion are in `studies/d2b_feet_rop/RESULTS.md`.

### Tests not run

- **Full MMWSpalling regression sweep:** deferred to the commit, per the
  packet's targeted-witness rule.
- **`test_feet` is check-only**, not a regression run. The 340 s 3-D
  configuration takes about 40 min with three jobs in parallel (over an hour
  alone).
- **The study runs used the binary built at 16:54**, before `foot_stall_time`
  replaced `foot_stalled_steps` (see takeaways). The change is diagnostic
  only: the witness hashes and `robin_feet` (a) and (g) are unchanged. The
  study tables compute the longest hold from `nozzle_z` instead.

## Implementation takeaways

- **Input keys** (all `surface_patch.`):
  - `nozzle_descent = prescribed | feet`, `foot_r_inner`, `foot_r_outer`,
    `foot_standoff`, `foot_rule = pads | max | mean` (default pads),
    `foot_npads` (default 3), `foot_pad_quantile` (default 0.9),
    `nozzle_feed_max` (default +inf);
  - `jet_T_ent_mode = fixed | exhaust` (exhaust requires
    `jet_entrained_mass = 1`).
  - Under `feet`, `nozzle_z0` and `nozzle_feed` abort, and
    `nozzle_collision_radius` defaults to `foot_r_inner`.
- **Not checkpointed:** `nozzle_z`, the feet init flag and T_rec. A restart
  re-lands the burner on the feet and restarts T_rec from T_ent.
- **Stall counter replaced (user direction during review).** The packet's
  `foot_stalled_steps` increments on almost every step at 2 mm, because pad
  heights move in whole cells (18,710 at 75 s; 37,399 of about 37,500 steps
  in Stage A).
  - It is replaced by `foot_stall_time`: simulated time since `nozzle_z`
    last decreased, 0 on a descending step.
  - `robin_feet` (e) checks it equals (steps − 1)·dt for a held burner and
    resets on descent.
- **The packet's freeze depth is off by the stand-off.** d(r) = 50·v/(v −
  v(r)) is the feet depth when column r freezes. The wall sits at the
  column's own depth, 50·v(r)/(v − v(r)). The hand model and analysis use
  the latter.
- **The centre pit runs away because of the 12 D stand-off clip in
  `h_expr`,** not because of exhaust recirculation alone.
  - Clipped, the centre keeps 1.86× the ring's h beyond 90 mm. The hand
    model with the unclipped Martin shape finds a J-M centre equilibrium at
    s_c ≈ 0.40 m (1.93 m/h).
  - Every exhaust run hit the bottom watchdog: J-M at 342 s, J-5 at 150 s,
    J-10 at 103 s, D at about 290 s. The 1436 K run was the only one to
    reach 600 s.
  - As a result **no run has a steady window** (`jet_s_c` drift 0.17–1.4
    mm/s against the 0.02 limit).
  - The J-M ROP is still rising slowly (1.20 → 1.40 m/h) as face power
    falls (12.4 → 6.4 kW ×4) and T_rec rises (1247 → 1561 K), toward the
    hand model's deep-pit limit of 1.82 m/h at 1679 K.
  - A planner who wants a steady scored window needs a taller domain or an
    h shape that decays beyond 12 D. Either is a model decision, not
    tuning.
- **Scoring decisions**, recorded in the RESULTS §0 addendum at 18:38,
  before any Stage B result, following review feedback:
  - volume rate and rock-side power are taken inside Ø 93 (Meier's numbers
    are ROP × nominal area; whole-top removal is 5.68 against 3.27 inside);
  - hole Ø is scored over the 0.5 m block, as the packet's target table
    says: columns above the nozzle keep their depth, active columns are
    projected with the run's own v(r) until the nozzle reaches them;
  - the depth-mean over the drilled depth alone (≈ 0.12 m, mostly mouth
    crater) is reported but not scored.
- **Stage A.**
  - The needle definition was kept (deeper than *every* 4-neighbour by
    more than 2 dz). No run has any, whether depth comes from k_top or
    summed h_applied.
  - A review count of 36 matches "deeper than some neighbour": wall steps
    steeper than 63°, which A2 has (74 at 150 s, 6 in the ring). Recorded as
    a caveat.
  - All three runs failed criterion 3 on never-fired columns outside the
    ring (71–89). A1 also had whole-hole T_top p99 887 K (ring maximum
    862 K). Rule 5 chose A2.
  - **A3 (`spall.surface_normal = 1`) is inert under the pinned closure**:
    the cos θ normal thickness is divided back out by the applied
    recession, so A3 = A2.
- **The fixed-293 K stall came out as predicted.** Both mass treatments
  stalled (at 162 s and 118 s) after 12 and 6 mm of burner descent. Read as
  "the jet alone does not clear the ring here".
- **Mesh (Stage D, trimmed to r ≤ 60 mm).**
  - Ring −4.8 %, burner −6.4 %, T_fire 0.00 % from 2 to 1 mm.
  - Pinned share 0.82 at 2 mm and 0.72 at 1 mm; unset face 0.02 and 0.01.
  - The trimmed domain recedes all the way to its 60 mm edge, so its hole
    shape is not meaningful.
- **Cost.** A 600 s run alone takes about 1.2 h (J-10 at 2 ms: about
  2.4 h). In practice the bottom watchdog ends runs at 25–42 min with three
  jobs in parallel. The 1 mm trimmed run took about 80 min to 295 s.
- **Hand model vs simulation.**
  - The hand model is a deep-pit limit. J-M is 25 % below it (1.36 vs
    1.82), 1436 K is 31 % below, and J-5 / J-10 are 44–52 % below (short
    runs, far from the limit).
  - A review check at 75 s: with the run's own ring gas temperature, the
    closed form gives 1.30 m/h against a 1.31 m/h ring recession. The
    closure is right; the gap is in the gas that reaches the ring.
- **`analyze.py` notes.**
  - `thermo.dat` is written every 0.1 s, not every step.
  - The t = 0 row has zeros in the D2b columns (it is written before the
    feet are placed) and is excluded from the plots.
  - `mesh_result()` is shared with `test_feet`.
- **Relinking `bin/` while study runs use it is unsafe on macOS.** The
  rebuild waited until all runs had finished; a syntax-only compile checked
  the source before that.

## Review findings

Verdict: accepted

Reviewer, 2026-09-16. No simulation or build was rerun. One analysis-only
script was run, for the reason given under Tests.

### What was checked

- **Scope (deliverables 1-9).** All are present:
  - `FeetDescent`, the exhaust mode, the collision default and the thermo
    columns;
  - the §0 hand prediction and Stage A rule;
  - `unit/robin_feet` with checks (a)-(j) and the aborts;
  - the study harness with Stages A, B, D and C, run in that order
    (`stageA.log` / `stageBDC.log`: B 18:23, D 19:05, C 20:23, done
    21:27);
  - `input_feet`, `test_feet` and the README note.
  - `input_drilling` and `test` are unchanged (not in `git status`).
  - One deviation, `foot_stalled_steps` → `foot_stall_time`, was
    user-directed and is documented.
- **Source diff.** I diffed against the implementer's saved
  `scratchpad/d2b/MMWSpalling.H.preD2b` (262 changed lines). The diff
  touches only:
  - the thermo registration;
  - the `BuildFlameColumns` signature and its `z_n` line;
  - the new `FeetDescent`;
  - the `exhaust` branch of `JetEnthalpyMarch` (T_mix, P_recirc, P_budget
    and the tolerance);
  - members and parse.

  `SurfaceCellFlux`, `PinRule`, `BuildPinEffective`, `EnergyLedger`,
  `PatchDiagnostics`, `bit.*` and `Removal.H` (mtime 09-15) are untouched.
  Under `fixed`, T_mix = T_ent and P_budget = P_nozzle, so the default path
  is the same arithmetic.
- **Correctness.**
  - `FeetDescent` (`MMWSpalling.H:1569`) matches the packet:
    - inclusive annulus;
    - sectors over the observed atan2 range;
    - nearest-rank quantile ceil(qn)−1;
    - mean over pads;
    - never-rise min/max descent with the feed cap;
    - runtime aborts;
    - it runs before `s` is computed (`:1494`).
  - The exhaust budget checks out by algebra: with m = mdot/φ,
    m·cp·(T_stag − T_ent) = mdot·cp·(T_n − T_ent) + (m − mdot)·cp·(T_mix − T_ent),
    so `jet_P_decay` = 0 holds.
  - The march floor stays T_ent, and T_rec is the previous step's T_gas
    (`:1780`, `:1853`).
  - Parse aborts are at `:4775` and `:4804-4860`; the collision default is
    at `:4875`.
- **Tests ran against the final code.**
  - `robin_feet/output/test_pass.log` (21:38): 20 PASS, 0 FAIL. Its outputs
    (21:29-21:32) are newer than the test edit (18:35) and the binary
    (21:28:19). The run metadata says it was compiled at 21:27:57.
  - The witness hashes `scratchpad/d2b/pre.md5` vs `post.md5` (17:00) and
    vs `post2.md5` (21:38) show 2501/2501 identical. They cover all 7
    listed units, Meier `dev2d` and the S1 smoke. `final_checks.log` shows
    every witness at rc 0.
- **Study binary.** The study used the pre-rename binary. I diffed the
  `diff.patch` embedded in `B_JM_A2` (identical in A_A1, D_1mm and C_1436)
  against the one in `robin_feet/output/F`. In `MMWSpalling.H` they differ
  only in the stall-counter lines, so the study numbers stand for the final
  code.
- **Numbers.** I recomputed from `B_JM_A2/thermo.dat` over 171-342 s:
  - burner ROP 1.363 m/h;
  - T_rec 1484 K;
  - face power ×4 7.87 kW;
  - ledger 1.8e-14;
  - gap to foot_standoff 4e-17, with no lift-off.

  D_1mm and D_2mm (1.334 / 1.426 m/h, so −6.4 %) and C_1436 (0.666 m/h)
  also match `RESULTS.md`.
- **Tests: one analysis-only run.** I ran `test_feet` because no passing
  log of it exists. It reads existing output only and ran no simulation.
  Result: PASS, rc 0, with every outcome as recorded (4 expected-fails).
- **Guardrails and caveats.**
  - No tuning: the anchor, mode, feet and quantile are as stated, and Stage
    A was decided by rule 5.
  - pinned_idle_cycles is 2.
  - The mask top is used for z_face.
  - There is no flux at s ≤ 0.
  - There are no edits to `ext/`, `bin/` and so on.
  - The non-checkpointed state is documented.
- **Takeaways.** They are thorough and usable. The freeze-depth correction,
  the 12 D clip cause, the fact that A3 is inert under pinned, the stall
  counter and the cost figures are all recorded.

### Findings for the planner (non-blocking)

- **The trimmed mesh check is not like for like under exhaust mode.** The
  packet assumed rock beyond 60 mm cannot change the ring gas. Under
  `exhaust` it does: T_rec is the end-of-march T_gas, fed back into T_stag.
  - The trimmed 2 mm run has T_rec 1579 K; the full-domain J-M run has
    1484 K.
  - The 2 vs 1 mm comparison is like for like (both trimmed), but "1 mm J-M
    ≈ 1.27 m/h" is an extrapolation to the full domain.
- **The minimum-Ø pass is a 0.5 m-block artefact.**
  - Over the drilled depth the minimum is 78.2 mm, below the 80 mm burner.
  - The hand model narrows toward 2·r_q = 78 mm.
  - The model has no burner-body collision; the collision radius is only
    28 mm. A deeper hole would fail ≥ 80 mm. D3 should say so.
- **Scoring changed mid-study.** The volume rate and face power were moved
  inside Ø 93 by the 18:38 addendum.
  - It was written after Stage A (which showed 2.99 inside Ø 93, near the
    band), while Stage B was already running.
  - It is disclosed and review-directed, and the packet's own
    wall-profile rate (5.67) fails too, so the verdict does not change.
  - Hole-Ø block extrapolation uses the run's own v(r) instead of the hand
    profile. Both fail (104.3 vs 112.1).
- **The ROP pass is a transient-window score.**
  - `jet_s_c` drifts 0.39 mm/s, and the rate rises from 1.20 to 1.40 m/h
    toward the 1.82 deep-pit limit.
  - A converged J-M could leave the band (1.92 is the band top).
  - The next packet should resolve the 12 D h-clip or domain height before
    treating 1.36 m/h as the model's prediction.
- **The mesh line is soft.** T_fire at +0.00 % is fixed by construction
  under the pinned closure, so the only real mesh number is ring −4.8 %
  (burner −6.4 %).
- **Stage A picked A2 by the tie-break, not a clean pass.** It carries:
  - 71 never-fired columns in the hole;
  - 74 wall steps steeper than 63° (6 in the ring).
- **Martin validity.** 88 % of J-M in-hole columns are outside the
  correlation's validity range at t_end, so the h shape there is
  extrapolated.

