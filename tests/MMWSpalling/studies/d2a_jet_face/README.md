# D2a: impinging-jet face source study

**Purpose.** D2a makes the C1 pinned Robin closure spatial. The flame
coefficient and gas temperature become input expressions of the radius r from
the jet axis and the column stand-off s = z_nozzle(t) − z_face(col):
`surface_patch.h_expr`, `surface_patch.T_flame_expr`, plus `nozzle_z0`,
`nozzle_feed` and `nozzle_collision_radius`. The nozzle follows a prescribed
path. The study measures the emergent face shape, hole radius and stand-off in
2-D (machinery) and in a 3-D quarter domain at three feeds.

This is a study, not a regression, so there is no `test`. The regression is
`tests/MMWSpalling/unit/robin_jet`. Results are in [RESULTS.md](RESULTS.md).
**No parameter was adjusted toward Meier's ROP or hole diameter.**

## How to run

```bash
VENV=/Users/tzetze20/Desktop/code/.venv/bin/python
$VENV jet.py                       # correlation arithmetic + h_loc verification
$VENV run.py --parser-check        # compiled expressions vs numpy (10 one-step 3-D runs)
$VENV run.py --part 1 --jobs 3     # 2-D slab machinery, 7 runs, 4 ranks each -> output/p1_*
$VENV run.py --probe               # 3-D cost probe (10 s) with wall-time estimates
$VENV run.py --part 2 --jobs 3     # 3-D quarter domain, 5 runs, 4 ranks each -> output/p2_*
$VENV analyze.py                   # RESULTS.md tables above <!-- DISCUSSION --> + PNGs
```

Everything is a CLI override of `validation/meier/sp_meier_pilot/input_2d_dev`
(Part 1) or `input_drilling` (Part 2); no input is copied or edited.

## Correlations

### Martin (1977), single round nozzle

H. Martin, *Heat and Mass Transfer between Impinging Gas Jets and Solid
Surfaces*, Adv. Heat Transfer 13 (1977) 1–60. The formula is area-averaged
over a disk of radius r:

- `Nu_avg = h_avg D / k = Pr^0.42 · G(r/D, H/D) · F(Re)`
- `G = (D/r)(1 − 1.1 D/r) / (1 + 0.1 (H/D − 6) D/r)`
- `F = 2 Re^0.5 (1 + 0.005 Re^0.55)^0.5`
- valid for 2000 ≤ Re ≤ 4×10⁵, 2.5 ≤ r/D ≤ 7.5, 2 ≤ H/D ≤ 12.

**Verification of the planner's constants** (they were quoted from memory):

- **G and the ranges** match Zuckerman & Lior, *Jet impingement heat
  transfer: physics, correlations, and numerical modeling*, Adv. Heat Transfer
  39 (2006), table "Single Round Nozzle, Source: Martin [2]". That table
  additionally lists 0.004 ≤ f ≤ 0.04, the relative nozzle area, which applies
  to arrays.
- **F: Zuckerman & Lior quote a piecewise power-law form**, e.g.
  F = 0.54 Re^0.667 for 3×10⁴ < Re < 1.2×10⁵. At Re = 5.98×10⁴ that is 4 %
  below the continuous form.
- **The continuous F** is the form the packet gives, and a secondary summary
  of Martin's correlation (web search) quotes it identically. It is used here.
  Incropera §7.7 was not checked directly.
- The 4 % difference is well inside the Martin-versus-bracket conflict studied
  here.

### Local coefficient

`h_loc(r) = (1/2r) d(r² h_avg)/dr`, with c = 0.1 (H/D − 6):

`h_loc = Pr^0.42 F (k/D) · (D/2r)(r² + 2cDr − 1.1cD²)/(r + cD)²`

This matches the planner's derivation. `jet.verify_hloc()` rebuilds G from the
Martin disk average at 2.5 D plus the integral of 2 r h_loc over 2.5 D → r. It
reproduces G to 1.0×10⁻¹² for 2.5–7.5 D at H/D = 2, 7 and 12.

### Properties and arithmetic (`jet.py`)

- Air is used for the lean CH₄/air combustion products (~73 % N₂ by mass);
  this is the stated approximation.
- Film temperature: T_film = (T_ref + T_fire)/2 = (1500 + 821)/2 = 1160.5 K.
- Sutherland's law (White, *Viscous Fluid Flow*):

| quantity | law / value at 1160.5 K |
|---|---|
| μ | 1.716e-5 (T/273.15)^1.5 (383.55)/(T + 110.4) = 4.535e-5 Pa·s |
| k | 0.0241 (T/273.15)^1.5 (467.15)/(T + 194) = 0.0728 W/m·K |
| c_p | 1175 J/kg·K |
| Pr | c_p μ / k = 0.732 |

- ṁ = (52 + 2.47)/3600 = 0.01513 kg/s and D = 7.1 mm, so
  Re = 4ṁ/(πDμ) = 5.98×10⁴, inside the valid range.
- F = 864.1, G(2.5 D, 7 D) = 0.2154, Pr^0.42 = 0.877.
- Nu_avg = 163.3, so h_avg(2.5 D) = 1674 W/m²K. The planner's estimate was
  1.4–1.7 kW/m²K.
- The local value is h_loc(2.5 D, 7 D) = 1527 W/m²K. This is the **J-M
  h_ref**, local because the profile below is normalised by h_loc at 2.5 D.

## Profiles (parser expressions built in `jet.py`)

- **h(r, s)** = h_ref · h_loc(max(r, 2.5 D), clamp(s, 2 D, 12 D)) / h_loc(2.5 D, 7 D).
  - It is flat inside 2.5 D, with the 1/r-like Martin decay beyond.
  - The patch radius is the domain extent, so there is no sharp edge inside the
    hole.
  - The stagnation sensitivity ("rise") multiplies by
    1 + 0.5·max(0, 1 − r/2.5 D), so h reaches 1.5 × h_loc(2.5 D) at r = 0.
- **T_gas(s)** = T_ent + (T_ref − T_ent)·min(7/5, 7D/max(s, 2D)).
  - It is flat in the 5 D potential core and falls as 1/s beyond, equal to
    T_ref at 7 D, and radially uniform.
  - Default T_ent = 293.15 K; the sensitivity uses 600 K.
  - The core value is T_ent + 1.4 (T_ref − T_ent), which is 1983 K for J-M.
    This exceeds Meier's ≈ 1900 K adiabatic flame; it is reported, not
    changed.
  - The max(s, 2D) only affects columns outside the collision radius once the
    nozzle has passed their top (s ≤ 0).
- **Anchors (h_ref, T_ref):**
  - J-M = (1527, 1500 K);
  - J-5 = (5000, 1200 K), the bracket low end;
  - J-10 = (10000, 1000 K), Meier Ch. 7.

### AMReX parser gotcha (must read before writing expressions)

AMReX 25.12 (`ext/amrex/Src/Base/Parser/AMReX_Parser_Y.cpp`,
`parser_ast_optimize`, PARSER_DIV) rewrites `f / F2(x, number)` as
`f * F2(x, −number)`. That identity is intended only for `pow`, but it fires
for any two-argument function whose **second** argument is a number:

- `1527/max(r, 0.5)` evaluates to 0;
- `1527/min(r + 0.5, 1)` evaluates to −1527.

Found in the first probe: the flame was effectively off, with 0.1 W on the
face. `ext/` is do-not-touch, so every `max`/`min` here takes its numeric
argument **first**, e.g. `max(0.01775, r)`.

`run.py --parser-check` verifies the compiled expressions inside ALAMO. It
takes one step on the 3-D geometry (every column on the face form at
T = 293.15 K) and compares the domain-summed `patch_P_robin` with a numpy
face-balance replica. The check covers s0 = 1.5, 3, 5, 7, 10 and 13 D, the
rise and T_ent = 600 K sensitivities, and J-5 and J-10. All ten agree to
≤ 1×10⁻⁶.

## Settings

- **Common:**
  - `robin_form = pinned`, `pinned_idle_cycles = 2` (J-10 idle 1 is the stated
    sensitivity);
  - `surface.follow_mask = 1`, `energy_ledger.enabled = 1`,
    `spall.removal_events_csv = 1`, `beam.P0 = 0`;
  - **`weibull.V0 = V_cell`** (user decision 2026-09-15): dz³ = 8e-9 at 2 mm,
    1e-9 at 1 mm.
- **Time step:** the overshoot rule q·dt/(ρCp·dz) ≤ 5 K.
  - q is the maximum stagnation closed-form flux: r = 0, s over [2, 12] D,
    T_fire = 821 K.
  - dt = 4 ms, halved until the rule holds; every run lands on 4 ms (1.6–4.6 K
    per step).
  - C1 used a 1 ms cap. The 4 ms cap is justified by Part 1 `JM_dt1ms`: end
    centre depth −0.3 %, T_fire −0.01 %, R_h unchanged. The explicit
    conduction limit dz²/(6α) ≈ 1 s is far away.
- **Part 1 (2-D slab):**
  - `input_2d_dev` (70×4×60 over 0.14×0.008×0.12 m), x0 = 0.07, y0 = 0.004,
    radius = 0.07 (the whole top), `nozzle_z0 = 0.120 + 7 D`, feed 0, 60 s,
    4 ranks.
  - The slab is 4 cells thick with y0 on its mid-plane, so r includes a 1 or
    3 mm y offset; the radial bins use that r.
- **Part 2 (3-D quarter):** `input_drilling` overridden to:
  - 60×60×100 over 0.12×0.12×0.20 m;
  - axis at x0 = y0 = 0, the xlo/ylo corner, where Neumann-0 acts as symmetry;
  - radius = 0.2 (the whole top), `nozzle_z0 = 0.20 + 7 D`, `bit.enabled = 0`.
  - Stop at 300 s or when the face would reach 40 mm above the bottom,
    estimated from the closed-form stand-off lead. That gives 190 s at
    3 m/h and 300 s otherwise.
  - Cost probe (`--probe`, 10 s of J-M at 1.5 m/h): 2.3×10⁷ cell-steps/s, so
    ≈ 20 min per 300 s run alone on 4 ranks. No run exceeded the 45 min
    guideline, so none was shortened.
- **`surface_patch.nozzle_collision_radius = 2.5 D`** in Part 2 (user decision
  2026-09-15).
  - With a prescribed feed the nozzle passes the original surface z = 0.20 at
    119 s (1.5 m/h) or 60 s (3 m/h).
  - After that, every unrecessed outer column would have s ≤ 0, and the
    packet's "any in-patch column" abort would stop the run.
  - The key restricts the abort to columns with r ≤ 2.5 D. Outside it, s ≤ 0
    enters the expressions, which clamp s at 2 D.
  - The default is the patch radius, which gives exactly the packet's abort.
  - `patch_min_standoff` is the minimum over the collision radius.
- **Column-top patch only.** There is no Robin on side-wall faces, and the
  per-column flux uses the horizontal area dx·dy with no 1/cosθ factor. The
  max face slope is reported as a known limitation.
- **Not checkpointed:** the pin state and the ledger.
