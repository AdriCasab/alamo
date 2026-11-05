> INSTALLED as `ACTIVE_STEP.md` on 2026-10-13 (status line, column-kind guardrail, c6 volume basis, descent-mode label, clamp sentence added). Historical draft.

# Active Step: P2-3 — the floor coefficient to Kant, the wall on, the gas marched up the bore

Status: DRAFT (advisor-written 2026-10-13; becomes `ACTIVE_STEP.md` after `/plan` archives P2-2f)

## Context

P2-2f closed the rock half: the Grimsel properties transfer to Meier unrefitted (firing T 852.6 K vs
843.8 on Kant; price 1.27 GJ/m³ = ×1.37 on 0.930), and the residual is **delivery**, three ways at
once: 6.0 kW leaves as flakes vs Meier's 3.5; volume 1.85×; rate 3.47 vs 1.5 m/h. The energy partition
of `ctrl8/m13` (log 2026-10-13 §1) says why and where: only **25 %** of the jet enthalpy reaches any
rock face, the gas over the whole floor sits within 250 K of the nozzle, and delivery is
**coefficient-limited**. The model's axis coefficient **1,449 W/m²K** is 1.4–1.9× the only measured
values on a same-lab methane-air burner (Kant & von Rohr 2016 Table 4: 766 / 1,009 / 873 at 50 / 80 /
110 mm) while its gas temperature at the rock (1,760 K) matches the same table (1,683–1,723 K). The
user's decisions (2026-10-13): **scale the floor map to Kant's coefficient; turn the wall on; scale the
wall law by the same factor as a variant; cool the gas as it rises by marching it up the bore the way
it is marched out across the floor.** Entrainment, the cuttings sink and the body-cooling loss stay
parked; so do 1 mm and the constant-property legs.

**Why the bore is thin, and what "wall on" must not repeat.** `side_face_factor_above = 0.0`
(carried from CTRL5) gives every above-plane side face zero flux: 4,126 wall columns per quarter at
`side_P_face_above = 0`; the bore is 79.3 mm = 2·`foot_r_outer`, the tool stance. The wall coefficient
in use (P2-1 derived law, p50 372 W/m²K) is reasonable; the gas temperature every side face sees is
**one scalar**, `side_T_eff` = `jet_T_exhaust` = 1,504 K at every height (`MMWSpalling.H:2481`).
That uniform hot gas is the mechanism that ballooned the mouth in D2v's M2 (Ø 137 vs Meier 96 at
10 mm). The wall march replaces the scalar with a height-resolved temperature that the wall itself
depletes. ⚠ Two traps found while reading the code: (i) the wall law applies only to faces with
r ≥ `side_face_h_r_min` = 0.040 m, and the bore wall sits at r ≈ 39.7 mm — **inside the gate** —
so "wall on" without moving the gate heats nothing (D2ab's "empties the wall band" in a new form);
(ii) `side_T_eff` has **five** read sites and the criterion's (`Removal.H:605`) is the one that gets
missed (P2-0 lesson).

⚠ Carried: D2w measured that side heating on the floor cone flanks (below-plane, **already on** in
CTRL8) doubled the rate on top of a 2× floor delivery; D2x found no wall-band treatment closed the
bore without giving rate back — both under the old delivery. The order here is floor first, wall
second, for that reason. ⚠ D2ad: a bare level cut under prescribed descent starved the pad ring
before the standoff closed; under feet the standoff is pinned, but the pads at r 28–40 mm will see
the map's own radial decay (G ≈ 0.54 at 35 mm) × 0.68 ⇒ ~0.35 MW/m², marginal against Kant's
threshold class. **The pads descending on the 60 s leg is the gate that protects everything after.**

## Sources

- `tests/MMWSpalling/validation/kant/kant_von_rohr_validation_test_for_mmwspalling_plan.txt:66-80` —
  Table 4: SOD 50/80/110/130/150 mm → h_fl 766/1009/873/747/633 W/m²K, T_flame 1450/1410/1360/1250/
  1100 °C, q_sens 720/866/755/609/450 kW/m² (sensor plate, ±17 %). Linear interpolation at the m13
  median standoff **0.0853 m** gives **h = 986 W/m²K ⇒ scale 986/1449 = 0.681**.
- `studies/p2_2c_meier/feet_argv.txt` — the CTRL8 argv (76 keys). The floor map is the `h_expr`
  string; its peak literal is `1449.0*` (first factor inside `if(s>0.0, ...)`); the ring plateau is
  `142.0`, the beyond-76 mm value `0.001`. Scaling the **peak only** (1449 → 985) keeps the ring at
  142 and flattens the profile, which is the direction D2ad's measured floor (flat 0.48–0.60 MW/m² to
  35 mm) asks for; scaling the whole string would keep the 2:1 peak. **Deliverable 1 scales the peak
  only** and records the alternative.
- `studies/p2_2c_meier/ctrl8/m13/run/thermo.dat`, `run_jet_profile.csv` (150–350 s, quarter):
  `jet_P_cap` 10,739 W = fresh 7,598 + `jet_P_recirc` 3,141; `jet_P_face` 1,870; `jet_P_exhaust`
  8,868 (net out 5,727); `jet_T_exhaust` 1,504 K; floor bins 21, T_gas 1,760 → 1,504 K; `side_P_pads`
  629, `side_P_bore` 132, `side_P_wall` 26, `side_P_face_above` 0; `side_h_used_p50` 372, max 388;
  march mass m = 1.55·ṁ ⇒ **m c_p = 7.3 W/K**.
- `src/Integrator/MMWSpalling.H:3066` `JetEnthalpyMarch` (design comment `:3012-3058`: budget identity
  `jet_P_face + jet_P_exhaust == jet_P_cap` to 1e-9, abort); `:3588` `jet_T_exhaust = T_gas` after the
  last floor bin; `:3584` `jet_T_rec_next`. `:2481` `side_T_eff` resolution under `side_face_T = live`.
  Above-plane classification: a column with `flame_s[c] <= 0` (`:2250`, `:2329`, `:8160`).
  `side_face_factor_above` sentinel −1 = inherit (`:9269-9275`); kernel use `:2322-2330`.
- **The five readers of `side_T_eff`:** `Removal.H:605` (lateral criterion), `MMWSpalling.H:2047`
  (kernel), `:3419` (march), `:4374` (ledger), `:5350` (wall-law gas properties). **The four callers of
  `SideFaceH(h, r2, col, k)`:** `Removal.H:696`, `MMWSpalling.H:2285`, `:3427`, `:4378`; definition
  `:5127`; law comment `:5020-5034` (wall jet below z_fill, duct above).
- `tests/MMWSpalling/unit/jet_d2c/test:377-381` — the budget-identity assertions to mirror for the
  wall march. `unit/lateral_flux_lock/test` — the site-list lock to extend to five.
- `studies/p2_2c_meier/reduce_price.py`, `bore_profile.py`, `run_ctrl8.py` — reductions to reuse
  (price from `n_voided`, both paths, pad clearing excluded; Ø(z) from the void mask with geometry
  parsed from the argv).
- Meier targets (`validation/meier/meier-2017-pilot-660kg-validation-reference.md` §7 and
  `meier_fig8_8_hole_profile.csv`): rate 1.5 m/h ±20 %; Ø 82–87 mm at 150–450 mm (visible width,
  lower bound, ±5); 96 mm at 10 mm; 1.270 L in 500 s (×4 of the quarter).
- D2v M2 (`ARCHIVE_DONE.md` ~L12618-12630): above-plane f = 1 with uniform hot gas ⇒ Ø 137/135/132 at
  10/25/50 mm — the ballooning the negative control must reproduce.

## Goal

1. **OFFLINE — the scale, derived and printed** (`studies/p2_3_wall/derive_scale.py` →
   `SCALE.txt`): Kant Table 4 interpolated at the m13 median standoff; the scaled `h_expr` written
   out verbatim (peak literal 1449.0 → 985.0, nothing else touched) and parsed back to prove it is
   the same string with one literal changed; the value of the map at r = 0, 20, 35, 40 mm and
   s = 85 mm before/after; the pad-ring flux estimate h·(T_g − T_s) at r = 35 mm with T_g 1,504 K
   and T_s 850 K, against Kant's threshold-class 0.45–0.6 MW/m².
2. **SOURCE — the wall march** (`MMWSpalling.H`, after `JetEnthalpyMarch`): `WallEnthalpyMarch`
   under new key `surface_patch.jet_wall_march` (default 0 ⇒ byte-identical). Takes the floor
   march's exit state (`T_gas` after the last bin, `m`), enumerates the **above-plane** side faces
   (the existing void-aware enumeration and wetted-area law), bins them by cell height k from the
   lowest above-plane face upward, and per bin: `P_k = Σ h_face·max(0, T_gas − T_face)·A_face`,
   `T_gas ← max(T_ent, T_gas − P_k/(m c_p))`, writing `wall_T_gas[k]`. `jet_T_exhaust` becomes the
   temperature after the **last wall bin** (the mouth); `jet_T_rec_next` keeps the **floor-exit**
   value (entrainment is parked — state this in a comment). New thermo columns `jet_P_wall`,
   `jet_T_mouth`, `wall_bins`, `wall_L_fit` (the cooling length fitted to `wall_T_gas[k]`), and a
   `<plot_file>_wall_profile.csv` (time, k, z, n_faces, T_gas, P_k) at `jet_profile_interval`.
   Budget identity extended: `jet_P_face + jet_P_wall + jet_P_exhaust == jet_P_cap` to 1e-9, abort.
3. **SOURCE — one resolver for the wall-gas temperature:** `SideFaceTgas(col, k)` returning
   `wall_T_gas[k]` for an above-plane face when the march is on, `side_T_eff` otherwise. Replace
   the scalar at **all five** read sites (Sources). ⚠ The criterion site in `Removal.H` is the one
   to prove, not assume.
4. **SOURCE — the radius gate:** the derived wall law must reach the bore wall. Either lower
   `side_face_h_r_min` in the argv to 0.036 or derive the gate from `foot_r_outer − dx`; whichever,
   the leg's `side_h_used_n` on above-plane faces must be > 0 (witness). Record which.
5. **SOURCE — `surface_patch.side_face_h_scale`** (default 1.0, byte-identical): multiplies the
   derived law's output at the resolver, nowhere else.
6. **Lock extension:** `unit/lateral_flux_lock/test` gains a clause enumerating the five
   `side_T_eff` readers and asserting every one goes through `SideFaceTgas`; RED pre-change.
7. **Fixture** `unit/wall_march/test` (mirror `unit/jet_d2c`): a short feet leg with the march on —
   (a) identity `face + wall + exhaust == cap` to 1e-9 every step; (b) `wall_T_gas[k]` monotone
   non-increasing in k; (c) **negative case**: march on with `side_face_factor_above = 0` ⇒
   `jet_P_wall == 0` exactly and `jet_T_mouth == jet_T_exhaust(floor)`; (d) rank-independence of
   `jet_P_wall` (1 vs 4 ranks); (e) with the march **off** and factor_above = 1, `side_T_eff` is
   what every face sees (the old behaviour reproduces).
8. **Tree gate** on the rebuilt binary with the a2-fixed scorer and the `d2aa` positive control;
   declared diff set written first: **EMPTY** (both new keys default-off; the radius gate is an argv
   value, not a default).
9. **Three 60 s feet legs** (4 min each), all with floor peak 985 and `side_face_factor_above = 1`:
   `W1` march on, `side_face_h_scale` 1.0; `W068` march on, scale 0.68; `NC` **march off** (uniform
   1,504 K on every wall face) — the negative control that must reproduce M2-style widening.
10. **Two 500 s legs:** `W1` and `W068` (m = 13, the P2-2f configuration otherwise). `t_stop` from
    the 60 s descent rate as in P2-2f b2.
11. **Reductions:** `reduce_price.py` and `bore_profile.py` reused unchanged; new `wall_profile.py`
    (T_gas(z) and P_k per bin, cooling length by least squares on ln(T − T_ent) vs z, mouth gas T);
    an **isotropy witness** (bore by azimuth sector, spread) since removal geometry is touched.
12. Verdict.

## Criteria

| # | run order | test | bar | gate / report / verdict | what a failure of this row means |
|---|---|---|---|---|---|
| a1 | 1 | scale derived from Table 4 at the m13 standoff; scaled `h_expr` round-trips with ONE literal changed | **0.68 ± 0.01**; parse-back identical except `1449.0 → 985.0`; pad-ring flux at r = 35 mm **stated in MW/m²** | **gate** | the map was edited in more than one place, or the scale is not what Kant's table gives |
| a2 | 2 | fixture (a)–(e) | (a) identity ≤ 1e-9 of the cap; (b) monotone; (c) `jet_P_wall` **exactly 0**; (d) 1-rank == 4-rank to 1e-12; (e) old behaviour reproduces | **gate** | the march leaks enthalpy, heats where it must not, or is rank-dependent — every leg below is void |
| a3 | 3 | lock: five `side_T_eff` readers routed through `SideFaceTgas` | RED pre-change, green after, **each site named** | **gate** | a reader still sees the scalar; the criterion and the kernel would disagree on the gas temperature |
| a4 | 4 | tree gate, rebuilt binary, declared set EMPTY | **0 changed paths of 2639**, control 26 | **gate** | a default moved; nothing after is attributable |
| a5 | 4 | above-plane faces reached by the wall law | `side_h_used_n` (above-plane) **> 0** on the 60 s `W1` leg; the r-gate value recorded | **gate** | the bore wall is inside `side_face_h_r_min` and "wall on" heats nothing (D2ab's trap) |
| b1 | 5 | 60 s `W1`, `W068`: pads descend, rc 0 | descent `v` **> 0.8 m/h**, `foot_stall_time` < 10 s, rc 0 | **gate** | the scaled map starves the pad ring (D2ad) and the burner jams — the **finding** is then the map SHAPE, not the level; do not revert the scale |
| b2 | 5 | 60 s legs: ledger | `ledger_err_solid` ≤ 5e-3 on all three | **gate** | the wall march or the resolver broke the energy accounting |
| b3 | 5 | 60 s `NC` vs `W1`: the negative control balloons and the march does not | Ø at 10 mm depth: `NC` − `W1` **≥ 10 mm**; `W1` ≤ Meier 96 + 10 | **gate** | if `NC` does not widen, the wall heating is inert and the march is protecting nothing; if `W1` widens too, the march is not cooling the gas |
| b4 | 5 | `t_stop` from the 60 s rate | deepest column ≥ 60 mm above the guard at `t_stop` | **gate** | the leg bottoms out |
| c1 | 6 | 500 s `W1`, `W068`: cooling length and mouth gas T from `wall_profile.py` | L_fit **0.25–0.55 m**; mouth T **850–1,150 K**; `jet_P_wall` **0.8–2.5 kW** (quarter) | **report** | outside ⇒ the wall is taking far more or less than the derived law × area implies; check areas (wetted law) before physics |
| c2 | 6 | bore Ø(z), 150–450 mm, both legs, with sector spread | **82–87 mm** on at least one leg; mouth (10 mm) ≤ 106; sector spread ≤ 2 mm | **report** — headline | below 82 on both ⇒ wall flux below threshold (then the derived law, not the gate, is the question); above 87 or mouth > 106 ⇒ the march under-cools (check `m c_p` and the face areas) |
| c3 | 6 | rate, depth-matched 150–350 s | **2.0–2.7 m/h** (prediction); Meier 1.2–1.8 quoted beside | **report** | inside Meier's band would be a FINDING (delivery closed more than the coefficient alone predicts), not a pass; provisional at 2 mm |
| c4 | 6 | volume at 500 s, ×4 | **1.5–1.9 L** (prediction); Meier 1.27 quoted | **report** | with c3 inside their bands, a miss here is the wall's volume share — report it split floor/wall |
| c5 | 6 | ledger on the 500 s legs | `ledger_err_solid` ≤ 5e-3 | **report** | a late-run leak the 60 s gate could not see |
| c6 | 6 | firing T (regime-1 `T_surf`, population stated) + fraction above 873.15 K | **845 ± 10 K**; clamp fraction **12–25 %** | **report** | a shift > 10 K means the delivery change moved the criterion, which it must not; the clamp fraction should FALL with lower flux |
| c7 | 6 | price, both paths, pad clearing excluded | **1.20–1.30 GJ/m³** | **report** | outside ⇒ numerator (`E_flake` terms) before physics; the price is arithmetic in c6 |
| c8 | 6 | pad-ring flux and recession profile across r 28–40 mm | reported with the ring's `h`, `T_gas`, `T_s` in stated units | **report** | the number that decides P2-4 (map shape) if b1 was marginal |
| d1 | 7 | verdict | — | **verdict** | states how much of the 1.7× delivery excess the coefficient closed, whether the wall widens the bore to Meier without a funnel, and which of the parked levers (entrainment, cuttings, body loss) the residual now points at |

## Pre-registered predictions

- **P1 — the floor scale alone:** flake power 6.0 → **~4.1 kW** (Meier 3.5), rate 3.47 → **~2.4 m/h**,
  `jet_T_exhaust` (floor exit) up ~50 K; firing T and price unchanged within their bars.
- **P2 — the wall march cools the gas on its own:** with m c_p = 7.3 W/K and a quarter perimeter of
  0.0628 m, L = m c_p/(h P) = **0.46 m at h ≈ 250** (W068) and **0.31 m at 372** (W1) ⇒ mouth gas
  **~1,040 / ~920 K** for a 0.4 m hole; wall flux at the floor edge 0.2–0.3 MW/m². No body-loss
  parameter is needed for the mouth to stay below threshold.
- **P3 — the negative control balloons:** `NC` (uniform 1,504 K) widens the mouth by ≥ 10 mm over
  `W1` within 60 s (M2 reached Ø 137 at 10 mm).
- **P4 — the bore:** `W1` lands **82–87 mm** in 150–450 mm; `W068` lands 79–84 (marginal flux);
  both show a flat wall from 125 mm down (sector spread ≤ 2 mm).
- **P5 — the ring risk:** pad-ring flux at r = 35 mm ≈ 0.35 MW/m² after scaling — **marginal**.
  If b1 fails, the finding is the map's SHAPE (D2ad's flat-profile requirement), queued as P2-4;
  the scale is not reverted.
- **P6 — this step moves no criterion:** c6 inside 845 ± 10 on both legs; clamp fraction falls.

## What would make this inert — witness each

- **The wall law not reaching the bore wall** — a5's `side_h_used_n` on above-plane faces.
- **The resolver not reaching the criterion** — a3's five-site lock; a2 (e).
- **The march on but heating nothing** — a2 (c) must be exactly zero ONLY under factor_above = 0;
  on `W1`, `jet_P_wall` > 0 (c1's lower bound).
- **The negative control not ballooning** — b3's two-sided bar.
- **The scale edited in two places** — a1's parse-back.
- **A 60 s descent rate read from the t = 0 row** — skip t = 0 (P2-2f takeaway 5).

## Guardrails

- ⚠⚠ **Floor first, wall second, in the same binary but read in that order:** P1 (floor-only
  numbers) is reported from the `NC` leg's floor rows before any wall number is quoted.
- ⚠⚠ **The scale is Kant's, not a fit:** 0.68 comes from Table 4 at the leg's own standoff and is
  not adjusted to hit a Meier number. If b1 fails, the next step is the map shape (P2-4), not a
  different scale.
- ⚠ **Entrainment, cuttings sink, body-cooling loss, 1 mm, constant properties: PARKED** (user,
  2026-10-13). `jet_T_rec_next` keeps the floor-exit value; the wall march does not feed it.
- ⚠ **Cell mode only; both-path charged volume; no quarter-domain factor on a ratio** (P2-2f).
- ⚠ **`ledger_err` (total) is not a gate on any flake leg** — `ledger_err_solid` is.
- ⚠ **Carried argv must round-trip**; `h_expr` keeps its quoting; multi-value keys intact.
- **Wait on a PID or a sentinel file, never on a `pgrep` of your own command line.**
- **Every launch under `caffeinate -dis`; `amr.thermo.plot_int` explicit; archive output dirs before
  launching; `pkill -f` forbidden.**
- **Do not touch** `h_expr` beyond the one literal, the Grimsel laws, `R`, the flaw, the criterion,
  `side_face_factor` (below-plane), or the floor march's binning.

## Commands

```bash
# ord 1 -- OFFLINE
/Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/studies/p2_3_wall/derive_scale.py \
  > tests/MMWSpalling/studies/p2_3_wall/SCALE.txt

# ord 2-4 -- build, fixture, lock, tree gate (note the witness/ path)
EIGEN=$PWD/ext CPLUS_INCLUDE_PATH=/opt/homebrew/include \
  LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/Cellar/gcc/15.2.0_1/lib/gcc/current make -j8
/Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/unit/wall_march/test
/Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/unit/lateral_flux_lock/test
bash tests/MMWSpalling/studies/tree_sweep_0922/witness/guarded_post_all.sh p2_3

# ord 5 -- three 60 s feet legs (harness writes argv.sh per leg, sentinel .done on exit)
caffeinate -dis /Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/studies/p2_3_wall/run_p2_3.py short W1 W068 NC
# ord 6 -- two 500 s legs
caffeinate -dis /Users/tzetze20/Desktop/code/.venv/bin/python \
  tests/MMWSpalling/studies/p2_3_wall/run_p2_3.py full W1 W068
# ord 7 -- reductions
/Users/tzetze20/Desktop/code/.venv/bin/python tests/MMWSpalling/studies/p2_3_wall/run_p2_3.py reduce

# Cost: 3.77-4.61 wall-s per simulated s uncontended (P2-2f) => 4 min per 60 s leg, 31-38 min per 500 s leg.
# Planned total ~1.5 h: 12 min short + ~70 min full + reductions.
```
