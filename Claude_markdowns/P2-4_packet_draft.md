# Active Step: P2-4 — the cold body in the annulus: a second sink in the wall march

> DRAFT rev 2, 2026-10-13, by the reviewer agent, for installation via `/plan` once `/verify` has
> ruled on P2-3. Not installed. Rev 2 folds in: the depth-margin gate (the other agent: a 1,000 s leg
> bottoms out under branch B), S0 first, the c4 bar restated on a statistic that resolves it, and the
> user's decisions of 2026-10-13: ONE long leg carries both the gates and the finding (B1 is the first
> 500 s of it), W1L is NOT run, runs stay serial. Supersedes the reviewer's reply that predicted a
> mouth "frozen at ~108": it assumed 300 mm of immersed body; a 500 s run has 52 mm on average.

Status: draft (not ready for implementation until P2-3 is archived — `ARCHIVE_DONE.md` *"P2-3"*)

## Context

P2-3 turned the wall on with the derived law and marched the floor's exhaust up the bore. The march
is built, witnessed and byte-identical off, and it cools the gas by 53 K over the bore: the rock wall
is a near-equilibrium surface, not a sink. The gas reaches the mouth at ~1,490 K, the wall above the
nozzle plane spalls for as long as it is exposed, and the bore above the plane runs away (W1 mouth
117.8 mm vs Meier 96; 42 % of removed volume is wall; power into rock rose 35 % over P2-2f because
the wall and the outer floor absorbed what the floor released and the exhaust re-entrained the rest).
Below the plane the bore is right (85.3 at 200 mm depth vs Meier 85).

The one cold surface in Meier's hole is the burner's 67 mm water-cooled mantle (160 L/h), standing
in an 85–95 mm bore from the nozzle flange upward. The model has no such surface. This step adds it
to the existing march as a second, one-sided debit per bin, under keys that default to shipped
behaviour, and measures what it does to the gas, the wall and the bore. ⚠ Sized with the duct
coefficient the wall law already produces (83 W/m²K window mean, 64 at 500 s) and the immersion a run
actually has, the sink is WEAK at 500 s (100–130 K at the mouth) and moderate at 800 s (190–270 K).
The long leg is 800 s, the longest that keeps ≥ 60 mm of rock under both rate branches; it pre-registers
that the near-plane wall-jet region keeps spalling, and it is the first bore to reach seven of Meier's
eight calibrated depths.

Carried caveats: the body does not feed the floor — the gas the floor jet re-entrains is the pit gas
below the plane and `jet_T_rec_next` keeps the floor-exit value (MMWSpalling.H ~3700, "Entrainment is
PARKED"). The pre-crossing lip (~6 mm over the P2-2f control, cut under the free-surface plateau while
the burner was still above the rim) is not this step's. The cap rule, cuttings, ambient air at the
mouth and 1 mm stay parked (user, 2026-10-13). ⚠ With no matched long control (user's decision), the
second-half comparison rests on W1's realised 300–500 s growth extrapolated; the stated weakness is
that a lip fed by pre-crossing preheat could slow on its own. If the result is positive and contested,
`W1L` (800 s) is the follow-up, not a prerequisite.

## Sources

- **Body geometry.** `Claude_markdowns/meier_thesis_notes.md` L24 (chamber shield / cooling mantle),
  L43–45 (cooling water **160 L/h = 186 W per kelvin of rise; ΔT not reported**); memory
  `meier-burner-drawings.md` ("Mantle above is a 67 × 2.5 tube (Ø 67)"; Ø 80 tube 50 mm below the
  flange = the standoff section, uncooled, with the three exhaust slots). The wall law already encodes
  the mantle: `MMWSpalling.H:5376` `WJ_R_SHIELD = 0.0335` (duct gap = r_bore − 0.0335 at :5395, :5763).
  ⇒ **r_body = WJ_R_SHIELD; no new geometry number.**
- **The march.** `MMWSpalling.H:3708–3800` (`if (jet_wall_march && side_on && side_total > 0)`; pass 2
  bin loop from :3748; per-bin `T_in`, `q_sum`, `P_raw`, `T_out = max(T_ent, T_in − P_raw/mcp)` at
  :3793–3795; `T_gas = T_out`). `mcp = m_jet·jet_cp` at :3416 (jet mass incl. entrained; updated :3607).
  Closure :3847–3858: `jet_T_exhaust = T_gas + P_wall/mcp` (floor-exit meaning kept), `jet_P_exhaust =
  mcp(T_gas − T_ent)` at the mouth; invariants from :3860 (`face + wall + exhaust == cap`, 1e-9 abort).
  Sector fraction `wall_sector` (:5682, `PatchAzimuthFraction`, 0.25 on the quarter deck); duct
  coefficient `wall_h_duct` (per-step scalar, :5750–5790); readiness `wall_law_ready` (:5788; needs
  `r_bore > WJ_R_SHIELD + 2dx = 0.0375`, :5704). Resolver `SideFaceTgas` :5473; `SideFaceH` :5434.
  Column registration order is load-bearing for the tree gate (CODE_CENSUS F1): append new columns
  LAST (after `side_h_used_n_above`). Parse pattern for the P2-3 keys at :10109–10120.
- **Measured on W1 (150–350 s unless stated).** `mcp` effective 6.1 W/K (325.5 W / 53.1 K);
  fresh ṁ·c_p = 3.7826e-3 × 1250 = **4.73 W/K** (the deck's `jet_mdot` is the quarter's share of
  15.1 g/s: CH₄ 2.47 + air 52 kg/h). `wall_h_duct` **83.2** (window), 73.5 (350–500 s), 63.5 at
  500 s; `wall_z_fill` 0.055 m (window), 0.090 at 500 s; bore above the plane (0.5 − `nozzle_z`)
  **0.052 m** window mean, **0.165 m at 500 s**; `wall_T_mouth` 1,494.8 (window), 1,477 at 500 s;
  deepest void **299 mm at 500 s** (`bore_profile`); top-layer area-equivalent Ø (Ø_eq, 0–10 mm,
  §12 of the 2026-10-13 log) 117.8 at 300 s → **134.7 at 500 s = 0.085 mm/s**; wall share of removed
  volume 42 % (W1), 31 % (W074); mouth Ø (mask, sector mean) 117.8 / 104.5 / 122.4 (W1/W074/NC) with
  sector spread **10.6 mm** on W1; rate 1.633 / 1.933 / 1.609 m/h; `side_h_used_p50` 262 → 167
  (composite: wall-jet below `z_fill`, duct above — **the body sees the duct value only**; P2-1's 372
  was the composite p50 at another datum and is not the body's coefficient).
- **Rock depth.** `input_drilling`: `prob_lo z = −0.1`, rock top 0.500 ⇒ 600 mm of rock; usable
  depth with the 40 mm bottom guard and the ≥ 60 mm margin P2-2f/P2-3 enforced: **500 mm**.
  Projection from 299 mm at 500 s: branch A (1.63 m/h) 435 mm at 800 s, 526 at 1,000; branch B
  (2.15 m/h) 478 at 800, 598 at 1,000 ⇒ **800 s is the longest safe stop under both branches.**
- **Body perimeter, quarter:** P_q = 2π · 0.0335 · 0.25 = **0.05262 m**; area per 2 mm bin
  **1.052e-4 m²**. Cooling length L = mcp/(h · P_q): **1.40 m** (h 83, jet mass), 1.81 m (h 64),
  0.83 m (h 140, the value the law gives at a stance-width bore); with fresh mass ×0.78.
- **Machine.** Apple M5 Pro, 5 performance + 10 efficiency cores, 26 GB; 4 MPI ranks per run;
  3.80 wall-s per simulated s uncontended (P2-3 takeaway 12) ⇒ **500 s = 32 min, 800 s = 51 min**.
  Two concurrent 4-rank runs would share five fast cores and slow each other; stay serial. The grid is
  5 boxes of 60³ on 4 ranks (one rank carries two), so `-np 5` balances it: an optional 3-minute
  identity check (W1_s at −np 5 byte-identical to W1_s) would license it for the long leg.
- **Meier targets (standing):** Fig. 8.8 Ø(z) 96/92/88/87/87/87/87/86/85/85/84/84/83/82/80/78 at
  10…480 mm (±5 mm, a sawn width); rate 1.2–1.8 m/h; volume 1.27 L / 500 s; price 1.27 GJ/m³;
  drilling time for the 0.5 m hole ≈ 1,300 s.

## Goal

1. **Keys** (all default to shipped behaviour; each aborts if named without its feature):
   `surface_patch.jet_body` (0; needs `jet_wall_march`), `surface_patch.jet_body_T` (320.0 K; needs
   `jet_body`), `surface_patch.jet_body_h_scale` (1.0; needs `jet_body`), `surface_patch.jet_wall_mdot`
   = `jet` | `fresh` (`jet`; needs `jet_wall_march`).
2. **The body debit, inside the existing bin loop** (:3748–3800), evaluated at the bin's inlet `T_in`
   beside the rock faces' `q_sum`, and only when `wall_law_ready`:
   `P_body_k = (jet_body_h_scale · wall_h_duct) · (2π · WJ_R_SHIELD · wall_sector · DX[2]) · max(0, T_in − jet_body_T)`;
   `T_out = max(T_ent, T_in − (P_raw + P_body_k) / mcp_wall)`. One-sided like the rock faces. `mcp_wall`
   = `mcp` under `jet`, `jet_mdot·jet_cp` under `fresh`. ⚠ The body's h is the law's duct scalar via
   ONE read; no typed fallback — when the law is not ready the term is zero and that step is counted.
3. **Closure.** `jet_T_exhaust = T_gas + (P_wall + P_body)/mcp_wall` (floor-exit meaning kept);
   under `jet`: `face + wall + body + exhaust == cap`; under `fresh`: the entrained share never leaves
   — `jet_P_retained = (m_jet − jet_mdot)·c_p·(T_floor_exit − T_ent)` and
   `face + wall + body + exhaust_mouth + retained == cap`. Both forms under the binary's 1e-9 abort.
4. **Columns, appended last:** `jet_P_body` [W], `body_h_used` [W/m²K, the scalar used], `body_z_imm`
   [m, = wall_bins·dz], `body_L` [m, = mcp_wall/(h·P_q)], `body_noLaw_n` (steps this interval with
   bins but no ready law), `jet_P_retained` (gated on `fresh`), and `ledger_E_flake_wall` [J,
   cumulative flake energy booked on regime-6 removals — the per-regime price has been owed since
   P2-3]. The profile CSV gains a `P_body_k` column.
5. **Fixture** `unit/wall_march/test` rows (f)–(k): (f) body on, rock faces off (`factor_above 0`):
   profile exponential — fitted L vs `body_L` within 2 %, and `jet_P_body == mcp_wall·(T_in,0 − T_mouth)`
   to print precision; (g) `jet_body_T = 5000` ⇒ `jet_P_body` exactly 0 and the P2-3 fixture state
   byte-identical; (h) 1-rank == 4-rank; (i) key off ⇒ columns absent; (j) `fresh`: `exhaust_mouth +
   retained` equals the `jet` form's exhaust to 1e-9 on step 1 (bookkeeping only); (k) body on before
   any above-plane face exists ⇒ `jet_P_body` exactly 0 with `wall_bins` 0.
6. **Lock** `unit/lateral_flux_lock/test`: the body reads its coefficient through the law's duct
   scalar and nowhere else (no private `h`); `SideFaceH` role count and cap-site count UNCHANGED (the
   body is not a face). Text negative control as in P2-3.
7. **Tree gate**, declared set EMPTY, pin `p2_3.md5`.
8. **Study** `studies/p2_4_body/`: `run_p2_4.py` (base argv = `p2_3_wall/legs/W1/argv.sh`, 105 lines,
   with the P2-3 keys as run; **a `--gate-at 500` mode that reads the live `thermo.dat` and event log
   at the 500 s mark, scores b1–b5, and kills the run on any failure**); `reduce_p2_4.py` importing
   `reduce_p2_3`, `wall_profile`, `bore_profile` with windows 150–350 and 500–800 s, a per-regime price
   from `ledger_E_flake_wall`, the top-layer Ø_eq(t) (area of voided cells in the 0–10 mm layer,
   quarter ×4 — the statistic c4 is scored on), and the **cap share by radius and regime** (fraction of
   firings with `h_f` or `h_f_lat` at `flake_h_max`, bands 0–20/20–28/28–40/40–48/48–60 mm; a
   diagnostic carried for the parked cap step, no bar).
9. **Legs** (literal key lines; everything else = W1's argv), SERIAL, in this order:
   - `S0`  : W1's argv with `surface_patch.side_face_factor_above=0.0 surface_patch.jet_wall_march=0
     stop_time=500` — the floor scale ALONE, wall off: the control P2-3 never ran; it settles P5
     before anything else is interpreted. 32 min.
   - `B1L` : W1's argv with `surface_patch.jet_body=1 surface_patch.jet_body_T=320.0 stop_time=800`
     — ONE run: its first 500 s are the mechanism leg "B1" (b1–b5, c1 scored at the 500 s mark from
     the live logs; killed on a failure), its 500–800 s are the finding (c2–c9). 51 min.
   - `B1F` (optional, deferred): as B1L with `surface_patch.jet_wall_mdot=fresh stop_time=500`.
   - `W1L` (NOT run; the follow-up if the result is positive and contested): W1's argv, `stop_time=800`.
   Caffeinate every launch; nothing concurrent.
10. **Reductions and verdict** per the table.

## Criteria

| # | run order | test | bar | gate / report / verdict | what a failure of this row means |
|---|---|---|---|---|---|
| a1 | 1 | fixture (f)–(k) | (f) L within 2 %, identity to print precision; (g) **exactly 0** and byte-identical; (h) 1e-12; (i) absent; (j) 1e-9; (k) exactly 0 | **gate** | the body leaks, heats when it must not, reads a wrong area (a missing sector fraction reads as 4× cooling in (f)), or is rank-dependent — every leg below is void |
| a2 | 2 | lock | duct scalar read once; role/cap counts unchanged; control RED | **gate** | a private coefficient entered beside the law |
| a3 | 3 | tree gate | **0 of 2639**, control 26 | **gate** | a default moved |
| a4 | 4 | `S0` ran, rc 0, ledger ≤ 5e-3, pads descend | report its rate (P5) | **gate** on validity only | S0 void ⇒ P5 cannot be decided; proceed with B1L regardless |
| b1 | 5 (500 s mark) | `B1L` pre-crossing identity | event log byte-identical to W1 through the first above-plane step (≈148 s) | **gate** | the sink acted before an annulus existed — it is reading something other than the march's bins |
| b2 | 5 | `B1L` reach | `body_h_used` > 0 on ≥ 95 % of steps with `wall_bins` > 0; `body_noLaw_n` fraction reported | **gate** | the law's readiness gate (r_bore > 37.5 mm) starves the body — the leg measures nothing (D2ab's trap) |
| b3 | 5 | `B1L` ledger at 500 s | `ledger_err_solid` ≤ 5e-3 | **gate** | the closure broke |
| b4 | 5 | `B1L` pads descend | v > 0.8 m/h over 150–350 s, `foot_stall_time` < 10 s | **gate** | the burner jammed (not expected: the body is above the plane) |
| b5 | 5 | **depth margin** | deepest void at 500 s + v(350–500 s) × 300 s **≤ 500 mm**; and at the end, deepest void ≤ 500 mm | **gate** (kill at 500 s if projected over) | the hole reaches the bottom guard and everything after the breach is void — the P2-2f/P2-3 margin rule |
| c1 | 5 | `B1L` mouth gas at the 500 s mark vs W1 at 500 s | `wall_T_mouth` lower than W1's 1,477 K by **100–130 K**; window 150–350 s lower by 30–50 K | **report** | less ⇒ h or area under-read (check `body_h_used`, `body_z_imm`); more ⇒ the sector fraction or mcp is wrong |
| c2 | 6 | `B1L` mouth gas at 800 s | **1,190–1,270 K** (immersion 0.30–0.35 m, L 1.4–1.8 m) | **report** | the exponential with the logged h and immersion is the bar; a miss is an accounting error before it is physics |
| c3 | 6 | body power and water rise | at 500 s **0.65–0.85 kW** quarter (2.6–3.4 full, water 14–18 K); at 800 s **1.15–1.65 kW** quarter (4.6–6.6 full, 25–35 K) | **report** | above ~60 K of water rise the sink is implausibly strong and the coefficient, not the physics, is the question |
| c4 | 6 | top-layer **Ø_eq** (area statistic, resolution ±1–2 mm — NOT the sector-mean mask, whose spread is 10.6 mm) | at 500 s **128–132** (W1 134.7); growth over 500–800 s **≤ 0.042 mm/s**, half of W1's realised 0.085 over 300–500 s (⇒ ≤ +13 mm vs +25 extrapolated) | **report** — headline | growth at W1's slope ⇒ the lip is cut by preheated rock the gas temperature no longer governs — the cap/preheat, not the sink, is the lever; ⚠ a slowdown is the sink's OR preheat exhaustion's — the body power and the lip gas T logged in the same run are the tie-break, and `W1L` the follow-up |
| c5 | 6 | bore Ø(z) 150–450 mm, per sector, **seven depths** (450 falls in the floor pad) | reported with Meier beside; prediction: still above Meier at every depth above the plane, closer than W1 in the 25–150 mm section | **report** | at Meier ⇒ the sink alone closes it (not expected); at W1's values ⇒ see c4 |
| c6 | 6 | wall share of removed volume, 150–350 and 500–800 s | 500–800 s **20–30 %** (W1 42 % in 150–350 s; its 350–500 s value reported beside); the lowest 55–90 mm above the plane still widening | **report** | share unchanged ⇒ c4's failure branch |
| c7 | 6 | rate, 150–350 and 500–800 s, legs S0 and B1L | two pre-registered branches (P5) | **report** | neither branch ⇒ the pads, not the wall, moved — check `foot_carry_cols`, standoff |
| c8 | 6 | firing T, clamp fraction by volume | 845 ± 10 K; 12–25 % | **report** | the step moved the criterion, which it must not |
| c9 | 6 | price per regime (floor / wall) from `ledger_E_flake_wall`; cap share by radius/regime | reported, no bar (diagnostics for parked steps) | **report** | — |
| d1 | 7 | verdict | — | **verdict** | how much of the upper-bore excess the cold body explains at the law's coefficient over the depth a run reaches, which rate branch held, and which parked lever the residual points at |

## Pre-registered predictions

- **P1 — the sink is weak at 500 s and moderate at 800 s.** ΔT_mouth = (T_plane − T_body)(1 − e^{−z/L})
  with z the immersed length and L 1.4–1.8 m: 36–45 K (window), 106–134 K (500 s), 190–270 K (800 s,
  z 0.30–0.35 m). ⚠ The reviewer's "mouth ~1,190 K, frozen at ~108" of 2026-10-13 assumed z = 0.3 m
  throughout and is withdrawn.
- **P2 — the cooling arrives where the wall is already quiet.** The body's cooling accumulates with
  height; the lowest 55–90 mm above the plane (below `wall_z_fill`, wall-jet h ~250) sees gas cooled
  by < 50 K and keeps spalling. The duct region above (h 65–85) was not reaching 850 K within 300 s of
  exposure anyway (Robin semi-infinite, θ ≈ 0.4 ⇒ ~700–760 K) — its firings in W1 are preheated rock
  (regime-6 `ell_lat` 6–10 mm). So the sink's geometric effect is at the LIP and in the upper bore late
  in the run, not near the plane.
- **P3 — the mouth.** Ø_eq 128–132 at 500 s; the lip's growth halves once its gas falls below ~1,250 K
  (after ~600 s). The pre-crossing lip is identical to W1's (b1).
- **P4 — energy.** `jet_P_wall` falls with the gas temperature (500–800 s: 150–250 W vs W1's ~400 at
  the end of its run); `jet_P_body` as c3; floor rows unchanged.
- **P5 — rate, two branches, S0 decides.** (A, the reviewer's reading from the floor bins: the ring
  flux is identical between W1 and W074 and the rate is pad-ring limited) S0 within ±0.2 of W1's
  1.63 and B1L 1.6 ± 0.3 in both windows. (B, the other agent's reading: the in-band rate is geometric,
  the wall's sideways removal slows the descent) S0 ≈ 2.0–2.3 and B1L rises toward it in 500–800 s.
  **Either branch is a finding. Under B the rate leaves Meier's band WHILE the bore improves, and that
  combination is the result: the wall is NOT turned back on and no coefficient is re-tuned to recover
  the rate.** b5's projection is computed under whichever branch S0 selects.
- **P6 — T_body carries no magnitude:** 320 ± 30 K moves `jet_P_body` by ±2.5 %. Immersion is measured
  (`body_z_imm`), not assumed.
- **P7 — criterion untouched:** c8 inside its bars on every leg.

## What would make this inert — witness each

- **No annulus yet** (before the crossing): `wall_bins` = 0 ⇒ term 0 — physically right; b1.
- **Law not ready** (r_bore ≤ 37.5 mm, marginal at the crossing: 37.8): `body_noLaw_n`; b2.
- **Sector fraction dropped or doubled**: fixture (f)'s L and identity.
- **Faces handed the un-cooled gas**: the resolver reads `wall_T_gas[b]` = each bin's `T_in`, and bin
  k+1's inlet is bin k's `T_out` including the body — fixture (f) monotone profile, (g) byte-identity.
- **Scale key silently 0**: abort on `jet_body_h_scale` ≤ 0.
- **The lip slowing for its own reasons** (preheat exhaustion): c4's tie-break — the lip's gas T and
  the body power are logged in the same run; a slowdown that precedes the gas falling below ~1,250 K is
  not the sink's.

## Guardrails

- No edit to the floor map, the criterion, `SideFaceH`/`SideFaceTgas` signatures, or the wall law's
  branches (the duct scalar is READ, not changed). `jet_T_exhaust` keeps its floor-exit meaning;
  `jet_T_rec_next` keeps the floor-exit value (the body does not feed the floor — stated, not fixed).
- Keys default to shipped behaviour; new columns appended last; no typed body coefficient.
- Parked: cuttings, ambient entrainment at the mouth, 1 mm, the cap rule, the free-surface plateau.
- ⚠ **The rate guard (P5):** a rate that rises out of Meier's band while the bore closes is the
  pre-registered finding. Do not restore wall heating or re-tune any coefficient to recover it.
- ⚠ **The depth guard (b5):** no leg runs past the point where the projected deepest void exceeds
  500 mm; `stop_time` is 800 s and is not raised without re-doing the projection.
- Runs are serial (five fast cores, four ranks each). `-np 5` only after the 60 s identity check.
- ⚠ Per-step columns are MEANS over a window, never sums (standing rule 7(i)); `ledger_err` (total)
  is never a bar on a flake leg; the 1e-9 identity is the binary's abort, not a `thermo.dat` read.

## Commands

```bash
# build (CLAUDE.md env), then:
python3 tests/MMWSpalling/unit/wall_march/test            # a1 (f)–(k)
python3 tests/MMWSpalling/unit/lateral_flux_lock/test     # a2
tests/MMWSpalling/studies/tree_sweep_0922/witness/guarded_post_all.sh p2_4   # a3, vs p2_3.md5
caffeinate -dis python3 tests/MMWSpalling/studies/p2_4_body/run_p2_4.py S0          # 32 min, a4 / P5
caffeinate -dis python3 tests/MMWSpalling/studies/p2_4_body/run_p2_4.py B1L --gate-at 500   # 51 min; b1–b5, c1 at the mark, killed on failure
python3 tests/MMWSpalling/studies/p2_4_body/reduce_p2_4.py S0 B1L                   # c2–c9, d1
```
