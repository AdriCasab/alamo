# Active Step: P2-17 — Buckling on in production, the flake cap shown dead as a bound, and the Höser uptake ladder

Status: **accepted by review** (2026-11-03; all six gating rows reproduce independently, no diff-review row tagged `bug`). ⚠⚠ The declared 12.0 mm span needs **0.432 MW/m²** and the onset deck delivers **0.401** ⇒ with B on, the deck the onset level was calibrated on removes NOTHING; the two calibrations are now resolved on different configurations by construction (declared onset-only; DR-6).

## Context

1. **P2-16 (accepted, archived)** retired the constraint depth. Kant's two granites cannot separate "one span for
   all rocks" (R0) from "span ∝ 1/UCS" (R1): R1-over-R0 holds only if Bethel's onset equals Gotthard's (review
   DR-1). Grimsel's span is 11.7–13.2 mm under either rule. On a flat deck buckling is a net-flux condition
   (L_req ≈ const·μ·c_f·k/q_net).
2. **The plan of record switches B on BEFORE the uptake factor is fitted** (user, 2026-11-03), so the factor is
   calibrated on the final form once. The form is settled before the merge slices (Meier's feet → MMW), so they port
   it once.
3. **The declared span is 12.0 mm**: one span for both Kant granites, inside R0's window under both onset readings
   (Δ 0: [11.992, 12.010) mm; Δ −30 K: [11.683, 12.010) mm). It is chosen on parsimony: P2-16's formal verdict
   favoured R1 (Grimsel 13.20 mm), which is run here as the alternative span. The span is a rock property, so it goes
   in `rock_grimsel.inputs`.
4. **The 1.5 mm flake cap is dead as a flake bound under B.** Both act on the normal flake thickness μ·c_f·ℓ, and B
   withdraws a line whose flake exceeds L/C_L(T_s) < 1.5 mm below a crossover T_s* (z1). ⚠ It is NOT a no-op to
   switch off: the same value sets every line's gather and remap length (ledger L4). P2-17 keeps the value, proves the
   bound inert with a counter, and relabels it. Decoupling the window is a later source change.
5. **The calibration deck is onset-only** (P2-16): it keeps B off through a declared override.
6. **The ladder asks whether B moves the uptake factor.** It runs the Höser face-list deck (P2-16 `HOS_R` = `input_lab`
   + G1's argv) with the flame coefficient scaled by S. S is a study lever in argv only, never a model setting. It
   compares the pre-P2-17 form, the production form and the alternative span against Höser's 0 MPa volume. The
   direction is known in advance (B only removes firings, so S* can only rise); the size is the finding.

Carry-forwards: P2-16 DR-4 → z2; DR-12 → the rock file's comment and the model card; DR-11 → no width bar here.

## Sources (opened 2026-11-03)

- **B:**
  - `FaceRemoval.H:488` sets `bk_on = sp_buckle_span > 0`, and :607 evaluates it only for a fired line.
  - :617–620 computes x = `L.mu * sp_flake_c_f * ell_eff` (through lines: `L.mu * L_av`).
  - :315–325 gives `BuckleRatio` bk = C_L·x/L with C_L = π/√(6(1+ν)ε_th); B withdraws the line iff !(bk ≤ 1).
- **The cap:**
  - :505 `ReconFromX(…, spall_h_spall_max / L.mu, sp_flake_h_max / L.mu)`;
  - fallback :536–538 `h_bnd = sp_flake_h_max / L.mu; if (h_bnd > 0 && hf > h_bnd)`;
  - counter `fr_cap += 1` for a fired non-through line with `h_f_capped || h_f_bounded` (:704, thermo `fr_cap`).
- **The cap's window role:**
  - :368–380 `hcap = flake_h_max > 0 ? min(flake_h_max, h_spall_max) : h_spall_max`;
  - `L.Dmax = ⌊h_spall_max/(μ dx)⌋`; `L.Lg = Dmax + 2 + ⌈hcap/(μ dx)⌉`;
  - consumers of `L.Lg`: `o.run` (:497, logged), `L_av` (:695), `RecedeAlongLine(o.run, …)` remap (:716), the gather
    (:438);
  - `nprof = min(Dmax, run − 1)` (:552), so the depth scan itself does not depend on Lg.
- **C_L(T_s)**, from `studies/p2_15_buckle/derive_span.py` `C_L()`: 12.75 / 12.26 / 10.26 / 9.23 at 850 / 873 / 1000 /
  1100 K.
- **On-disk firings** (event logs in `studies/p2_16_span_rule/legs/`):
  - max T_surf: CAL 962.5 K (p99 913.9), HOS_R 1021.4 K (p99 919.0);
  - normal flake μ·c_f·ℓ_eff on HOS_R: median 0.557, p99 0.681, max 0.926 mm;
  - cap binding 0.00 % on both legs.
- **Column-path what-if** (`studies/whatif_uptake/legs/*/run_removal_events.csv`; B off, cap 1.5). Flake c_f·ℓ
  median / p90 in mm:

  | S | median / p90 | capped |
  |---|---|---|
  | 0.15 | 0.66 / 0.85 | — |
  | 0.10 | 0.84 / 1.04 | — |
  | 0.087 | 0.91 / 1.12 | 4.1 % |
  | 0.058 | 1.17 / 1.70 | 13.7 % |

- **The heat scale** (`studies/p2_11_heat/run_p2_11.py` `h_expr(S)`):
  - the G1 formula with leading constant C = S × 7272.3958, formatted `%.6f`;
  - G1 = S 0.1636 (C 1189.763953, = HOS_R's);
  - `SSTAR.txt`: G1_S14 (S 0.14) V(f) 30 s **24.922 cm³**; G1_S20 47.430 cm³;
  - S* rule: log V vs log S between the bracketing legs;
  - volume `reduce_p2_11.py` `fr_volume`: (cells − Σ fr_f)·dx²·dz·4 [cm³].
- **Höser 0 MPa** (`validation/hoeser/hoeser_2016_drilling_series.csv`, row `load,38,8,7.9,0,30,7.4,26,5.12`; read off
  Fig. 5.12 to ±5 % of axis span): 26 cm³/min ⇒ **13.0 cm³ at 30 s**.
- **P2-16:**
  - CAL (h 500) needs 12.84 mm at its earliest top firing (t 6.184 s, T_s 819.61 K);
  - span-off CAL flux legs `studies/p2_14a_cref_face/legs/{C0_h150, CAL, C0_h1000, CAL_h2000}` need
    49.4 / 12.84 / 6.50 / 3.61 mm;
  - R1's Grimsel band is 10.30–17.77 mm.
- **CAL deck keys:** `surface_patch.h_expr = 500.0`, `T_flame_expr = 1673.15`, `losses.T_amb = 293.15`,
  `losses.h_conv = 10.0`; ε = 0.8 (rock file).
- **AMReX `query` reads the LAST occurrence** (`ext/amrex/Src/Base/AMReX_ParmParse.cpp:1397`).
- **`check_shared.py`:** (i) accepts a deck-set shared key only via `--allow` (:104–111); (ii) resolves an allowed key
  to the deck's value (:117–119).
- **Runtime readers of the shared files:** only `unit/face_heat/test` reads `input_lab` at run time (its argv sets
  `spall.flake_h_max=0.0015`, no span). The other unit inputs are self-contained.

## Goal

1. **`validation/shared/rock_grimsel.inputs`:** a new block, "the buckling span (P2-15 condition B)", containing
   `spall.buckle_span = 0.012`. Its comment carries:
   - the class DERIVED-ON-KANT;
   - "one span for both Kant granites", with both R0 windows from Context 3;
   - the alternative 13.2 mm (∝ 1/UCS) and the band 10.30–17.77 mm;
   - "⚠ one calibration with `spall.flake_c_f`: B sees only L/c_f (SR 1)";
   - the minimum net flux from z2, and its side of h 500.
2. **`validation/shared/model_form_flame.inputs`:**
   - delete the `spall.buckle_span` line and its comment (moved to the rock file);
   - keep `spall.flake_h_max = 0.0015`, with the comment rewritten: "DEAD as a flake bound under B (P2-17 e0); still
     sets each line's gather/remap window (`FaceRemoval.H:368–380`): a NUM role, to decouple and retire with the
     column path".
3. **`validation/kant/grimsel_calibration/input`:**
   - directly after the two FILE lines: `# check_shared: allow spall.buckle_span` and
     `spall.buckle_span = 0.0   # ONSET-ONLY: B off by definition (P2-16); overrides the rock file`;
   - the header's onset-only paragraph and `README.md` say the override is declared.
4. **`validation/shared/check_shared.py`:** a deck line matching `^#\s*check_shared:\s*allow\s+(\S+)` adds that key to
   `allow`. Read it by its own scan, since `parse()` strips comments. Update the docstring.
5. **Study `tests/MMWSpalling/studies/p2_17_buckle_on/`:**
   - `run.py` (the legs below; per leg: `argv.sh`, `ALLOW.txt`, `BINARY.txt`, `.done`);
   - `offline.py` (z1, z2);
   - `chk_directive.py` (z0);
   - `reduce.py` (c-rows, e0, e1) and `chain.sh`;
   - `inputs_pre/` (cp -p of the three edited inputs and `check_shared.py`).
6. **`docs/project/MODEL_CARD.md`:**
   - §6 row 21: the cap is dead as a bound under B; its remaining role is the NUM window; no longer TUNED in effect;
   - row 28: 12.0 mm declared, the band, the alternative, q_net,min and the side of h 500;
   - §9: the onset deck's declared override;
   - the uptake-gap entry: the ladder's S*(O) and S*(N) and the gap on both S and absorbed energy.

## Reachability ledger (step zero resolves every row; stop on a failed row)

| # | key / quantity | site and guard | units |
|---|---|---|---|
| L1 | `spall.buckle_span = 0.012` in the rock file | Parse → `sp_buckle_span`; `bk_on` :488; evaluated `if (bk_on && o.fire)`. The column path refuses > 0. **Step zero greps that only `grimsel_calibration/input` and `hoeser/input_lab` include the rock file**, and that no tree input does | m |
| L2 | CAL's `spall.buckle_span = 0.0` after the FILE lines | AMReX LAST occurrence (:1397) ⇒ `bk_on` false | m |
| L3 | **THE CONJUNCTION: is the cap dead under B?** | the cap compares along-line h_f with `flake_h_max / L.mu`, i.e. the normal μ·h_f with 1.5 mm. B compares C_L(T_s)·μ·c_f·ℓ_eff with L. Both act on μ·c_f·ℓ, so obliquity cancels. The cap can act on a passing line only if L/C_L(T_s) > 1.5 mm, i.e. T_s > T_s*(L). Witness: thermo `fr_cap` (:704) | m; K |
| L4 | ⚠ the cap's window role | `hcap` → `L.Lg` → `o.run` / `L_av` / remap length (:368–380, :497, :695, :716). Off ⇒ every line goes from 19 to 32 cells at 1 mm and the deep enthalpy remap moves. **Kept at 0.0015**; the step zero records each consumer | cells |
| L5 | the heat scale | only `surface_patch.h_expr`'s leading constant, C = S × 7272.3958 (`h_expr(S)`); S 0.1636 reproduces the base's 1189.763953 | W/m²K |
| L6 | B's counters | thermo `bk_n_eval`, `bk_n_block`, `_floor/_wall/_roof`, `bk_max_pass`; present only under `bk_on`. **Step zero states whether they and `fr_cap` are per-step or cumulative** (SR 7(i): a fraction is two sums) | counts |
| L7 | volume and target | `fr_volume` at the 30 s plotfile, quarter domain ×4; target 13.0 cm³ | cm³ |
| L8 | the directive | `# check_shared: allow KEY` in the DECK file ⇒ KEY ∈ allow | — |

## Legs (exact argv tokens)

**Base (B0):** `studies/p2_16_span_rule/legs/HOS_R/argv.sh` token for token, with these changes:
- `stop_time=30.0`;
- `plot_file=…/p2_17_buckle_on/legs/<leg>/run`;
- `surface_patch.h_expr=<h_expr(S)>` replacing the base's (the only heat change).

`amr.plot_int=1875` stays, giving plotfiles at 15 s and 30 s.

| form | added tokens | `check_shared` |
|---|---|---|
| **O** (the pre-P2-17 form) | `spall.buckle_span=0.0` | `--allow spall.buckle_span` |
| **N** (production) | none | — |
| **A** (alternative span) | `spall.buckle_span=0.0132` | `--allow spall.buckle_span` |

**S grid:**
- O and N at S ∈ {0.20, 0.1636, 0.14, 0.12, 0.10, 0.087, 0.07};
- A at S ∈ {0.14, 0.12, 0.10, 0.087, 0.07};
- names `<form>_S<round(1000·S)>`: O_S200, O_S164, O_S140, O_S120, O_S100, O_S087, O_S070, and likewise for N and A
  (19 legs).

**Other runs:**
- **CAL1:** `studies/p2_16_span_rule/legs/CAL/argv.sh` verbatim (only `plot_file` re-pointed), on the edited deck.
- **face_heat:** `tests/MMWSpalling/unit/face_heat/test` as is.

## Criteria

| # | run order | test | bar | gate / report / verdict | what a failure of this row means |
|---|---|---|---|---|---|
| z0 | 1 | `chk_directive.py`, two-sided: the edited CAL deck; a scratch copy without the directive; a copy whose directive names an unset key | CAL (i) PASS, 1 declared; no-directive copy (i) FAIL naming `spall.buckle_span (deck)`; unset-key copy (i) FAIL | **gate** | the override is not declared, so the onset deck's B-off status is unchecked |
| z1 | 1 | cap-dead crossover T_s*(L), where C_L(T_s*) = L / 1.5 mm, for L = 12.0, 13.2 and 17.77 mm; beside it the on-disk max firing T_s | report | report | P1 |
| z2 | 1 | **DR-4, minimum flux at a span:** h_min(L) by log-log interpolation of L_req vs h between the two span-off CAL flux legs bracketing L (L_req recomputed by P2-15's discriminant as P2-16 `reduce.py` imports it, not copied); q_net,min = h_min·(T_g − T_s,e) − εσ(T_s,e⁴ − T_amb⁴) at that leg's earliest firing [MW/m²] (step zero confirms against `SurfaceCellFlux`). For L = 10.30, 12.0, 13.2, 17.77 mm; which side of the CAL deck's own q_net | report | report | P2 |
| a1 | 2 | CAL1 vs P2-16 `legs/CAL` (= P2-13) | event log byte-identical; kept plotfiles identical; thermo non-timer columns identical; `check_shared` PASS with 1 declared; metadata `spall.buckle_span` 0 | **gate** | the shared-file edit moved the onset deck; nothing below is attributable |
| a2 | 2 | face_heat on the edited files | rc 0, verdict PASS (a moved number is declared together with that run's Σ `bk_n_block`) | **gate** | B reaches the face-heat accounting; stop and declare |
| a3 | 2 | O_S164's event-log rows with time ≤ 15.0 vs P2-16 HOS_R's whole event log; its 15 s plotfile vs HOS_R's | identical | **gate** | the O tokens do not restore the pre-P2-17 form |
| a4 | 2 | O_S140 V(f) at 30 s vs P2-11 G1_S14 | 24.922 ± 0.001 cm³ | **gate** | the heat-scale path or the volume instrument moved since P2-11; S* would not compare |
| a5 | 2 | N_S070 reaches B: `check_shared` (ii) resolves `spall.buckle_span` 0.012; Σ `bk_n_eval` > 0 | yes | **gate** | the span never reached the binary; the ladder would be inert |
| c1 | 3 | blocked fraction Σ `bk_n_block` / Σ `bk_n_eval` per N and A leg, with n_eval beside it (counts only if n_eval < 100), and the floor/wall/roof split | report | report | P4 |
| c2 | 3 | **S\*(O), S\*(N), S\*(A)** (V(f) 30 s vs 13.0 cm³, P2-11's rule); S\*(N)/S\*(O), S\*(A)/S\*(N); the gap 1/S\*; `ledger_E_src` ×4 at 30 s, interpolated to S\* the same way [kJ] | report | report | P3, P4, P5 |
| c3 | 3 | flux exponent n = Δln V / Δln S on the pair bracketing each S\*, O vs N | report | report | P6 |
| c4 | 3 | isotropy witness: removed volume per 11.25° sector (8 sectors, quarter domain, r ≤ 25 mm) at 30 s, divided by each sector's cell census; spread = max/min − 1; N vs O at S 0.10 and 0.07 | report | report | N's spread above O's by > 0.05 ⇒ B introduces a grid-direction bias |
| e0 | 4 | **the cap is dead as a bound:** Σ `fr_cap` = 0 on every N and A leg; positive control Σ `fr_cap` > 0 on O_S070 and O_S087; max firing T_s on every N/A leg vs T_s*(its span) | all three | **verdict** | a cap count on N/A ⇒ the bound still acts and the relabel is wrong. A zero positive control ⇒ the cap never acted on Höser, and its deadness rests on z1 alone |
| e1 | 4 | **verdict** | — | **verdict** | B on is the production form (a1–a5); the cap's status (e0); P4's branch; the span choice's effect on S\* (P5); the declared span's q_net,min and side of h 500 (z2) |

## Pre-registered predictions

- **P0 (a rows, z0):** every identity holds; a5 reads `bk_n_eval` > 0.
- **P1 (z1):** T_s* = 1273 K (12.0 mm), 1152 K (13.2 mm), **893 K (17.77 mm)**. The on-disk maxima are 962.5 / 1021.4 K,
  so the bound is dead at 12.0 and 13.2 mm, but not across R1's whole band.
- **P2 (z2):**

  | span | h_min | q_net,min |
  |---|---|---|
  | 12.0 mm | ≈ 536 W/m²K | **≈ 0.44 MW/m²**, above the CAL deck's ≈ 0.41 |
  | 13.2 mm | ≈ 488 W/m²K | ≈ 0.40 MW/m², below |
  | 10.30 / 17.77 mm | — | ≈ 0.51 / 0.30 MW/m² |

  Bars are ±5 %. ⇒ At the declared span, Grimsel's minimum flux sits ABOVE h 500: the onset deck would not spall with B
  on, which is why it is onset-only.
- **P3 (c2, O):** S*(O) = 0.097 (0.088–0.106). This extrapolates G1_S14's 24.922 cm³ with n 1.5–2.2; P2-11's
  0.14/0.20 pair gives n 1.80.
- **P4 (c1, c2, N vs O):**
  - **Direction:** S*(N) ≥ S*(O), because B is a pure withdrawal (P2-15 DR-3). S*(N) < S*(O) by more than 1 % is a
    finding: B's withdrawals raised removal elsewhere.
  - **Blocked fraction:** < 1 % at S 0.20 and 0.1636 (P2-15 HOS_S: 0.12 % at 11 mm); 2–15 % at 0.10; 25–75 % at 0.07.
    The face list's flakes run ≈ 0.85× the column path's at equal S (HOS_R 0.557 vs what-if 0.66 mm).
  - **Size:** S*(N)/S*(O) **1.01–1.25, central ≈ 1.05**; the gap 1/S* narrows from ≈ 10× to ≈ 9×.
  - **Branches, pre-written:**
    - **(a)** 1.02–1.5: B trims the uptake gap but leaves it open (> 6×). The open problem stands; item 5 fits S with
      B on.
    - **(b)** < 1.02: B is inert at Höser's uptake level, and B's role on Höser is the minimum flux only.
    - **(c)** > 1.5: B carries a large part of the uptake story, or the declared span is too short for Höser. Before
      item 5, check the span against Höser's own flux.
- **P5 (c2, A vs N):** S*(A)/S*(N) ∈ [0.95, 1.00]. A longer span blocks less; for Höser's calibration the span choice
  is immaterial.
- **P6 (c3):** n(N) ≥ n(O) on the S* bracket, with n(O) 1.5–2.2.
- **Refuting evidence on disk:** none. B has never run on the face list below S 0.1636; P2-15 HOS_S is consistent with
  P4's < 1 %.

## Guardrails

- **No `src/` change, no build.** The tree is not re-run (L1's grep).
- **`spall.flake_h_max` keeps 0.0015** (L4). Only its comment changes, and no code is deleted.
- **S lives in leg argv only.** S* is not written into any input. The span is never tuned on Höser; A is the declared
  alternative, not a fit.
- **Edits are limited to Goal 1–4, Goal 6 and the study.** `input_lab` is untouched. `inputs_pre/` is copied before the
  first edit.
- **Plotfiles:** keep 15 s and 30 s on every ladder leg; CAL1 keeps P2-16 CAL's set.

## Commands

```bash
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; S=tests/MMWSpalling/studies/p2_17_buckle_on
# before editing: cp -p rock_grimsel.inputs model_form_flame.inputs grimsel_calibration/input check_shared.py -> $S/inputs_pre/
# ord 1 (zero compute): z0, z1, z2
$PY $S/chk_directive.py && $PY $S/offline.py
# ord 2: CAL1, face_heat, O_S164, O_S140, N_S070 -> a1-a5 (stop on a failed gate)
# ord 3: the remaining 16 ladder legs (serial, -np 4) | ord 4: reduce.py -> RESULTS.txt (c1-c4, e0, e1)
caffeinate -dis bash $S/chain.sh
```

## Claude completion notes

**Files changed (no `src/` change, no build):**
- `validation/shared/rock_grimsel.inputs`: a new block "the buckling span (P2-15 condition B)", setting
  `spall.buckle_span = 0.012` (Goal 1, plus DR-12 and z2's numbers).
- `validation/shared/model_form_flame.inputs`: the span line removed; `spall.flake_h_max` relabelled NUM (dead as a
  bound; the window role kept).
- `validation/kant/grimsel_calibration/input` and `README.md`: the directive and the `spall.buckle_span = 0.0`
  override after the FILE lines, declared in both.
- `validation/shared/check_shared.py`: the `# check_shared: allow KEY` directive, scanned on its own, deck file only;
  docstring updated.
- `docs/project/MODEL_CARD.md`: §6 rows 21 and 28, §9's B row, and §11 item 1 (the uptake gap). The old card is
  archived as `MODEL_CARD.pre_p2_17.md`.
- Study `tests/MMWSpalling/studies/p2_17_buckle_on/`: `run.py`, `offline.py`, `chk_directive.py`, `reduce.py`,
  `chain.sh`, `chain_resume.sh`, `inputs_pre/`, `Z0–Z2`/`A1–A5`/`A3.first`/`RESULTS`/`RESULTS.first.txt`,
  `face_heat_output_pre/` (SR 9), and `legs/` (20 legs, 4.9 GB with the 15 s and 30 s plotfiles).

| # | measured | verdict | output |
|---|---|---|---|
| z0 | D1 the edited CAL deck: (i) PASS, 1 declared. D2 no directive: (i) FAIL, `UNDECLARED: spall.buckle_span (deck)`. D3 an extra directive for the unset key `weibull.m`: (i) FAIL, `nobody sets: weibull.m` | **PASS** | `Z0.txt` |
| z1 | T_s* 1272.7 K (12.0 mm; C_L at the table's end 7.997 against 8.000), 1151.9 K (13.2 mm), 894.0 K (17.77 mm). On-disk max firing T_s: CAL 962.5 K, HOS_R 1021.4 K (p99 913.9 / 918.8) | report | `Z1.txt` |
| z2 | 12.0 mm: h_min 535.5 W/m²K, **q_net,min 0.4316 MW/m², ABOVE** the CAL deck's 0.4014. 13.2 mm: 487.7, 0.3909, below. 10.30 mm: 0.5081. 17.77 mm: 0.2940. All with the code's h_conv term; the packet's form reads +0.005 MW/m² | report | `Z2.txt` |
| a1 | CAL1 vs P2-16 CAL: event log identical; 98 plotfile files, 0 differ; thermo 7500 rows × 77, 0 differ; `check_shared` PASS (1 declared, all 64 keys resolved); metadata `spall.buckle_span` 0.0 | **PASS** | `A1.txt` |
| a2 | face_heat rc 0, PASS; log identical to UNIT_p2_16 (masked). B evaluated 0 lines on its one-step legs (Σ `bk_n_block` 0 / Σ `bk_n_eval` 0) | **PASS** | `A2.txt` |
| a3 | **first scoring FAIL** (5070 vs 5066 rows with time ≤ 15.0; `A3.first.txt`). Restated: the rows of the steps HOS_R ran (start < 15.0 − dt/2) are 5066 / 5066, byte-identical; the 4 excluded rows belong to the step starting at t = 15 s, stamped 14.999999999999234. 15 s plotfile: 14 files, 0 differ | **PASS** (restated) | `A3.txt` |
| a4 | O_S140 V(f) 30 s **24.9222 cm³** against 24.922 | **PASS** | `A4.txt` |
| a5 | N_S070: `check_shared` PASS, metadata `spall.buckle_span` 0.012, Σ `bk_n_eval` 27,740 (sampled) | **PASS** | `A5.txt` |
| c1 | blocked fraction (evaluation-weighted), N: 0.0070 / 0.0316 / 0.178 / 0.609 / 0.789 / 0.956 / 0.993 at S 0.20 … 0.07; A: 0.035 / 0.282 / 0.472 / 0.906 / 0.985 at S 0.14 … 0.07. **Per firing (removals N/O):** 1.003 / 0.998 / 0.997 / 0.991 / 0.971 / 0.924 / 0.673. Roofs 0 everywhere; walls take most blocks at S 0.12–0.087, floors at 0.07 | report | `RESULTS.txt` c1 |
| c2 | **S\*(O) 0.09944** (bracket O_S087–O_S100), **S\*(N) 0.10135** (N_S100–N_S120), **S\*(A) 0.10028**. S\*(N)/S\*(O) **1.0192**, S\*(A)/S\*(N) 0.9894. Gap 1/S\* 10.06× → 9.87×. `ledger_E_src` ×4 at S\*: 24.71 / 25.18 / 24.89 kJ (N/O 1.019) | report | `RESULTS.txt` c2 |
| c3 | n(O) 1.832 (on 0.087–0.10), n(N) 2.046 (on 0.10–0.12), so n(N) − n(O) = +0.214. ⚠ The two brackets are different pairs | report | `RESULTS.txt` c3 |
| c4 | spread O / N: 0.1067 / 0.1447 at S 0.10, 0.2177 / 0.2551 at S 0.07; N − O +0.038 / +0.037 (under the 0.05 bar) | report | `RESULTS.txt` c4 |
| e0 | Σ `fr_cap` = 0 and event-log cap flags = 0 on all 12 N/A legs. Positive control: O_S087 3 and O_S070 1 cap events (thermo, sampled: 0 and 1); O_S100 also 2. Every N/A max T_s (≤ 1057.9 K) is below its T_s* | **verdict: the cap is DEAD as a bound under B** | `RESULTS.txt` e0 |
| e1 | B on is the production form (z0, a1–a5 PASS). Cap dead (e0). **P4 branch (b): 1.0192 < 1.02, by 0.0008, so at the (a)/(b) edge.** The span choice is immaterial (P5). Declared span: q_net,min 0.432 MW/m², above h 500 | **verdict** | `RESULTS.txt` e1 |

**Prediction scores.**
- **P0** held. a3 needed one comparator restatement; the identity itself held byte-for-byte.
- **P1** held: 1272.7 / 1151.9 / 894.0 K against 1273 / 1152 / 893.
- **P2** held: 535.5 / 0.432 (against ≈ 536 / 0.44), 487.7 / 0.391, 0.508 / 0.294, all inside ±5 %.
- **P3** held: S\*(O) 0.0994 in [0.088, 0.106].
- **P4:**
  - direction held (1.019 ≥ 1);
  - size in range (1.01–1.25) but under the 1.05 central;
  - **blocked-fraction bands MISSED** (3.2 % at 0.1636 against < 1 %; 79 % at 0.10 against 2–15 %; 99.3 % at 0.07
    against 25–75 %), because they were written per firing and the counter is per evaluation (takeaways);
  - branch (b) by 0.0008.
- **P5** held: 0.989 in [0.95, 1.00].
- **P6** held: n(N) 2.046 ≥ n(O) 1.832, and n(O) is in 1.5–2.2.

## Implementation takeaways

- **Deviation: a3's filter.** The packet's "time ≤ 15.0" also admits the step that STARTS at t = 15 s. The event log
  stamps a step's start (SR 7(iii)), and the time accumulated in floating point to 14.999999999999234, while HOS_R
  (stop 15.0) never ran that step. Restated by step (start < 15.0 − dt/2); first scoring kept in `A3.first.txt`. The
  chain resumed from a4 (`chain_resume.sh`); O_S164 was not re-run. ⇒ **A packet that truncates a longer run to a
  shorter one's horizon must cut by step index (or start < stop − dt/2), never by `time ≤ stop`.**
- **`bk_n_eval` / `bk_n_block` count evaluations, not lines.** A fired line is evaluated once per step. A line whose
  onset keeps passing while B blocks it is re-counted every step, so the blocked fraction is evaluation-weighted. At
  S 0.10, 79 % of evaluations are blocked but only 2.9 % of removals are lost. The per-firing measure is the removal
  count N/O (c1's second block). Any later bar on "how much B blocks" should name which of the two it means.
- **Every B and cap thermo column is PER-STEP, and the Höser deck samples thermo 1 in 6** (`amr.thermo.plot_int = 6`;
  `FaceRemoval.H:341` resets `fr_*` each step, as P2-15 does `bk_*`). Ratios are fine; counts are samples. e0 needed
  the event log: O_S087's sampled Σ `fr_cap` read 0 while its event log carries 3 cap events (rc_flag bit 16).
- **The positive control is tiny:** 1–3 cap events in about 1,700–3,800 removals on O_S070–O_S100, and 0 above
  S 0.10. On Höser the 1.5 mm cap barely ever bound even with B off; its "deadness as a bound" is real but cheap. Its
  window role (L4) is the one that matters, and it is untouched.
- **B at S 0.20 adds removals** (N/O 1.003, V +0.28 %): withdrawals reshape the face, and later firings elsewhere
  slightly outnumber them. That is below P4's "> 1 %" finding threshold, but it is the sign P2-15 DR-3's "pure
  withdrawal" does not guarantee on a drilling deck.
- **c3's exponents come from different S pairs.** S\*(N) moved across the S 0.10 grid point, so n(O) is read on
  0.087–0.10 and n(N) on 0.10–0.12. Part of the +0.214 is the curvature of V(S), not B.
- **z2 uses the code's balance:** `SurfaceCellFlux` subtracts h_conv (T_s − T_amb) with h_conv 10 W/m²K, which the
  packet's formula omits (−0.005 MW/m², about 1.2 %). T_s at h_min is linear in ln h between the bracketing legs
  (819.6–821.6 K; immaterial).
- **The directive is read from the DECK file only.** An included file cannot declare its own override; that is
  deliberate, so a shared file cannot launder a change to itself.
- **`face_heat` stays rc 0 with the span in `input_lab`.** Its one-step legs never fire a line (Σ `bk_n_eval` 0).
  `R_nolist` (face_list = 0) still aborts naming `face_heat`, because the `surface_patch` parser's check
  (`MMWSpalling.H:7914`, parsed at :190) runs before B's face-list abort (:950). ⚠ Any future deck that includes
  `input_lab` or the rock file with `face_list = 0` and no `face_heat` key will now abort on `spall.buckle_span > 0`.
- Nothing is committed. `legs/` is 4.9 GB (15 s and 30 s plotfiles on every ladder leg, per the guardrail).


## Review findings

Reviewer, 2026-11-03. Independent reduction throughout: my own `cmp` / `diff -rq` / per-column byte comparisons, my
own volume, thermo, event-log and sector reducers, my own log-log solves, my own 1 mK scan of `C_L`, my own mutation
copies of the CAL deck, and the source read of the cap-vs-B conjunction. No test, simulation or build was re-run —
every artefact cited below was already on disk with the binary's hash (`38c8623885c5b8b6…`) beside it, and I
re-verified that hash against the live binary and that **no file under `src/` is newer than it**. The only zero-compute
scripts I executed were `check_shared.py` (on five scratch copies of the CAL deck) and my own reducers.

| # | test | bar | **derived independently** | **read from** | row |
|---|---|---|---|---|---|
| z0 | the directive, two-sided | D1 PASS/1 declared; D2 FAIL naming the key; D3 FAIL naming the unset key | D1 `(i) PASS … (1 declared)`; D2 `(i) FAIL … UNDECLARED: spall.buckle_span (deck)` **and (ii) FAIL** (the directive is load-bearing for the resolution too); D3 `(i) FAIL … nobody sets: weibull.m`. **My two extra mutations:** directive kept, override line deleted → `(i) FAIL … nobody sets: spall.buckle_span`; directive with collapsed whitespace → PASS | my own 5 scratch copies of `validation/kant/grimsel_calibration/input` through `validation/shared/check_shared.py` | **reproduces (gate)** |
| z1 | T_s\*(L) with C_L(T_s\*) = L/1.5 mm | report | **1272.7 / 1151.9 / 894.0 K** at 12.0 / 13.2 / 17.77 mm; C_L at the table's end **7.9973** and **flat above it** (C_L(1500 K) = 7.9973); at R1's 10.30 mm band edge **T_s\* does not exist**. On disk max T_surf 962.5 K / p99 913.9 (CAL, 10,804 firings), 1021.4 / 918.8 (HOS_R, 5,066) | my own 1 mK monotone scan of `studies/p2_15_buckle/derive_span.py` `C_L`; my own `T_surf` reduction of `p2_16_span_rule/legs/{CAL,HOS_R}/run_removal_events.csv` | reproduces |
| z2 | h_min(L) and q_net,min | report | 12.0 mm: **h_min 535.53, q_net,min 0.43156 MW/m², ABOVE** CAL's **0.40137**; 13.2 mm: 487.78 / 0.39092, below; 10.30 mm 0.5081; 17.77 mm 0.2940. ⚠ The verdict reduces to **L vs CAL's own L_req 12.837 mm** (6.5 % margin at 12.0 mm) — the interpolation lands at f = 0.099 and f = 0.979, i.e. next to the CAL leg both times | my own log-log interpolation and my own hand evaluation of the face balance over Z2's L_req / T_s table | reproduces |
| a1 | CAL1 vs P2-16 CAL | event log, plotfiles, non-timer thermo identical; `check_shared` PASS 1 declared; metadata 0 | event log **byte-identical, 10,805 lines** (`cmp`); **98 plotfile files, 0 differ** (103 files − metadata/thermo/diff.patch/2 visit); thermo **7,500 rows × 79 cols, 0 of the 77 non-timer columns differ**; metadata differs **only** in `plot_file`, one added line `spall.buckle_span = 0.012`, `HASH` and 3 timestamps — the resolved last value is 0.0 in **both** runs, which is why the bytes match; `check_shared` PASS, 1 declared, 64/64 keys | my own `cmp`, `diff -rq`, column-by-column string diff, `diff` of `metadata` | **reproduces (gate)** |
| a2 | face_heat on the edited files | rc 0, PASS, Σ `bk_n_block` declared | rc 0; my masked **line-by-line** diff against `p2_16_span_rule/UNIT_p2_16/face_heat.log`: **52 lines, 0 differ**; my own Σ over all six running legs: **bk_n_eval = bk_n_block = 0** — and the `bk_*` columns **are present**, so the span did reach the binary and the gate is not vacuous | `unit/face_heat/output/*/run/thermo.dat`, `studies/p2_17_buckle_on/face_heat.{log,rc}` | **reproduces (gate)** |
| a3 | O_S164 vs HOS_R to 15 s | identical | **5,066 rows byte-identical** to HOS_R's whole log under **my own step-index filter** `round(t/dt) < 1875` (not the implementer's `start < 15 − dt/2`); the packet's `time ≤ 15.0` admits **exactly 4** extra rows, all of step 1875 stamped `1.4999999999999234e+01`; 15 s **cell + node** plotfiles identical (14 files). Argv differs from HOS_R in **exactly three tokens** (`plot_file`, `stop_time`, `+spall.buckle_span=0.0`): `h_expr` is **byte-equal**, and 0.1636 × 7272.3958 = 1189.763953 reproduces it | my own filter + `diff -rq` + token-by-token argv comparison | **reproduces (gate)** |
| a4 | O_S140 V(f) 30 s | 24.922 ± 0.001 | **24.9222 cm³** | my own `(size − Σ fr_f)·dx²·dz·4` over `legs/O_S140/run/03750cell`, geometry from that leg's own argv | **reproduces (gate)** |
| a5 | N_S070 reaches B | metadata 0.012, Σ `bk_n_eval` > 0 | metadata's **last** `spall.buckle_span` = **0.012**; Σ `bk_n_eval` = **27,740**; `check_shared` PASS 64/64. My `bk_span` thermo read confirms **every N leg ran 0.012 and every A leg 0.0132**, and O/CAL1 carry **no** `bk_*` columns at all | my own metadata grep, thermo sums and `bk_span` census | **reproduces (gate)** |
| c1 | blocked fraction, split | report | N **0.0070 / 0.0316 / 0.1782 / 0.6086 / 0.7888 / 0.9563 / 0.9930**; A **0.0349 / 0.2816 / 0.4722 / 0.9062 / 0.9853**; **roofs 0 on every leg**; walls the majority at S 0.12–0.087 (1670/2495 at 0.10), floors **85 %** at 0.07. Per firing N/O **1.0028 / 0.9978 / 0.9971 / 0.9906 / 0.9712 / 0.9244 / 0.6734** | my own thermo sums and event-log row counts | reproduces |
| c2 | S\*, the gap, E_src | report | S\*(O) **0.099435**, S\*(N) **0.101354**, S\*(A) **0.100285**; N/O **1.0192**, A/N **0.9894**; gap **10.057× → 9.866×**; E_src ×4 at S\* **24.706 / 25.180 / 24.886 kJ** (N/O 1.019) | my own log-log solve over my own V(f) and the last thermo row of `ledger_E_src` | reproduces |
| c3 | n on the bracket | report | n(O) **1.8321**, n(N) **2.0483**, Δ **+0.216**. ⚠ The implementer's caveat ("part of the +0.214 is curvature") resolves **in their favour**: on **matched** brackets Δ is **+0.149** (0.10–0.12) and **+0.457** (0.087–0.10), both positive, so P6 holds without the bracket mismatch | my own paired fits | reproduces |
| c4 | isotropy spread | report | **0.1067 / 0.1447** at S 0.10, **0.2177 / 0.2551** at S 0.07; N − O **+0.0380 / +0.0374**, both under 0.05 | my own 11.25° sector census (Σ(1−fr_f) over r ≤ 25 mm per sector ÷ that sector's cell count) | reproduces |
| e0 | the cap is dead as a bound | all three | Σ `fr_cap` = 0 **and** event-log bit 16 = 0 on **all 12** N/A legs; controls **O_S087 3, O_S100 2, O_S070 1** bit-16 events (thermo sampled Σ `fr_cap`: 0 and 1); **bit 8 (the 15 mm search depth) = 0 everywhere**; max N/A T_s **1057.9 K** (N_S120) against T_s\* 1272.7 (N) / 1151.9 (A) | my own `rc_flag` bit census over all 19 ladder event logs + my own thermo sums | **reproduces (verdict)** |
| e1 | verdict | — | gates z0, a1–a5 all reproduce; the cap is dead as a bound; P4 **branch (b)**, which I confirm is **robust to the estimator** (see DR-2); P5 immaterial; q_net,min 0.432 MW/m² above h 500 | — | **reproduces** |

⚠ The packet marks **nothing** unquotable (no "unquotable" / "do not quote" string anywhere in it), so there is
nothing to violate. Checked against both `ROADMAP.md` `## Standing rules` and `## Known Stale/Important Notes`: the
notes engaged here are SR 1 (the rock file's comment names the `c_f` pairing itself), SR 3 (V and E_src both ×4, and
the ratio uses both sides), SR 5 (the reducer takes geometry from each leg's own argv — no `LZ`), SR 7(i) (every
`bk_*`/`fr_cap` column is per-step and is read as a ratio or declared "sampled"; `ledger_E_src` is cumulative and read
at a row, not summed), SR 7(iii) (this is exactly what a3 hit), SR 8 (argv + deck + allow + binary hash per leg),
SR 9 (`face_heat_output_pre/`, `inputs_pre/`, `MODEL_CARD.pre_p2_17.md`; no leg was re-run — no `.old*` exists), SR 10
(c4 is per sector), SR 11 (`caffeinate -dis -w $$`), SR 14 (L3 is a source assertion, verified below), SR 20 (L3 names
the μ on both sides). The column-path caveats in the notes are untouched: no `src/` change, and no column-path deck was
edited. SR 18's share is unchanged — I measure **20.7 %** of CAL's and **21.5 %** of O_S164/N_S164's firings above
873.15 K, in line with the rule's own figures.

### Diff review

Scope: the four edited input/script files (diffed against `inputs_pre/`, which `cp -p` preserved — Nov 2 10:03 and
Nov 3 09:33/09:34 mtimes), the model card (against `MODEL_CARD.pre_p2_17.md`), the five study scripts, and the two
source sites the relabel rests on. All four edited files are **untracked**, so `git diff` shows nothing — the
`inputs_pre/` copies are the only diff basis, and their mtimes support the claim that they predate the edits.

| # | hunk / call / test | class | finding | tag |
|---|---|---|---|---|
| 1 | `rock_grimsel.inputs` +16, `spall.buckle_span = 0.012` | 7 duplication/readers | reader at `MMWSpalling.H:9160` (`query_default`, default **0.0**), `bk_on` at `FaceRemoval.H:488`, evaluated under `bk_on && o.fire` at :607. Metadata on all 19 legs resolves it as expected. Duplication: **none** — one declaration, one reader | ok |
| 2 | `model_form_flame.inputs` − the `spall.buckle_span = 0.0` line | 9 removed check | authorised by Goal 2. Behaviour-neutral: the parse default is 0.0, so a deck including only this file is unchanged. Verified that **exactly two live decks** include either shared file (`validation/hoeser/input_lab`, `validation/kant/grimsel_calibration/input`) and **no tree-sweep input does** — L1 holds, so the 2,639-file gate genuinely does not need re-running | ok |
| 3 | `spall.flake_h_max` comment: TUNED → "NUM … DEAD as a flake bound under B" | 2 consistency | **The "DEAD" label is a rock-dependent fact written in the rock-independent file.** Deadness holds because T_s\*(L) = 1272.7 K at the *rock file's* 12.0 mm; at R1's 13.2 mm T_s\* falls to 1151.9 K and at 17.77 mm to 894 K, **below** HOS_R's own 1021.4 K max. A second rock with a shorter span revives the cap, and `model_form_flame.inputs` would then be lying. The condition is stated in the rock file's comment but not in the form file's | **risk** |
| 4 | CAL deck: directive + `spall.buckle_span = 0.0` after the two FILE lines | 2 consistency | AMReX last-wins (`AMReX_ParmParse.cpp:1397`) — and `a1`'s `metadata` is the proof, carrying **both** occurrences in order with 0.0 last, in both CAL1 and the P2-16 reference. `reduce.py`'s own metadata read also takes the **last** match, mirroring the parser | ok |
| 5 | `check_shared.py`: `directives()` + `allow \|= set(declared)` | 1 duplication / 2 consistency | Duplication: **none** — a separate scan is *required* because `parse()` strips comments, and the docstring says so. Consistency: `declared` is merged into `allow` **before** anything downstream, so the directive and `--allow` take the identical path through both (i) and (ii). I verified this by running O_S100 (`--allow`) and the CAL deck (directive) and getting the same resolution behaviour | ok |
| 6 | the directive is read from the **deck file only** | 4 scope vs source | Deliberate and documented: an included file cannot launder a change to itself. The one place it could bite — a future shared-file-level override — would fail loudly on (i) | ok |
| 7 | `run.py` `leg_argv`: token replacement with `assert seen == set(rep)` | 2 consistency | The assert is the right guard: a base argv missing `stop_time`/`plot_file`/`h_expr` aborts instead of silently shipping the base's value. Proven live — O_S164's argv differs from HOS_R in exactly the three intended tokens | ok |
| 8 | `run.py` `FORMS` ↔ `ALLOW.txt` coupling | 7 readers | `["spall.buckle_span"] if FORMS[form] else []` — O and A declare, N does not, which is exactly what `check_shared` needs. Verified by running it myself on N_S100 (0 declared), A_S100 and O_S100 (1 declared), all PASS, 64/64 | ok |
| 9 | `reduce.py` a3: `float(t) < 15.0 - 0.004` | 5 hard-coded | `dt` is in each leg's argv (`timestep=0.008`) and in `metadata`; the comparator hard-codes half of it. A leg at another `dt` would silently shift the cut. My own filter derives the step index from `dt` and lands on the same 5,066 rows | nit |
| 10 | `reduce.py` `E_src`: docstring "at the last thermo row (t = 30 s)" | 11 false comment | The last thermo row is at **29.952 s** (625 rows at `thermo.plot_int = 6`), so c2 pairs a 30.000 s volume with a 29.952 s energy — 0.16 % and identical on both forms, so the N/O ratio is unaffected, but the label is wrong | nit |
| 11 | `reduce.py` a2: `moved = [x for x in nl if x not in rl]` | 10 test strength | One-sided set membership, order-insensitive — it would not see a line **removed** from the reference, and it is not part of a2's `o`. The packet's bar only requires a moved number to be *declared*, so this is in scope; my own masked line-by-line diff gives 0 differing over 52 lines | nit |
| 12 | `reduce.py` `SPAN = {"N": 0.012, "A": 0.0132}` feeding e0's T_s\* | 5 hard-coded | The span is logged — `bk_span` is a thermo column and `metadata` carries the resolved value — yet e0 computes T_s\* from a module constant, so a leg that ran a different span would be scored against the wrong T_s\*. a1/a5 check metadata for CAL1 and N only; the **A legs' 0.0132 is never checked by the script**. I closed it by reading `bk_span` on all 19 legs: every N leg 0.012, every A leg 0.0132 | **risk** |
| 13 | `reduce.py` e1: gate roll-up re-reads `Z0.txt`/`A*.txt` from disk | 10 test strength | A stale `A*.txt` would pass the roll-up silently. Mitigated in fact: `CHAIN.log` records the order, the a3 FAIL, the resume reason, and `A3.first.txt` is kept; `RESULTS.first.txt` → `RESULTS.txt` differs **only** by the added per-firing block, so nothing was re-scored | nit |
| 14 | `offline.py` `t_star` at L = 12.0 mm | 6 boundary / 12 discriminator | `want = 8.000` against `C_L(table end) = 7.9973`: T_s\* = **1272.7 K** is a **0.45 K sliver** at the end of the strain table, above which C_L is flat forever. The honest statement at the declared span is "the cap can never bind", not "above 1272.7 K"; the implementer does flag the 7.997-vs-8.000, and the conclusion is **strengthened**, not weakened. At 10.30 mm `t_star` correctly returns `None` and e0 treats it as dead | **risk** |
| 15 | the L3 conjunction, read in source | 14 / 20 one criterion one kernel | **Verified and correct.** The cap compares along-line `h_f` with `sp_flake_h_max / L.mu` (`FaceRemoval.H:504` through `ReconFromX`, and the fallback at :536–538) ⟺ normal μ·c_f·ℓ vs 1.5 mm; B compares `x = L.mu * sp_flake_c_f * ell_eff` (:617–620) with L/C_L. Same normal quantity, μ cancels. **And it is not circular:** `ReconFromX` caps only `r.h_f` and **never** `r.ell_eff` (`MMWSpalling.H:9020`+), so B sees the **uncapped** flake. `fr_cap` increments only for a fired **non-through** line (:704) and resets every step (:341); `rc_flag` bit 16 is `h_f_bounded` (:893) | ok |
| 16 | e0's evidence scope | 12 discriminator | Σ `fr_cap` / bit 16 count only **removed** flakes, so the claim proven is "no removed flake was capped" — which **is** the claim the relabel makes. A line B blocks is never logged and never counted; that is correctly out of scope. Separately, max firing T_s is partly self-limited (SR 15), but here B pushes it **up** (N_S120 1057.9 K against HOS_R's 1021.4 K), so the T_s margin is conservative, not flattering | ok |
| 17 | N/O = **1.0028** and V +0.28 % at S 0.20 | 12 carried threshold | P2-15 DR-3's "B is a pure withdrawal" was established on a **flat** deck (576/576 survivors identical). On a drilling deck it is false: withdrawals reshape the face and later firings elsewhere **outnumber** them. The implementer says so plainly. The standing claim now needs the qualifier "on a flat deck", or the next packet will reuse it to argue a direction it does not license | **risk** |
| 18 | `unit/face_list/input`, `unit/face_removal/input` | 7 readers | Both are **generated** by `make_input.py`, which flattens `input_lab` **plus the shared files** into a self-contained input. Neither carries `spall.buckle_span` today (both predate it, Oct 27–28), so the takeaway's "only `face_heat/test` reads `input_lab` at run time" is correct — but **regenerating either now switches B on** in two unit tests whose assertions have never seen it, and nothing in the tree gate covers them | **risk** |
| 19 | the RUNNABLE registry | 7 readers | Four of its entries run the shared files (`HOSER input_lab`, `HOSER_47 G1`, `HOSER_0 HG`, `HOSER_0_S s100`) and now resolve span 0.012 — with B **on** — while `KANT_cal` resolves 0. The pass was not re-run (outside Goal 6), and at its `stop_time = 1e-9` B cannot fire, so `RUNNABLE.txt` is almost certainly still accurate. ⚠ Any bare re-run of G1/HG/s100 no longer reproduces its banked numbers; the `O` token is what restores them, and a3 proves it does so byte-for-byte | nit |
| 20 | the production form vs the onset deck | 2 consistency | q_net,min(12.0 mm) = 0.432 MW/m² against the onset deck's 0.401 ⇒ **with B on that deck removes nothing**, so the onset level (850.43 K) and the span are resolved on **mutually exclusive configurations**. This is declared (deck header, README, rock file, model card §9, pre-registered as P2) and is physically coherent — B is what supplies Kant's minimum flux, and the fixture deliberately sits below it. But it is the fact the next packet most needs in front of it | **risk** |
| 21 | MODEL_CARD rows 21 / 28, §9, §11 item 1 | 5 hard-coded / 3 units | Checked against my own numbers: row 21's `L.Lg = Dmax + 2 + ⌈min(cap, h_spall_max)/(μ dx)⌉` matches the source (it omits only the harmless `min(span+1, ·)` clamp), and 19 → 32 cells at μ = 1, dx = 1 mm is right. §11's "2.9 % of firings at S 0.10 and 0.2 % at S 0.1636" matches my event-log counts (2.885 %, 0.220 %). `TARGET_V = 13.0` traces to `hoeser_2016_drilling_series.csv:12` (`load,38,8,7.9,0,30,7.4,26,5.12`): 26 × 30/60 = 13.0 | ok |
| 22 | c4's 0.05 bar | 12 discriminator | Both forms' spread roughly **doubles** from S 0.10 to 0.07 (0.107 → 0.218 for O), so a fixed 0.05 bar on N − O sits on a moving baseline; the measured Δ is flat (+0.038, +0.037) across that doubling, which is the reassuring reading and the one worth keeping | nit |

**Tally: 22 rows — 0 bug, 6 risk, 6 nit, 10 ok.** No row touches a gating row or a standing rule as a `bug`.

### What I credit

- **DR-1 — the retirement of `d` holds and the directive is a real mechanism, not a label.** z0's two sides plus my
  two extra mutations (directive without the override; whitespace-collapsed directive) show the check discriminates in
  both directions, and D2 fails **(ii) as well as (i)** — so the directive is load-bearing for the *resolution*, not
  only the declaration. a1 then proves the whole shared-file edit moved the onset deck by **zero bytes** over 10,805
  event-log lines, 98 plotfile files and 77 thermo columns, with `metadata` carrying both span occurrences and
  resolving the last.
- **DR-2 — branch (b) is robust to the estimator, which the packet did not claim.** The implementer flagged that
  S\*(O) and S\*(N) come from **different** brackets and that the 1.0192 sits 0.0008 under the 1.02 cut. I tested
  whether the branch survives that: taking both forms on 0.10–0.12 gives **1.0190**, both on 0.087–0.10 gives
  **1.0178**, and a 3-point log-log quadratic through each form gives **1.0186**. Every estimator lands **below 1.02**.
  So branch (b) is the right call, and the same check turns the c3 caveat around — n(N) > n(O) on **matched** brackets
  too (+0.149 and +0.457).
- **DR-3 — the one honest negative.** `N/O = 1.0028` at S 0.20: B **adds** removals. The implementer could have
  rounded it into the noise and instead named it as the sign that P2-15's "pure withdrawal" does not survive a
  drilling deck. That is the single most useful line in the takeaways.
- **DR-4 — the counter's meaning is pinned before any bar is read off it.** `bk_n_eval` is per **evaluation**, so
  79 % blocked at S 0.10 coexists with 2.9 % of firings lost. The notes say which measure each number is, and c1
  reports both. P4's blocked-fraction bands are scored **MISSED** rather than quietly re-interpreted.
- **DR-5 — the a3 stop was taken rather than argued around.** `CHAIN.log` shows the chain stopping at a3, the
  restatement, the resume point and that O_S164 was **not** re-run; `A3.first.txt` keeps the failed scoring. My own
  step-index filter independently lands on the same 5,066 rows and the same 4 excluded rows.
- **The source argument for the relabel is the right kind of argument** (SR 14/20): a stated conjunction on one normal
  quantity with μ named on both sides, which I could check at the three sites and which turns out to be non-circular
  because `ReconFromX` caps `h_f` but never `ell_eff`.

### For the planner

- **DR-6 (row 20) is the carry-forward.** The production form cannot reproduce the deck the onset was calibrated on.
  It is declared everywhere it should be, but the next packet should state which *measured* minimum flux 0.432 MW/m²
  is being held to, since that is now the production form's sharpest prediction.
- **Row 3** — move the cap's "DEAD" claim, or its condition, out of the rock-independent file.
- **Row 17** — qualify "B is a pure withdrawal" to flat decks before it is reused.
- **Row 18** — `unit/face_list` and `unit/face_removal` are generated from the shared files and are now stale with
  respect to them.
- **Row 12** — e0 should read `bk_span` instead of a module constant; the A legs' span is otherwise unchecked.

## Verdict: accepted

- **Gating rows.** z0, a1, a2, a3, a4, a5 all reproduce from my own reductions (table above). No diff-review row is
  tagged `bug`.
- **Guardrails.** No `src/` change (binary sha256 unchanged from P2-16's recorded `38c86238…`, and nothing under
  `src/` is newer than it); `spall.flake_h_max` still `0.0015` with only its comment changed; S appears in leg argv
  only and in no input; `input_lab` untouched; `inputs_pre/` predates the edits by mtime; 15 s and 30 s plotfiles
  present on every ladder leg; nothing in `ext/`, `bin/`, `obj/`, `build/`, `compile_commands.json`, `configure` or
  `LICENSE` touched; nothing committed.
- **Reruns.** None. No test, simulation or build was re-run — every number above came from artefacts already on disk.
  The only things I executed were `check_shared.py` on five scratch copies of the CAL deck and my own reducers, all
  zero-compute.
- **Files inspected.** `ACTIVE_STEP.md`; `ROADMAP.md` (`## Standing rules`, `## Known Stale/Important Notes`);
  `validation/shared/{rock_grimsel.inputs, model_form_flame.inputs, check_shared.py}` against
  `studies/p2_17_buckle_on/inputs_pre/`; `validation/kant/grimsel_calibration/{input, README.md}`;
  `docs/project/MODEL_CARD.md` against `MODEL_CARD.pre_p2_17.md`; `studies/p2_17_buckle_on/{run.py, offline.py,
  chk_directive.py, reduce.py, chain.sh, chain_resume.sh, CHAIN.log, Z0–Z2, A1–A5, A3.first, RESULTS.first, RESULTS}`;
  all 20 `legs/*/{argv.sh, ALLOW.txt, DECK.txt, BINARY.txt, rc, .done, run/metadata, run/thermo.dat,
  run_removal_events.csv, run/0{1875,3750}cell}`; `p2_16_span_rule/legs/{CAL,HOS_R}` and `UNIT_p2_16/face_heat.log`;
  `p2_15_buckle/derive_span.py`; `p2_8_face_list/reduce_p2_8.py` and `p2_11_heat/reduce_p2_11.py` (the imported
  instruments); `unit/face_heat/output/*/run/thermo.dat`; `unit/face_{list,removal}/{input, make_input.py}`;
  `docs/project/model_card/runnable.py`; `src/Integrator/MMWSpalling/FaceRemoval.H` (:336–384, :484–545, :600–632,
  :698–712, :876–893) and `src/Integrator/MMWSpalling.H` (:947–955, :9020+, :9160–9163).
- The user should run `/plan` next; the planner will archive this step.
