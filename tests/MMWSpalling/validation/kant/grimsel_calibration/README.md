# Kant 2017 — the Grimsel calibration deck (P2-13)

**What it is.** The one place the Grimsel rock's onset level is calibrated, in the **production form**: the same
`validation/shared/rock_grimsel.inputs` and `model_form_flame.inputs` every flame deck includes, on today's removal
rule (the face list). Same rock ⇒ same properties, calibrated once (`docs/project/MODEL_CARD.md` §1, §12).

**The deck.** A 24 × 24 × 64 mm Grimsel block at 1 mm, 20 °C, Neumann on every face, under a uniform flame over its
whole top: T_flame 1673.15 K (Kant 2017's ≈ 1400 °C), h 500 W/m²K (Kant does not report h), `jet_closure none`. Every
column is the same 1-D onset experiment. `weibull.V0` = V_cell (moves with the mesh).

**An onset-only deck (P2-16, 2026-11-03).** The deck calibrates the onset rule A with the buckling condition B
(`spall.buckle_span`) **off by definition**, and nothing else. B is a pure withdrawal (a firing that passes it is
identical to the span-off firing, P2-15 DR-3, 576/576), so the onset measured here **is** the production onset
wherever production fires. The flame's h 500 W/m²K is a **convention** (Kant 2016's minimum coefficient, chosen in
P2-13): Kant 2017 states no standoff, coefficient or flux, so the deck carries no evidence about B or the minimum flux,
and a span that would block it is not in conflict with it (P2-15 DR-1 closed). With B off the onset is flat in flux
(P2-14a c1: 847.7 / 850.4 / 851.4 / 853.0 K at h 150 / 500 / 1000 / 2000), so the calibration does not depend on the
convention.

**The declared override (P2-17).** Since P2-17 the rock file carries the buckling span (`spall.buckle_span = 0.012`).
The deck sets `spall.buckle_span = 0.0` directly after its two `FILE =` lines (AMReX reads the last occurrence) and
declares it with the directive `# check_shared: allow spall.buckle_span`, so `check_shared.py` accepts exactly that one
override and still refuses any other.

**The statistic.** Per column, the first event-log row with `dir = 5` (top face) and `ts_fb = 0`; its `T_surf` (the
reconstructed surface temperature the criterion fired on). The level is the median over columns.

**The decision rule (c1, d1).**
- The **Grimsel onset band** is the overlap of Kant 2017's measured onset at 0 MPa (553–694 °C) and Höser's IR onset
  on Grimsel (430–600 °C): **553–600 °C = 826.15–873.15 K**.
- CAL's median **inside** the band ⇒ the banked constants (R 0.531294, flaw scale 8 µm, m 13) **stand**; nothing moves.
- **Outside** ⇒ the flaw scale is recalibrated **once**, here: RECAL_6 and RECAL_10 (a0 6 / 10 µm) must bracket 850 K;
  a0 is log-log interpolated to 850 K (4 significant figures), RECAL_x must land within ±3 K, and the value is written
  to `rock_grimsel.inputs` (both flaw keys). G1's S\* must then be re-derived. If the bracket fails, stop: one flaw
  scale cannot reach the level.
- ⚠ The band is a target for the **level** only. The flux (c4), confinement (c5) and mesh (c6) rows are reports and
  are never tuned to.

**Files.**
- `input` — the deck. `run.py LEG ...` — the legs (serial, `caffeinate`, argv + binary per leg in `legs/<leg>/`);
  `run.py --argv LEG` prints a leg's command; `run.py UNITS` runs a4.
- `offline_onset.py` — the offline twin (uniform stress, Tada weight, `Grimsel.H` via `KIC_CPP.txt`) → `A0.txt`.
- `reduce.py [A1|A2|A3|B1|D1]` — the rows; no argument writes `RESULTS.txt`.
- `chain_ord2.sh` — the packet's run order from the gates to the results. `pre/` — `input_lab` before P2-13.

**Run** from the repo root (the `FILE =` paths are relative to the cwd), on 4 ranks:
```bash
PY=/Users/tzetze20/Desktop/code/.venv/bin/python; K=tests/MMWSpalling/validation/kant/grimsel_calibration
$PY $K/run.py CAL && $PY tests/MMWSpalling/validation/shared/check_shared.py $K/input $K/legs/CAL/argv.sh
```
