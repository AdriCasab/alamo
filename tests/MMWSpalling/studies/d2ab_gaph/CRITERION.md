# D2ab — scoring criterion, FIXED BEFORE THE FIRST LEG

Hashed with **zero outputs on disk**. `PREFLIGHT.md` (sha
`9840fa2bd373d7b8674b6b6ac2e7dd18697e9cb03104a77dc22b1c3eae35f0fd`) and
`preflight.py` (`322791fd988d9e8c35d0d3267f074f4d7e0458e3aacf10170d9e8041c30cc9e6`)
were frozen first; every number below that came from a derivation came from there.

**No re-windowing and no re-choosing a metric after a result is seen.**

## Legs

| leg | `side_face_h` | `side_face_T` | `side_face_h_r_min` | base |
|---|---|---|---|---|
| `GL`  | 38.4  | live | 0.040 | SW |
| `GH`  | 169.3 | live | 0.040 | SW |
| `WLr` | 81.0  | live | 0.040 | SW |
| `WLv` | 81.0  | live | 0.0   | SW |
| `SWc` | —     | —    | —     | SW (bridge/identity) |
| `S0`  | —     | —    | —     | RB_D (identity) |

`RT` is **DROPPED**. The packet makes it conditional: "either item 0 names a real
source with flux-rate or flux-temperature pairs, or `RT` is dropped". Item 0
found none — the repo holds Rauenzahn's **surface-temperature** compilation
(Fig. 2.12, `meier_thesis_notes.md:123`), against which the model's 528 K already
sits inside the band (`:222`), and every other figure on record is **fuel-basis**.
A prescribed-flux flat face returns `q / (rho c_p dT_spall)` by construction, so
it would have tested nothing the constant 1.32 J/mm³ has not already shown.
Carried as a one-line citation in `RESULTS.md`, not as a leg.

## Window and rate definition

Rate is the **burner DESCENT slope over the 100–150 mm stretch**, the
depth-based definition imported unchanged from `d2x_areasplit/analyze.py`.
Window-independent by design. Power/temperature means are over **150–400 s**.

---

## (a) SCORED — the derived bracket against its OWN pre-registered prediction

From `PREFLIGHT.md` §4:

| leg | h [W/m²K] | predicted rate [m/h] | 15 % bar |
|---|---|---|---|
| `GL` | 38.4  | **2.87** | 2.44 – 3.30 |
| `GH` | 169.3 | **3.16** | 2.69 – 3.64 |

**Bar: measured rate within 15 % of the predicted rate.** Ø@100, Ø@125, pit lead
and volume are **reported** against the §4 interpolations, not scored.

⚠ This scores **the derivation, not Meier**. A leg outside Meier's 1.3–1.6 m/h
band that matches its own prediction is a **PASS**.

⚠ `PREFLIGHT.md` records a **disagreement with the packet's P2** (2.3–2.9 m/h)
and the reason for it. Scoring is against item 0, as the packet's (a) specifies;
P2 is reported alongside and scored separately as a prediction.

## (b) SCORED — the attribution pair, one edge each

`WL81` (frozen: const T 1463.5, all faces) → `WLv` (live T, all faces) → `WLr`
(live T, wall only). Report rate, side power and the `side_P_{bore,pads,wall}`
split for all three.

* **Bar 1: `WLv` rate in 1.85–2.00 m/h.**
* **Bar 2: `WLr` rate strictly above `WLv`.**

⚠ `WL81` carries no `side_face_h_r_min`, so its split is over ALL faces; `WLv`
reproduces that population with only the temperature changed, which is what
makes the edge single.

## (c) SCORED — the face assignment is demonstrated, not asserted

⚠ **The packet's literal form is unscoreable and this is the substitute, fixed
here before launch.** The packet asks that `side_P_bore` and `side_P_pads` be
"unchanged from the inherited-map legs to within 5 %". The only all-inherited
leg is `SW`, and **`SW` was frozen before the region columns existed** — checked
on disk: 0 of `side_P_{bore,pads,wall}` present in its thermo header. That is
D2aa's open item-0(a) debt. `SWc` is run to supply the missing baseline, and (c)
is scored in three parts:

* **(c1) against `SWc`, the packet's own test.** On `GL`, `GH` and `WLr`,
  `side_P_bore` and `side_P_pads` within **5 %** of `SWc`'s.
* **(c2) the gate is insensitive to `h`.** Across `GL`/`GH`/`WLr` the *same*
  two bands must agree within **5 %** of each other, although `h` spans
  38.4 → 169.3, a factor of **4.4**. Only `side_P_wall` may move.
* **(c3) the positive control — the gate is not a no-op.** `WLr` vs `WLv` at the
  **same** `h = 81` and the same live T: `side_P_bore` and `side_P_pads` must
  **differ by more than 5 %**. If they do not, the radius gate did nothing and
  (c1)/(c2) would pass vacuously.

## (d) SCORED — the power identity

`|jet_P_face − (tops + sides + side_loss)| ≤ 0.02 W` on **every** Meier leg.
This is the one check that catches a kernel/march/ledger split, which is why the
live temperature is resolved once per step in `BuildFlameColumns` and all three
sites read that member.

## (e) SCORED — identity

* `S0` vs `d2z1_delete/output/RB_D`: `_removal_events.csv` and `_jet_profile.csv`
  byte-identical, **thermo header delta EMPTY**, 0 columns differing.
  ⚠ `metadata` and `diff.patch` are run provenance (compile date, plot_file path,
  the ParmParse table, a start stamp) and are excluded, the same rule the
  2639-file tree list applies.
* `SWc` vs frozen `d2z2_multicell/output/SW`: `_removal_events.csv` byte-identical
  and every **shared** thermo column identical. ⚠ `SWc` legitimately gains no new
  column here — `side_T_wall` is registered only when `side_face_h` is given, and
  `SWc` gives none — so the header delta must also be EMPTY.
* **Tree gate 2639/2639** against `witness/d2z1.md5`.

## (f) REPORTED — Meier, residuals stated not hidden

Rate **profile by depth in 50 mm stretches** for every leg (deliverable 3: the
rate is not steady — `A` falls 2.35 → 1.92 → 1.70 m/h across 100–200/200–300/
300–400 s, and a single number hides that). Ø(z), volume, ring table per ring,
never averaged.

⚠ **Meier Fig. 8.8 is a visible width = a LOWER BOUND on diameter**, band
85–95 ± 5 mm. **No agreement claim tighter than ± 5 mm.** Meier's rate is
**one-sided** (operator advance).

---

## Pre-registered predictions

* **P1** `GL`/`GH` within 15 % of `PREFLIGHT.md` §4 (2.87 / 3.16 m/h).
* **P2** (the packet's) `GL`/`GH` land 2.3–2.9 m/h. ⚠ Item 0 derived 2.87–3.16
  instead and said so before launch; both are scored.
* **P3** `WLv` 1.85–2.00 m/h; `WLr` above `WLv`.
* **P4** the derived bracket **contains** 81, so `GL`/`GH` bracket `WLr`.
  Already true on paper at pre-flight (38.4 < 81 < 169.3); the leg test is that
  the **rates** bracket too.
* **P5** the charged riser area keeps **growing** as `h` falls (68 → 133 → 143 cm²
  on the D2aa legs). Refuted if it saturates.
* **P6** `jet_T_exhaust` is **steady at the thermo cadence, not oscillating**,
  on every live leg. The feedback is real (less heat to the wall ⇒ hotter gas ⇒
  more heat to the wall). **Refuted ⇒ the loop needs damping before any leg is
  scored.** Scored as: the 150–400 s standard deviation of `jet_T_exhaust` is
  below 2 % of its mean, and no growing oscillation in the trace.
* **P7** **no leg reaches Meier's 1.3–1.6 m/h band.** Held ⇒ the residual is the
  bore and pad-ring lateral faces and the successor is **lateral recession**.

## Guardrails in force

* `h` came from geometry; the rate is the prediction. **No point value inside the
  bracket.** Nothing here may be re-chosen to land a rate.
* `spall.per_cell_removal` OFF on every leg. `side_face_h_max` never set.
* `142` is the annulus derivation's D2aa value and is not varied here.
* The pad ring stays the **open-floor** map (`PD` jammed for an architectural
  reason: undercut is impossible in this code). Stated as a limitation.
* Frozen studies are read-only, imported by path: `WL`, `WL81`, `PD`, `BO`, `SW`,
  `SWs`, `SW0`, `A`, `D`, `RB_D`.
