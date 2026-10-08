# D2v criteria — fixed and hashed BEFORE the first campaign leg launched

## 0. Disclosure

**No campaign leg (`M0`, `M1`, `M2`) has been launched when this file is
hashed.** What has run: the two ~10 s trace legs `TF`/`TP` of item 3, whose own
criteria were hashed separately in `HMAP_TRACE.md` (`90d84f2f8b2172c7…`) before
any temperature was read. They **passed T1–T10**: the compiled map is the
intended one on **both** branches (floor inverts to 1407.5 against a known 1449;
wall reads **142.00 flat** where a blend would read 559–695), and `max dT` beyond
`r = 76 mm` is **0.013 K** over 2463 columns. **N1 is excluded before launch.**

## 1. The metric — ONE definition, applied identically to every leg

Final `k_top` per column from `_removal_events.csv` (`k_top = k_top(0) −
Σ n_voided`; columns with no event keep the initial top) ⇒
`z_face = (k_top + 1)·dz`. Binned on **fixed physical 2 mm radial bins**
(asserted in the analyzer; a mesh-dependent binning would make the metric
self-referential). `z(r)` = **mean** face height per bin, population always
reported. **`Ø(depth)` = twice the crossing radius, LINEARLY INTERPOLATED
between bracketing bins** — D2t's `2·r_centre` snapped to a 4 mm lattice, and
4 mm is 4.7 % of an 85 mm hole and equals the reference's own ±5 mm.

**Scored at t = 500 s on every leg.** ⚠ The frozen table was measured at 572 s;
at 500 s its 200 mm entry reads 80.4 mm rather than 86.6 because the hole was
still deepening. Shallow depths are unaffected. All cross-leg comparisons are
at 500 s.

⚠ **My `M0` reference, measured with this metric:** Ø = **144.3 / 129.1 / 114.3
/ 105.3 / 99.5 / 92.6 / 86.6 mm** at 10 / 25 / 50 / 75 / 100 / 150 / 200 mm
(572 s), against the packet's 141 / 128 / 113 / 105 / 99 / 92 / 86. Five
binning conventions failed to reproduce the packet's values; see `PREFLIGHT.md`
§1. **Criterion (a)'s `M0` is therefore 37.2 mm, not 35.3.** The bar is
unchanged and the comparison is internally consistent because one metric is used
for all legs.

The reference is Meier's **visible width** — a **LOWER BOUND** on Ø at **±5 mm**.
Magnitude is scored against the **85–95 ± 5 mm band** (thesis D = 85 mm;
3.42 L / 0.5 m implies a mean of 93 mm); the csv is used for **SHAPE only**.
**Never quote agreement tighter than ±5 mm.**

## 2. Criteria

* **(a) THE SCORED ONE — THE FUNNEL.** Mean `|Ø_model − Ø_Meier|` over depths
  **10 / 25 / 50 mm** on `M1`. **Bar: <= 20 mm.** `M0` measures **37.2 mm**
  (excess 48 / 37 / 26). Report the whole ladder and `M2` too.
* **(b) THE WALL.** `M2 − M1` on Ø at **100 / 150 / 200 mm**, and against the
  85–95 ± 5 band. **Reported, not gated** — one mesh cannot settle it.
* **(c) ITEM 0's WITNESS.** `side_P_face_above` **> 0 on `M2`, exactly 0 on
  `M1`**. **Zero on `M2` ⇒ inert ⇒ STOP and report.** ⚠ On `M0` the column does
  not exist at all (the diagnostic is registered only when
  `side_face_factor_above >= 0`, and `M0` sets no side-face key), so `M0` is
  recorded as **absent**, not 0 — the same deviation D2r-1 recorded.
* **(d) NO BLANKET — THIS ONE CAN KILL THE PACKET.** **Bar: mean recession
  beyond `r = 76 mm` <= 2 mm on EVERY leg.** Score the **mean**, not the max:
  the frozen baseline itself has a max of 6.00 mm (three single-cell events in
  one column) against a mean of 0.216 mm, so a max-based bar would fail on the
  control. Report the max and the beyond-90 mm figure alongside, plus the count
  of above-plane columns with `h > 10 W/m²K` (must be confined to the annulus
  band), plus what each leg strips from the r = 50–76 mm ring (`M0`: 32.5 mm).
* **(e) VOLUME.** Total excavated × 4 against Meier's 3.42 L / 0.5 m pro-rated
  to the depth reached. **Reported.** D2l found the model removes ~1.9× too much.
* **(f) THE LATCH.** Above-plane column count and plane crossings per 100 s
  block. **Reported.**

## 3. Pre-registration, with numbers

* **P1 — `M0` reproduces the frozen `d2i_gate/LI0_2mm` event log BYTE FOR BYTE
  to 500 s.** **Any difference ⇒ STOP**: the baseline this packet attributes
  against would not be what it claims. ⚠ Byte identity to 500 s does **not**
  reproduce the 572 s table's 200 mm entry — the hole deepens after 500 s. That
  is a statement about the *time*, not about identity.
* **P2 — `M1` reduces the collar.** Direction only. **Refuted if Ø @ 10 mm on
  `M1` exceeds `M0`'s.**
* **P3 — `M2` is wider than `M1` at 100–200 mm depth.** **Refuted if
  `M2 <= M1` there.**
* **P4 — (d) HOLDS: no blanket.** **The one that kills the packet.** Refuted if
  any leg strips more than 2 mm mean beyond r = 76 mm. ⚠ Note the mechanism is
  real and unbounded: the march skips above-plane columns, so `M2`'s above-plane
  heat is **not** capped by the jet's enthalpy (`PREFLIGHT.md` §0).
* **P5 — (a)'s MAGNITUDE IS DECLARED OPEN.** No number predicted. Declaring it
  open is the discipline D2t honoured and D2s failed.
* **P6 — THE NUMERICAL FLOOR APPLIES AT SHALLOW DEPTH ONLY.** The 16 mm w1
  blend hands over at a centre of 48 mm, so this configuration **cannot cut
  narrower than Ø ≈ 96–100 mm at 2 mm** where the mouth taper and early
  formation set Ø, i.e. **10–50 mm depth on `M1`/`M2`**. **Refuted if `M1` or
  `M2` cuts below 92 mm at 10–50 mm depth.** ⚠ It does **not** hold at depth and
  is **not** scored there: `M0` reads 86.6 mm at 200 mm with no blend at all, so
  a blanket "below 92 at depth" clause would fire on the control. **Deep Ø is
  REPORTED, not gated.**

**A pass on (a) is "the funnel is removed down to the numerical floor", NOT
"the model reproduces Meier's diameter".**

## 4. Guardrails binding on the scoring

* **Score shape and volume, never rate** (D2l: rate agreement was a coincidence
  of two errors).
* `h_ann = 142`, `side_face_factor` 0.2 / 1.0, and `jet_T_above = 1690` are
  **declared and frozen**; the taper radii come from geometry (`foot_r_outer`,
  the hole wall, the measured mouth ring). **Nothing is adjusted to move Ø(z).**
* **Both blends are 16 mm = 8 coarse cells.** Taper width is a **numerical**
  parameter, narrowed only with the mesh, **never to move a result**.
* **DOMAIN WATCH:** assert the deepest column stays **>= 20 mm** above the
  domain floor; the depth result is void if not.
* **NO 1 mm RUN, NO CONVERGENCE CLAIM, NO SHAPE VERDICT FROM 2 mm ALONE.** D2t
  showed the side-face lever's effect appears as a *mesh gap*, so a single mesh
  cannot see it. This is a **screening** run.
* **No Meier validation claim.** The reference is a lower bound at ±5 mm.
