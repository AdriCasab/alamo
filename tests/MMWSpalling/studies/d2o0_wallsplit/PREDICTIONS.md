# D2o-0 — predictions and decision rule
# (AMENDED 14:53:01, before any result was seen; never edit again)

## AMENDMENT RECORD (read first)

**First hash** `9bbf85bc3e1683916787d548d3e1a73e17d56b811f1851241d16a5541854735f`, 14:47:16.
Amended because **my P2 definition did not reproduce the packet's own control reference**, found
by checking it against that reference — not by looking at any result. The four runs were killed
~2 minutes in (t ≈ 2–75 s of 450), `output/` was deleted, **no `.done` was written and no D2o-0
output was inspected**; the killed log is kept as `campaign_killed_defn.log`.

**What was wrong.** §5 defined "radial gain in the band" as `(min Ø below 10 cm − 80)/2`. On the
control that gives **−0.5 mm/side**, while ACTIVE_STEP states the control is **≈ 33 mm/side**. The
packet's number is the **mouth** gain: the control's Ø at z = 10 mm is 144.9, so
(144.9 − 80)/2 = **+32.4** — an exact match. The corrected definitions are in §5, and all three
gain measures are now reported with the control's value for each, so nothing is hidden.

**Two things that checking the definition exposed, both stated before any result:**
- **P1 on the shaft mean is already satisfied by the control.** Control shaft (z ≥ 100 mm):
  min 79.0, **mean 90.0**, max 98.9 — mean inside 85–95. So "shaft Ø 85–95" only discriminates if
  it is required of the **whole shaft profile**, not its mean.
- **P2 at the mouth is largely a start-up comparison.** At t = 0 the nozzle sits 50 mm above a flat
  surface, so **s = 50 mm > 45 everywhere and the entire surface is briefly in the FLOOR zone**,
  carving a wide shallow crater before the hole forms. Above the nozzle plane the control is not
  heated either (its `s ≤ 0` rule), so the control's 134 mm at z = 20 mm falling to 82 mm at
  z = 220 mm is mostly that transient. **The control's *steady* band gain — fully-passed rock
  (Ø 89 at z = 180) against freshly cut rock (Ø 79 at the bottom) — is +5.0 mm/side, which is
  already inside P2's 2.5–7.5 target.** Both measures are therefore reported, and the steady one is
  the physically meaningful one.

---


Key-only. No source change, no build, no `make`. Binary `bin/mmwspalling-3d-g++` =
`62dea451e6b1d1908e0ea12f730cdc3f7c01790f19d6d63b8c2b88befe0f67fa`, to be unchanged at the end.

**The question.** D2m showed the model drills a **cone with no edge** (v/v_c = 0.62 at 39 mm, never
reaching 0.1·v_c inside 70 mm) while Meier drilled a **cylinder**. This packet prescribes the
floor/wall split by hand — impingement on the floor, grazing on the wall — and asks whether a
cylindrical hole appears, **before** any closure is written.

## 1. Three numbers the packet does not state, resolved from the named control

ACTIVE_STEP specifies the wall coefficient but not its gas temperature, and `jet_closure = none`
removes the march that supplied *any* gas temperature. Each is taken from the named control
(`d2i_gate/output/LI0_2mm`), not invented:

1. **Wall gas temperature = 1610.1 K** — the control's own `jet_T_mouth`, mean over 150–450 s.
2. **Floor gas temperature = the control's own radial profile**, quadratic fit to its
   `_jet_profile.csv` T_gas(r) averaged over 150–450 s:
   `T = 1742.48 − 435.739 r − 42654 r²` (r in m), reproducing it to **3.4 K** over the floor zone.
   Clamped below at 1400 K so the unused branch can never go ≤ 0.
3. **The `s ≤ 0` rule must be written by hand.** Under `jet_closure = enthalpy` the code zeroes
   flux above the nozzle plane itself; under `none` it does **not**, and it strictly requires
   h > 0 and T_gas > 0 at every in-patch column (`MMWSpalling.H:1690`). Without the guard the whole
   undisturbed surface out to the 0.2 m patch radius would be heated. Emulated as **h = 1e-3**, the
   `d2j0b_plane` convention.

**Pre-flight (`run.py --list`, run before this file was hashed):** no empty zone, no silent zero.
At t = 150 / 250 / 350 s the floor zone holds 332 / 325 / 318 columns at r = 1–41 mm with
h = 668–1040 W/m²K, and the wall zone 598 / 293 / 171 columns at r = 40–69 mm with h = 40 exactly.
**The s = 45 mm corner lands at r ≈ 40 mm — the skirt radius — which is the packet's geometric
claim, confirmed on the control's own geometry.**

## 2. The correlation, named and justified in advance

Annulus D_h = 2 × gap = 0.020 m, u ≈ 23 m/s, exhaust at ≈ 1600 K (ν ≈ 2.4e-4 m²/s, k ≈ 0.10 W/mK,
Pr ≈ 0.72) ⇒ **Re = u·D_h/ν ≈ 1.9e3**, below the 2300 transition and thermally developing over the
whole half metre.

| branch | Nu | h [W/m²K] |
|---|---|---|
| laminar, fully developed, one side heated | ≈ 4.0 | **20** |
| laminar, thermally developing (Hausen, D/L = 0.04) | ≈ 6.1 | 30 |
| turbulent, Dittus–Boelter 0.023 Re^0.8 Pr^0.4 | ≈ 9.2 | 46 |

**ACTIVE_STEP states 40 W/m²K as the duct estimate; that is the PRIMARY and is used as given, not
re-derived.** The pre-registered sensitivity `Wsens_2mm` is **the other branch, laminar fully
developed at 20** — a factor 2, the spread the packet names. **Neither value is ever tuned.**

## 3. The packet's targets

- **P1** — shaft Ø **85–95 mm** below 10 cm depth (Table 8.2 gives 85; the rim reading of 95 is
  biased outward and is not quoted alone).
- **P2** — radial gain within the band **2.5–7.5 mm per side**, against a control of ≈ 33 mm/side.
- **P3** — end-of-run wall temperature **600–800 K** (altered, not spalled).
- **P4** — wall heat **1.5–3.5 kW**, integrated over the prescribed wall field. With
  `jet_closure = none` there is no march, so it is computed in `analyze.py` from the field and the
  wall area — stated here so it is not mistaken for a code diagnostic.
- **P5, conditional** — if the wall lands **within ~50 K of the 821 K threshold**, expect mesh
  divergence and treat any single-mesh result as **void**.

## 4. What I actually expect — and it is that P1–P4 FAIL, low

These are derived from the field, not fitted to the targets. **If they are right the packet's
decision rule sends us to "report what a prescribed field could not produce", and the field must
not be adjusted to rescue it.**

**The wall never reaches the spallation threshold.** A semi-infinite solid under a step Robin
condition has T_s = T_0 + (T_g − T_0)·[1 − exp(β²)erfc(β)], β = h√(αt)/k, with
α = k/ρc_p = 6.90e-7 m²/s and k = 1.5 W/mK. A column's residence in the 45 mm band at the control's
1.491 m/h is **109 s**, giving √(αt) = 8.7 mm, **β = 0.231**, and

    T_wall ≈ 293 + 1317 × 0.216 ≈ **578 K.**

- **Q1 — the wall does not spall.** 578 K is **243 K below** the 821 K firing temperature, so the
  wall cannot recede at all. **P5 does not fire** (it needs within 50 K), and the knife-edge
  divergence it anticipates should be absent.
- **Q2 — P2: mouth gain falls far below the control's +32.4, to roughly +5 to +15 mm/side**, carried almost entirely by the flat-start transient rather than by wall erosion; **steady band gain 0–2 mm/side**, below the control's +5.0. Not exactly zero only because
  of a start-up transient: at t = 0 the nozzle is 50 mm above a flat surface, so **s = 50 mm > 45
  everywhere and the whole surface is briefly in the FLOOR zone** — about **12 s** at 0.414 mm/s
  before s drops through 45 mm. At floor h the top cell heats at ≈ 200 K/s, so it fires within
  ≈ 3 s. That transient widens the top of the hole and **must not be read as wall erosion**; it is
  why P2 and P1 are scored **below 10 cm depth**, where no rock ever saw the flat-start transient.
- **Q3 — P1 fails low: shaft Ø ≈ 80 mm**, i.e. exactly the stance (D2m: min Ø = 2·`foot_r_outer` to
  within 2 mm), not 85–95.
- **Q4 — P3 is marginal and probably just below: 550–650 K** at the band, against a target of
  600–800. Note 578 K straddles the ≈ 570 K alteration threshold, so **"altered but not spalled" is
  reproduced even though the number may miss the band's lower edge** — these are different claims
  and are scored separately.
- **Q5 — P4 fails low: wall heat ≈ 0.4–0.9 kW**, not 1.5–3.5. At T_w ≈ 600 K,
  q = 40 × (1610 − 600) ≈ 40 kW/m²; the band's true (conical) area from r = 40 → 69 mm over a
  45 mm height is ≈ 0.018 m² for the full hole, giving ≈ 0.7 kW. Reaching 1.5 kW would need either
  ~2.5× the area or ~2.5× the coefficient.
- **Q6 — the cone does not become a cylinder.** With the wall inert, v(r) should fall off *faster*
  than the control, not flatter: the hole becomes a **narrower** cone of the stance diameter, not a
  cylinder. **The diagnostic quantity is the shape of v(r), not the rate.**

**If Q1–Q6 hold, the honest conclusion is that a single grazing coefficient on the wall cannot
produce Meier's cylinder** — it produces the stance-width hole D2m already measured. That bounds
any closure built on the same idea, which is exactly what this test is for.

**What would refute me:** any radial gain above 2.5 mm/side below 10 cm depth, a wall above 700 K,
or a flatter v(r) than the control.

## 5. Scoring — fixed here, before any result

- **Ø(z)**: 2 × azimuth-mean `r_wall`, 9 sectors (the `d2c_steady/score.py` definition reused
  through `d2l_scan/analyze.py`), reported as a **full profile**, never one depth, and compared at
  **matched feet depth**, never matched time.
- **Radial gain, three measures, all reported, with the control's value for each so the
  definition is auditable:**
  - **mouth gain** = (Ø at z = 10 mm − 80)/2. **Control +32.4**, which is ACTIVE_STEP's "≈ 33
    mm/side", so this is the packet's measure. It is dominated by the flat-start transient.
  - **steady band gain** = (Ø of fully-passed rock, at z = feet depth − 55 mm − Ø at the bottom
    of the profile) / 2. **Control +5.0** — already inside P2's 2.5–7.5 band. This is the
    physically meaningful measure of what the band does.
  - **shaft gain** = (mean Ø over z ≥ 10 cm − 80)/2. **Control +5.0.**
  **P2 is scored on the packet's mouth measure, with the steady measure reported beside it and
  given equal weight in the discussion.**
- **v(r)**: D2m's estimators — 2 mm annuli to 70 mm, from the removal log, **reported absolute,
  never normalised** (three metric traps so far, all from normalisers that moved).
- **Wall temperature**: surface temperature of columns in the band at end of run, by depth.
- **Wall heat**: Σ over band columns of h·(T_gas − T_s)·dA × 4, from the prescribed field.
- **Mesh**: 2 mm / 1 mm gap by 100 s block from 150 s. **A single mesh is not quotable.**
- Nothing re-windowed, re-binned or re-chosen after results are seen.

## 6. The decision

| outcome | reading |
|---|---|
| **P1 and P2 hold** | the floor/wall split is the mechanism → write the D2o closure, gas-path following, annulus correlation as named here, **no free wall parameter**; P3 and P4 become its acceptance targets |
| **P1 or P2 fails** | **report what a prescribed field could not produce**, which bounds any closure built on the same idea. **Do not tune the field to make it pass.** |
| wall within 50 K of 821 K | **P5 fires** — expect mesh divergence, treat single-mesh results as void |
