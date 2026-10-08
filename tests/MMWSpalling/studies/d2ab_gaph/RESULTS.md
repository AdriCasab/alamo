# D2ab — the gap-derived wall coefficient, and the face assignment

**The derived coefficient is INERT, and the reason is geometric, not numerical.**
Restricting the annulus pair to `r > foot_r_outer` removes *all* wall heating,
the hole then stops widening at the pad edge, and the wall band **empties
itself** — so `h = 38.4` and `h = 169.3`, a span of 4.4×, produce the same run to
0.7 %. The binding term is the bore and pad-ring lateral faces, exactly as the
packet predicted, and the successor is **lateral recession**.

## 0. Provenance

* **`PREFLIGHT.md`** sha `9840fa2bd373d7b8674b6b6ac2e7dd18697e9cb03104a77dc22b1c3eae35f0fd`,
  **`preflight.py`** `322791fd988d9e8c35d0d3267f074f4d7e0458e3aacf10170d9e8041c30cc9e6`
  — frozen chmod 444 **before `CRITERION.md`**, both with **zero outputs on disk**.
* **`CRITERION.md`** sha `19feec7e23be0155e1bf21c8df5f1cc0447daaeb2f3bb962790d9706c15a67ad`,
  chmod 444, **zero outputs on disk** at the time of hashing.
* **Campaign binary** `e3a7f1601fac897a55bbb42fe9074ac3ec5ac4185dc47d9ac73a8b3d123b1476`
  (`BINARY.sha256`). All six scored legs ran against it, no rebuild mid-sweep.
* **Post-campaign binary** `ea833f2f03855145a93eb98508e04a92cbc1aff80bfda40b1cc32bf96776c471`
  — adds only the `live`-without-march parse guard. `S0R` re-runs `S0` against it
  to **prove** the guard is parse-only rather than asserting it.
* Nothing committed. Frozen studies imported by path, never re-run.

## 1. Item 0 — the derivation, and the D2l debt

**The blocking D2l factor-of-four audit is DISCHARGED.** D2l *did* apply it:
`vdot()` returns `4.0 * …`; `tse_rock = 4.0*p_face/vmm3` puts the 4 on both sides
so it cancels; `tse_hhv` pairs the full 38 kW burner with the full-hole volume.
Reproduced independently from `d2l_scan/output/M_f10/thermo.dat`: quarter-hole
`jet_P_face` 3204.2 W × 4 = **12.82 kW**, matching `RESULTS.md:190` exactly.
**"17.3 % of HHV" and "needs 9.1 %" stand as written. Nothing to restate.**

Derived bracket, from Meier's **measured** Ø85–93 against the Ø80 skirt (the
leg's own diameter is a model *output* and was not used), the **nozzle** mass
flow 15.13 g/s, and gas properties at the legs' own exhaust temperature:

| | Ø85 | Ø93 |
|---|---|---|
| `D_h = D_o − D_i` | 5.0 mm | 13.0 mm |
| Re | 2134 | 2035 |
| laminar, Nu = 5.385 | 99.7 | **38.4** |
| Dittus-Boelter | **169.3** | 62.7 |

`Re = 2ṁ/(π(R_o+R_i)μ)`, so **the hole diameter very nearly cancels** — the
85–93 mm uncertainty moves Re by < 5 %. What moves Re is the **mass-flow basis**,
by 2.3×. The nozzle feed is primary because under `jet_T_ent_mode = exhaust` the
march's entrained gas *is* the hole's own recirculated exhaust, so counting it in
the annulus double-counts mass that never crosses the exit plane.

**Bracket: h_low 38.4, h_high 169.3 W/m²K.** No point value inside it was taken.
It contains von Rohr's 81, so P4 was already true on paper.

**Two disagreements with the packet, reported not replaced.** The gap top end is
**6.5 mm**, not 6 — `(93−80)/2`. And item 0 predicted **2.87 / 3.16 m/h** against
the packet's P2 band of 2.3–2.9, the difference being the bore+pads term; that
was written down and hashed *before* launch.

## 2. Scorecard

| | result | |
|---|---|---|
| **(a)** derived bracket vs its own prediction | **0/2 FAIL** | GL +35.8 %, GH +24.2 % |
| **(b)** attribution pair | **2/2 PASS** | both bars |
| **(c)** face assignment | **3/6** | c1 FAIL, c2 PASS, c3 PASS |
| **(d)** power identity | **6/6 PASS** | ≤ 1.33e-2 W, `ledger_err` 1.24e-14 |
| **(e)** identity | **4/5** | one criterion defect of mine, below |
| **(f)** Meier | reported | |

## 3. The headline — the wall band empties itself

| leg | h | r_min | bore | pads | **wall** | D@100 | ROP |
|---|---|---|---|---|---|---|---|
| `SWc` | inherit | — | 402.8 | 2063.9 | **1074.2** | 87.6 | 3.694 |
| `GL` | 38.4 | 0.040 | 394.1 | 2199.4 | **0.0** | 79.6 | 3.897 |
| `GH` | 169.3 | 0.040 | 395.1 | 2206.9 | **0.0** | 79.7 | 3.925 |
| `WLr` | 81.0 | 0.040 | 393.8 | 2181.2 | **0.0** | 79.6 | 3.858 |
| `WLv` | 81.0 | 0.0 | 159.3 | 493.4 | 329.1 | 85.3 | 1.887 |

Powers W, full hole. `side_A_wall` is **0.00 cm²** on the gated legs — not zeroed
power, **no faces at all**. The ring table shows the cliff: recession 244 mm at
the 38 mm ring and **exactly 0.0 at 40 mm and beyond**, where `SWc` still cuts
247.8 at 40 and 62.1 at 42.

**The mechanism.** `foot_r_outer` is 40 mm and the burner is Ø80. Gating the pair
to `r ≥ 40 mm` means the only heat available beyond the pad edge is the annulus
coefficient itself. At 38–169 W/m²K against an inherited impinging map of many
hundreds, that is far too weak to *initiate* widening — so the hole never crosses
40 mm, the wall band stays empty, and `h` never acts on anything. It is
self-reinforcing: no wall heat ⇒ no widening ⇒ no wall faces ⇒ no wall heat.

**The hole ends up 79.6 mm — narrower than the Ø80 burner that drilled it.** That
is only possible because undercut is impossible in this code; a real burner could
not descend. It is the same architectural limit that jammed `PD` in D2aa, showing
up from the other side.

⚠ **`GL` ≡ `GH` to 0.7 %**, so **P4 is refuted in a degenerate way**: the rates
cannot bracket `WLr` because all three are the same run. The bracket *contains*
81 as predicted, but the containment is vacuous.

## 4. The attribution pair — one edge each

| leg | T | faces | ROP | side P | pads |
|---|---|---|---|---|---|
| `WL81` (frozen) | 1463.5 const | all | 1.780 | 825.5 | 416.2 |
| `WLv` | **live** (1622 K) | all | **1.887** | 981.8 | 493.4 |
| `WLr` | live (1531 K) | **r>40** | **3.858** | 2574.9 | 2181.2 |

* **The frozen-T blocker was real but small: +6.0 %** (1.780 → 1.887). D2aa's wall
  numbers were self-inconsistent, and correcting them moves the rate by 6 %, not
  by the 22–27 % the ΔT error suggested. `WLv` lands inside its pre-registered
  1.85–2.00 band (**P3 HELD**).
* **The face assignment is worth +104 %** (1.887 → 3.858). Pads go 493 → 2181 W
  because they revert to the inherited impinging map instead of being throttled
  to 81 W/m²K.

So of D2aa's wall treatment, **two thirds of the effect was misapplied heat**, as
the packet argued from the region split — confirmed here by direct measurement.

## 5. (c) — what passed and what did not

* **(c2) PASS.** Across `GL`/`GH`/`WLr`, bore spreads **0.35 %** and pads
  **1.17 %** while `h` spans 4.4×. The gate does not leak with `h`.
* **(c3) PASS — the positive control.** `WLr` vs `WLv` at the *same* h = 81 and
  same live T: bore **+147 %**, pads **+342 %**. The gate is emphatically not a
  no-op. **Without this leg, (c1)/(c2) would have passed vacuously on a run where
  the mechanism did nothing** — it is the check that earns the others.
* **(c1) FAIL, for a physical reason.** Bore sits within 2.2 % of `SWc` but pads
  run **+5.7 to +6.9 %** against a 5 % bar. This is *not* a gate leak: it is the
  march's enthalpy coupling. `SWc` spends 1074 W heating the wall; the gated legs
  spend none, so their gas stays **68 K hotter** (1531 vs 1463 K) and the *same*
  inherited-map pad faces receive more power. A 5 % bar on a band that is coupled
  to the others through the gas temperature was too tight. **Scored as written
  (FAIL); the mechanism is identified rather than explained away.**

## 6. Why (a) failed — and the correction it forces

The **power** prediction was excellent: item 0 said 6.86–7.30 kW, measured
7.00 kW (**+2 %**). The **rate** prediction missed by 32 %. So what failed is the
`rate ∝ P^1.59` scaling inherited from D2aa, not the derivation.

**Specific energy is invariant at 1.33–1.42 J/mm³ on every leg**, including the
new ones — D2aa's finding reproduced. That fixes *volume* rate to power. But
**descent** rate also depends on how much rock a metre of descent costs:

| leg | P kW | V L | D@100 | P/V | ROP |
|---|---|---|---|---|---|
| `GL` | 7.00 | 0.91 | 79.6 | 7.69 | 3.897 |
| `GH` | 7.01 | 0.93 | 79.7 | 7.54 | 3.925 |
| `WLr` | 6.98 | 0.92 | 79.6 | 7.59 | 3.858 |
| `SWc` | 8.26 | 1.08 | 87.6 | 7.65 | 3.694 |
| `WLv` | 5.25 | 1.21 | 85.3 | 4.34 | 1.887 |
| `WL` | 5.65 | 1.23 | 86.7 | 4.59 | 2.073 |
| `WL81` | 5.09 | 1.22 | 85.2 | 4.17 | 1.780 |
| `A` | 4.18 | 1.27 | 83.8 | 3.29 | 1.251 |

**ROP correlates with P/V at r = 0.9974, against 0.9255 with P alone**, over
eight legs from three studies. `SWc` has the *most* power (8.26 kW) and is
*slower* than `GL` (7.00 kW), because it cuts a wider hole.

> **D2aa's "rate follows TOTAL power; only shape follows the distribution" is too
> strong.** It held there because every D2aa leg cut a similar diameter. The
> correct statement is: **specific energy is invariant, so VOLUME rate follows
> power; DESCENT rate follows power divided by the volume excavated per metre.**
> Shape and rate are *not* separable, and any future rate prediction must carry a
> diameter prediction with it.

## 7. Meier

**No leg reaches the 1.3–1.6 m/h band (P7 HELD).** Rates 1.887 to 3.925. The
gated legs are ~2.5× fast *and* cut a hole 79.6 mm wide — below Meier's 85–95 mm.

⚠ **Fig. 8.8 is a visible width = a LOWER BOUND on diameter**, band 85–95 ± 5 mm;
no agreement claim tighter than ± 5 mm is made. Meier's rate is **one-sided**
(operator advance).

**Rate is not steady** (deliverable 3), reported in 50 mm stretches:

| leg | ROP | 50–100 | 100–150 | 150–200 | 200–250 | 250–300 | 300–350 |
|---|---|---|---|---|---|---|---|
| `GL` | 3.897 | 3.87 | 3.71 | 3.67 | 3.76 | 3.95 | 3.59 |
| `GH` | 3.925 | 3.74 | 3.74 | 3.81 | 3.73 | 3.93 | 3.58 |
| `WLr` | 3.858 | 3.74 | 3.83 | 3.72 | 3.73 | 3.82 | 3.76 |
| `WLv` | 1.887 | 1.76 | 1.93 | 1.96 | — | — | — |
| `SWc` | 3.694 | 3.65 | 3.58 | 3.64 | 3.77 | 3.83 | 3.67 |

The gated legs are **flat** to ±5 % over 300 mm — a consequence of the hole being
pinned at the burner radius, so the geometry stops evolving.

## 8. Identity and the deferred guard

* **`S0` vs `RB_D`: PASS.** `_removal_events.csv` and `_jet_profile.csv`
  byte-identical, thermo **header delta EMPTY**, 0 columns differ.
* **`SWc` vs frozen `SW`: 0 SHARED columns differ**, removal log byte-identical.
  ⚠ **Criterion defect of mine, recorded not reinterpreted.** `CRITERION.md` (e)
  says `SWc`'s header delta "must also be EMPTY". That sentence is **wrong** and
  contradicts the check made minutes earlier that put `SWc` in the campaign at
  all: frozen `SW` predates D2aa, so `SWc` *necessarily* gains D2aa's seven
  region columns. **Scored as written (FAIL);** the substantive result — `SWc`
  reproduces frozen `SW` exactly on every column both carry — is what
  re-validating the baseline required, and it holds.
* **Tree gate vs `witness/d2z1.md5`: 2639/2639, ZERO differing hashes**, no path
  added or removed. With the new keys unset not one file moves by a byte. The
  guard confirms the binary did not move mid-sweep. (`beam_void_closure` rc=1 is
  D2z-1's documented T3 debt, the only red, unchanged here.)
* **`S0R` vs `S0`: `_removal_events.csv`, `_jet_profile.csv` and `thermo.dat` all
  IDENTICAL.** The `live`-without-march guard was written during the campaign and
  deliberately **not** built (never rebuild mid-sweep), which left the six scored
  legs on the pre-guard binary. Re-running the identity leg against the
  post-guard binary **proves** the guard is parse-only rather than asserting it,
  so the scored legs remain valid against the final source. Detail in `GATE.txt`.
* **Unit test `unit/side_face_step_height`: 46/46 PASS** (38 inherited + 6 new
  D2ab checks + 5 new aborts).

## 9. RT — dropped, as the packet's own condition requires

Item 0 found **no flux-rate or flux-temperature source** in the repo. What exists
is Rauenzahn's **surface-temperature** compilation (Fig. 2.12,
`meier_thesis_notes.md:123`), against which the model's 528 K already sits inside
the band (`:222`); every other figure on record is **fuel-basis** (Meier's TSE
15.4 J/mm³ on HHV, Browning and Fleming 8.7–60) against the model's
delivered-to-rock 1.33. A prescribed-flux flat face returns
`q/(ρ c_p ΔT_spall)` by construction, so it would have tested only the spall
temperature. **Carried as a citation, not a leg.**

## 10. What this packet does and does not establish

**Does:** the annulus coefficient is derived, not chosen; the face assignment is
demonstrated by a positive control; the frozen-T blocker is fixed and is worth
+6 %; the misapplied-heat share of D2aa's wall treatment is two thirds; specific
energy is invariant across eight legs; and descent rate follows **P/V**, not P.

**Does not:** reach Meier on rate or diameter; produce a usable wall coefficient,
because under this configuration **no face is in the annulus**; or say anything
about the bore and pad-ring faces, which are in impinging flow and now carry
**93 %** of the charged side power (2593 of 2593 + 0 W).

**The successor is lateral recession / undercut.** It is now the critical path
for three independent reasons: it is the only way the wall band can ever become
non-empty; it is what jammed `PD`; and it is what would stop the model drilling a
hole narrower than its own burner.

## 11. Named debt

* **Lateral recession — critical path.**
* **`beam_void_closure` T3** — still the only red in the tree, D2z-1's documented
  debt.
* **`A` at 1 mm** — user-parked.
* **The 5 % bar in (c1)** was set without allowing for the march's enthalpy
  coupling between the region bands; a successor should bar the *population*
  (area, face count), not the power, when the bands are thermally coupled.
