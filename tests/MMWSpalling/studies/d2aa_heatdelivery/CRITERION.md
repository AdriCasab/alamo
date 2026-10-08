# D2aa CRITERION — hashed before the first leg, zero outputs on disk

Separate the two physical levers for the 3.54 kW side-face surplus: the floor
under the pads, and the hole wall. Base is **frozen `SW`**, one edge per leg.

Import the equal-burner-depth machinery (`shape`, `rop`, `rings`, both volume
bases) **unchanged** from `studies/d2x_areasplit/analyze.py`. Frozen studies are
read by path and never re-run: `SW`, `SWs`, `SW0`, `A`, `D`, `M1`, `RB_D`.

⚠ **All Meier thermo powers are quarter-hole sums** (`x0 = y0 = 0` on a
0.12 x 0.12 m domain). Full-hole = 4x. `d2x_areasplit/analyze.py:104` already
multiplies volume by 4.

## Legs — one edge each from `SW`

| leg | edge | keys |
|---|---|---|
| `WL` | von Rohr grazing pair, **annulus derivation** | `side_face_h=142`, `side_face_T=1463.5` |
| `WL81` | same, **Dittus-Boelter end of the band** | `side_face_h=81`, `side_face_T=1463.5` |
| `PD` | pad ring fully shielded | `pad_shield_open=0.0` |
| `BO` | both primaries | `side_face_h=142`, `side_face_T=1463.5`, `pad_shield_open=0.0` |
| `S0` | `RB_D`'s argv, all new keys absent | — |

⚠⚠ **`WL` is run at BOTH ends of the band and neither end is chosen.** The packet
writes `side_face_h=<von Rohr band>` with no number, and P3's 0.4–0.5 kW implies
h ≈ 81–97 — but the standing constraint is that **142 is the annulus derivation's
value, never varied**. Picking 81 would land P3 and would be choosing `h` to hit a
rate, which is the guardrail this packet exists to respect. **142 is the primary**;
81 is the band's other end. `side_face_T = 1463.5 K` is `jet_T_exhaust` **measured
on `SW`** over 150–400 s (1463.46), not chosen.

⚠⚠ **`PD` is BOUNDED, not fitted.** The slot width is on record nowhere
(`meier_thesis_notes.md:40` says only "three exhaust slots"), so the ring's open
fraction cannot be derived. `pad_shield_open = 0` is the packet's own physical
statement ("under them there is essentially no gas") and is **one end of the
bracket; the other end is frozen `SW` itself** (open = 1). The true value lies
between them and closing it needs drawing sheet 15. **No intermediate value is
invented.** Per risk 3 the shield is **adiabatic-except-at-the-slots**, NOT a
sink: a sink needs a contact conductance and a body temperature, neither known.

---

## (a) SCORED — rate, shape and volume together, per leg
Mean |ΔØ| at equal burner descent (depths capped 50 mm above the burner, scored
only on rock the nozzle plane has passed), Ø at 100 and 125 mm, burner descent,
pit lead, volume on **both** bases, and the 2 mm ring table 30–50 mm **per ring,
never averaged**.
⚠ **Rate definition: burner descent over the 100–150 mm stretch** (`A` reads
1.251 there, 1.296 over 150–400 s).
⚠ **A rate match bought with a deep pit must be visible** — pit lead is scored
alongside, not reported after.
⚠ **Pit-lead ceiling (item 0(c)):** a leg whose pit lead exceeds **150 mm** is
reported as disqualified from the rate comparison. See (g) for why the clip is
nevertheless not a discriminator.

## (b) SCORED — THE FALSE-PASS TEST
**An in-band rate is a PASS only if all three hold:**
1. the 30–40 mm rings do **not** fall the way `A`'s did (125 → 90 → 69 mm);
2. `wall_step_f_multi_pads` does **not** rise across the window;
3. no 1 mm-equivalent claim is made.

**If the rate is in band and (1) or (2) fails, the verdict is "the throttle
returned", not "the leg passed."**

## (c) SCORED — the power identity
`jet_P_face − (patch_P_robin + side_P_face + side_P_loss)` on every leg.
**Bar: ≤ 0.02 W** (6-significant-figure floor 0.015 W).
⚠ The only check that catches a kernel↔march split. Both new mechanisms were
therefore applied where the **march** sees them — `pad_shield_open` in
`BuildFlameColumns` and the `(h, T)` pair at the march's own side term — not at
the Robin sites alone.

## (d) SCORED — the region split closes
`side_P_bore + side_P_pads + side_P_wall` equals `side_P_face` to **±0.01 %**, and
the three `side_A_*` sum to the charged lateral area. Accumulated from the same
per-face term, so any gap is a binning bug, not round-off.

## (e) SCORED — identity
`S0` byte-identical to `d2z1_delete/output/RB_D` (thermo, `_removal_events.csv`,
`_jet_profile.csv`), **header delta empty**. Tree gate against `witness/d2z1.md5`.
⚠ `metadata` and `diff.patch` are run provenance and excluded, as the tree list
excludes them.

## (f) REPORTED — specific energy
Report power-to-rock / volumetric rate per leg against the model's own
**1.35 J/mm³**. ⚠ If rock really takes 4–6 J/mm³ the model's spall energetics are
wrong by ~3x — a larger, separate finding. **Not a licence to keep the power**,
and not acted on in this packet.

## (g) REPORTED — item 0's three audits
Lateral-area decomposition, the D2l factor-of-four audit, and the decay-clip
status.

⚠⚠ **The packet's decay-clip figure is CORRECTED here, before any leg runs.**
`SW` and `A` both set `jet_core_length = 8.0` and `jet_D = 0.0075`, so the core
length is **8 D = 60 mm, not "twelve nozzle diameters = 90 mm"**. Measured, the
stand-off exceeds it in **100 % of rows on every leg** (`A` 157.3 mm mean, `SW`
87.8, `D` 81.5) — the decay is always engaged, never clipped at phi = 1. It is
nevertheless **not a discriminator**: under `jet_T_ent_mode = exhaust` the mixing
temperature is already 1463 K, and measured `T_stag` is **1776–1780 K across all
three legs** despite `s_c` spanning 81–157 mm. The pit-lead ceiling in (a) is
registered as the packet asks, but it is not load-bearing.

---

## Pre-registered predictions

* **P1 — reviewer's reading, PRIMARY.** `WL` rate **1.8–2.2 m/h**, wall band
  falling back to **84–85 mm**. Basis: `A`'s 1.251 was the step drain, D2x
  measured the un-throttled pad-ring rate at **2.09–2.11 m/h**.
* **P2 — planner's alternative.** `WL` rate **1.3–1.6 m/h**, band near **87 mm**.
* **P3** `WL` side path lands **0.4–0.5 kW** full-hole, down from 3.54 kW. Held
  inside **0.3–0.7 kW**, refuted outside. Scored on the **h = 142 primary**;
  `WL81` reported as the band's other end.
* **P3b** `PD` reduces the pad-ring top-face power by **≥ 5x** and the descent
  with it. Held / refuted on the measured ring power, not the rate.
* **P4** `BO` lands **below** both single edges on rate. ⚠ May fail for a reason
  that is not a defect: the march **conserves gas enthalpy**, so heat not
  deposited under the pads arrives at the slots and the wall. **The two levers
  interact by construction.**
* **P5 — the false-pass guard.** On any leg landing in band,
  `wall_step_f_multi_pads` is **flat or falling** and the 30–40 mm rings hold.
  Refuted ⇒ the throttle returned.
* **P6** more than half the charged side area is **under the pads** (packet:
  bore 1 / pads 10 / wall 7 cm² per quarter), so `PD` moves the side term too and
  `WL` alone cannot account for all of the 3.54 kW. Held / refuted on
  `side_P_pads` vs `side_P_wall`.
* **P7 — the pads-leg jam signature, registered in advance.** A rate collapsing
  below **~0.5 m/h** with the 28–40 mm rings frozen is **the model's
  lateral-recession limit** (undercut is impossible in this code), **not
  physics** — read as D2o-0's 0.11 m/h jam, not as a failed leg.

## Refuting evidence on disk

* **`A` at 2.27 mm mean |ΔØ| is still the best shape on record.**
* **D2x side power `C` 223.2 / `D` 1234.7 / `E` 453.8 W** (quarter).
* **`C` 2.09 / `E` 2.11 m/h** is the un-throttled pad-ring rate. `WL` near 2.1
  confirms P1; near 1.4 confirms P2; **outside both** means something changed
  more than one thing.
* **No Meier agreement claim tighter than ±5 mm** — Fig. 8.8 is a visible width,
  i.e. a LOWER BOUND on diameter; band 85–95 ± 5 mm, rate band 1.3–1.6 m/h.
