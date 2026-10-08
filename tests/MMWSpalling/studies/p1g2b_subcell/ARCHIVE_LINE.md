# P1g-2b `(z)` — the archive line

Deliverable 7. Lift verbatim into `docs/project/ARCHIVE_DONE.md`. Every number is
from `SCORE_p1g2b.txt`, `SUBCELL_PRE.md`, `WALL_FINDINGS.md` and
`MECHANISM_PRE.txt`. **No source was edited**; the binary hashed `7908bd6b…`
before and after both legs.

---

**P1g-2b — the lateral sub-cell mode at Meier scale (2026-10-06).** Two Meier legs
on the shipping binary settled three things that had been conflated, and retired a
number that three readers had quoted.

**What banking is worth: exactly nothing.** `CTRL3` (`lateral_subcell=0`) against
`CTRL2` (`=1`) — **one key, one binary** — gives a bore difference of **+0.00 mm in
all five scored depth bands**, with the firing temperature, the time-to-fire
(116 → 115 s) and the regime-1 count all unmoved. The mechanism is visible in one
pair of numbers: lateral **firings** rise **×1.35** (29,407 → 39,795) while lateral
**voids** stay **×1.00** (29,407 → 29,370). Banking changes how a firing is
*counted*, not how much rock *leaves*. **P1 confirmed**, and more sharply than
predicted — the effect is identically zero, not merely small.

**What the thermal half is worth: −3.0 mm, and it is the first measurement of it at
scale.** `CTRL2` against `CTRL4` (`lateral_subcell=1 lateral_thermal=1`) — again one
key, one binary — **narrows every band**: −3.74 / −2.31 / −2.53 / −3.02 / −3.51 mm,
mean **−3.02 mm**, carrying the depth-matched error from **+19.4 to +16.4 mm**.
Lateral firings fall ×0.86 and voids ×0.85, and the time-to-fire *lengthens*
(115 → 119 s) — i.e. moving enthalpy with the receding face leaves less heat behind
to drive the next firing, exactly the mechanism **P3** predicted, at the small
magnitude it predicted. ⚠⚠ **This figure is REPORTED, NOT ACCEPTED: `(b3)` FAILED on
`CTRL4`** — see the defect below.

**What the drift is worth: +4.95 mm, and it is the binary, not the mode.** The
**+4.2 to +6.6 mm** between the archived `Leg B` and `CTRL2` has been quoted three
times as the sub-cell mode's effect. It is not: those legs differ by one key **and by
the binary** (`e8d8db0a…` vs `7908bd6b…`), with P1g-0b/1/2/3 of source between.
`CTRL3` shares Leg B's keys and CTRL2's binary, so it splits the difference cleanly,
and the split is decisive — the lateral-firing rise factorises as **binary ×18.57**
against **banking ×1.35**, product **×25.12**, against an observed **×25.12**. Since
banking's bore effect is exactly zero, **the entire +4.95 mm bore drift is four steps
of source work**, and the likeliest single cause remains **P1g-0b's coefficient fix**.
**P2 confirmed.** The binary owns the time-to-fire too (139 → 116 s).

**⚠⚠ A DEFECT THE STEP FOUND: the thermal shift does not close the energy books at
scale.** `(b3)` bars `max |ledger_err_solid|` at 5e-3. `Leg B`, `CTRL2` and `CTRL3`
all read **7.78e-4** — identical, and an order under the bar. `CTRL4` reads
**8.94e-3**, over the bar on **67.1 % of steps**, first crossing at **t = 164.7 s**
and growing monotonically with the cumulative shifted enthalpy (6.2e-4 at 153 J
shifted → 7.9e-3 at 3,965 J → 8.9e-3 at 8,273 J). The residual is only **~1e-6 of the
enthalpy the shift relocates**, `lat_overlap` is **0**, and the shift's independent
discrete estimate agrees to **5.3 %** (`lat_shift_E_est/lat_shift_E = 0.9472`, better
than the 14.2 % P1g-2 measured on its fixture). ⇒ this reads as an **accumulating
bookkeeping residual proportional to relocated enthalpy**, not a gross leak — but the
bar is absolute and **the gate is not waived here**. ⚠ For the planner: `CTRL4`
relocates **8.3 kJ**, orders more than any leg the 5e-3 bar was set against, so
whether the bar is mis-specified for this regime or the shift term is genuinely short
is the open question — and it is a planner's call, not an implementer's.

**Wall heating is NOT the width lever, and the wall-jet law is deferred for the
*brake*, not the lever.** P2-0's case was a feedback: a duct coefficient climbs
steeply as the hole narrows. In a radial wall jet momentum is conserved over an
expanding annulus, so `v_m ∝ 1/R` at roughly constant `Re` and the brake is **linear**
— **2.8–3× weaker than the duct's on a common diameter pair** (4.04× vs 1.45× over
Ø84→Ø121.8; 1.75× vs 1.19× over Ø102→Ø121.8). ⚠ That ratio is the one wall-jet
quantity that is settled, because a `c_f` convention change is a constant multiplier
and every ratio survives it. The duct value **50.3 W/m²K** reproduces as the
far-field limit, correctly — it is simply not where the faces are. ⚠⚠ **Do not quote
an absolute wall-jet `h`, "169.3 is right in the mean", or "−3.9 mm" as the lever:**
the `c_f` definition is open (corner **268 or ~536**; the no-½ convention lands on
Martin proper, 529, to 1.3 %), so the void-weighted lever spans **0.82× to 1.65×** of
the constant — **a sign change**. **Open item, named: Launder & Rodi 1981's `c_f`
definition.** ⚠ Rev 1 declared it settled by matching 268 against 251 W/m²K read out
of Leg B's own `h_expr` — but that 251 is Martin **×** a hand-set `0.09/s` decay, a
**fitted map**. A physics estimate is validated against the correlation, never
against a fitted map.

**P4 confirmed, and it names the successor.** Neither half of the sub-cell mode is the
width lever: banking is 0.00 mm and the thermal half is −3.0 mm against a
depth-matched error of **+16.8 mm** (⚠ **not +10.7** — that figure scored every band
against Meier's single 96 mm *mouth* value, understating the bore error by ~40 %).
**+16.4 mm survives both.** ⇒ **P2-1 is the floor taper at 40–56 mm.** ⚠ Meier's
profile is a *visible width*, a lower bound on Ø at ±5 mm, so every error above is an
**upper** bound on the model's excess.

⚠ Also measured and worth carrying: the model cuts the wall at a median **90 mm above
the advancing floor**, with **0 of 27,270** voids within 25 mm of it — and **that zero
is SCOPE, not physics**. The freshly cut bore within 25 mm of the front is Ø27.9–70.0
mm, i.e. radii of 13.9–35.0 mm, **inside the 40 mm `side_face_h_r_min` gate**, so no
face there is in lateral scope at all. A wall-jet law changes nothing at the corner
unless the gate is lifted, and lifting it puts the law onto the pad-ring risers that
feed the floor rate (D2w/D2x/D2o-0) — a rate-containment decision, not a side effect.
