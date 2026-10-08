# D2ad — calibrate the floor flux on Fig. 8.8's standoff: THE WINDOW IS EMPTY

**No `(F_q, T_s)` setting in this two-lever family reproduces Fig. 8.8's
geometry.** The two halves of the window are satisfied by **disjoint** sets of
legs: the standoff comes down only when the floor flux is cut, and cutting the
flux starves the **pad-ring annulus** before it brings the standoff inside
100 mm. There is no overlap, and the failure is **structural, not a resolution
problem** — it is a radial threshold the two-lever model cannot represent.

## 0. Provenance

* `PREFLIGHT.md` `427d727be0548d38029331d6fcf876d8e365228624ca64fd974d184256a73eeb`,
  `preflight.py` `a684615cc66d30087625a111196fd4936f5eeb0278299ddeca77735f145d2197`
  — frozen chmod 444 after item 0, before `CRITERION.md`.
* `CRITERION.md` `4a1ec1c80f6c530061711c235977896a5f7ba4c9938b5d34f6590cb07955b09e`,
  chmod 444, hashed before the first ladder leg.
* **Binary `ea833f2f…` before AND after** — **no source change**, as required.
* Tree gate **2639/2639, zero differing** (`--force`). Nothing committed.

## 1. Item 0 overturned the ladder, through the packet's own mandates

| item | result |
|---|---|
| **(b)** `kic.scale` → T_s | **measured** 821.3 / 864.6 / 895.5 / 919.1 / 950.1 K at scale 1.0–2.2; firing window **10.5–18.8 K** (bar 30) ⇒ the lever moves temperature, not the criterion. 870 K → **1.3523**, 920 K → **1.9088** |
| **(c)** back-reaction | **FAILED** the 40 K bar at **53.5 K**. C1's linear `1900 − 235.3·F_q` replaced by a quadratic through the three measured points |
| **(d)** `G1600` floor T_g | **1465.8 K measured**, not the estimated 1410 |
| **(e)** `R90` | **4.237 m/h**, not 3.3 — **+28.4 %**, over the packet's own 10 % trigger |

Both recompute mandates fired, so the ladder was recomputed: **`s* = 254.2·F`,
not `198·F`**, contact at **`F ≤ 0.354`** not 0.4545, window `s*` ∈ (90, 100].
That **inverted the pre-registered predictions** — the packet's predicted pass
`Q5T870` moved to 134 mm, and `Q35T870`, a leg the packet *dropped* as "contact
under C4" on the refuted `R90` = 3.3, became a predicted pass. I kept `Q5T870`
so P1 could be tested directly rather than quietly dropped.

**Two clean witnesses.** The `h_expr` prefix is exact — floor power ratios
**1.000 / 0.520 / 0.350** at F_q = 1.00 / 0.50 / 0.35. And `SMOKE_Q10`
reproduced frozen `FH` to **−0.0 K**, a self-check on the whole lifting chain.
That second one only appeared after fixing a bug in my own reader: a 60 s run's
last profile row is t = 50.016, so I had been comparing smokes at 50 s against
`FH` at 60 s — exactly the like-for-like error item 0(c) exists to prevent.

## 2. The result — (a), both halves

| leg | pred `s*` | `s_end` | slope | (i) ≤100 | 1st contact | worst | (ii) | verdict |
|---|---|---|---|---|---|---|---|---|
| `C0` | 254 | 180.0 | 7.41 | FAIL | never | −3 mm | **PASS** | fail (i) |
| `G1600` | 194 | 144.7 | 2.48 | FAIL | 62 s | +151 | FAIL | fail both |
| `Q5` | 153 | 110.7 | 1.36 | FAIL | 62 s | +183 | FAIL | fail both |
| `Q5T870` | 134 | **96.7** | 1.89 | **PASS** | 62 s | +275 | FAIL | fail (ii) |
| `Q35` | 111 | **54.7** | 0.43 | **PASS** | 62 s | +289 | FAIL | fail (ii) |
| `Q35T870` | 97 | 28.0 | −3.56 | FAIL | 61 s | +198 | FAIL | **ABORT** |
| `Q30` | 96 | 26.0 | −3.83 | FAIL | 61 s | +194 | FAIL | **ABORT** |
| `Q30T870` | 84 | 20.0 | −9.28 | FAIL | 61 s | +106 | FAIL | **ABORT** |

**(i) passes on `Q5T870` and `Q35`. (ii) passes only on `C0`. The intersection is
empty.**

Three legs **aborted on the bore collision guard** at column (13,3), **r = 27.9
mm** — the nozzle itself reached rock, which is stronger than pad contact. Their
traces are valid to the abort and scored there.

⚠ Contact is **persistence, not onset**. The feet start exactly on the flat
surface by construction, and the question is whether the pad ring ever pulls
away. Measured gap [mm], positive = ploughing through rock:

| leg | t=10 | 30 | 60 | 100 | 200 | 400 |
|---|---|---|---|---|---|---|
| `C0` | +2.2 | +2.5 | +1.0 | **−20.3** | −58.7 | −85.3 |
| `G1600` | +4.2 | +6.5 | +9.0 | +11.7 | −8.7 | −15.3 |
| `Q5` | +4.2 | +8.5 | +13.0 | +15.7 | +15.3 | +16.7 |
| `Q35` | +4.2 | +10.5 | +19.0 | +23.7 | +39.3 | +122.7 |

**`C0` is the only leg that clears and stays clear.** `G1600` ploughs to 150 s
then clears — under the frozen criterion ("no contact after t = 60 s") it fails
as written, and that refinement is reported, not substituted.

## 3. Why every `s*` predictor failed — the mechanism

Recession by radius at each leg's last time [mm]:

| leg | 0–8 | 8–16 | 16–24 | 24–32 | 32–40 | 40–48 | **shoulder lag** |
|---|---|---|---|---|---|---|---|
| `C0` | 402.0 | 401.3 | 398.1 | 392.2 | 383.2 | 238.1 | **9.8** |
| `Q5` | 394.0 | 392.4 | 389.0 | 383.8 | 355.0 | 46.6 | **10.2** |
| `Q5T870` | 380.0 | 378.3 | 375.0 | 368.8 | 275.8 | 25.7 | **11.2** |
| `Q35` | 338.0 | 336.7 | 332.6 | 323.3 | 136.9 | 23.8 | **14.7** |
| `Q35T870` | 201.2 | 199.3 | 195.5 | 159.6 | 57.3 | 13.9 | **41.6** |
| `Q30` | 200.0 | 197.7 | 194.0 | 160.3 | 62.8 | 17.6 | **39.7** |

**The floor does not slow uniformly as the flux falls.** The outer bore drops
below the firing threshold and *stops*, while the axis keeps going. The shoulder
lag (how far r = 24–32 mm trails the axis) sits at 10–11 mm on the high-flux
legs and jumps to **40 mm** on the aborted ones. The burner descends on a
schedule, lands on that shoulder, and the guard fires at r = 27.9 mm.

> **This is a THRESHOLD effect.** The `F` model is linear in `(T_g − T_s)` and
> cannot represent a *radius at which firing stops*. That is why **both** `s*`
> predictors — the packet's and item 0's — overpredicted the equilibrium
> standoff on **every** leg, and why the two "in-window" legs collapsed into
> contact instead of settling at 96–97 mm.

✅ **P2 HELD**: the ordering is monotone, `C0` 180 > `G1600` 145 > `Q5` 111 >
`Q5T870` 97 > `Q35` 55. **The ladder ranks correctly even though every absolute
prediction is wrong** — the levers work; the equilibrium model does not.

## 4. The flux gate (b)

| leg | F_q | T_s | T_g measured | q100 | vs Kant 0.45–0.87 | vs Meier 0.5–0.6 |
|---|---|---|---|---|---|---|
| `C0` | 1.00 | 823.3 | 1684.3 | **1.05** | **ABOVE** | above |
| `G1600` | 1.00 | 823.3 | 1452.8 | 0.77 | inside | above |
| `Q5` | 0.50 | 823.3 | 1776.7 | 0.58 | inside | inside |
| `Q5T870` | 0.50 | 870.0 | 1786.3 | 0.56 | inside | inside |
| `Q35` | 0.35 | 823.3 | 1848.6 | 0.44 | below | BELOW |
| `Q35T870` | 0.35 | 870.0 | 1862.1 | 0.42 | below | BELOW |

**P3 first clause HELD** — `C0` at 1.05 sits above Kant's measured band, the
known error restated on a correct basis. The third clause is **vacuous**: there
is no passing leg.

⚠ The legs that satisfy **(a)(i)** — `Q5T870` 0.56 and `Q35` 0.44 — straddle
Meier's required band. So the standoff half of the window and the flux gate are
*not* in conflict; the conflict is entirely with **(a)(ii)**, the pad ring.

## 5. Scorecard

| | result |
|---|---|
| **(a)** | **NONE — the window is EMPTY** |
| **(b)** | reported; `C0` above Kant (P3 clause 1 held) |
| **(c)** Ø 120–250 mm vs 85±5 | **0/8** — every leg outside |
| **(d)** volume ≤ 3.42 L | **8/8 PASS** |
| **(e)** firing window ≤ 30 K | **8/8 PASS** |
| **(f)** power identity ≤ 0.02 W | **9/9 PASS** |
| **(g)** identity + tree gate | **3/3 PASS**, 2639/2639 zero differing |

⚠ **(c) and (d) are VOID from sustained contact on every leg but `C0`**, and are
read only up to first contact. Their numbers are therefore weak evidence, and
(c)'s 0/8 in particular should not be read as a shape result.

## 6. What this establishes, and what it kills

**Establishes.** Both levers are live and calibrated: `F_q` scales floor power
exactly (1.000/0.520/0.350) and `kic.scale` moves T_s 821 → 950 K with a firing
window of 10–19 K. The ladder orders correctly. `R90` = 4.237 m/h and the
back-reaction are now measured, not assumed.

**Kills.** The two-lever calibration premise. Reducing floor flux does **not**
translate the whole floor downward; it shrinks the *fired radius*, so the
pad-ring annulus starves first. Holding Fig. 8.8's standoff and keeping the
burner off the rock are **incompatible** in this family — and no amount of
resolution in `F_q` or `T_s` closes a gap whose cause is radial, not scalar.

**⚠ The equilibrium `s*` model is retired.** Both versions overpredicted every
leg. Any successor must predict the *floor profile*, not a single standoff.

## 7. Named debt

* **Lateral recession** — still the critical path, and this packet sharpens why:
  the shoulder at r ≈ 28 mm cannot recede laterally, so the burner lands on it.
* **The radial fall-off of the map is now the binding term**, not its level.
  A successor should test the map's *shape* (how fast `h` falls with r across
  28–40 mm), which is what starves the pad ring.
* ⚠ **`prescribed` has no force balance**: the pads pass through rock by up to
  289 mm on `Q35`. Any future prescribed leg must treat (c)/(d) as void from
  contact, as done here.
* ⚠ Tree-gate cache still keyed on the binary alone ⇒ `--force` mandatory on an
  inputs-only packet.
* `beam_void_closure` **T3** — unchanged, the only red.
* `A` at 1 mm — user-parked. AMR parked. No third mesh.
