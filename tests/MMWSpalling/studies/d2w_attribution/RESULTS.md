# D2w — the 2×2: is the near-cylinder the MAP or the SIDE FACES?

**The attribution is clean and it corrects D2v.** The **confined map owns the
funnel**; the **side-face path owns the rate doubling**. D2v's headline
"`M1 − M0` is the funnel fix" conflated them.

| cell | map | side path | (a) mean \|ΔØ\| | ROP [m/h] | median step |
|---|---|---|---|---|---|
| `M0` | no | no | 37.23 mm | **1.502** (in band) | 6.90 mm |
| **`A`** | **yes** | no | **2.46 mm** | **1.296** | 3.92 mm |
| `B` | no | **yes** | 28.20 mm | **3.087** | 3.79 mm |
| `M1` | yes | yes | 3.82 mm | **3.123** | 1.33 mm |

Scored at the common time **401.6 s**; metric imported unchanged from
`d2v_hmap/analyze.py`. `M0`/`M1` imported by path, never re-run.

## 1. Who owns what

**The funnel is the MAP.** `A` alone gives **2.46 mm — better than `M1`'s
3.82** — while `B` alone only moves 37.23 → 28.20. `|A − M1| = 1.36 mm`, so
**P3 holds** and the map owns it. The sharpest single number: **at 10 mm depth
`M1` and `A` are both 100.4 mm, identical to a tenth of a millimetre.** The side
path adds *nothing* at the collar.

**The rate is the SIDE PATH.** `B` = 3.087 against `M1`'s 3.123 —
`|B − M1| = 0.04` versus `|B − M0| = 1.58`, so the bar passes and **P2 holds**.
The side-face path reproduces essentially the whole doubling **on its own**.

**The flat floor needs BOTH.** `M1`'s median radial step is 1.33 mm; `A` (3.92)
and `B` (3.79) are both roughly half-way from `M0`'s 6.90. **P4 holds** — `B` is
marginally closer — but the honest reading is that neither mechanism alone
flattens the floor.

## 2. ⚠ The headline P5 declared open: `A` gives the best shape at the band edge

`A` is the **best shape of all four cells** (2.46 mm) at **1.296 m/h** against
Meier's 1.3–1.6 band — **0.004 m/h, or 0.3 %, below the lower edge**. Its
diameter tracks the reference from 100.4 mm at the collar to 83.8 at 100 mm,
where Meier reads 96 and 87.

So the answer to P5 is: **very nearly yes.** One mechanism, the confined map,
delivers the funnel removal *and* a rate at the bottom edge of the band —
a simpler and stronger result than D2v claimed, and it does **not** need the
side-face path that causes the doubling.

⚠ Still **one mesh, no convergence claim, no validation claim**: the reference
is a **lower bound** on Ø at ±5 mm and the burner is methane/air at λ = 1.2.

## 3. ⚠ P1 fired on the number and is REFUTED on its own witnesses

P1's lower branch ("`A` < 1.3 m/h ⇒ the map alone costs rate") fired at
**1.2963 m/h — by 0.0037 m/h, a 0.3 % margin on a pre-registered threshold.**
The *direction* is solid (a **13.7 %** drop from `M0`); the *branch label* rests
on a hair and should not be leaned on.

**More importantly, the mechanism P1 hypothesised is contradicted by all four of
its own required witnesses:**

| witness | `M0` | `A` | reading |
|---|---|---|---|
| pad-ring recession, r = 28–40 mm | 1.69 m/h | **1.60 m/h** | 5 % slower — **not stalled** |
| `foot_stall_time` | 1.28 s | 1.44 s | same order |
| `foot_cols` | 162 | **162** | identical |
| lateral ΔT at the wall | 312.3 K | **288.5 K** | **LOWER than `M0`** |

There is no drain and no stall. **`A` is slower because it absorbs less power**:
tops go **1427.5 → 1044.9 W, −27 %**, which is the confined map deliberately no
longer heating everything outside the mouth. A 27 % power cut buying a 13.7 %
descent cut is an energy-budget effect, not a drain.

**⇒ The packet's conditional follow-on should NOT be taken.** It says that if
`A` drops below `M0`, the side path was *"compensating a drain the confined map
creates"*, so the successor is a **wall BC** rather than **deriving `f`**. But
**`B` runs at 3.087 m/h with the stale map and no manufactured cold wall at all**
— so the side path is not compensating anything the map creates; it adds rate
independently. **The supported successor is "derive `f`".**

## 4. (d) Additivity — rate yes, shape no

**Rate is nearly additive:** `M1 − (A + B − M0)` = **+0.242 m/h, 7.7 % of `M1`**.
The side path adds ~1.59 m/h and the map subtracts ~0.21, largely independently.

**Shape is not**, and the packet's (a) form overstates it. (a) is a mean of
absolute deviations and is bounded below by zero, so an interaction there is
partly an artefact (+10.38 mm). On **raw diameters** the interaction is
**+9 to +12 mm at 10–125 mm depth and ≈ 0 below 150 mm** — because **the map
saturates the collar**: `M1` = `A` = 100.4 mm at 10 mm, so there is nothing left
for the side path to add there.

## 5. (e) Power — three columns, and the identity holds exactly

Means over 150–400 s [W]:

| cell | tops | sides | side loss | **TOTAL** | identity |
|---|---|---|---|---|---|
| `M0` | 1427.5 | — | — | 1427.5 | +0.00 |
| `A` | **1044.9** | — | — | 1044.9 | +0.00 |
| `B` | 1148.5 | **697.7** | 34.3 | 1880.5 | +0.00 |
| `M1` | 1080.5 | 648.9 | 31.0 | 1760.3 | −0.00 |

`jet_P_face` is the **TOTAL** drawn from the stream, not the tops — the march
gathers side faces too. The identity
`jet_P_face = patch_P_robin + side_P_face + side_P_loss` holds to **±0.00 W on
all four cells**. `B` carries **more** side power than `M1` (697.7 vs 648.9),
consistent with its unconfined map exposing more wall.

## 6. (f) Volume — the basis changes the verdict, as the packet warned

| cell | volume [L] | pit depth | ×(pit) | burner descent | ×(descent) |
|---|---|---|---|---|---|
| `M0` | 1.861 | 280.0 mm | 0.97 | 168.7 mm | **1.61** |
| `A` | 1.263 | 276.0 mm | 0.67 | 148.7 mm | 1.24 |
| `B` | 2.377 | 358.0 mm | 0.97 | 300.7 mm | 1.16 |
| `M1` | 2.042 | 362.0 mm | **0.82** | 306.7 mm | **0.97** |

D2v quoted 0.82× for `M1` on the pit basis alone; on the descent basis it is
0.97×. **Always print both.**

## 7. Witnesses and gates

* **N1** — `A` registers **no** side-face column at all (keys absent, not
  zeroed: `:6099` aborts on `side_face_factor` without `side_face_flux`).
* **N2** — `B`'s `h_expr` is byte-equal to `M0`'s. **N5** — `A`'s is byte-equal
  to `M1`'s 751-character token, copied from `M1.log`, never rebuilt.
* **N4** — `B`'s `side_P_face_above` is **+0.00 max and +0.00 min**: it neither
  heats nor cools its above-plane faces.
* **N3** — end times differ (`A` 500.0 s ok; `B` 408.4 s and `M1` 401.6 s on the
  frozen harness's bottom watchdog). Scored at the common 401.6 s; pit depths
  reported.
* `ledger_err` ≤ **1.24e-14** on all four. Regressions `spall_event`,
  `side_face_flux`, `convective_patch`, `jet_enthalpy` **PASS**. **No `src/`
  change**, so D2v's 2639/2639 identity still stands.

## 8. What this corrects in D2v

1. **The attribution.** "`M1 − M0` is the funnel fix" conflated two mechanisms.
   Correct: **map → funnel, side path → rate**.
2. **The rate was never reported.** `M1` runs at 3.12 m/h, twice Meier's band,
   on the configuration D2v presented as a success. D2v's packet said "score
   shape and volume, never rate" and I complied — but a guardrail against
   *scoring* a quantity is not a licence to leave it **unreported**.
3. **`jet_P_face` was quoted as the top-face power.** It is the total.
