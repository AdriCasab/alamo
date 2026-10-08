# D2w criteria — fixed and hashed BEFORE the first campaign leg launched

## 0. Disclosure

**Neither new leg (`A`, `B`) has been launched when this file is hashed.** What
has run: a 20 s smoke of each, recorded in `PREFLIGHT.md` §0, which established
only that `A` registers **no** side-face columns and `B`'s
`side_P_face_above` is exactly 0.00 in both directions. No shape, rate or power
result has been read from either.

`M0` and `M1` are **imported by path and never re-run**; their numbers were
reproduced exactly in `PREFLIGHT.md` §1.

## 1. The factorial

| cell | confined map | side-face path |
|---|---|---|
| `M0` | no | no |
| `A` | **yes** | no |
| `B` | no | **yes** |
| `M1` | yes | yes |

**This is an attribution packet: every criterion scores WHICH CELL OWNS AN
EFFECT, not agreement with Meier.**

## 2. The metric

Imported unchanged from `d2v_hmap/analyze.py`: final `k_top` from
`_removal_events.csv`, **fixed physical 2 mm** radial bins, mean face height per
bin, crossing radius **linearly interpolated** between bracketing bins. One
metric across all four cells or the factorial is meaningless.

**ROP = the slope of `nozzle_z` over 150–400 s.** Never pit depth — the burner
rides the feet, not the deepest point. Quoting pit depth as ROP is the error D2v
made.

**(a) is scored at the common time** (the smallest `t_end` over the four cells),
with each leg's end time and depth reported. D2v established (a) is
window-insensitive: `M0` reads 37.2 mm at both 401.6 and 500 s.

## 3. Criteria

* **(a) SCORED — WHO OWNS THE FUNNEL.** Mean |Ø_model − Ø_Meier| over
  10 / 25 / 50 mm, on all four cells. Baselines `M0` **37.23**, `M1` **3.82**.
  **Pre-registered: `A` within 5 mm of `M1`'s 3.82 ⇒ the map owns it.** If `A`
  is materially worse, the side faces were doing part of the shape work and
  D2v's framing is wrong in a second way.
* **(b) SCORED — WHO OWNS THE RATE.** Burner descent on all four, against `M0`
  **1.50** and `M1` **3.12** m/h. **Bar: |`B` − `M1`| < |`B` − `M0`|.** The
  finding is which single key group reproduces the doubling.
* **(c) WHO OWNS THE FLAT FLOOR.** Median radial step in `z(r)` inside
  r = 76 mm, against `M0` **6.9 mm** and `M1` **1.3 mm**.
* **(d) ADDITIVITY.** Interaction `M1 − (A + B − M0)` on (a) and (b). **A large
  interaction is itself a result** and would mean the mechanisms are not
  separable.
* **(e) POWER — THREE COLUMNS, NOT TWO.** `patch_P_robin` (tops),
  `side_P_face` (sides), `side_P_loss`, `jet_P_face` (**TOTAL**), mean over
  150–400 s, on all four, plus the identity
  `jet_P_face = patch_P_robin + side_P_face + side_P_loss` (verified to +0.00 W
  on `M1`). ⚠ `jet_P_face` is the total drawn from the stream, **not** the
  top-face term — the march gathers side faces too. `B` should carry side power
  close to `M1`'s 624.8 W.
* **(f) VOLUME ON BOTH BASES.** Pro-rate by **pit depth AND by burner descent**,
  printed together. D2v reported 0.82× on the first basis alone; on the second
  `M1` is ~0.97× and `M0` ~1.31×. **The basis changes the verdict.**

## 4. Pre-registration, with numbers

* **P1 — `A`'s ROP is TWO-SIDED, and the low side is the interesting one.** `A`
  puts the confined map's cold wall at 40–56 mm right beside the pads (28–40 mm)
  with **no side heating at all** — the configuration D2j-0b and D2o-0b showed
  draining neighbouring columns.
  * **Upper branch:** "map is rate-neutral" refuted if `A` > **2.0 m/h**.
  * **Lower branch:** if `A` < **1.3 m/h**, the map alone **costs** rate.
  * **Witnesses, all required:** pad-ring recession rate over r = 28–40 mm;
    `foot_stall_time`, `foot_carry_cols`, `foot_cols`; and the lateral ΔT at the
    wall (`M0` 312 K, `M1` 357 K).
  * ⚠ **The reading matters more than the number.** If `A` drops well below
    `M0`, the side-face path in `M1` was **compensating a drain the confined map
    creates**, not simply adding power — and then `f = 0.2` stops being a fudge
    and becomes **the wall boundary condition the model lacks**. The successor
    is then *"give the wall a derived BC"*, **not** *"derive `f`"*.
* **P2 — `B` reproduces most of the doubling, ROP > 2.5 m/h**, consistent with
  D2k's 2.7× at `f = 1`. **Refuted if `B` < 2.0 m/h.**
* **P3 — `A` removes most of the funnel:** (a) on `A` within 5 mm of 3.82.
  **Refuted if (a) on `A` > 15 mm.**
* **P4 — the flat floor belongs to `B`:** `B`'s median radial step closer to
  1.3 than to 6.9 mm. **Refuted if `A`'s is the closer.**
* **P5 — DECLARED OPEN.** Whether `A` delivers shape **and** in-band rate
  together. **No number predicted** — the discipline D2t and D2v honoured and
  D2s failed.

## 5. Guardrails binding on the scoring

* **`f = 0.2` is frozen and NOT scanned.** ⚠ Its status has changed: it is now
  implicated in a **factor-of-two rate error**, not a quiet placeholder. **Do
  not tune it to recover a rate.**
* **Report ROP as burner descent only**; report volume on **both** bases,
  always together.
* **`M0`/`M1` imported, never re-run.**
* **No 1 mm run, no convergence claim, no Meier validation claim.** The
  reference is a **lower bound** on Ø at ±5 mm; the burner is methane/air at
  λ = 1.2.
* **`M2`, `jet_above_h`, `jet_T_above` and the enthalpy march are OUT OF SCOPE.**
  The above-plane accounting defect (`:2138`, `:2224`) is real and is a later
  packet.
* **No `src/` change.** If one seems necessary, STOP and report.
