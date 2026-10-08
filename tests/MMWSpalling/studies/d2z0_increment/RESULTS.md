# D2z-0 RESULTS — retire the half-cell removal increment

Criteria: `CRITERION.md`, sha256 `ce81f371962f6588…`, hashed read-only with zero
scored-leg output on disk. Metric imported **unchanged** from
`studies/d2j0b_plane/analyze.py`. Frozen studies read by path, never re-run.

Legs: `H2` (A_f02's argv, key absent, 1.8 min), `F2` (+ `full_cell`, 2 mm,
1.8 min), `F1` (+ `full_cell`, 1 mm, 15.7 min). Full output in `analysis.txt`.

---

## Verdict

| criterion | result |
|---|---|
| **(a)** rate change +0.2…+0.6 % per 50 s block | **FAIL** |
| **(b)** mesh-gap shift < 0.2 pts per block | **FAIL** |
| **(c)** plumbing: `firings_per_void` → 1.000, one row per void | **PASS** |
| **(d)** identity: `H2` ≡ frozen `A_f02`; tree 2639/2639 | **PASS** |

| prediction | outcome |
|---|---|
| **P1** rate change +0.2…+0.6 % at both meshes | **REFUTED** (−1.73 … +4.35) |
| **P2** mesh-gap shift < 0.2 pts in every block | **REFUTED** (max 3.42) |
| **P3** `firings_per_void` → exactly 1.000 | **HELD** |
| **P4** `H2` byte-identical to frozen `A_f02` | **HELD** |

**The mechanism works exactly as designed. Two of the four criteria are written
at a time resolution the estimator cannot support, and they fail for that
reason, not because `full_cell` does something unexpected.**

---

## (c) and (d) — the mechanism and the default, both clean

| leg | firings | voids | firings/void | `n_voided` | `h/dz` | thermo `fpv` | thermo `ric` |
|---|---|---|---|---|---|---|---|
| `A_f02` (frozen) | 18,713 | 9,562 | 1.9570 | {0, 1} | 0.50851 | absent | absent |
| **`F2`** | 9,606 | 9,606 | **1.0000** | **{1}** | **1.00000** | 1.00000 | 1.00000 |
| `A_f02_1mm` (frozen) | 135,859 | 69,240 | 1.9621 | {0, 1} | 0.50846 | absent | absent |
| **`F1`** | 69,702 | 69,702 | **1.0000** | **{1}** | **1.00000** | 1.00000 | 1.00000 |

* The **phi tie holds**: `n_voided == 1` on every one of 79,308 `full_cell`
  firings across both meshes. `h = dz` lands exactly on a cell boundary, which
  `Removal.H:145-148` warns can "miss or add one"; the in-code assertion never
  fired and the logs confirm it independently.
* `H2` is **byte-identical** to frozen `d2l_scan/output/A_f02` — 24 `Cell_D`
  files, `thermo.dat` and the removal log. Both ran the full 250 s, so the
  comparison is the whole run, not a prefix.
* **Tree identity 2639/2639** with the key unset, run after the last source edit
  under `guarded_post_all.sh`, which confirms the binary did not move mid-sweep.
* New regression `unit/removal_increment`: **14/14 PASS**.

⚠ **Baseline correction.** `ACTIVE_STEP.md` gives the origin as **1.977**. That
is **Meier `A`'s** value, reproduced here exactly. On the **arbiter** this packet
actually runs it is **1.957** (2 mm) and **1.962** (1 mm). The bar is the
destination, so nothing changes; the number is corrected, not substituted.

## (a) and (b) — why they fail, and what the run actually did

Scored, per 50 s block:

| mesh | 50–100 | 100–150 | 150–200 | bar |
|---|---|---|---|---|
| 2 mm | −1.73 % | **+4.35 %** | −1.38 % | +0.2…+0.6 % |
| 1 mm | +0.81 % | +0.88 % | −0.37 % | +0.2…+0.6 % |

Mesh-gap shift: **+2.59 / −3.42 / +1.00 pts** against a < 0.2 pt bar.

**The same estimator on one 150 s window instead of three 50 s blocks:**

| mesh | half_cell | full_cell | change | cells voided |
|---|---|---|---|---|
| 2 mm | 1.60215 | 1.60780 | **+0.352 %** | 4,024 → 4,053 (+0.72 %) |
| 1 mm | 1.60685 | 1.61388 | **+0.438 %** | 33,113 → 33,219 (+0.32 %) |

Packet point prediction: **+0.36 %**. Aggregate mesh gap: +0.29 % → +0.38 %,
a shift of **+0.085 pts** — inside (b)'s bar.

⚠ **These aggregate figures were computed AFTER the block result was known.**
They are diagnostic, not a re-score. (a) and (b) fail as hashed and P1/P2 are
refuted as hashed. The aggregate explains the failure; it does not overturn it.

**Why the blocks scatter.** `full_cell` changes *when* each cell leaves, so the
removal sequence is reshuffled across block boundaries even though almost the
same number of cells leaves overall. Cells voided per block at 2 mm are
**+0.08 / +2.03 / +0.07 %** — already less extreme than the rate metric's
−1.73 / +4.35 / −1.38 %, because `col_rate` fits per-column recession inside the
window and is sensitive to which side of a boundary a flip falls on. On 25 s
sub-blocks the scatter is larger still and sign-inconsistent:

```
  50-75  +1.82 %    100-125  +0.73 %    150-175  -3.30 %
  75-100 -1.64 %    125-150  +3.40 %    175-200  +3.60 %
```

Any implementation of `full_cell` shifts removal timing — that is the entire
point of the change — so this reshuffling is unavoidable and the ±0.2–0.6 %
per-block band is **not achievable by any correct implementation** on this
arbiter. The bar states a true physical expectation at a window length where
run-to-run phase noise is an order of magnitude larger than the effect.

Note the 1 mm blocks are much tighter than the 2 mm ones (+0.81 / +0.88 / −0.37
against −1.73 / +4.35 / −1.38): four times as many columns in the disc and twice
the flip rate average the phase noise down, as one would expect if it is noise.

---

## An observation the packet's mechanism story does not cover

The packet excludes `prescribed_T` because "16d's fix already applies and the
flat half-cell does not arise". In the new regression's `prescribed_T` control
the surface is genuinely hot (≈1100 K) against a 655 K top cell, so the lerp is
a real gradient — and the scan **still** returns `h_scan = 0.98437 mm ≈ dz/2 −
a_f`, the same landing point as computed-T.

That points at the **gate location** rather than the flat lerp: `Sp_top >= 1` is
tested at the top cell's **centre** (`z = dz/2`) while removal depth is measured
from the **face**, so a marginal firing lands at ≈ `dz/2` whichever mode supplies
`T_face`. Reported only — the guardrail forbids touching the scan, and
`prescribed_T` is excluded as specified. It belongs to whoever revisits 16d.

## Stage C — not run

Re-baselining the four D2z references under `full_cell` (~65 min) is
**the user's call** and was not taken. Every frozen comparator stands, and the
default is unchanged: `half_cell` is byte-identical to the key unset.

## What this does and does not buy D2z

**Does:** a column now fires once per cell, so D2z's void-event filter and its
half-cell warning both disappear — one of its three estimator traps is gone.
`removal_increment_cells` and `firings_per_void` give it a live check.

**Does not:** it changes nothing physical. The aggregate rate moves +0.35 % to
+0.44 % and the mesh gap moves by under a tenth of a point. This was always a
hygiene packet with a tiny predicted effect, and that is what it delivered.
