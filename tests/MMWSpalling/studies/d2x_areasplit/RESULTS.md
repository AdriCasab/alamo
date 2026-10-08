# D2x — close `A`'s wall band without losing its rate

**Answer: no. `A` stands.** Every side-face treatment closes the wall band and
**overshoots it**, and every one costs rate. No leg meets the four bars.

| | `M0` | **`A`** | `B` | `M1` | `C` | `D` | `E` |
|---|---|---|---|---|---|---|---|
| (a) @ descent 150 | 28.40 | **2.27** | 19.93 | 2.73 | 5.51 | 7.56 | 14.23 |
| (a) @ descent 180 | 25.12 | **2.51** | 17.36 | 2.55 | 5.38 | 7.86 | 14.37 |
| Ø @100/125 mm | 99.5/95.3 | **84.4/82.7** | 93.2/90.1 | 87.9/86.4 | 91.9/90.7 | 95.1/94.0 | 101.7/100.2 |
| descent [m/h] | 1.50 | **1.25** | 3.09 | 3.14 | 2.09 | 3.80 | 2.11 |
| pit lead @150 | 108 | 128 | 58 | 58 | **82** | **32** | 78 |
| volume, descent basis | 1.68× | 1.24× | 1.41× | 1.11× | 1.16× | **1.13×** | 1.34× |

## 1. The decisive result

The packet asked whether the wall band's last 3–6 mm could be closed **without
giving back** `A`'s shape, rate, pit lead and volume. **It cannot, by any of the
three treatments.**

`A` is short of the 87 mm lower bound by **2.6 and 4.3 mm**. The treatments do
not close that gap — they **sail past it**: `C` reaches 91.9/90.7, `D`
95.1/94.0, `E` 101.7/100.2. The wall band is not delicately restored; it is
over-driven. And the cost is the rate: `A` 1.25 → `C` 2.09, `E` 2.11, `D` 3.80,
none inside 1.3–1.6.

**All three legs FAIL the four bars.** `C` and `D` pass bar 4 only.

## 2. ⚠ But two of `A`'s four defects ARE fixed — this is a trade, not a null

`C` and `D` **improve** the two defects a wall treatment was never expected to
touch:

* **pit lead** 128 mm → **82** (`C`) and **32** (`D`);
* **volume** 1.24× → **1.16×** (`C`) and **1.13×** (`D`) on the descent basis.

So the side-face path buys pit lead and volume at the price of shape and rate.
**That is a real coupling, and it says the four defects are not independent.**

## 3. P2 HELD — the derived number that replaces `f`

**`side_w_mean` = 0.6670 on `D` and 0.6782 on `E`**, inside the pre-registered
0.55–0.75, and consistent with risers of n ≈ 2–3 cells. **The honest per-edge
area is 3.3× the `f = 0.2` it replaces**, and `f = 0.2` is below the 0.4142 that
a one-cell drop already forces. `f` was not merely mis-tuned — it was outside
the range the geometry permits.

⚠ **And the honest area does not help.** `D` (derived area, all faces) is the
*worst* leg on rate (3.80 m/h) and second-worst on shape. **Deriving `f` was one
of the two live readings from D2w's review; this packet refutes it as a route to
Meier's hole.** The other reading — that `f` was masking a missing wall boundary
condition — is not supported either, since `C` (the band-only fudge) also
overshoots.

## 4. The ring witness, per ring and never averaged

`A`'s collapse in 40–50 mm is real (−29 / −59 / −87 / −96 % of `M0`). Every
treatment **over-corrects the whole flank, not just the collapsed rings**:

| ring | `A` | `C` | `D` | `E` |
|---|---|---|---|---|
| 30–32 | −3.6 % | **+14.5** | **+53.6** | **+14.7** |
| 36–38 | −4.0 % | **+34.6** | **+83.2** | **+37.0** |
| 42–44 | −59.1 % | **+56.0** | **+115.0** | **+59.2** |
| 44–46 | −86.6 % | +24.8 | **+128.3** | **+68.9** |

**The side term cannot be aimed at the collapsed annulus alone.** Even `C` and
`E`, gated to r ≥ 38 mm, lift the 30–38 mm rings by 14–37 % — because the gate
selects *columns*, and a column inside the gate still has risers facing inward.
That is why the wall band overshoots.

## 5. Gates and witnesses

* **`C0` PASSES**: with both new keys at their defaults, **62,501 event rows and
  2,083 thermo rows byte-identical to `M1`**. ⚠ Compared at 199.9 s, not 200.0:
  `C0`'s `stop_time` lands exactly on a step, so its final step at
  t = 199.99999999997618 is never produced while `M1` (which ran to 401 s) has
  it. `C0` is a strict **prefix** of `M1` — zero differing rows anywhere both
  have data. **This witness is mandatory and the tree gate cannot replace it**:
  no input in the 2639-file set turns `side_face_flux` on.
* **Tree identity 2639/2639**, zero differing lines, run **twice** (before and
  after the diagnostic fix), binary pinned each time.
* **`side_P_face_above` exactly 0.000e+00 at max AND min on `C`, `D`, `E`** — no
  above-plane leakage.
* **`ledger_err` ≤ 1.24e-14** on every leg; the identity
  `jet_P_face = tops + sides + side_loss` holds to **±0.000 W on all seven cells**.
* `r_min` reaches the block: `C` has **−33 %** faces vs `M1` (bar −30 %).
* ⚠ `area_exact` witness **marginal**: `D`'s `side_P_face` is **1.90×** `M1`'s
  against a "≥ 2×" bar. Not inert — `side_w_mean` = 0.667 on 5858 faces proves
  the law is live — but the bar as written is missed, and it is reported as
  missed rather than reinterpreted.

## 6. ⚠ A bar that its own control cannot meet

**Bar 3 requires descent in 1.3–1.6 m/h, and `A` reads 1.251.** The bar was set
from `A`'s **1.296**, which is the slope over the **time** window 150–400 s;
criterion (c) defines the rate as the slope over the **100–150 mm descent
stretch**, which gives **1.251**. Same leg, same data, two definitions,
0.045 m/h apart — and the band edge falls between them. **The control fails a
bar written from the control.** Reported, not adjusted: no leg was near 1.3–1.6
anyway, so nothing turns on it.

## 7. Pre-registration outcomes

| | prediction | outcome |
|---|---|---|
| **P1** | `D` rate > 4.0 (refuted < 3.5) | **3.80 — neither** |
| **P2** | `side_w_mean` on `D` in 0.55–0.75 | **0.6670 — HELD** |
| **P3** | `C` restores 42–46 mm within 20 % of `M0` | **+56.0 / +24.8 % — partial (overshoots)** |
| **P4** | `C` closes the band to 87 ± 2 | **91.9/90.7 — partial: gains +7.5/+8.0 mm where 2.6/4.3 were needed** |
| **P5** | `C` rate in 1.3–1.9 (refuted > 2.4) | **2.09 — neither** |
| **P6** | `E` between `C` and `D` | **HELD** (2.11 between 2.09 and 3.80, though `E` ≈ `C`) |

**P4 and P5 together were the packet.** P4 lands "partial-overshoot" and P5
"neither", so the clean conclusion the packet hoped for — *the wall band was a
drain, `f` is retired, shape and rate both from derived quantities* — **does not
arrive**. Per the packet's own follow-on: **`A` stands as the reference
configuration, the side-face path is off the critical path, and the open items
are the pit lead (128 mm) and the volume (1.24×)** — except that `C`/`D` show
those two *are* reachable, at a cost.

## 8. What is NOT claimed

* **No Meier rate, diameter, volume or energy claim.** Fig. 8.8 is a **lower
  bound** at ±5 mm; fuel is methane/air at λ = 1.2.
* **One mesh, no convergence claim.**
* `f = 0.2` was **not scanned**; `D`/`E` do not set it at all.
