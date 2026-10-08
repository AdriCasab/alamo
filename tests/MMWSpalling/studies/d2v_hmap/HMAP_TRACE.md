# D2v Goal item 3 — the compiled h map, BOTH branches: criteria, pre-registered

**This is a STOP gate.** `flame_h` is a plain `std::vector` (`:5825`), written
to no plotfile, no `thermo.dat` and no csv, so the compiled `h(r)` has no direct
readout. A numpy re-evaluation of the expression is only a filter: the D2o-0 bug
lived **inside the AMReX parser**. The temperature field is the only witness.

## 0. Disclosure

`TF` and `TP` had **already been run** when this file was written and hashed
(~10 s each). Inspected before hashing: both exited 0; both wrote a `00030cell`
plotfile at t = 0.48 s; both `_removal_events.csv` have **zero** data rows;
`TF_jet_profile.csv` has 2550 rows and **`TP_jet_profile.csv` has only a
header** — i.e. no column marches in `TP`, which is the design (every column is
above the plane) and a precondition of the method, not a result. **No
temperature has been read.**

## 1. Why two legs

**A feet-descent run can never execute the wall branch.** The nozzle rides
`foot_standoff = 50 mm` above the rock, so at t = 0 every column has
`s = +50 mm > 0`. A short feet run therefore tests the floor branch **only**
(N7).

* **`TF`** — feet descent, `M1` keys, `s = +50 mm` everywhere ⇒ the **floor**
  branch, `h = [h_imp(r,s)·w1 + 142·(1−w1)]·w2 + 1e-3·(1−w2)`.
* **`TP`** — prescribed descent with the plane at `z = 0.39 m`, **below** the
  original surface at 0.40 m, so `s = −10 mm` everywhere ⇒ the **wall** branch,
  `h = 142·w2 + 1e-3·(1−w2)`. Needs `jet_above_h = 1`, all feet-only keys
  stripped (they are `query_required` only under feet, and ParmParse
  strict-checks unused keys), `nozzle_collision_radius = 0` (it defaults to
  `radius` = 0.2 m without feet and would abort immediately), and
  **`jet_free_surface` stripped explicitly** — it does not abort under
  prescribed descent, it silently treats every column as free surface, because
  its test is `r > foot_r_outer` and `foot_r_outer` defaults to 0.0 and is never
  parsed outside the feet branch.

## 2. Method

30 steps of 16 ms from a uniform 293.15 K, stopped at 0.48 s — long before the
first spall, which the frozen log puts at **3.184 s**. Before any removal every
column has `pin_T = 0`, so `PinRule` returns 0 and **every column is on the face
form**: the flux is the Newton surface balance in `h`, not a pin.

Measured per column: `E = rho·c·dz·sum_k (T − 293.15)`, the column enthalpy
gain. Inverted to `h` by a 1-D discrete reference solved in Python with the
**same** `dz`, `dt`, face-form balance and explicit update as the kernel, and
bisected on `h`. **This re-implements the thermal response, never the map.**

**`T_gas` differs between the two legs and this is the correction D2t's recipe
needed:**

* `TF`: the march makes `T_gas` **radius-dependent**, so the raw rise is not
  `h(r)`. Normalise by `T_gas(r) − T_0` read **per bin** from
  `TF_jet_profile.csv` (`r_lo`, `T_gas`), at the profile row nearest the
  plotfile time.
* `TP`: no column marches, so `T_gas` is the constant `jet_T_above = 1690 K`.

Binned in **2 mm fixed physical radial bins**.

## 3. Criteria — all evaluated before any campaign leg launches

Peak reference: `h_imp` is clamped to `r >= 18.75 mm`, so the floor branch is
flat at **1449 W/m²K** inside that radius at `s = 50 mm` — a **known value**,
which is what makes T1 a self-calibration rather than an assumption.

| # | leg | test | bar |
|---|---|---|---|
| **T1** | TF | **self-calibration**: inverted `h` over `r < 18.75 mm`, intended value known exactly | **1449 ± 5 %**. Failing T1 voids the **estimator**, not the map — re-derive, do not reinterpret. |
| **T2** | TF | monotone decreasing in `r` out to the domain corner | no bin-to-bin rise above `0.01 ×` peak |
| **T3** | TF | **is it the intended map** | `max |h_meas(r) − h_intended(r)| <= 0.05 ×` peak over every bin |
| **T4** | TF | the annulus plateau is reached | `h_meas` over `56 <= r <= 60 mm` within **±15 %** of `h_ann = 142` |
| **T5** | TF | at the floor by `r_out` | `h_meas <= 0.01 ×` peak for every bin with `r >= 76 mm` |
| **T6** | TF | **THE FUNNEL WITNESS — no heating beyond the mouth** | `max dT <= 0.05 K` over all columns with `r > 76 mm`. The **stale** map gives `h = 362 W/m²K` there and would show tens of K. This is the single number that says the funnel's cause is removed. |
| **T7** | TP | **the wall branch is FLAT `h_ann`** | `h_meas` over `r < 56 mm` within **±15 %** of 142 |
| **T8** | TP | **the wall branch is NOT blended** — the rejected form | `h_meas` at `40 <= r <= 56 mm` must be `<= 1.3 × 142`. A blend would read **356–553 W/m²K** there. |
| **T9** | TP | wall branch reaches the floor by `r_out`, no rise beyond | as T5/T6 |
| **T10** | both | zero spall events | confirmed in §0 |

**T2–T10 are STOP conditions.** Any failure stops the packet and is reported as
N1 (map ≠ intended) or N3 (blanket).

## 4. What this cannot show

The trace sees `h(r)` on a flat, unremoved surface at t ≈ 0. It says nothing
about the funnel's *removal*, the wall, or the latch. Those are the campaign's
job and are scored in `CRITERION.md`.
