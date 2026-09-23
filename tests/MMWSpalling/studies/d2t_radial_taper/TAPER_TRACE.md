# D2t Goal item 2 — the compiled-taper thermal trace: criteria, pre-registered

**This is a STOP gate.** If the compiled profile is not the intended one, the
packet stops here. D2o-0 lost a whole packet to a parser expression that did
not mean what it read like (`a/max(x, num)` evaluates to 0), and a numpy
re-evaluation of the formula cannot catch that, because the bug lives *inside
the AMReX parser*. `flame_h` is a plain `std::vector` (`:5825`) and is never
written to a plotfile, `thermo.dat` or any csv, so the compiled `h(r)` cannot
be read out directly. The temperature field is the only witness.

## 0. Disclosure

`TR_below` and `TR_above` had **already been run** when this file was written
and hashed (they take ~10 s each). What was inspected before hashing: that both
exited 0, that a `00010cell` plotfile exists, and that both
`_removal_events.csv` contain **only a header row** — i.e. zero spall events,
which is a precondition of the method, not a result. **No temperature has been
read.** The estimator below was written before any temperature was loaded.

## 1. Method

Both legs are 10 steps of 16 ms at 2 mm, `stop_time = 0.16 s`, stopped long
before the first spall (the frozen logs put that at 4.58 s). Before any
removal every column has `pin_T = 0`, so `PinRule` returns 0 and **every
column is on the face form** — the flux is the Newton surface balance in `h`,
not a pin. Every column starts at 293.15 K, so the temperature rise after
0.16 s is a direct readout of `h(r)`.

Measured per column: `E = rho*c*dz * sum_k (T(i,j,k) - 293.15)` [J/m²], the
column's enthalpy gain, which equals the time-integral of `q_in`. Lateral
conduction over 0.16 s reaches `sqrt(alpha*t) = 0.33 mm`, far below the
taper's 16 mm width, so columns are independent to well under a percent.

`E` is inverted to `h` by a 1-D discrete reference solved in Python with the
**same** `dz`, `dt`, face-form balance (`G = 2k/dz`, radiative + ambient
losses) and explicit update as the kernel, and bisected on `h`. **This
re-implements the thermal response, never the taper**: the taper enters only
through the compiled run. The inversion carries its own control — see T1.

Binned in **2 mm fixed physical radial bins**, the same convention as the
shape metric.

## 2. Criteria — all evaluated before any campaign leg launches

| # | test | bar |
|---|---|---|
| **T1** | **self-calibration**: inverted `h` averaged over `r < r_flat`, where the intended value is known exactly | `700 ± 3 %`. **Failing T1 voids the estimator, not the taper** — re-derive the inversion, do not proceed and do not reinterpret. |
| **T2** | flat inside `r_flat` | `max abs(h/700 - 1) <= 0.02` over `r < 18 mm` |
| **T3** | monotone in `r` | `h` non-increasing bin to bin out to `radius`, tolerance `0.005 x` peak |
| **T4** | **is it the intended profile** | `max abs(h(r)/700 - taper(r)) <= 0.02` over every bin with `r <= radius` |
| **T5** | at the floor by `radius` | `h <= 0.01 x` peak for every bin with `r >= r_taper` |
| **T6** | **no rise beyond `radius`** (N2: the cliff moved, not removed) | `max dT <= 1e-3 K` over columns with `r > radius` |
| **T7** | **branch association**: `TR_above/TR_below` in `h` | `142/700 = 0.202857 ± 2 %` at every bin with `r <= r_flat`. This is what proves the taper multiplies the whole `if(s>0,700,142)` and is not bound to one branch. |
| **T8** | zero spall events in both legs | both csvs header-only (**already confirmed, see §0**) |

**T2–T8 are STOP conditions.** Any failure stops the packet and is reported as
N1 or N2 from the packet's inert table.

## 3. What this cannot show

The trace sees `h(r)` on a flat, unremoved surface at `t = 0`. It says nothing
about the taper's *consequences* — the lagging band, the latch, or shape
convergence. Those are the campaign's job and are scored in `CRITERION.md`.
