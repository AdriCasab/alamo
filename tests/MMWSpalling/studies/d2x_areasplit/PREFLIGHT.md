# D2x pre-flight — written BEFORE any campaign leg launched

## 0. Reachability ledger — every leg, every call site

| leg | key diff from `A` | sites that read it | reached because |
|---|---|---|---|
| `C0` | 5 keys, **both new ones at their defaults** | kernel + ledger + march | `side_face_flux=1` opens the side branch; the new keys take the pre-D2x path |
| `C` | `C0` + `r_min=0.038`, no `area_exact` | all three | new radius test inside `col_in_patch` |
| `D` | flux, `above=0.0`, `area_exact=1` | all three | `area_exact && !above` branch |
| `E` | `C`'s band + `area_exact=1`, **no** `side_face_factor` | all three | both new guards |

### ⚠ The packet names TWO call sites. There are THREE.

Goal 5 says "both call sites changed identically — `:1629-1657` and
`:2838-2872`". A third site applies the side factor: **`JetEnthalpyMarch`**,
which bills the jet for what the side faces draw (`side_f_m * a_rob * sN[idx]`).
**Criterion (e)'s own identity is what proves it must change with them**: with
the kernel depositing at weight `w` and the march billing at `f`,
`jet_P_face = patch_P_robin + side_P_face + side_P_loss` would break on `D`/`E`
— a plumbing artefact that would look like a physics result. So the march
carries the same per-edge weights (the weight is folded into its `sN`, and its
factor becomes 1), and the radius gate too. **Reported, not silently expanded.**

### Empirical demonstration — not a grep

**Five guards, each made to fire on the real Meier key set:**

| set | result |
|---|---|
| `side_face_r_min` without `side_face_flux` | **ABORT** "requires surface_patch.side_face_flux = 1" |
| `side_face_area_exact` without `side_face_flux` | **ABORT** same |
| `area_exact=1` **and** an explicit `side_face_factor` | **ABORT** "do not set both" |
| `side_face_r_min=-1` | **ABORT** "must be finite and >= 0" |
| `side_face_area_exact=2` | **ABORT** "must be 0 or 1" |
| `area_exact=1` alone / `r_min` alone | **run, rc = 0** |

⚠ **A bug my own test caught before launch.** The first version read
`pp.contains("side_face_factor")` *after* `query_default` had already named that
key — and `query_default` **adds** the key to the table (the D2c lesson, already
a comment in this file). So `area_exact=1` alone aborted as if `f` had been set.
Fixed by capturing "was it given" before the first `query_default`.

### The regression proves the keys are inert at their defaults

`tests/MMWSpalling/unit/side_face_area/` (new, on the D2k staircase geometry):
**`DEF` is byte-identical to `ON`** across `thermo.dat`, the removal log and all
12 `Level_0/Cell_D` files with both new keys explicitly at their defaults; the
radius gate cuts 130 faces to 6 at `r_min = 5 mm` and to 0 at 50 mm;
`side_w_mean` lands in **[0.4142, 0.8608]** over 8205 rows; and
`side_P_face(AEX) == side_w_mean · side_P_face(ON)` to **1.2e-6**.

**Tree identity with both keys unset: 2639/2639, zero differing lines**, binary
pinned across the sweep.

## 1. The weight, checked against the packet's table

`w(n) = (sqrt(dx² + (n·dz)²) − dx)/(n·dz)`; cubic cells give

| n | 1 | 2 | 3 | 4 | 100 |
|---|---|---|---|---|---|
| computed | **0.4142** | **0.6180** | **0.7208** | **0.7808** | 0.9900 |
| packet | 0.4142 | 0.6180 | 0.7208 | 0.7808 | → 1 |

Exact. **`w >= 0.4142` by construction**, so the `f = 0.2` it replaces is below
the minimum the geometry permits.

## 2. Reference reproductions — all exact

**(a) at equal burner descent**, depths capped 50 mm above the **pads** (which
start at the rock surface, so the cap is `descent − 50 mm`, giving exactly the
packet's depth lists):

| | `M0` | `A` | `B` | `M1` |
|---|---|---|---|---|
| descent 150 mm (10–100) | **28.40** | **2.27** | **19.93** | **2.73** |
| descent 180 mm (10–125) | **25.12** | **2.51** | **17.36** | **2.55** |

All eight reproduce to two decimals. Sampling times match too: `A` reaches
150 mm at **405.3 s** (packet 405.3) and tops out at **183.3 mm** (packet 183.3).

`A`'s defects at descent 180 mm: Ø **84.4** at 100 mm and **82.7** at 125 mm —
both exact. The ring collapse, `A` vs `M0`, 150–400 s:
**−3.6 / −1.9 / −2.4 / −4.0** out to 38 mm, then
**−13.1 / −29.1 / −59.1 / −86.6 / −95.7 %** — every figure exact.

At the **primary** target (descent 150 mm): pit lead **128.0 mm** (packet 127)
and volume **1.24×** on the descent basis (packet 1.24, exact). ⚠ These two move
with the sampling point — at descent 180 mm they read 136 mm and 1.19× — so
this packet reports them at the primary target and says so.

**No disagreements.**

## 3. What this corrects in my own D2w write-up

D2w scored (a) at a **common time**. At 401.6 s `A`'s burner had descended only
148.7 mm while `M1`'s had gone 306.7 mm, so `A`'s 125–200 mm entries were
**unfinished pit compared against `M1`'s finished wall**. At equal descent the
ordering reverses: **`A` 2.27 against `M1` 2.73** — `A` is the better shape, and
gets there at less than half the rate. **Equal-time scoring across legs with
different rates is not a comparison.**
