# D2t — retire the mask edge with a smooth radial taper, and score SHAPE

**Verdict: (a) PASSES.** The project has a mesh-converged shape metric for the
first time. The cliff is **removed**, not heated: under the sharp mask the
innermost never-spalling column sat where `h` was at **full strength**; under
the taper every column under strong heating spalls.

| gate | result |
|---|---|
| **(a) SCORED — shape convergence, `P_1mm` vs `P_2mm`** | core **2.00 mm** (bar 4) **PASS**; band **1 bin** (bar 1 bin) **PASS** |
| **(b) smoothness witness** | **PASS** — zero never-spalled columns inside 26 mm on all four legs |
| **(c) wall convergence, `PA_1mm` vs `PA_2mm`** | core **2.02 mm**, band **1 bin** — **PASS** |
| **(d) core not broken** | **PASS** — −3.2 % to **−5.0 %** of the recomputed reference (bar 5 %) |
| **(e) the latch** | **reported, not gated** — confirmed, see §5 |
| item 2 compiled-taper trace | **PASS T1–T8**, self-calibration **699.99 vs a known 700** |
| tree identity | **2639/2639 byte-identical**, zero differing lines, binary pinned across the sweep |
| regressions | 8 unit + `jet_d2c` + `dev2d` + `s1_smoke` via the sweep, all rc=0; plus **`side_face_flux`, `convective_patch`, `seed_void`, `sp_lefm` PASS** |
| `ledger_err` | ≤ **8.1e-15** on all four legs |
| domain watch | clearance **40–42 mm** (bar 10) |
| no `src/` change | none made — the packet expected none and none was needed |

**Frozen, imported by path, never re-run:** `d2r1_aboveplane_sides/output/`,
`d2j0b_plane/` (arbiter keys), `d2l_scan/output/A_f02*`.

## 1. The cliff is gone — the single clearest number in this packet

| leg | in-patch cols | never-spalled | innermost never-spalled | `h` there |
|---|---|---|---|---|
| frozen `B_1mm` (sharp mask) | 707 | 11 | **29.71 mm** | **1.000 × peak** |
| frozen `T_1mm` (sharp + above-plane lever) | 707 | 11 | 29.71 mm | 1.000 × peak |
| **`P_1mm` (taper)** | 1965 | 1359 | **27.72 mm** | **0.334 × peak** |
| **`PA_1mm` (taper + lever)** | 1965 | 1359 | 27.72 mm | 0.334 × peak |

Under the sharp mask, rock sitting under the **full-strength jet** never spalled
at all. That is the artefact four packets chased. Under the taper, every column
where `h >= 0.5 ×` peak spalls; the columns that never fire are the ones that
are genuinely barely heated, which is physics, not quantisation. The large
never-spalled count is simply the floor region now inside a 50 mm patch.

## 2. (a) — the scored metric

**Core**, `max |z_1mm(r) − z_2mm(r)|` over `r < 18 mm`: **2.00 mm**, one coarse
cell, against a 4 mm bar. The core profile is flat at 108 mm (2 mm) / 110 mm
(1 mm) recession out to r ≈ 11 mm and the 2 mm difference is uniform — it is a
**cell-quantisation offset, not a divergence**.

**Band**, `D(z)` on the fixed 10–100 mm ladder: **adjacent bin at worst**, at 6
of 10 ladder depths identical and at 4 of 10 one bin apart.

## 3. ⚠ The two caveats that qualify the pass

**(i) The `D(z)` metric's resolution EQUALS its bar.** `D = 2·r_centre` lives on
an exact 4 mm lattice, so every nonzero difference is a multiple of one bin. The
band test therefore has only two passing outcomes — identical bin, or adjacent
bin — and **"1 bin" means adjacent bin, not "4.0 mm measured"**. It is a real
pass against a pre-registered bar, but it is the **weakest** part of this
result. *(This also produced a spurious FAIL on the first scoring run: comparing
in millimetres, one bin evaluates to 4.0000000000000036 > 4.0. Fixed by
comparing in integer bins. The bar was not changed — `CRITERION.md` is
chmod 444 and hashed.)*

**(ii) `z(r)` in the band does NOT converge, and `D(z)` hides that by design.**
Max `|z_1mm − z_2mm|` in 18–34 mm is **34.0 mm** for `P` and **16.4 mm** for
`PA`, at r = 23 mm. `D(z)` still passes because the wall there is near-vertical,
so a large *vertical* disagreement maps to a small *radial* one. The packet
chose `D(z)` for exactly this reason and said so in advance — but it must not be
read as "the band converged". **The band converged in diameter, not in depth.**

## 4. (c) — what annulus heating does to a wall, the first measurement

`PA − P` on `D(z)` is **+1 bin, and never negative**, at 7 of 10 ladder depths
at 1 mm. P2 is **held**, though only at the metric's resolution limit. The rate
decomposition is far more informative:

| band 18–34 mm, mesh gap per 50 s block | 50–100 | 100–150 | 150–200 | 200–250 |
|---|---|---|---|---|
| `P` (lever off) | −7.0 % | −24.2 % | **−34.7 %** | **−48.9 %** |
| `PA` (lever on) | −6.4 % | −12.9 % | −11.8 % | **−6.3 %** |

**Without the lever the band's mesh gap grows without bound; with it the gap
roughly halves and stops growing.** `side_P_face` more than doubles (121 → 282 W
at 2 mm, 134 → 299 W at 1 mm). This independently re-confirms D2r-1's mechanism
on a population ~50× larger than D2r-1's 13-column rim, which is what the packet
wanted from this leg.

The **core** is untouched by the lever (`−0.8 / +2.3 / −0.3 / +0.2 %` for `P`
versus `−0.8 / +2.4 / −0.0 / +0.4 %` for `PA`), confirming P3's mechanism: the
side-face branch fires on a solid cell with a **void** lateral neighbour, and a
core column is the *deeper* one, so its neighbours at its own `k_top` are solid.
The extra flux goes to the band column, not the core.

## 5. (e) — the latch, confirmed

Plane crossings per block collapse while the above-plane population saturates:

| leg | 50–100 | 100–150 | 150–200 | 200–250 |
|---|---|---|---|---|
| `P_1mm` | 1525 above, **146** crossings | 1573, 46 | 1605, 32 | 1626, **21** |
| `PA_1mm` | 1503 above, **124** crossings | 1529, 24 | 1534, 5 | 1539, **5** |

Columns cross the plane early and then **stop crossing** — the two states are
separated by the feed rate, so each column latches. This is D2o-0b's
sharp-corner pathology moved from `r` into `s`, and the taper made it
load-bearing by manufacturing a large population at the crossing. **The lever
reduces crossings ~4× by the last block**, i.e. it stabilises the wall rather
than un-latching it.

## 6. Pre-registration outcomes

| | prediction | outcome |
|---|---|---|
| **P1** | above-plane wall, ~159 / 629 columns in 19–34 mm; refuted under 100 at 1 mm | **HELD** — 124/159 at 2 mm, **571/629** at 1 mm |
| **P2** | `PA − P` on `D(z)` is LARGE; refuted under 1 bin | **HELD**, but only *at* 1 bin — the metric cannot show more |
| **P3** | (d) passes, core untouched | **HELD** |
| **P4** | (b) passes | **HELD** |
| **P5** | **declared open — no prediction made** | (a) passes |
| **P6** | the latch may prevent band convergence | **not triggered** — the band converged in `D(z)` |

P5 was recorded as "unknown" on purpose. Inventing a number there is what killed
D2s, and this is the first time the project has run a shape metric at two meshes.

## 7. Two corrections to my own pre-flight

1. **The packet's band count 159 / 629 is exactly right.** My `PREFLIGHT.md`
   §1 recorded **156 / 624** and labelled it a "direct grid count"; it was an
   *area* approximation, `(π/4)(34² − 19²)/dz²`. A direct grid count reproduces
   the packet exactly. **The error was mine, not the packet's.** `PREFLIGHT.md`
   is left as written rather than amended after the fact.
2. The other three PREFLIGHT disagreements stand: recession is **108–110 mm**,
   not ~104; the front sits 21.9–23.6 mm below the plane, not 18–21; and
   `B_1mm` has **110** above-plane columns from 27.50 mm, not 111 from 27.3.

## 8. ⚠ What is NOT claimed

* **No Meier number, no validation claim, no claim the shape is right.** The
  arbiter is a flat plate under a descending disk. This packet shows shape is
  **numerically convergent**; it cannot show it is **correct**.
* **The band is converged in diameter only.** See §3(ii).
* **(d) passes with little margin at 2 mm** (`P_2mm` at **−4.97 %** against a
  5 % bar). The core is genuinely ~3–5 % slower than the frozen sharp-edge
  reference. The plausible mechanism is the drain: the newly-built wall at
  r ≈ 19–21 mm is cold relative to the core's outer bins and conducts heat out
  of them. **It is inside the pre-registered bar, but it is not nothing**, and a
  successor that widens the taper would push on it.
* `side_face_factor = 0.2` remains **underived**. It was not scanned here.
