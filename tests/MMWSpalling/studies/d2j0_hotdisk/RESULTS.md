# D2j-0 results — the hot-disk edge test (2026-09-21)

**Verdict: P2 FAILS, so per `PREDICTIONS.md` the answer is stop.** A stalled annulus does form
at the patch edge at every mesh, but **`f_stall` does not grow as dx falls** (7.0 / 6.7 / 7.4 %
at 4 / 2 / 1 mm) — it is the same 28–30 mm ring at every resolution. `w_edge` is **converging**
(+119 % from 4 to 2 mm, then +15.3 % from 2 to 1 mm), so `S05` does not run.

**The bare lateral edge does not reproduce the D2i failure.** It produces a mesh-dependent
*offset* — the disk-mean recession falls 10.5 % per halving at the finest step — but D2i's gap
**grows in time without bound** (−2.1 → −4.0 → −16.3 → −36.1 % over successive 100 s blocks).
An open loop with fixed h and a static patch has no mechanism that grows with time, and none
appeared. The edge supplies a standing offset; something else supplies the runaway.

**The next packet must isolate the stand-off feedback**, not the removal law: the D2h open-loop
pair, whose `nozzle_collision_radius = 0.028` fix is already in `d2h_openloop/run.py` and which
was never rerun. **Do not start the continuous-front redesign on this evidence.**

## 0. Predictions and binary

- `PREDICTIONS.md`, sha256 `db007f1bf9a397c7b1245b373d96dd3b9f446e9d85f841cb19f11301bff16084`,
  mtime **2026-09-21 09:35:52**, read-only, written **before any run** (the output directory did
  not exist at that moment). Not amended; every metric, window and bin below is as fixed there.
- Binary `bin/mmwspalling-3d-g++` sha256
  `6d047b50bdcde619952aff907d88f89ed0cae8f9d3ed645b43ba0ac2cd647005`, verified at the start and
  at the end. **No source edit, no build, no relink.**
- Standalone study: `run.py` and `analyze.py` import none of the d2c→d2i chain.

## 1. Runs

| run | set | dz | cells | dt | K/step | V0 [m³] | t_end | wall | max abs ledger_err |
|---|---|---|---|---|---|---|---|---|---|
| `S4` | sharp | 4 mm | 20²×25 | 4 ms | 0.24 | 6.4e-08 | 149.9 s | 1.7 min | 6.6e-15 |
| `S2` | sharp | 2 mm | 40²×50 | 4 ms | 0.48 | 8e-09 | 149.9 s | 20.9 min | 7.9e-15 |
| `S1` | sharp | 1 mm | 80²×100 | 4 ms | 0.96 | 1e-09 | 149.9 s | 90.7 min | 8.7e-15 |
| `T4` | tapered | 4 mm | 20²×25 | 4 ms | 0.24 | 6.4e-08 | 149.9 s | 0.6 min | 7.0e-15 |
| `T2` | tapered | 2 mm | 40²×50 | 4 ms | 0.48 | 8e-09 | 149.9 s | 3.2 min | 9.4e-15 |
| `T1` | tapered | 1 mm | 80²×100 | 4 ms | 0.96 | 1e-09 | 149.9 s | 14.3 min | 9.6e-15 |
| `S05` | — | — | — | — | — | — | **not run** | — | — |

`dt` is 4 ms at every mesh, far below both the explicit limit (0.24 s at 1 mm) and the overshoot
rule. `weibull.V0 = V_cell = dz³`, the scored convention, so it moves with the mesh. The energy
ledger closes to ≤ 9.6e-15 everywhere. Wall times are under `caffeinate`; the three jobs of each
set shared the machine, so `S1`'s 90.7 min is contention, not cost.

**Two deviations from the packet's literal text, both recorded before the runs were scored:**

1. **`amr.blocking_factor = 1` at every mesh.** The packet's 0.10 m box is 25 cells deep at
   4 mm, which is not divisible by the input's `blocking_factor = 2`, and `S4` aborted on it.
   Setting 1 uniformly was preferred to varying the key across a mesh sweep or moving the
   packet's geometry. **Verified harmless:** `S2` at `blocking_factor = 1` is **byte-identical**
   to the same run at 2 (`thermo.dat` and removal log both; kept as `output/S2_bf2`).
2. **Tapered patch radius 0.0325 m, not 0.030.** The parser aborts on h ≤ 0 inside the patch, so
   the patch edge is placed exactly where the cosine ramp reaches zero; otherwise the ramp would
   be truncated at 30 mm and the "tapered" set would carry a half-height jump. Minimum in-patch
   h is 51.8 / 3.34 / 0.230 W/m²K at 4 / 2 / 1 mm — positive at every mesh, checked by
   `run.py --list` before launch.

## 2. Sharp set (h = 700 inside r = 30 mm; h drops to 0 across one cell)

| run | dz | v_c [m/h] | r_half [mm] | r(0.9) | r(0.1) | **w_edge** [mm] | w_edge/dz | f_stall | r_stall_inner 50/100/150 s |
|---|---|---|---|---|---|---|---|---|---|
| `S4` | 4 mm | 1.582 | 27.87 | 23.16 | 28.86 | **5.70** | 1.43 | 7.0 % | 28 / 28 / 28 mm |
| `S2` | 2 mm | 1.580 | 26.78 | 16.77 | 29.29 | **12.51** | 6.26 | 6.7 % | none / none / 28 mm |
| `S1` | 1 mm | 1.584 | 25.02 | 14.48 | 28.91 | **14.43** | 14.43 | 7.4 % | none / none / 28 mm |

## 3. Tapered set (raised cosine, 700 at 27.5 mm to 0 at 32.5 mm)

| run | dz | v_c [m/h] | r_half [mm] | r(0.9) | r(0.1) | **w_edge** [mm] | w_edge/dz | f_stall | r_stall_inner 50/100/150 s |
|---|---|---|---|---|---|---|---|---|---|
| `T4` | 4 mm | 1.581 | 28.15 | 21.56 | 30.20 | **8.64** | 2.16 | 4.7 % | 28 / none / none |
| `T2` | 2 mm | 1.581 | 26.27 | 16.40 | 28.99 | **12.59** | 6.29 | 7.3 % | none / none / 28 mm |
| `T1` | 1 mm | 1.583 | 24.66 | 14.38 | 28.87 | **14.49** | 14.49 | 8.2 % | none / none / 28 mm |

`v(r)/v_c` by 2 mm annulus is in `analysis.md`; `d2j0_vprofile.png` plots both panels. Empty
bins at 4 mm are real — a 2 mm annulus can contain no 4 mm cell centre.

## 4. Predictions scored

- **P1 — HOLDS.** `v_c` = 1.582 / 1.580 / 1.584 m/h across 4 / 2 / 1 mm, spread **0.3 %**
  (limit 5 %), against a 1.632 m/h hand estimate. The 1-D energy balance is mesh-independent,
  the configuration is sound, and P2–P4 are scoreable. (Caveat: at 4 mm only one column has
  r < 6 mm, so `v_c` there rests on a single column — it agrees anyway.)
- **P2 — FAILS.** A stalled annulus exists at every mesh, but `f_stall` is 7.0 / 6.7 / 7.4 %
  (sharp) and 4.7 / 7.3 / 8.2 % (tapered): **it does not grow as dx falls.** It is the outermost
  one or two bins, 28–30 mm, at every resolution — the same physical ring, not a growing front.
  This is the pre-registered fail condition, and it decides the packet.
- **P3 — FAILS.** `w_edge` 5.70 → 12.51 → 14.43 mm: +119 % from 4 to 2 mm, then **+15.3 %** from
  2 to 1 mm, below the 20 % threshold. The transition width is converging on ≈ 15 mm, and 4 mm
  is simply outside the asymptotic range (as `PREDICTIONS.md` anticipated: the ablation layer
  α/v = 1.52 mm spans 0.38 cells at 4 mm). **`S05` therefore does not run.**
- **P4 — nominally holds, but the evidence is too weak to lean on, and it is reported that way.**
  The coded test scores "appeared where there was none" as inward motion, which is what happens
  at 1 mm (none at 100 s → 28 mm at 150 s) against a static 28 mm at 4 mm. **The raw numbers
  show nothing walking inward at any mesh**: the ring sits at 28–30 mm wherever it appears, and
  at `T4` it even disappears between 50 and 150 s. `r_stall_inner` as defined is too coarse at
  these bin widths to decide P4. Do not cite P4 in either direction.

## 5. The sharp-vs-tapered discriminator — neither branch

The discriminator anticipated that the sharp set might scale with dx while the tapered set did
not, which would have pointed at the gas closure (`jet_r_reach`). **What happened is that the
two sets are nearly identical:**

| | 4 mm | 2 mm | 1 mm | 2 → 1 mm |
|---|---|---|---|---|
| sharp `w_edge` | 5.70 | 12.51 | 14.43 mm | +15.3 % |
| tapered `w_edge` | 8.64 | 12.59 | 14.49 mm | +15.1 % |

At 2 and 1 mm the profiles agree to within ~0.005 in `v(r)/v_c` at every radius. **Smearing the
h edge over 5 mm changes essentially nothing**, so the recession edge is not set by the
sharpness of h. It is set by lateral conduction into the cold rock outside the disk, which is
physical and is resolved alike at 2 and 1 mm. That removes the gas-closure branch as well as the
h-edge branch: neither is the D2i driver.

## 6. The decision, and the one number worth carrying

**P2 fails → stop**, exactly as pre-registered. The next packet isolates the **stand-off
feedback** — rerun the D2h open-loop pair (`OL_2mm` / `OL_1mm`), which aborted on a collision
guard and was never rerun; the fix is already in its `run.py`. The continuous-front redesign in
`Claude_markdowns/2026-09-21b.md` is **not** to be started on this evidence.

**But the edge is not innocent, and this is the number D2j-a or its successor should carry.**
Mean recession over the disk (r < 30 mm), 50–150 s:

| set | 4 mm | 2 mm | 1 mm | 4 → 2 mm | 2 → 1 mm |
|---|---|---|---|---|---|
| sharp | 1.3126 | 1.1331 | 1.0142 m/h | −13.7 % | **−10.5 %** |
| tapered | 1.3428 | 1.1082 | 0.9986 m/h | −17.5 % | **−9.9 %** |

So the open-loop edge **does** cost about 10 % of the disk-mean rate per halving at the finest
step, in the **same direction** as D2i's ring lag, and it is not converged. What it does not do
is **grow with time** — with h fixed and the patch static, the offset is standing. D2i's gap
grows block on block. **A static ~10 % mesh offset cannot produce a −36 % and rising gap**,
which is why the edge is not the cause even though it is a real defect worth fixing later.

Read together with D2e (clip), D2h (idle clock) and D2i: three mesh seeds found and removed, a
fourth measured here and shown to be bounded, and the runaway still unexplained. The remaining
candidate that has never been tested in isolation is the feedback loop itself.

---
**Reviewer addendum (2026-09-21, evening, `studies/d2j0b_plane/RESULTS.md`).** Re-run to 250 s
with a descending plane, the *same* sharp disk (C0, no exclusion) does **not** give a standing
offset: the 1 mm/2 mm disk-mean gap grows −10.4 → −15.3 → −20.1 % over 50 s blocks and the
stalled ring walks inward after 200 s (26–28 mm stops at 216 s / 239 s). The "bounded" reading in
§6 above was a property of the 150 s window. The D2i cascade is reproduced on the disk with an
emulated exclusion and no feet or jet; heating the exposed wall does not change it. See D2j-0b.
