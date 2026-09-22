# AMR surface-band gate — results

Binary `bin/mmwspalling-3d-g++-amr` (POSTFIX=3d-g++-amr build of the working tree on
2026-09-22; the D2k campaign binary was not touched). Case: D2j-0b C1 hot disk.
Criterion fixed before the first run in `CRITERION.md`.

## 0. Single-level check (code is key-off exact)

`U2_250` (uniform 2 mm, max_level 0, new binary) `_removal_events.csv` is **byte-identical**
to `d2j0b_plane/output/C1_2mm_removal_events.csv` (old binary, 11 148 rows).

## 1. G1 identity — PASS

`F1_60` (2 mm base, level 1 covering the whole domain) vs `U1_60` (uniform 1 mm), 60 s:
`_removal_events.csv` **byte-identical** (26 842 rows, total h_applied 13.49879 m over all
columns). Wall 6.9 min vs 5.6 min (the coarse level is extra work here). No regrid repair
messages (the level never changes shape).

## 2. G2 band — PASS

`B1_250` (2 mm base, level 1 = 20 mm below / 4 mm above the surface, regrid every 2 s,
117 regrid repairs of ~1.3–1.4 k cells each, no coverage abort) vs `U1_250` (uniform 1 mm),
same binary, 250 s. D2j-0b metrics unchanged (`d2j0b_plane/analyze.py`):

| block [s] | U1_250 [m/h] | B1_250 [m/h] | gap |
|---|---|---|---|
| 50–100 | 0.996 | 0.996 | +0.0 % |
| 100–150 | 0.690 | 0.686 | −0.5 % |
| 150–200 | 0.414 | 0.408 | −1.4 % |

v_c 1.459 vs 1.459 m/h. Band t_stop / t_cross (median, s): 18–20 mm 162/177 vs 160/174;
20–22 139/153 vs 138/153; 22–24 118/130 vs 118/130; 24–30 mm identical. r_front 26/22/18 mm
at 100/150/200 s in both. The residual sits in the outermost bands late in the run, where
the 2 mm level carries the deep conduction — consistent with the band-edge approximation.

Cells: U1_250 960 000 (16 boxes); B1_250 120 000 on level 0 + 174 160 on level 1 (339
boxes) = 3.3× fewer. **Wall time 19.9 vs 20.0 min (ratio 1.00)** with both runs sharing the
host with the D2k campaign: the cell saving is eaten by per-step column gathers
(ColumnTopSolid, CheckSurfaceCovered, the removal pass are per-column, not per-cell),
two-level FillPatch, and 339 small level-1 boxes. Speed is the next item, not correctness:
larger level-1 `max_grid_size` / `blocking_factor`, `grid_eff`, and moving the coverage
check to regrid steps only.

## 3. Status

AMR (max_level = 1) on the thermal + sp_weibull removal path is correct on this case:
identity gate byte-identical, band gate within 1.4 %. Not yet exercised under AMR: feet
descent, jet enthalpy closure, foot body clearance, energy ledger (aborts by design),
mechanics-on paths, max_level ≥ 2.
