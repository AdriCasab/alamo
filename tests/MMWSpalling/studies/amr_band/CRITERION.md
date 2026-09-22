# AMR surface-band gate (fixed before the first run; never edit)

**What is tested.** `amr.max_level = 1` on the thermal + sp_weibull removal path with the
finest level owning the surface (`SurfaceOwner`), band tagging (`surface.refine_depth /
refine_height`), the regrid repair of the column-linear level set, and the owner gates on
every surface closure. Case: D2j-0b C1 hot disk (static disk, exclusion plane, no feet, no
jet). Binary: `bin/mmwspalling-3d-g++-amr`.

**G1 identity (F1_60 vs U1_60).** A 2 mm base with level 1 covering the whole domain must
reproduce the uniform 1 mm run: `_removal_events.csv` identical row for row (time, col_i,
col_j, regime, k_top, h_applied, n_voided exact; T_top, Sp_top, h_scan to 1e-10 relative).
PASS = identical. Any difference is a bug in the ownership/regrid plumbing, not a physics
result.

**G2 band (B1_250 vs U1_250).** 2 mm base + 1 mm band (20 mm below / 4 mm above the
surface, regrid every 2 s) against uniform 1 mm, same binary, D2j-0b metrics reused
unchanged (`d2j0b_plane/analyze.py`): disk-mean rate per 50 s block (r < 30 mm), band
t_stop / t_cross table, r_front. PASS requires the disk-mean rate within **5 %** of U1_250
in every block 50–100 / 100–150 / 150–200 s and no band stalling in one run but not the
other. The band edge is the only approximation (coarse conduction below 20 mm), so the
expected residual is set by the thermal penetration depth sqrt(kappa t) ~ 13 mm at 250 s.

**Reported, not deciding.** Wall time ratio B1_250 / U1_250; cells per level from the
plotfiles; regrid repair messages (cells rebuilt) in the log; any "surface left the refined
band" abort (which is a FAIL of the band sizing, reported as such).
