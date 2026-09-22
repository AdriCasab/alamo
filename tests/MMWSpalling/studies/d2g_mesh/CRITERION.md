# D2g Stage 1 mesh criterion (fixed before launch; never edit)

**Pair.** `L09_2mm` and `L09_1mm`: D2f Q09 keys (pads, q 0.9, no clip),
0.40 m full domain, stop cap 700 s, ended by the bottom watchdog.

**Burner ROP.** The least-squares slope of `nozzle_z` against time over the
thermo rows with a ≤ t ≤ b (and t > 0), in m/h. This is `d2e run.blocks()`
and `analyze_pairs.window_rop`.

**t_m** = min(t_end of the two runs), where t_end is each run's last thermo
row.

**Gap** = ROP(1 mm) / ROP(2 mm) − 1.

**PASS requires both:**
- **(a)** the matched-window gap over [150 s, t_m] is within **±5 %**;
- **(b)** every **100 s** block is within **±5 %**. The blocks start at
  150 s: [150, 250], [250, 350], …, and the last one is [150 + 100k, t_m].
  That last block is dropped if it is shorter than 50 s.

**Otherwise FAIL.** No other window, block length or offset is evaluated for
the verdict.

**Reported, not deciding:**
- the 50 s blocks from 0 s ([0, 50], [50, 100], …, the last ending at t_m if
  ≥ 25 s), for continuity with D2f;
- the growth diagnostic: the least-squares slope of the 100 s block gaps
  against their block mid-times, in percentage points per 100 s.

**Outcomes:**
- **PASS:** D2g Stage 2 (scored runs) is planned, quoting "mesh-checked to
  t_m at 0.40 m".
- **FAIL:** stop. The next packet diagnoses the late drift. The leads are
  named from the logs (rim spall rate, flank removal, T_rec, s_c,
  `jet_r_reach`, pinned share) and not tested here.
