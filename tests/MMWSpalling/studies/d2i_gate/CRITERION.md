# D2i mesh gate for the support rule with the idle clock off (fixed before launch; never edit)

**Pair.** `LI0_2mm` and `LI0_1mm`: D2f Q09 keys (feet descent, `foot_rule =
pads`, `foot_pad_quantile = 0.9`, `foot_body_clearance = 0`) with
`surface_patch.pinned_idle_cycles = 1.0e9` (the idle clock off), i.e. D2h's
I0 keys. 0.40 m full 0.12 × 0.12 quarter domain, stop cap 700 s, ended by the
bottom watchdog. Binary `6d047b50…`, no build.

**Burner ROP.** The least-squares slope of `nozzle_z` against time over the
thermo rows with a ≤ t ≤ b (and t > 0), in m/h. This is `d2e run.blocks()` /
`analyze_pairs.window_rop`, the same estimator as D2f, D2g-1 and D2h.

**t_m** = min(t_end of the two runs), where t_end is each run's last thermo
row.

**Gap** = ROP(1 mm) / ROP(2 mm) − 1.

**PASS requires both:**
- **(a)** the matched-window gap over [150 s, t_m] is within **5 %**;
- **(b)** every **100 s** block is within **5 %**. The blocks start at 150 s:
  [150, 250], [250, 350], …, and the last one is [150 + 100k, t_m]. That last
  block is dropped if it is shorter than 50 s.

**Otherwise FAIL.** No other window, block length or offset is evaluated for
the verdict.

**Reported, not deciding:**
- the 50 s blocks from 0 s ([0, 50], [50, 100], …, the last ending at t_m if
  it is ≥ 25 s), for continuity with D2g-1 and D2h;
- the growth diagnostic: the least-squares slope of the 100 s block gaps
  against their block mid-times, in percentage points per 100 s;
- the silent-fraction table by 4 mm radial band at 250 / 350 / 450 / 550 s
  with pad firings per column per 100 s, in the D2h §2 format. A column is
  silent at time t if the time since its last removal event exceeds
  2·t_cell, t_cell = ρc_p·dz·528 / 3e5;
- per 100 s block: rim spall by band, pinned share, idle columns, `T_rec`,
  `jet_s_c`, `jet_r_reach`, `jet_P_face`, flank (40–60 mm) removal, and the
  worst `ledger_err`.

**Silence above the nozzle plane does not count against the gate.** The
44–52 mm wall bands go silent under the `s ≤ 0` rule at both meshes; that is
the next physics packet, not this gate.

**Identity check (not part of the verdict, but reported):** each run's
`thermo.dat` rows with t ≤ 349.9 s must be byte-identical to the matching
D2h `I0` run — same binary, same keys, same plotfile spacing; only
`stop_time` differs. `I20_1mm` sets a different key, so it has no identity
partner.

**Outcomes:**
- **PASS:** the support rule is mesh-converged to t_m at 0.40 m. D2g Stage 2
  (the scored Meier runs) is planned, and must quote "mesh-checked to t_m at
  0.40 m; height alone moved ROP +4 % in D2d".
- **FAIL:** stop. Report the new front's position and timing, the pinned
  share, idle columns, `jet_s_c`, `jet_P_face` and flank removal per block,
  and name what appears next. Do not test it here, and do not draw a
  mechanism from the final block alone — it sits against the bottom stop.

**Idle sensitivity (reported, not part of the verdict).** Burner ROP over
250–350 s **at 1 mm** for `pinned_idle_cycles` = 2 (`d2g_mesh/output/L09_1mm`),
20 (`I20_1mm`, new, to 350 s) and 1e9 (`d2h_openloop/output/I0_1mm`), with the
idle-column count and the silent-fraction row by band at 350 s for each. 1 mm
is where the front appeared; at 2 mm the clock is already known to be
immaterial. The question is **whether any finite clock eventually cascades or
whether a long one is safe**, which is what must be known if a physical idle
timer is ever reintroduced. No value is chosen to change the rate.

---

**Amendment record.** A first version of this file was written and hashed
(sha256 `051e93a50497d0e1089e98e56f77d3d64db6d8b7750e5f7c10d32061803939e3`,
2026-09-20 19:01:50) against the packet as it then read. The user revised
`ACTIVE_STEP.md` afterwards, and this file was amended to match before any
run that this gate scores: the sensitivity run moved from 2 mm to 1 mm, and
the plotfile spacing from 100 s to 50 s at both meshes with the full D2h
field list. **The deciding text — the pair, the ROP estimator, t_m, the gap,
and the PASS conditions (a) and (b) — is byte-unchanged from that first
version.** The two runs started under the first version were killed at
164.8 s (2 mm) and 11.7 s (1 mm) and their output deleted; every run scored
here starts from t = 0 after this amendment.
