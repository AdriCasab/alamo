# Notes archive

Caveats that no longer share an identifier with the active packet or the standing rules.
**Nothing here is known to be wrong** — the liveness check found only 7 of 339 notes whose
named references are actually gone. These are parked for size, not retired for staleness, so
**grep here before concluding something was never known.**

Split from `ROADMAP.md` on 2026-10-02: 274 of 339 bullets, 2,373 lines.
Regenerate the index with `python3 docs/project/tools/notes_index.py`.

<!-- INDEX:BEGIN -->
- L288: THE EXPOSURE DEFINITION IS NOT GRID-ISOTROPIC
- L295: A GATE SCOPED TO A SUB-DOMAIN MUST STATE WHAT IT EXCLUDES AND WHY THAT IS SAFE
- L301: A MEDIAN IS THE WRONG STATISTIC FOR A DEFECT WITH ANGULAR STRUCTURE
- L305: RETRACTED, both mine
- L311: THE BORE MUST BE REPORTED PER AZIMUTH SECTOR, NEVER AS A SINGLE DIAMETER
- L313: The rising support share was a DENOMINATOR EFFECT
- L317: The rule LEVELS; the criterion REMOVES
- L321: spall.support_step_min = 6 mm is a population separator, not a dial
- L325: DOCUMENT HYGIENE IS OVERDUE AND IS NOW A KNOWN DEBT
- L337: THE FIRING CRITERION IS GRADIENT-BLIND BY CONSTRUCTION
- L365: ON A FLAT PROFILE THE RAW RECONSTRUCTION IS CORRECT — DO NOT HOLD, SMOOTH OR CLAMP
- L375: A FLAT TWO-CELL RATIO MEANS THE GRADIENT IS UNRESOLVED, NOT THAT THE FLAKE IS
- L396: A QUASI-STEADY, LATERALLY UNIFORM TEST COLUMN CANNOT ATTRIBUTE A CHANGE THAT ACTS
- L407: REMOVED VOLUME IS ENERGY-LIMITED; RATE IS NOT
- L412: A QUANTITY THAT GROWS WITH HOLE DEPTH MUST BE QUOTED OVER THE WINDOW IT IS
- L421: RUNS IN /tmp ARE NOT RESULTS — AND ABSENCE OF EVIDENCE IN THE STUDY DIR IS
- L431: A HARNESS BAR THAT WAS NEVER TIGHTENED WILL PRINT PASS THROUGH AN 18×
- L436: CHECK WHICH BINARY PRODUCED THE SCORED RUNS
- L442: SMOOTHING AN ESTIMATOR THAT FEEDS A THRESHOLD CAUSES PREMATURE FIRING
- L454: A THRESHOLD CRITERION'S OWN VARIABLE IS A VACUOUS WITNESS — AND IT HAS NOW
- L465: THE TWO-CELL ESTIMATOR NEEDS THE HEATED LAYER TO SPAN MORE THAN ONE CELL
- L472: SMOOTHING A MESH-SENSITIVE ESTIMATOR CAN WIDEN THE SPREAD, NOT NARROW IT
- L479: ONE RULE FOR THE SURFACE TEMPERATURE, TWO REGIMES: EXPONENTIAL WHERE THE FRONT
- L490: A DIRECTLY-REMAPPED TEMPERATURE AND A BISECTED ONE DIFFER AT THE INVERSION
- L500: SHARE OF POWER IS NOT SHARE OF ERROR
- L507: A NEW BOUNDARY CLOSURE NEEDS ITS OWN robin_form, NOT A BRANCH ON ANOTHER
- L513: MESH-INDEPENDENCE AND ENERGY CONSISTENCY ARE TWO DIFFERENT GATES
- L522: surface.follow_mask = 1 IS MANDATORY on any case with active removal and a
- L530: A THRESHOLD CRITERION'S OWN VARIABLE IS A VACUOUS MESH WITNESS
- L539: THE DEPTH SCAN RETURNS h_scan ≈ dz/2 − a_f REGARDLESS OF THE GRADIENT
- L548: WHY D2z-0/-1 RETIRED THE SUB-CELL INCREMENT — the record, not the folklore
- L560: A CORRECT SURFACE TEMPERATURE ALONE MAKES THE MODEL 80–130 % FASTER, AND THAT
- L571: A RECONSTRUCTED SURFACE TEMPERATURE CLOSES A POSITIVE FEEDBACK LOOP
- L584: weibull.V0 = 8e-09 IS V_cell AT 2 mm
- L589: THE FLOOR GAS TEMPERATURE IS MEASURED, NOT ASSUMED: 1664.7 K
- L601: THE ROCK-SIDE STORY IS QUANTITATIVELY REFUTED AS A ROUTE TO THE FACTOR OF
- L619: KANT'S MEASURED h_fl,min = 500 ± 100 W/m²K IS A HARD FLOOR ON ANY FLUX
- L626: THE EQUILIBRIUM-STANDOFF PREDICTOR IS CALIBRATED — BUT IT ONLY EXISTS ABOVE
- L644: THE STAIRCASE "WETTED" AREA IS THE EXACT RAMP sqrt(1+δ²) − 1, NOT δ²/2
- L653: DO NOT BAR PHYSICS WITH AN ESTIMATOR CHECK
- L660: A DEFAULT-OFF KEY HIDES THE BLAST RADIUS OF ITS OWN CONDITION
- L671: THE HALF-CELL IS GONE; beam_void_closure T3 IS THE DEBT IT LEFT
- L681: THE TREE IDENTITY REFERENCE IS NOW witness/d2z1.md5
- L687: A SENTINEL THAT SHARES A TYPE WITH A MEASUREMENT WILL BE READ AS ONE
- L694: A RATE CRITERION HAS A TIME RESOLUTION. PRE-REGISTER THE WINDOW LENGTH,
- L706: A COLUMN FIRES TWICE PER CELL — event counts are NOT removal counts
- L713: A ONE-CELL STEP IN THE MODEL IS A SUB-CELL STEP IN THE ROCK
- L723: A IS ALREADY CASCADE-THROTTLED AT 2 mm — "A has no f, therefore no
- L731: HEATING AREA AND CONDUCTION AREA ARE DIFFERENT QUANTITIES AT A SUB-CELL
- L740: A FOURTH SIDE-PATH SITE, :2354-2355, AND ITS DEFAULT HIDES IT
- L750: THE POWER IDENTITY IS THE ONLY CHECK THAT CATCHES A KERNEL↔MARCH SPLIT
- L756: THE KERNEL HAS NO EXPLICIT TOP-FACE AREA
- L763: THE MEIER RUN KEYS ARE ON THE ARGV, NOT IN input_drilling
- L772: THE ARBITER MESH GAP HAS TWO MASKS AND THEY DISAGREE BY AN ORDER OF
- L780: D2k HEALS STEPS BY ACCELERATING THE TALLER COLUMN, NOT BY FEEDING THE
- L790: SLOPE STATISTICS MUST BE QUOTED ON THE BELOW-PLANE FLOOR
- L796: A CLIFF TEST MUST BE SCALE-FREE AND n_min A FRACTION OF THE DISC
- L805: THE RADIUS GATE SELECTS COLUMNS, NOT FACES
- L811: DO NOT ZERO JITTER RISERS WHEN FIXING THE AREA OVER-COUNT
- L820: TWO AREA ERRORS THAT A PER-FACE WEIGHT CANNOT FIX
- L828: SMOOTHING STENCILS MUST BE A FIXED PHYSICAL LENGTH
- L836: THE SURFACE IS BIMODAL, SO ONE STENCIL CANNOT SERVE IT
- L843: A WITNESS OVER A POPULATION MUST MASK ON THAT POPULATION
- L851: LEAVE NO GAP BETWEEN A PREDICTION AND ITS REFUTATION THRESHOLD
- L857: SCORE Ø(z) AT EQUAL BURNER DEPTH, NEVER AT EQUAL TIME — and only on rock
- L870: A (the confined map alone) ALREADY MEETS MEIER ON SHAPE AND RATE
- L879: f IS ALGEBRAICALLY UNFALSIFIABLE AS WRITTEN
- L886: THE DERIVED SIDE-FACE WEIGHT HAS A HARD FLOOR OF 0.4142, SO f = 0.2 IS
- L897: MEIER Ch. 3's 0.4–1.6 kW/m²K ARE IMPINGING-RIG VALUES — check before
- L903: THE SKIRT IS NOT IN THE MODEL
- L912: jet_P_face IS THE TOTAL STREAM DRAW, NOT THE TOP-FACE POWER
- L920: A CONFINED h MAP WITHOUT SIDE HEATING MAY COST RATE, NOT JUST SHAPE
- L929: ROP ON THE MEIER GEOMETRY IS THE SLOPE OF nozzle_z, NEVER PIT DEPTH
- L933: THE BELOW-PLANE SIDE-FACE PATH IS NOT A PASSENGER ON MEIER
- L939: side_face_factor REQUIRES side_face_flux — setting it to 0 ABORTS
- L944: THE ABOVE-PLANE ENTHALPY DEFECT, QUANTIFIED
- L953: REPORT EXCAVATED VOLUME ON BOTH PRO-RATING BASES
- L956: REFERENCE ERRATUM
- L960: jet_free_surface = 1 DOES NOT ABORT UNDER PRESCRIBED DESCENT — it silently
- L968: THE TAPER'S DIAMETER FLOOR APPLIES AT SHALLOW DEPTH ONLY
- L976: THE FLOOR/WALL h MAP MUST BRANCH ON s, NOT BLEND ACROSS IT
- L983: THE INNER BLEND MUST START AT THE SKIRT RADIUS foot_r_outer = 0.040
- L987: AN 8-COARSE-CELL TAPER IMPOSES A DIAMETER FLOOR
- L992: input_feet_d2c IS NOT THE FROZEN D2i CONFIGURATION
- L997: A FEET-DESCENT RUN NEVER EXECUTES THE ABOVE-PLANE BRANCH
- L1002: Recession beyond r = 76 mm on the frozen run is mean 0.08 mm but MAX 4.00 mm
- L1006: THE DOMINANT MEIER SHAPE ERROR IS A FUNNEL — THE MODEL IS ~45 mm TOO WIDE
- L1020: A SCALAR ABOVE-PLANE h IS A BLANKET — 89 % OF COLUMNS ARE ABOVE THE
- L1027: side_face_factor_above MUST BE SET TO 0.0 FOR A CLEAN CONTROL
- L1031: input_feet IS STALE: 0.3 m deep at 4 ms
- L1035: THE MEIER MAGNITUDE BAND IS 85–95 ± 5 mm, and the csv is SHAPE ONLY
- L1038: side_face APPEARS ZERO TIMES IN input_feet
- L1040: INTERPOLATE THE CROSSING RADIUS BEFORE ANY Ø CLAIM
- L1044: THE MEIER SHAPE REFERENCE
- L1050: THE s-SWITCH IS NOW THE SHARPEST CORNER IN THE MODEL
- L1056: A SINGLE MESH CANNOT SEE THE LEVER'S EFFECT
- L1060: THE CORE'S MARGIN OVER THE NOZZLE FEED IS ONLY ~5 %, AND h STEPS 4.9x AT
- L1070: THE RESIDUAL MESH GAP IS THE MASK EDGE, NOT THE INTERIOR (2026-09-23,
- L1087: "ABOVE-PLANE" IS THE MASK RIM, NOT THE DRILLING FRONT
- L1096: D2r-1's CRITERION (e) WAS MISREAD — THE NOTE THAT STOOD HERE IS WITHDRAWN
- L1102: A STAIRCASE FACE NEEDS A >= 1-CELL DROP, so any column owning one has local
- L1109: f = 0.2 IS RETIRED AS A MESH-GAP QUESTION BUT NOT AS A PHYSICS ONE
- L1114: EVERY METRIC IN THIS PROJECT IS A DISK-MEAN RATE, AND THAT IS HOW THE ARTEFACT
- L1120: The heated side faces are the STEP FACES OF THE DRILLING FRONT, not a pit
- L1129: D2r-1's (e) ratio does reach +19.06 % in the unscored 200–250 s block
- L1134: D2r-1 REDUCED THE DRAIN BY 6.8 % IN WATTS, not by eliminating it
- L1139: Hash CRITERION.md BEFORE THE FIRST LEG, not before the scored ones
- L1144: The AZIMUTH is unused and deliberately deferred
- L1150: The 0.5 mm third mesh is PARKED by the user (2026-09-23)
- L1154: STANDING PRE-FLIGHT ITEM 0, bought by THREE consecutive dead packets
- L1163: ABOVE-PLANE TOP HEATING IS ALREADY KNOWN INERT
- L1170: THE IDENTITY SET DOES NOT EXERCISE side_face_flux
- L1176: f AND THE CAP ARE TWO FUDGES ON A RAW FACE COUNT
- L1186: CHECK A PREDICTION AGAINST MEASUREMENTS ALREADY ON DISK BEFORE HASHING IT
- L1193: THE WATER-COOLED SHIELD BIASES BOTH SCORED QUANTITIES, SO IT GOES IN
- L1199: DERIVE A BAND'S ENDPOINTS FROM THE RUN'S OWN STATE, NOT FROM AN INPUT KEY
- L1206: AN ANCHOR MUST NAME EVERY KEY THAT DIFFERS (D2r-0)
- L1212: THE ENTHALPY BOUND DOES NOT BIND INSIDE THE HOLE (D2r-0)
- L1219: CALIBRATE A WALL MODEL ON A MEASUREMENT THE PROJECT ALREADY HAS (D2r-0)
- L1226: THE REFRAMING (2026-09-23): warm the wall, do not remove it — and widening
- L1234: jet_P_exhaust IS NOT LOST, IT IS UNCONTACTED
- L1239: (d)'s k·ΔT/dx diagnostic is STILL OWED
- L1243: THE WATER-COOLED SHIELD IS A MISSING COMPETING SINK
- L1249: test_feet, test_feet_d2c and unit/amr_microstructure_regrid are
- L1253: CONTIGUITY TESTS ARE THIS LINE'S RECURRING BUG CLASS (D2q-2b, 2026-09-22)
- L1263: A probe needs a guard that the state under test actually exists (D2q-2b)
- L1269: The "don't rebuild during a sweep" rule needs harness enforcement, not
- L1274: surface.follow_mask = 1 is ALREADY SET everywhere this line scores
- L1282: The 2639-file identity set is SINGLE-LEVEL, so the regrid-repair site is
- L1287: A pointwise corner flux must not be a convergence gate once removal is
- L1293: A default-off key keeps byte-identity available as the plumbing oracle
- L1299: THE BINARY sha256 IS NOT AN IDENTITY ORACLE ON THIS MACHINE (D2q-1,
- L1309: A pre-check that detects one direction of a bug does not clear the other
- L1315: THERMAL LENGTH CORRECTED: κ/v ≈ 1.6 mm, not 3.2 (2026-09-22)
- L1324: Anchor a ceiling to an artefact-free reference (D2p-1, planner error)
- L1332: A criterion a fix satisfies by construction is not a criterion
- L1338: Corner regularisation points the WRONG WAY and must not be revisited
- L1345: The height function cannot express lateral recession — and the continuous
- L1360: Before a packet claims a mechanism "fits the existing representation", name
- L1365: G1-style field gates must be ONE-SIDED (D2o-0b review, 2026-09-22)
- L1374: The mesh-divergence trigger must watch k·ΔT/dx at the transition against
- L1379: Depth-mean Ø is not a single number for a cone (D2o-0b)
- L1385: The 85–95 mm width carries ~a quarter of the volume and is not clearance
- L1390: AMR is not a cheap parallel win, and will not buy back the 6 h fine legs
- L1395: Fig. 8.8 colour original: request it as a CHECK, not an anchor
- L1401: Ø(z) may only be scored on rock the nozzle plane has already passed (D2o-0,
- L1415: A prescribed field is not the control's field unless it is bounded in
- L1431: The floor/wall corner must be radial, not on stand-off (D2o-0, 2026-09-22)
- L1442: A shape criterion must bind on the whole profile, not its mean (D2o-0)
- L1450: Wall metrics must be defined on geometry, not on stand-off (D2o-0)
- L1459: The "do not rebuild" block from D2k/D2l is LIFTED (2026-09-22)
- L1467: Idle clock (pinned_idle_cycles) is mesh-dependent (D2h, 2026-09-20)
- L1474: D2e retroactive review (2026-09-19, accepted) — carry to the commit
- L1483: How Meier ran the pilot test vs the model (2026-09-18, full Ch. 8
- L1518: D2f planning decisions (2026-09-18, revised same day)
- L1526: D2e planning decisions (2026-09-18)
- L1538: D2d planning decisions (2026-09-17, planner, before any run)
- L1546: surface_mf is φ-windowed
- L1552: The S1 criterion finding (2026-09-15)
- L1567: Overshoot time-step rule for removal runs
- L1570: spall.h_col_events_csv is edge-triggered
- L1575: weibull.enabled = 0 + spall.enabled segfault — FIXED by A1
- L1578: Weibull V0 (user decision 2026-09-15): use weibull.V0 = V_cell in all
- L1587: Pinned Robin closure (C1) caveats
- L1598: Burner drawings (planner, 2026-09-16) — the D2b blocker is RESOLVED
- L1631: D2a2 lessons that constrain D2b and later jet work (2026-09-16)
- L1662: Process (planner decision 2026-09-16, from D2a2 takeaway 12)
- L1668: weibull.V0 = V_cell switches the size effect off, not just its mesh
- L1672: Jet h magnitude conflict (planner, 2026-09-15; D2a2 narrowed it, D2b settles)
- L1677: AMReX rank-count identity
- L1679: Energy ledger (energy_ledger.enabled)
- L1681: Kant p = 0 onset now prints 446.1 °C
- L1686: P4 cuttings/debris shielding — WITHDRAWN 2026-09-15 (design kept for
- L1698: Refactor track (R-series) is decoupling src/Integrator/MMWSpalling.H
- L1714: The long plan still contains an old Zhang/Oglesby note with omega_0 = 0.20 m
- L1716: Step numbering contains inserted/moved validation steps such as 4b, 6d, 7b,
- L1718: Homogeneous Hu reproductions should stay homogeneous. Voronoi heterogeneity is
- L1720: Step 8 changed AMR microstructure repair: passive D is still re-evaluated
- L1724: Step 8b is complete and must remain an indicator validation with
- L1727: Step 8b fixed a mechanics coupling issue: surface_patch prescribed
- L1730: Step 8b passes Hu surface f_b targets but relaxes L6 deep f_b to <0.25
- L1733: Quartz alpha-beta transformation strain is implemented only on the
- L1737: Step 9 CZM is default-off under cohesive.enabled. Its fields
- L1741: Step 11 removal is default-off under spall.enabled || vapor.enabled
- L1746: Step 11 void handling exists only on the microstructure heat path; the
- L1749: For Hu onset/removal work that needs damage or spall fields, use homogeneous
- L1752: Step 11b is a calibrated v1 onset validation, not a fully predictive damage
- L1766: AMR multi-level mechanics for hu_spall_onset is a known limitation, not a
- L1779: revised-model-approach.md (Stages 1–5 of the LEFM Sp + Weibull
- L1783: Step 12 ships the LEFM Sp + K_Ic infrastructure (toggle, Tada-weight
- L1792: Kant onset closed form is Kant Eq. 15, NOT Eq. 14 Sp_red
- L1802: Kant's Table-2 390–560 °C band was computed with frozen room-T
- L1812: The lateral roller box (el.bc.type = constant, all 26 regions
- L1823: The elastic MLMG solver fails on anisotropic cells
- L1827: Weibull first-fire is governed by order statistics, not the scale
- L1844: The Sp criterion stays sign-agnostic in src/Numeric/SpCriterion.H
- L1850: h_spall_field_mf and Rossi's measured observable are different
- L1862: Step 16 geometry was forced off-packet by MLMG conditioning + free-
- L1874: Sub-cell h_col events are scan-quantization-suppressed
- L1889: σ_xx(z) is local-per-z under 1-D confinement
- L1898: The crack_event_count_mf field is cell-centred and is
- L1911: Step 16c's option (L) σ_xx reconstruction has a sampling
- L1929: rossi-validation-diagnostic-design.md is the canonical
- L1951: Rossi R1 deferral checklist (do not re-litigate without
- L1961: Step 16d's σ_xx reconstruction is gated to surface_patch.mode =
- L1974: Step 17 per-cell v_n + regime_field + local-Q infrastructure
- L1996: §15c shipped DetachMode { PerFace, ConnectedCluster } +
- L2006: Step 14b patched PrincipalStressRatio (Step 10 production spall
- L2014: The Hu validation chain (6d → 7b/7d → 8b → 11b → 14) is closed on
- L2018: ARCHIVE_DONE.md and the Step 11b prose still mention
- L2023: src/BC/Operator/Elastic/ZloRoller321.H is untracked in git but is
- L2027: The kept working-tree changes to src/Solver/Nonlocal/Newton.H and
- L2036: Run test and validation simulations with 4 MPI ranks by default on this
- L2039: Step 6b fixed the old microstructure path that evolved Temp directly and
- L2041: Step 6c split mineral topology from grain topology: phase is a mineral ID,
- L2044: microstructure.grain_expr is optional in expression mode; avoid
- L2047: Old Step 4c is now Step 6d so Hu conduction validates the thermal driver
- L2049: The current Hu validation text has constants, BC trendlines, and qualitative
- L2053: Step 6d added surface_patch. on the explicit global-material path;
- L2057: surface_patch.T_expr must be one ParmParse token with no whitespace, e.g
- L2059: Step 15a wired surface_patch.mode = prescribed_T | convective_flame
- L2067: Step 7b is a homogeneous Hu reproduction and now passes single-level. It ships
- L2070: Do not tighten the Step 7b L3 stress ordering without revisiting the
- L2073: Step 7d is complete with shipped 20^3 + max_level=1 AMR Hu reruns. Keep the
- L2076: Exploratory 16^3 + max_level=2 Hu AMR runs showed static mechanics
- L2082: For yt/numpy regressions on this machine, use
- L2086: THE STEP-FACE MECHANISM — the cause of the D2d–D2i mesh failures
- L2130: D2k's outcome and the two live explanations (2026-09-22)
- L2162: Shape, not rate, is the honest score (D2l, 2026-09-22)
- L2189: How the feet reach the heat distribution (code path, for D2m)
- L2202: The burner is a body in the hole, not a point source (second reading of
- L2225: Nozzle temperature — what is and is not sourced (2026-09-22)
- L2254: A WITNESS MUST BE SCORED ON THE SAME POPULATION AS THE THING IT MEASURES
- L2270: A CRITERION BAR MUST NOT BE TIGHTER THAN THE INSTRUMENT (2026-09-28, D2z)
- L2278: side_face_area_law DEFAULTS TO step (:6819), and run.py's SHARED
- L2285: f IS GONE, AND THE RIM ARTEFACT WENT WITH IT (2026-09-28, D2z)
- L2296: delta_fallback_cols — BOTH READINGS ARE RIGHT, ON DIFFERENT WINDOWS
- L2303: EVERY MEIER THERMO POWER IS A QUARTER-HOLE SUM (2026-09-28)
- L2311: A's IN-BAND RATE WAS THE STEP-DRAIN THROTTLE — do not build on it
- L2317: side_faces COUNTS FACES THAT CARRY NO AREA — and it cost the planner a
- L2338: The model's specific energy is 1.35 J/mm³, and that is a SEPARATE finding if the
- L2347: A FALSE PASS IS AVAILABLE ON ANY SHIELDED-PAD LEG (2026-09-28)
- L2355: SPECIFIC ENERGY IS CONSTANT AT ~1.32 J/mm³, SO HEAT CANNOT BE MOVED TO BUY
- L2364: THE SIDE TERM DEFENDS ITSELF AGAINST A LOWER h (2026-09-28, D2aa)
- L2371: A FROZEN GAS TEMPERATURE IS SELF-INCONSISTENT — THE MARCH CONSERVES ENTHALPY
- L2380: THE RATE IS NOT STEADY — report it as a profile against depth (2026-09-28)
- L2386: PD JAMMED: shielding the pad ring is BLOCKED BY LATERAL RECESSION
- L2396: side_face_h REACHES EVERY LATERAL FACE, AND TWO-THIRDS OF THEM ARE NOT ON
- L2408: THE 4–6 J/mm³ THREAD IS CLOSED — it was a units mismatch (2026-09-28)
- L2422: MEIER'S RATE IS ONE-SIDED AND THE MODEL COMPUTES A DIFFERENT QUANTITY
- L2436: THE DECAY-CLIP WARNING WAS BACKWARDS AND IS RETIRED (2026-09-29)
- L2444: side_faces vs the WALL BAND: removing wall heat LATCHES it off (2026-09-29,
- L2451: A RATE RELATION THAT WORKS, AND ITS FAILURE MODE (2026-09-29, D2ab)
- L2458: FIG. 8.8 BOUNDS THE PIT LEAD, AND IT KILLS THE "CATEGORY ERROR" EXCUSE
- L2478: MEIER'S 3.42 L IS NOT A DIAMETER (2026-09-29)
- L2487: P0b (h) FAILS AT THE PERIPHERY ONLY, AND NO BOUND CAN FIX IT (2026-09-30)
- L2501: THE AXIS RATE IS NOT AN INDEPENDENT TEST UNDER PRESCRIBED FEED (2026-09-30)
- L2510: REMOVAL IS ENERGY-LIMITED — VOLUME IS THE WRONG INVARIANT (2026-09-30)
- L2517: QUARTER DOMAIN: THE ×4 IS REQUIRED FOR ANY MEIER VOLUME COMPARISON (2026-09-30)
- L2524: A DIAGNOSTIC THAT HAS NEVER BEEN NON-ZERO HAS NOT BEEN TESTED (2026-09-30)
- L2531: SIDE-POWER BINS DECOMPOSE GROSS, NOT NET (2026-09-30)
- L2539: DEEPEN A BLOCK DOWNWARD (2026-09-30)
- L2543: THE KANT ANCHOR SCORES A VACUOUS VARIABLE (measured 2026-10-01)
- L2561: ONLY 29 % OF REMOVED ROCK LIES BEYOND MEIER'S WALL (2026-10-01)
- L2569: THE 40–50 mm BAND IS FLOOR, NOT WALL — AND P1's GAIN THERE IS ARTEFACT-ASSISTED
- L2584: THE ESTIMATOR'S ROOT CAUSE: BOTH FITS ASSUME DECAY TO AMBIENT (2026-10-01)
- L2592: NO MEIER RATE, VOLUME, WIDTH OR STANDOFF IS QUOTABLE FROM P1 (2026-10-01)
- L2599: MEIER'S VOLUME PER mm OF AXIS IS 5.93 mL/mm, NOT 7.1 (2026-10-01)
- L2605: A CRITERION WHOSE DENOMINATOR THE TREATMENT REMOVES CANNOT BE SCORED (2026-10-01)
- L2612: (c)'s SHIPPED-ONSET CLAIM HOLDS ON ONE LEG OF THE PAIR (2026-10-01)
- L2619: A DIFFERENCE IS ONLY WELL DEFINED ON ONE PROFILE (2026-10-01)
- L2625: THE REFIT AND THE CRITERION VERY NEARLY CANCEL — AND THAT WAS NEVER PRE-REGISTERED
- L2634: THE EXPOSED-FACE CENSUS REPRODUCES EXCEPT ONE ROW (2026-10-01)
- L2642: a0 AND THE DRIVE FRACTION ARE ONE CALIBRATION KNOB, NOT TWO (measured 2026-10-01)
<!-- INDEX:END -->

- ⚠⚠ **THE EXPOSURE DEFINITION IS NOT GRID-ISOTROPIC.** This replaces my earlier
  *"the successor is the D2ad 25-50 mm wall corner"*, which is **retracted**: on a Cartesian grid a
  round inward step presents **two** inward-lower face neighbours on its **diagonals** (Δr = 0.71
  cell each) and **one** on its **axes** (Δr = 1 cell), so `support_min_sides = 2` drills the 45 deg
  directions of a smooth round wall and each firing creates the next. rev 7's criterion firings in
  that band are flat to **1.16x** max/min, so the flux map is **not** the cause. Any rule keyed on
  counting face neighbours inherits this and must be tested per azimuth.
- ⚠⚠ **A GATE SCOPED TO A SUB-DOMAIN MUST STATE WHAT IT EXCLUDES AND WHY THAT IS SAFE.** P1b's
  criterion (k) — *"the rule must be inert on a healthy floor"* — was exactly the right gate and was
  written with an `r < 40 mm` cut (`exposure.py:R_BURNER`). The artefact lives on the **40-50 mm
  wall**, where the same rule fires **53.4 %** of rev 7's healthy columns (21 % axes, **89-91 %**
  diagonals) and **26.3 %** of the 50-60 mm mouth. **Predictable at zero compute from a file already
  on disk, before a 500 s leg launched.** The exclusion, not the bar, was the defect.
- ⚠ **A MEDIAN IS THE WRONG STATISTIC FOR A DEFECT WITH ANGULAR STRUCTURE.** P1b's criterion (f)
  passed on a median over r < 40 mm while an **18.21x** azimuthal spread sat just outside it. The
  isotropy witness is a **per-15 deg-sector median** with a max/min ratio: rev 7 reads **1.05x** on
  the 44-50 mm ring and **1.04x** on 40-50 mm; P1b reads **18.21x** and **5.19x**.
- ⚠⚠ **RETRACTED, both mine.** (1) *"The hole came out cylindrical."* Straight-sided is real — axis
  radius 44.4 mm and diagonal 48.6 mm each constant from 100 to 250 mm depth, so P1's funnel IS gone
  — but the cross-section is **9-12 % out of round** (≈8 mm on diameter) against rev 7's 1-3 %. The
  max-radius-over-all-azimuths measure **picks the diagonal**. (2) *"The bore equals the heating edge
  to a millimetre."* Cause and effect inverted: the diagonal reaches the heating edge at 48-50 mm, the
  axes stop 4-5 mm inside it.
- ⚠⚠ **THE BORE MUST BE REPORTED PER AZIMUTH SECTOR, NEVER AS A SINGLE DIAMETER.** This is the direct
  consequence of the retraction above, and it applies to every future shape claim.
- ⚠ **The rising support share was a DENOMINATOR EFFECT.** 4.8 -> 25.4 % per 100 s, but support
  firings are 3276 -> 6295 -> 4990 -> 3989 -> **4250** (flat to falling) while criterion firings
  collapse **64,922 -> 12,485** (5.2x) as the standoff opens. "Becoming a larger share of removal"
  was misleading.
- ⚠ **The rule LEVELS; the criterion REMOVES.** P1b predicted the rule would miss the 22 % of frozen
  columns with one exposed side and that 30-60 would stay stuck. **Zero did.** Removing a laggard's
  proud neighbours lowers the step it sits in, its profile re-steepens, and the **ordinary** criterion
  takes it. ⇒ a more restrictive `min_sides` can still cure the stall, which is the premise of P1c.
- ⚠ **`spall.support_step_min` = 6 mm is a population separator, not a dial.** On a healthy floor the
  qualifying fraction collapses 32.9 -> 12.0 -> **0.0 %** from 2 to 6 mm; on a stalled floor it barely
  moves, 47.5 -> 43.4 -> 40.8 %. 6 mm is the gap **between** a healthy floor's one-to-two-cell phase
  lag and a stalled floor's four-plus-cell chasm.
- ⚠ **DOCUMENT HYGIENE IS OVERDUE AND IS NOW A KNOWN DEBT.** `CLAUDE.md` says keep `ROADMAP.md` short
  enough to read every session; it is **~4,780 lines / ~52,000 words**, with ~1,560 lines of completed
  history in `## Current Status` and ~3,100 lines / 325 bullets here. The planner spec's
  *"do not delete caveats unless the user says they are resolved"* is a **ratchet with no release**.
  Agreed plan: (i) mechanical half now-ish — a heading index for the 14k-line `ARCHIVE_DONE.md`,
  collapse `## Current Status` to one line per step, put **numbers** in `CLAUDE.md`'s budget rules,
  and add a `## Criteria` table with **run-order** and **what-a-failure-means** columns to the planner
  spec; (ii) the 325-bullet triage **with the user**, splitting out a hard-capped **standing-rules**
  block (unit conventions, scoring traps, coupled calibrations such as `a0`+drive-fraction) that is
  **never** archived. ⚠ Not yet done — this entry is the reminder.


- **⚠⚠ THE FIRING CRITERION IS GRADIENT-BLIND BY CONSTRUCTION.** It is a Tada edge crack
  in uniform compression over a flaw of **4–20 µm** — far smaller than a cell — so it
  samples an effectively uniform stress. **A flat slab at the firing temperature fires
  exactly as a steep skin does**; "a preheated slab should spall less" is **not in the
  model anywhere** and must be *added*, not restored. The only place a gradient enters
  today is the thickness law `h_f = c_f·ell_eff`, which is why an unresolved gradient
  produces an absurd *thickness* rather than suppressed firing. ⇒ **any thickness cap is a
  STOPGAP**; the gradient-aware criterion (buckling / energy release with a flake
  thickness) is the real item, and a `K_I(z)` scan on the reconstructed profile does **not**
  substitute for it — on a flat profile it returns a thick slab too.
  ⚠⚠ **THE GRAIN-SCALE-FLAW ROUTE IS THE WRONG SIGN AND IS WITHDRAWN (2026-09-30).** It was
  recommended here and in P0b revision 7; it is refuted. A larger flaw samples stress
  **deeper**, and a flat profile has *more* hot rock down there. Planner-verified through the
  **shipped** Tada integral at fixed surface temperature, `K_I` relative to the flat case:

  | flaw | front 0.42 mm | ell 1 mm | ell 2 mm | ell 7 mm | flat |
  |---|---|---|---|---|---|
  | 0.02 mm (today) | 0.97 | 0.99 | 0.99 | 1.00 | 1.00 |
  | **1.0 mm** | **0.31** | 0.57 | 0.75 | 0.92 | **1.00** |

  ⇒ the physical sharp front would get **3.2× less** drive than the unresolved wall. It
  **rewards exactly the columns that are the error**.
  ⚠ **`K_I` also cannot select a thickness**: `K_I ∝ √a` while the differential stress
  saturates for `a ≫ ell`, so maximising over `a` runs to the scan limit — unbounded
  behaviour in a new costume. Making the thickness an **output** needs an energy-release
  criterion and is a separate step from the firing test.
  ✅ **The route that works is DIFFERENTIAL** — see the next note.

- **⚠⚠ ON A FLAT PROFILE THE RAW RECONSTRUCTION IS CORRECT — DO NOT HOLD, SMOOTH OR CLAMP
  `x`.** The reconstruction's factor is `x/(1−e^{−x})`, which tends to **exactly 1** as the
  profile flattens, i.e. `T_surf = T_top` — which is the right answer for uniformly
  preheated rock. **Substituting a "last resolved" steeper `x` fires the column hundreds of
  degrees early**: planner-verified on a 702 K Meier cell, a held `x` = 2 reads **1239 K**
  and fires at `T_top` ≈ **494 K**, against the raw 758 K firing at ≈702 K. That is the same
  over-read that produced the h500 clamp failure. ⚠ And on a side-heated floor the flat
  state is the **persistent** state, not a transient to bridge. ⇒ **only the thickness law
  is wrong on a flat profile; bound `h_f` alone.**

- **⚠⚠ A FLAT TWO-CELL RATIO MEANS THE GRADIENT IS UNRESOLVED, NOT THAT THE FLAKE IS
  THICK — AND THE MODEL'S RESPONSE TO IT HAS THE WRONG SIGN.** The surface
  reconstruction infers a decay length from the ratio of the top two cell temperatures.
  On a **side-heated, cratered floor** the rock is preheated well below the surface, the
  ratio goes to ≈1, and `ell_eff` blows up: measured on Meier at 5.911 m/h the physical
  `κ/v` is **0.42 mm** while the model read **7.54 mm — 18×** — from a ratio of **1.3
  where a receding front demands 117**. Because `h_f = c_f·ell_eff`, the model then
  sheds slabs of **3–15 mm** (p50 3.02, p90 15.00, against 0.60 / 2.60 one revision
  earlier), bounded only by the 15 mm **search depth**.
  ⚠⚠ **Spallation is driven by the GRADIENT** — differential expansion between a hot
  skin and cold rock beneath — so **a uniformly preheated slab should fire LESS, not
  remove 15 mm.** This is an inverted response, not a tuning error.
  ⚠ A time-average over the estimate **masks** it (by averaging steeper moments in) at
  the price of premature firing elsewhere — so smoothing is not the fix.
  ⚠⚠ **A crack scan on the reconstructed profile does NOT fix it either** — run on a
  flat profile it returns a thick slab too. A **resolved-gradient guard** (hold the last
  resolved decay length) plus a **physical thickness bound** (Wilkinson's measured flake
  band tops near 1 mm) are prerequisites for any thickness-dependent law.
  ⚠ **Do not let a SEARCH DEPTH act as a physical bound**: `spall.h_spall_max` (15 mm)
  was bounding `h_f` by accident on 13.5 % of Meier firings.

- **⚠⚠ A QUASI-STEADY, LATERALLY UNIFORM TEST COLUMN CANNOT ATTRIBUTE A CHANGE THAT ACTS
  ON A CRATERED, SIDE-HEATED FLOOR.** "The 1-D column moved < 2 %, so the Meier change
  must be X" is not an attribution — the column is blind to lateral preheating, to
  unresolved gradients and to standoff feedback, which is where the difference lives.
  Attribute from the **event logs of the two runs being compared**, not from a
  configuration that lacks the mechanism.
  ✅ Worked example: Meier's apparent doubling was **not** the wall closure. Total
  removed volume rose **22.4 % on 12.7 % more delivered power** — energy-limited, and
  that is the closure's real size. The **rate** doubled because flakes got 5–25×
  thicker, which opened the standoff **134 → 205 mm** and *dropped* floor power a third.

- **⚠ REMOVED VOLUME IS ENERGY-LIMITED; RATE IS NOT.** Total volume tracks delivered
  power (+22 % on +13 %), while the descent rate depends on how that removal is parcelled
  into flakes and on the standoff feedback the map creates. ⇒ **quote volume when asking
  how much heat arrived, and rate only with the flake statistics beside it.**

- **⚠ A QUANTITY THAT GROWS WITH HOLE DEPTH MUST BE QUOTED OVER THE WINDOW IT IS
  USED IN.** Meier's side-face share of `P_robin` reads **0.363 at t = 60 s, 0.530 at
  160, 0.636 at 500** — it tracks annulus **wall area**, which grows as the hole
  deepens. The often-quoted **62 %** is the end-of-run figure; over the **60–160 s
  window the rate is actually quoted on it is 47 %**. ⇒ a large part of it is wall
  band **above** the floor: physically real in Meier, and largely **rate-inert**.
  **Share of power is not share of error** — bin by radius at the 40 mm gate and by
  whether the receiving cell later fires, over the rate window.

- **⚠⚠ RUNS IN `/tmp` ARE NOT RESULTS — AND ABSENCE OF EVIDENCE IN THE STUDY DIR IS
  NOT EVIDENCE OF ABSENCE.** P0b wrote every leg to `/private/tmp` and nothing under
  its study directory. A review then searched, found nothing, and concluded the Meier
  numbers had no run behind them; the run was in `/private/tmp/p0b_fh/` (63 MB event
  log, 188,785 rows, thermo to t = 500.0 s) and every quoted number reproduced from
  it. **Both lessons stand:** a packet must keep its runs, or a `MANIFEST.txt` with
  per-leg argv, output path, **binary sha256** and event-log sha256 — *and* a reviewer
  must check the scratch paths named in the packet's own commands before calling a
  result unbacked.

- **⚠⚠ A HARNESS BAR THAT WAS NEVER TIGHTENED WILL PRINT `PASS` THROUGH AN 18×
  REGRESSION.** P0b's h500 anchor moved **0.116 % → 2.089 %** against a legacy 3 %
  bar inherited from P0, and printed PASS. **Score each packet against the previous
  packet's MEASURED value, not the bar it was allowed.**

- **⚠⚠ CHECK WHICH BINARY PRODUCED THE SCORED RUNS.** P0b's 18 column legs and 12
  anchor legs ran at 14:58–15:20 against a **final binary built at 15:41:51**, with
  the flake-thickness cap — on the scored path — landing in between; one leg of
  eighteen was re-run after. **Record the binary sha per leg and assert it matches
  the shipped build.**

- **⚠ SMOOTHING AN ESTIMATOR THAT FEEDS A THRESHOLD CAUSES PREMATURE FIRING.** P0b's
  1 s EMA on `ell_eff` started from the first raw value (huge on a cold start) and sat
  on the `ell ≥ dz/4` floor, firing the h500 anchor columns at **0.30–0.81 s** where
  P0's unsmoothed estimator fired them at **1.44–1.81 s** with no clamp. The clamped
  legs reported a **471–481 K** surface rise where the true value was **208–316 K**.
  ⚠ The counter-argument that "without smoothing the raw estimate has no content at
  those instants" is a **circular** one: without smoothing the column does not fire
  at those instants — it keeps heating until the estimate has content.
  ⚠ The real residual: **`ell ≥ dz/4` is a MESH quantity setting a PHYSICAL one**, and
  when the heated layer is thinner than a cell the two-cell fit is mesh-dependent.
  That is a P1 item.

- **⚠⚠ A THRESHOLD CRITERION'S OWN VARIABLE IS A VACUOUS WITNESS — AND IT HAS NOW
  BITTEN TWICE.** P0 built `spall.surface_T = diagnostic` precisely so the *physical
  state at firing* could be scored instead of the criterion's own threshold
  variable; P0b then ran that control **only at `h_conv` = 150**, and the
  `h_conv` = 500 legs silently fired on the `ell ≥ dz/4` floor clamp — **two of three
  at exactly `x` = 4.00, the maximum amplification the clamp permits** — reporting a
  471–481 K surface rise where the true value was **208–316 K**. Three clamped legs
  agreed to 2.089 % and that agreement meant nothing.
  ⇒ **EVERY scored anchor needs its diagnostic control AT ITS OWN OPERATING POINT**,
  and a scored leg whose `rc_flag` shows the clamp is **invalid, not merely noisy**.

- **⚠ THE TWO-CELL ESTIMATOR NEEDS THE HEATED LAYER TO SPAN MORE THAN ONE CELL.**
  At `h_conv` = 500 the onset comes at t = 0.3–0.8 s, where thermal penetration
  `sqrt(4·alpha·t)` is **0.91–1.50 mm against cells of 2.0 / 1.0 / 0.5 mm**; at
  `h_conv` = 150 onset is at t ≈ 18–23 s and penetration is **7.15 mm, 3.6 cells**.
  ⇒ an anchor taken where the layer is of order one cell is **outside the
  reconstruction's domain of validity**, whatever the smoothing does.

- **⚠⚠ SMOOTHING A MESH-SENSITIVE ESTIMATOR CAN WIDEN THE SPREAD, NOT NARROW IT.**
  Measured: raw `ell_eff` spreads **9.6 %** across 2/1/0.5 mm, the 1 s EMA spreads
  **17.1 %** — and the rate spread scales with it (2.57 % smoothed, ≈1.44 % raw by
  linear scaling). A time filter whose constant is not tied to the event interval
  biases each mesh differently. ⇒ **tie any such constant to the physical interval,
  or do without.**

- **✅ ONE RULE FOR THE SURFACE TEMPERATURE, TWO REGIMES: EXPONENTIAL WHERE THE FRONT
  RECEDES, LINEAR WHERE IT DOES NOT.** A receding front has an advection–diffusion
  balance and an exponential profile (hence the two-cell reconstruction). A **hole
  wall that cannot recede laterally** has a purely **diffusive** profile whose
  penetration grows without bound: it reaches `dz` at **1.45 s** (2 mm cells) and is
  **8.3 cells by 100 s**, so the top-cell profile is nearly **linear** — exactly what
  the existing `face` half-cell Newton already assumes. **Only ≈0.3 % of wall cells
  are ever in the regime where the linear form is wrong** (a 2 mm cell is "fresh" for
  1.45 s, during which a 3.2 m/h hole deepens 1.3 mm). ⇒ **the hole wall needs no new
  physics**, only the `face` branch called laterally.

- **⚠ A DIRECTLY-REMAPPED TEMPERATURE AND A BISECTED ONE DIFFER AT THE INVERSION
  TOLERANCE.** `InvertHtoT` recovers T from H by **bisection with a finite tolerance**
  (`:3493`). Any scheme that remaps H **and** writes T directly (to keep the kernel's
  `T_old` in step) creates two versions of T that differ at round-off. **Mixing them —
  interior cells re-inverted from H, ghost cells filled from the directly-written T —
  is visible ONLY at a lateral box seam**, which is the signature of a decomposition
  failure that is round-off-sized, appears after the first remap, and only with a
  lateral split. ⚠ In a laterally uniform column a reduction-order explanation is
  hard to sustain, since every per-column buffer holds identical data.

- **⚠ SHARE OF POWER IS NOT SHARE OF ERROR.** Side faces carry **62 % of `P_robin`**
  on Meier, but heat into the annulus wall **above** the floor is physically real and
  barely touches the floor rate. Before sizing any lateral-face work, bin it: by
  radius at the `side_face_h_r_min` = 40 mm gate (hole wall vs floor riser), by
  **whether the receiving cell later fires**, and by multi-cell riser share (the
  owning-column proxy is exact for a one-cell step and wrong beyond it).

- **⚠ A NEW BOUNDARY CLOSURE NEEDS ITS OWN `robin_form`, NOT A BRANCH ON ANOTHER
  KEY.** Selecting a new flux form on an unrelated flag (e.g. `sp_surface_T_apply`)
  silently makes the other forms **identical by construction**, destroying the very
  control legs a gate depends on. A new form keeps `cell`/`face`/`pinned`
  byte-identical for free, and every frozen study from `d2v` to `d2ad` used `pinned`.

- **⚠⚠ MESH-INDEPENDENCE AND ENERGY CONSISTENCY ARE TWO DIFFERENT GATES.** "Do not
  score the level against the continuum ablation limit" (right — a finite flake
  legitimately sits above it) is **not** "do not score the energy balance". **Spread
  alone passes a mesh-independent but energetically wrong scheme.** Always carry a
  scheme-internal balance (`v·ρc·(T̄_flake − T0)` vs the ledger's `q_net`) and a
  domain closure (`ΔE_H = ∫(P_robin − P_loss)dt − ΣE_flake`) alongside any spread
  criterion. For a flake of `h_f = c_f·ell` the level multiple is
  `1/((1−e^{−c_f})/c_f)` = **1.16 / 1.21 / 1.27** at `c_f` = 0.3 / 0.4 / 0.5.

- **⚠⚠ `surface.follow_mask = 1` IS MANDATORY on any case with active removal and a
  `surface_patch`.** Without it the patch is applied at the **fixed domain top**, so
  heating stops dead at the first void: measured, the new top cell sat at
  **295.9 K after 59 s** of nominal 0.55 MW/m². It silently invalidated the first
  build of the P0 1-D column and **anyone building from the Kant template will hit
  it.** Note `robin_form = pinned` already requires it at parse (`:585`); `cell` and
  `face` do not.

- **⚠⚠ A THRESHOLD CRITERION'S OWN VARIABLE IS A VACUOUS MESH WITNESS.** "It fires
  at 462 K on every mesh" proves nothing — a cell-average criterion fires at 462 K
  of *cell average* just as flatly (spread 0.01 %) while the **true surface
  temperature at its own firing** moves **14.06 %** (557.6 / 508.8 / 484.9 K over
  2 / 1 / 0.5 mm). Score the **physical state at firing**, not the criterion's own
  threshold variable. `spall.surface_T = diagnostic` exists for exactly this: it
  computes and logs the reconstruction while leaving the criterion bit-for-bit
  pre-P0.

- **⚠⚠ THE DEPTH SCAN RETURNS `h_scan ≈ dz/2 − a_f` REGARDLESS OF THE GRADIENT —
  UNRESOLVED.** D2z-0 takeaway 7 measured this in a `prescribed_T` control with a
  **1100 K face against a 655 K top cell** (a real gradient, 16d's fix active). The
  cause is the **GATE LOCATION**: `Sp_top >= 1` is tested at the cell **centre**
  (`z = dz/2`) while depth is measured from the **face**. ⇒ **never take a physical
  flake thickness from the scan** — take it from the temperature profile's own 1/e
  depth. Whoever revisits Step 16d starts here, and `prescribed_T` probably carries
  the same artefact, so excluding it may be unnecessary.

- **⚠ WHY D2z-0/-1 RETIRED THE SUB-CELL INCREMENT — the record, not the folklore.**
  The reason was a **plumbing artefact**, not a pad-rim cycle: `firings_per_void`
  **1.957 → exactly 1.000**, `n_voided` a singleton over **79,308 firings**,
  `h/dz = 1.00000`; the rate effect was only **+0.352 % / +0.438 %**. ⇒ reinstating
  partial recession makes `firings_per_void > 1` and `h/dz < 1` **return by design**
  (firing twice to shed a 1 mm flake from a 2 mm cell is *correct*), so those
  assertions must be **restated** — exactly one void per zero-crossing, one or more
  firings per void — and gated on the key.
  ⚠ **D2z-1's real lesson was that its override was OVER-BROAD**: the blast radius
  was `validation/hu/hu_end_to_end` (damage-law, removes rock, **never run**).
  **Scope any sub-cell increment to the `sp_weibull` computed-T path only.**

- **⚠⚠ A CORRECT SURFACE TEMPERATURE ALONE MAKES THE MODEL 80–130 % FASTER, AND THAT
  IS THE RIGHT ANSWER.** With `weibull.m = 20` the criterion is **flux-blind**, so
  re-anchoring on Kant makes the floor fire at `T_surf ≈ 750–830 K` everywhere;
  against today's implied 1215 K surface the rate ratio is
  `(1684−754)/(754−293) ÷ (1684−1215)/(823−293) = 2.28`. **A Meier rerun after the
  surface-temperature fix alone will look dramatically WORSE, with the axis running
  to the block base** — that is what removing a rate-suppressing mesh artefact looks
  like, and it is why the flux-dependent flaking temperature must land before any
  Meier judgement. **+10–25 % is the figure for the two together** (`T_s ≈ 970 K` at
  the axis flux). ⚠ Do not write a stop rule that treats the +80–130 % as a failure.

- **⚠⚠ A RECONSTRUCTED SURFACE TEMPERATURE CLOSES A POSITIVE FEEDBACK LOOP.**
  `v ↑ → ℓ = κ/v ↓ → dz/ℓ ↑ → T_surf ↑ → σ ↑ → K_I ↑ → fires at a lower `T_top` →
  numerator `(T_g − T_surf,fire)` ↑ → v ↑`, and with the present removal energy
  `ρc(T̄_cell − T0)` the **denominator shrinks** as firing gets earlier, so the loop
  is unbounded as `ℓ → 0`. **The stabiliser is energy-consistent removal**: let the
  flake carry `ρc(T_s − T0)` and the rate becomes `v = q/(ρc(T_s − T0))` —
  self-limiting and `dz`-independent. ⇒ **never ship the surface-temperature
  reconstruction without it**, and make the **1-D steady-state check
  (`v = q/(ρc(T_s−T0))` ±5 % at three meshes) the runaway witness.**
  ⚠ Also: a correct surface temperature makes the model **FASTER by ~20 %**, not
  slower — it is a correctness and mesh-independence fix, and **the level must close
  on the gas side** (blowing, `cp(T)`, the unmeasured gas temperature).

- **⚠ `weibull.V0 = 8e-09` IS V_cell AT 2 mm**, so `(V_cell/V0)^{1/m} = 1` **only
  at 2 mm**. Any mesh change silently rescales the flaw draw unless `V0` is held
  fixed. Hold it and say so in the packet. (V0 = V_cell was the user's choice; it
  switches the Weibull SIZE EFFECT off rather than merely removing mesh dependence.)

- **⚠ THE FLOOR GAS TEMPERATURE IS MEASURED, NOT ASSUMED: 1664.7 K.**
  Power-weighted over all 23 bins at t = 490 s on `d2ac_feed/output/FH_jet_profile.csv`
  (1717.3 K at r = 0 falling to 1624.9 K at r = 44 mm; floor power 1408.9 W on the
  **quarter** domain, x4 = 5.64 kW). The 09-29a brief assumed **1707 K**; every
  floor-factor table built on that number is ≈5 % out. Nozzle is 1900 K, so the
  march gives up **235 K** between nozzle and floor.
  ⚠⚠ **THE MARCH BITES BACK: cutting `h` does NOT cut flux proportionally.** The
  march conserves enthalpy, so gas that gives up less heat stays hotter
  (`T_g(F_q) ≈ 1900 − 235.3·F_q`). A 0.5x cut in `h` yields a floor factor of
  **0.570, not 0.50**. Any flux ladder must be pushed lower to compensate, and the
  realised factor must be **measured per leg**, never assumed.

- **⚠⚠ THE ROCK-SIDE STORY IS QUANTITATIVELY REFUTED AS A ROUTE TO THE FACTOR OF
  TWO.** Meier's block is Grimsel granite; the site rock is Central Aar granite /
  Grimsel granodiorite, and it is measured. The model's firing temperature
  **823 K (550 °C rise)** sits in the **middle** of every direct measurement —
  Wilkinson & Tester 670–820 K, Kant & von Rohr mean ≈870 K, and **Kant's own model
  for Central Aare granite uses a 550 °C rise at the failure plane**. Raising `T_s`
  to the top of the measured range buys **15–25 %**, not 2x.
  ⚠ **The measurably wrong rock properties push the WRONG WAY:** measured
  **E = 47–53 GPa against the model's 30**, and α rises steeply above 400 °C — both
  raise the thermal stress and **lower** `T_s`. Measured ρc (1.96–2.08 vs the
  model's 2.17 MJ/m³K) also pushes the wrong way.
  ⚠ **Nasseri's 0.43 / 0.22 MPa√m at 650 / 850 °C are post-treatment ROOM-T tests.**
  There is **no at-temperature granite K_Ic data above 400 °C** at all.
  ⇒ **The factor of two must come mostly from the FLOOR HEAT DELIVERY** (model
  ≈1.2 MW/m² over Ø85 at Meier's standoff; Kant spalls granite at 0.45–0.87;
  Meier's balance needs ≈0.5–0.6). **A required `T_s` > 920 K is a refutation of
  the rock-side explanation, not a calibration.**

- **⚠ KANT'S MEASURED `h_fl,min` = 500 ± 100 W/m²K IS A HARD FLOOR ON ANY FLUX
  CALIBRATION.** A calibration that reproduces Fig. 8.8 by driving the floor
  coefficient below the value at which granite is **measured** to spall is not a
  calibration — it relocates the error into the **standoff decay law**
  `pow(0.09/max(0.09,s),1.0)`. Score the effective `F_q·h` at each leg's own
  equilibrium standoff against 500 ± 100 and against Kant's 630–1010 band.

- **✅ THE EQUILIBRIUM-STANDOFF PREDICTOR IS CALIBRATED — BUT IT ONLY EXISTS ABOVE
  90 mm.** Beyond 90 mm the map decays as `1/s`, so `s* = 0.09·R90·F/1.5` m, i.e.
  **`s* ≈ 198·F` mm** at `R90 = 3.3 m/h`. At F = 1 it predicts **198 mm**;
  D2ac's `FH` extrapolates by its own geometric increment decay (ratio ≈0.70) to
  **192 mm** — 3 % on a leg it was not fitted to. Carry ±5 %.
  ⚠⚠ **INSIDE 90 mm THERE IS NO EQUILIBRIUM AT ALL, ONLY CONTACT.** The map is
  **flat** there (`max(0.09,s)` clips), so the floor recedes at a constant `R90·F`
  with **no restoring force**: if `R90·F < feed` the burner simply catches the
  floor. ⇒ **contact threshold `F = feed/R90` = 0.4545** at 1.5 m/h, and any `s*`
  the formula returns **below 90 mm is not an equilibrium — it is a contact
  prediction**. The usable window is therefore only **`0.4545 < F ≤ 0.505`,
  i.e. 10 mm wide in standoff**, against a predictor carrying ±5 % (±5 mm) ⇒
  **a reading of 95–105 mm is NOT a decisive pass or fail.**
  ⚠ `2026-09-29a.md` quotes **two mutually inconsistent `R90` values** — 3.3 in §1's
  `s*` formula and 4.2 in its contact edge (thresholds 0.4545 vs 0.357). **3.3 is
  the validated one**; 4.2 would predict 252 mm for `FH` against its measured 192.
  **Measure `R90` from the control leg's early phase rather than quoting either.**

- **⚠ THE STAIRCASE "WETTED" AREA IS THE EXACT RAMP `sqrt(1+δ²) − 1`, NOT `δ²/2`.**
  `δ²/2` is the small-step limit and **over**charges where it matters: at δ = 2 it
  gives **2.0 footprints against a true 1.236**, at δ = 3 **4.5 against 2.162** —
  and **over half of `A`'s stepped pairs are ≥2 cells apart**. `RiserWeight` already
  is the ramp formula; the continuous law is that formula with `δ` in place of the
  cell count. ⚠ Consequence: on multi-cell risers a full-face law legitimately draws
  **more** power than D2x's ramp, so "more power than `D` means δ is inverted" is
  false — split such checks **per population**.

- **⚠ DO NOT BAR PHYSICS WITH AN ESTIMATOR CHECK.** Twice now a D2z draft has gated
  a live quantity against a band measured on a **different configuration**: first
  at-birth vs time-weighted, then closure-leg δ against the `f = 0.2` legs' band. A
  closure that heals better legitimately drops **below** that band, so correct code
  fails. **Gate the estimator against the offline estimator on the SAME leg's own
  log (~10 %); report the cross-configuration comparison as physics.**

- **⚠⚠ A DEFAULT-OFF KEY HIDES THE BLAST RADIUS OF ITS OWN CONDITION.** D2z-0's
  whole-cell override was scoped `reg_choice == 1 && !presc_T` and looked clean for
  a whole packet **because the key was off on every input but its own two legs**.
  The unconditional flip in D2z-1 revealed it also caught the **damage-law /
  `spall.prescribed_ratio`** path, where a multi-cell thickness is *intended*:
  `unit/spall_event` failed five assertions and
  `validation/hu/hu_end_to_end/input_granite2` (real removal, 120/260 cells, **never
  run in the campaign**) would have shifted silently. ⇒ **Scope a condition to its
  mechanism (`sp_weibull_spall`), scope its assertion identically, and do not treat
  a green campaign under a default-off key as evidence the condition is right.**

- **⚠ THE HALF-CELL IS GONE; `beam_void_closure` T3 IS THE DEBT IT LEFT.** Under
  `h = dz` the sub-cell remainder is single-valued (`dz/2 − O(1e-16)`, drifting
  ~1.6e-16 over a run) so T3's histogram cannot populate both half-cells **under any
  parameter choice**. `prescribed_T` is unavailable in that test
  (`surface_patch.enabled = 0`) and a vapour-driven variant removed nothing.
  ⚠ The pre-change split was **171 / 2229** — the **upper** half, where the E1 leak
  lived, was the **normal** operating point. **Restoration needs the initial
  remainder to VARY ACROSS COLUMNS (a stepped or offset initial surface) or a
  genuine sub-cell removal path.**

- **⚠ THE TREE IDENTITY REFERENCE IS NOW `witness/d2z1.md5`**, not
  `witness/tree0922.md5` (preserved untouched). See `witness/REFERENCE.md`.
  ⚠ **An md5-only reference cannot audit direction** — the sweep overwrites the
  outputs it hashes, so a criterion wanting "counts up, never down" needs the sweep
  to record counts, not just hashes.

- **⚠ A SENTINEL THAT SHARES A TYPE WITH A MEASUREMENT WILL BE READ AS ONE.**
  `fit_rop`'s `< 3 samples → 0.0` floor was implicitly calibrated on *firings*;
  whole-cell halved the rows for the same rock and the test printed
  "0.0000 m/h, +inf %", which reads exactly like a column that stopped receding. The
  real rate was **0.7279 → 0.7178 m/h**. Use a distinguishable sentinel (NaN, or a
  separate count column).

- **⚠⚠ A RATE CRITERION HAS A TIME RESOLUTION. PRE-REGISTER THE WINDOW LENGTH,
  NOT JUST THE BAR.** D2z-0 predicted +0.36 % and scored it on 50 s blocks, where
  phase reshuffling is ±2–4 pts with no consistent sign. Block-length
  sensitivity with one estimator: spread **6.90 → 1.95 → 0.85 → 0.00 pts** at
  2 mm (25 → 50 → 75 → 150 s) **with the mean invariant**; 1 mm is tighter because
  four times the columns average the noise down — itself evidence it is noise.
  **Two reasonable estimators disagreed on the aggregate by ~0.4 pts, the size of
  the effect.** ⇒ **Direction and structure are resolvable at sub-percent effect
  sizes; magnitude is not.** ⚠ And `col_rate` **amplifies** it: raw cells voided
  read +0.08 / +2.03 / +0.07 % where the rate metric read −1.73 / +4.35 /
  −1.38 %. **If a criterion is about how much rock left, COUNT CELLS.**

- **⚠ A COLUMN FIRES TWICE PER CELL — event counts are NOT removal counts.**
  98.6 % of consecutive firings alternate void / no-void, and the gap is
  **16 ms** (one timestep) after a non-voiding firing against **3536 ms** after a
  voiding one. **Anything computed from removal-event counts or consecutive
  firing intervals is out by ~2x**; filter on **`n_voided > 0.5`**. A phase built
  on all firings gives 0.716 where the void-only answer is 0.034.

- **⚠⚠ A ONE-CELL STEP IN THE MODEL IS A SUB-CELL STEP IN THE ROCK.** A column's
  removal interval varies by only 1–11 %, so recession is linear between flips
  and the **true step height** is recoverable as a phase lag. Measured on the
  arbiter, r < 25 mm: a **cured** floor (`f = 0.2`) shows steps of
  **0.05–0.09 cells** with **no** ≥2-cell pairs; a **cascading** floor (`f = 0`)
  shows **0.17 → 0.45 cells** with 1–14 % ≥2-cell pairs, worse at 1 mm.
  ⇒ **`f = 0.2` over-charges a healthy riser ~3x and D2x's derived `w = 0.4142`
  ~6x.** Long-lived near-full-cell steps are **the cascade itself**, not a
  separate phenomenon.

- **⚠⚠ `A` IS ALREADY CASCADE-THROTTLED AT 2 mm — "`A` has no `f`, therefore no
  problem" is REFUTED.** Independently measured on `A`'s central floor
  (r < 25 mm): **53.2 % of stepped pairs are ≥2 cells apart**, median true step
  height **0.55 cells**. `D` has **0.0 %** at 0.37. **`A` has the worst steps of
  any run.** Its 1 mm behaviour is unknown and the run is **~6 h, not 3** (16x
  the cells) unless stopped at the scoring window. **Parked by the user.** Expect
  the rate gap to appear; the **shape** is the new information.

- **⚠ HEATING AREA AND CONDUCTION AREA ARE DIFFERENT QUANTITIES AT A SUB-CELL
  STEP.** The lateral interface between the shorter column's hot top cell and the
  taller column's cell is **fully solid on both sides** — no void — so that
  correction is **thermal** (the cell is treated as buried interior while sitting
  `φ·dz` below the surface). The **wetted** area of a smooth sub-cell lag is
  `dx*dy*(√(1+φ²) − 1) ≈ dx*dy*φ²/2`, i.e. **0.0024 footprints at φ = 0.07, not
  0.07 — thirty times smaller**. Any closure charging `φ` to both must **report
  the ratio**; if they do not scale together, healing is far weaker than drain.

- **⚠⚠ A FOURTH SIDE-PATH SITE, `:2354-2355`, AND ITS DEFAULT HIDES IT.**
  `side_f_m = side_area_exact ? 1.0 : surface_patch.side_face_factor`, and
  `side_face_factor` **defaults to 1.0** (`:6199`). So any *new* area key that is
  not `side_area_exact` falls through to a number that is **accidentally
  correct**, and since **neither `A`'s nor `D`'s argv sets `side_face_factor`
  (grep count 0)**, a whole campaign can pass without exposing the mismatch. It
  turns into a silent kernel↔march split the moment a run sets the factor.
  **The side-face area now has FOUR sites: `:1675`, `:2218`, `:2354`, `:2958`.**
  Each generation of this work has found one more than the last.

- **⚠ THE POWER IDENTITY IS THE ONLY CHECK THAT CATCHES A KERNEL↔MARCH SPLIT.**
  `jet_P_face = patch_P_robin + side_P_face + side_P_loss`. A per-column area or
  partition assertion **cannot substitute** — it stays true no matter what the
  march bills. Keep the identity as a **scored criterion**, not an inert-table
  witness, in any packet that touches side-face area or weighting.

- **⚠ THE KERNEL HAS NO EXPLICIT TOP-FACE AREA.** `:1606` is
  `div_kgrad += q_in * inv_dz`, i.e. `q * dx*dy / (dx*dy*dz)` — **the area
  cancels against the cell volume**. Any scheme giving the top face an area other
  than `dx*dy` must introduce a **new explicit multiplier**; there is no existing
  factor to edit. The companion sites are `:1723` (loss), `:2917-2918` (ledger)
  and `:2337` (march `dA`, applied `:2415-2416`, `:2437-2438`).

- **⚠ THE MEIER RUN KEYS ARE ON THE ARGV, NOT IN `input_drilling`.**
  `validation/meier/sp_meier_pilot/input_drilling` contains **none** of
  `surface_patch.radius`, `foot_r_inner`, `foot_r_outer`, `foot_standoff`, and
  **contradicts** three more: it has `timestep = 0.001` (runs use 0.016),
  `amr.n_cell = 70 70 60` (runs use 60 60 200) and
  `geometry.prob_hi = 0.14 0.14 0.120` (runs use 0.12 0.12 0.4). The real values
  live on `d2w_attribution/output/A.log`'s command line. **Cite the argv**
  (D2v takeaway 12; this recurred in D2y's first draft).

- **⚠⚠ THE ARBITER MESH GAP HAS TWO MASKS AND THEY DISAGREE BY AN ORDER OF
  MAGNITUDE.** `studies/d2r1_aboveplane_sides/analyze.py:39` masks on
  **r < 0.030**, where the `f = 0.2` baseline reads **−6.4 / −8.0 / −14.0 %**.
  The **interior, r < 25 mm**, reads **−0.2 / +3.0 / −1.8 %**. The difference is
  the **25–30 mm rim stall** — the mask-edge artefact that D2r-0/0b/0c and D2s
  all chased. **Never quote one mask's numbers against the other's bar**, and say
  which mask a mesh-gap figure uses.

- **⚠⚠ D2k HEALS STEPS BY ACCELERATING THE TALLER COLUMN, NOT BY FEEDING THE
  RECEIVING CELL.** For a one-cell step the exposed riser face sits on the
  **taller column's own top cell** (the neighbour's void is at that level), so
  that cell gets **1 + f = 1.414** footprints and catches up. The buried cell one
  level down — the one the shorter neighbour actually drains into — has **no
  exposed face under any scheme**. ⇒ **Any strictly area-conserving scheme
  removes the healing surplus and should be expected to reopen the mesh gap.**
  The remedy is then a **numerical closure separate from area**, or sub-cell
  continuous recession — **not a factor**.

- **⚠ SLOPE STATISTICS MUST BE QUOTED ON THE BELOW-PLANE FLOOR.** Over r < 76 mm
  including untouched block top the median `|∇z|` reads **0.00**; on the
  **below-plane floor** it is **2.24 on `D` and 5.70 on `A`**, with a quarter of
  `A`'s above 12.75. Single-cell risers there are **44.8 % of `D`'s, carrying
  32.7 %** of the weighted area (not the 64.6 % / 54.5 % of the wider mask).

- **⚠ A CLIFF TEST MUST BE SCALE-FREE AND `n_min` A FRACTION OF THE DISC.** With
  `L = 6 mm` and an **absolute** slope cap of 2, **42.5 % of `A`'s below-plane
  floor columns** fail a fixed six-neighbour test and fall back to the raw
  treatment (cap 3 → 23.2 %, cap 5 → 7.4 %). `L = 12 mm` drops it to 3.4 % **only
  because 6 is fixed while the disc grows 28 → 112 neighbours**, so an `L`-scan
  moves an insensitivity criterion for a bookkeeping reason and it reads as
  physics. Use a multiple of the disc's **median** implied slope, and set the
  minimum as a **fraction** of the disc population.

- **⚠ THE RADIUS GATE SELECTS COLUMNS, NOT FACES.** D2x's `C`/`E` were gated to
  r ≥ 38 mm and still lifted the **30–38 mm** rings by 14–37 %, because a column
  inside the gate has risers facing inward. **The side term cannot be aimed at an
  annulus.** Any band-local treatment must gate **per face**, on the neighbour's
  radius.

- **⚠⚠ DO NOT ZERO JITTER RISERS WHEN FIXING THE AREA OVER-COUNT.** Single-cell
  risers are exactly where the **step-face drain** lives: a firing column's top
  cell drains laterally into the taller neighbour's buried cell, and heating that
  face is what **D2k** cured. The D2j-0b interior converged at `f` between 0.1
  and 0.2 and **did not converge at zero**. A scheme that conserves area by
  giving the top face its footprint and the risers "the remainder" starves them
  and **reopens the 1/dx cascade**. Distribute proportionally to raw face area
  over the top face **and** all exposed risers instead.

- **⚠ TWO AREA ERRORS THAT A PER-FACE WEIGHT CANNOT FIX.** **(1) Jitter** — on
  D2x's leg `D`, **64.6 % of floor risers are single-cell** carrying **54.5 %**
  of the weighted area, on a surface whose **median |∇z| is 0.00**; wall
  (r ≥ 38 mm) is 20.2 % / 11.9 %, i.e. real slope. **(2) The 2-D decomposition** —
  the code sums `(√(1+z_x²)−1) + (√(1+z_y²)−1)` where the metric is
  `√(1+z_x²+z_y²)−1`: **+5.0 % at slope 0.5, +13.2 % at 1, +23.6 % at 2**. Exact,
  mesh-independent, unaffected by smoothing. A riser exists only at slope ≥ 1.

- **⚠ SMOOTHING STENCILS MUST BE A FIXED PHYSICAL LENGTH.** With quantisation
  noise `σ ≈ dz` over `N = L/dx` points, OLS slope error is `σ√12/(L√N)`:
  **fixed physical `L` gives `(dx/L)^1.5` and converges; a fixed cell count gives
  `1/N^1.5` and never does.** Area error is third order **only where the true
  slope is ≈ 0**; on a real slope `d/dg[√(1+g²)−1] = O(1)` so it falls as
  `(dx/L)^1.5`. Removal events are spatially correlated through conduction, so
  the OLS rate is an **upper bound** — measure the order, never assume it.

- **⚠ THE SURFACE IS BIMODAL, SO ONE STENCIL CANNOT SERVE IT.** `|∇z|` over the
  core: **p50 = 0.00, p75 = 2.55, p90 = 15.8, p95 = 35.2, p99 = 45.0**; 14.4 % of
  columns exceed slope 5 (79°). A plane fit straddling a cliff is **fitting a
  discontinuity** and returns whatever the stencil width dictates. Exclude
  cliff neighbours and hand over to the discrete wall treatment, where a vertical
  face of height H carries exactly `H·dy`.

- **⚠ A WITNESS OVER A POPULATION MUST MASK ON THAT POPULATION.** D2x's
  `side_w_mean` first read **0.1471**, below a floor the packet called impossible,
  because **79 % of faces late in a Meier run are above-plane** under
  `side_face_factor_above` — a different law — and were counted in the divisor at
  weight zero. The flux was always right. Its regression missed it because the
  `AEX` variant never set `side_face_factor_above`. **Mask on
  `side_faces − side_cols_above > 0`**, and test a witness with the exclusions on.

- **⚠ LEAVE NO GAP BETWEEN A PREDICTION AND ITS REFUTATION THRESHOLD.** Two of
  D2x's six landed in dead zones (`D` 3.80 against ">4.0 / refuted <3.5"; `C`
  2.09 against "1.3–1.9 / refuted >2.4"). ⚠ And **state which definition a rate
  bar uses**: `A` reads **1.251 m/h** as the slope over the 100–150 mm *descent*
  stretch and **1.296** over the *time* window 150–400 s — band edge in between.

- **⚠⚠ SCORE `Ø(z)` AT EQUAL BURNER DEPTH, NEVER AT EQUAL TIME — and only on rock
  the nozzle plane has passed.** The plane rides **50 mm above the pads**, so a
  wall at depth `d` is unfinished until the burner has descended `d + 50` mm.
  Equal-time scoring **structurally rewards whichever leg drilled fast enough to
  finish the ladder**, which is exactly the doubled-rate legs. This invalidated
  D2v's and D2w's deep diameter entries and produced the false claim that `A`
  "pinches to 65.7 mm at 200 mm depth" — at 401.6 s `A`'s plane had only reached
  ~99 mm. **The tell was already on disk:** D2v's own `analysis.txt` prints `M0`
  at 401.6 s and 500 s as **identical at 10–100 mm** and different only at
  125–200. Re-scored at equal burner depth: at descent 150 mm `M0` 28.40 /
  **`A` 2.27** / `B` 19.93 / `M1` 2.73 mm. Restates D2o-0's rule (2026-09-22),
  which had been violated. **Measure rate over an equal DEPTH interval too.**

- **⚠ `A` (the confined map alone) ALREADY MEETS MEIER ON SHAPE AND RATE.**
  2.27 mm at descent 150 mm — better than `M1`'s 2.73 — at **1.296 m/h** against
  `M1`'s 3.123, with **no side-face path and no free parameter**. The claim "the
  model cannot get shape and rate right at the same time" is **FALSE**. `A`'s
  real open items are the wall band (**84.4 / 82.7 mm at 100 / 125 mm depth**
  against 87), the collapsed rings at **40–48 mm**, a **pit lead of 127 mm**, and
  **volume 1.24x** per metre of descent. ⚠ `A` is **0.67x on the pit basis and
  1.24x on the descent basis** — always print both, the basis flips the verdict.

- **⚠ `f` IS ALGEBRAICALLY UNFALSIFIABLE AS WRITTEN.** `side_face_factor`
  multiplies an **area** inside a Robin term `h*A*(T_gas − T)`, where an area
  factor and an `h` factor are **indistinguishable**. So `f = 0.2` is a single
  number absorbing two unrelated errors — a wrong area *and* a wrong
  heat-transfer coefficient — and **no model output can separate them**. It must
  be replaced by a derived area, not scanned.

- **⚠ THE DERIVED SIDE-FACE WEIGHT HAS A HARD FLOOR OF 0.4142, SO `f = 0.2` IS
  BELOW WHAT THE GEOMETRY PERMITS.** A riser face exists only where the surface
  drops at least one full cell, so the local slope is ≥ 1. The honest per-face
  weight is `w = (sqrt(dx² + Δz²) − dx)/Δz` = `(sqrt(1+n²) − 1)/n` for an
  `n`-cell drop: **0.4142 / 0.6180 / 0.7208 / 0.7808 → 1**. (Derivation: the true
  surface area over one cell is `dy*sqrt(dx² + Δz²)`; the top face already
  carries `dx*dy`; the riser carries the remainder, spread over its `n` cells.)
  **The D2s value of 0.18 was wrong**; the measured values were 0.595 and 0.683,
  consistent with n ≈ 2–3. **`w` is PER DIRECTION** — each lateral neighbour has
  its own drop — so the existing `nfx`/`nfy` lumping must be split into four.

- **⚠ MEIER Ch. 3's 0.4–1.6 kW/m²K ARE IMPINGING-RIG VALUES — check before
  citing them as annulus numbers.** They come from the calibration rig with a
  solid-stream nozzle, not confined annular flow. The von Rohr group's **measured
  confined-annular-jet** coefficients are the right class of number for a wall
  heat-transfer law and are a separate citation.

- **⚠ THE SKIRT IS NOT IN THE MODEL.** Meier's burner has a **Ø56 mm bore with
  three exhaust slots**, so the pad ring at **r = 28–40 mm sits under 12 mm of
  steel** and is heated only through the slots
  (`Claude_markdowns/meier_thesis_notes.md:40`). The model treats that ring as
  **open impingement**, and it is exactly where the burner rides and where the
  rate is set. Prime suspect if the cone flank proves over-driven.
  ⚠ Do not confuse the **Ø56 bore (radius 28)** with the confined map's blend
  radius of **56 mm**, which is just `foot_r_outer` 40 + eight coarse cells.

- **⚠ `jet_P_face` IS THE TOTAL STREAM DRAW, NOT THE TOP-FACE POWER.** The
  enthalpy march gathers **side faces as well as tops** (its bin loop pushes
  cells with void lateral neighbours into `hits`). Verified identity at t = 400 s
  on D2v's `M1`: **`jet_P_face` 1717.5 = `patch_P_robin` 1063.8 (tops) +
  `side_P_face` 624.8 (sides) + `side_P_loss` 28.9**. So the side-face path at
  `f = 0.2` is a **+59 % addition to the top-face power**. Always report the three
  terms separately; quoting `jet_P_face` as "the tops" misreads every power
  comparison.
- **⚠ A CONFINED `h` MAP WITHOUT SIDE HEATING MAY COST RATE, NOT JUST SHAPE.** It
  places a cold wall immediately outside the pad ring — D2j-0b's and D2o-0b's
  drain configuration, whose only established remedy is a wide taper. **If a
  map-only leg falls below the no-map baseline, the side-face path was
  COMPENSATING that drain rather than adding power, and `f = 0.2` has been
  standing in for the wall boundary condition the model lacks.** That is a
  different successor from "derive `f`". Witnesses: pad-ring recession over
  r = 28–40 mm, `foot_stall_time`, `foot_carry_cols`, lateral ΔT at the wall.

- **⚠⚠ ROP ON THE MEIER GEOMETRY IS THE SLOPE OF `nozzle_z`, NEVER PIT DEPTH.**
  The burner rides on the feet, not the deepest point, so the two differ by a
  factor of ~2. Measured 150–400 s: **`M0` 1.50 m/h (IN Meier's 1.3–1.6 band)**,
  `M1` 3.12, `M2` 3.10. Quoting pit depth as ROP is the error D2v made.
- **⚠ THE BELOW-PLANE SIDE-FACE PATH IS NOT A PASSENGER ON MEIER.** At
  `f = 0.2` it adds **`side_P_face` = 649 W** against **`jet_P_face` = 1760 W**
  through the tops, and it is the obvious cause of the 2x rate — consistent with
  D2k's independently measured **2.7x at `f = 1`** on this geometry.
  **`f = 0.2` is therefore implicated in a factor-of-two rate error and is no
  longer a quiet placeholder.** Derive it; never tune it to recover a rate.
- **⚠ `side_face_factor` REQUIRES `side_face_flux` — setting it to 0 ABORTS**
  (`:6099`). A "side-faces-off" control is the keys being **absent**, not zeroed.
  And `side_face_factor_above` must be set to **0.0 explicitly** whenever the
  side path is on, or the sentinel falls back to `side_face_factor` and the leg
  **cools** its above-plane faces (`h = 0` against 293 K ⇒ `q_in < 0`).
- **THE ABOVE-PLANE ENTHALPY DEFECT, QUANTIFIED.** `JetEnthalpyMarch` skips
  `!(flame_s[c] > 0.0)` at **`:2138`** and **`:2224`**. `M2` draws **1462 W**
  above-plane yet its exhaust falls only **12 K** vs `M1` (1528 -> 1516 K), where
  a counted draw at `m·cp` ~ 12.2 W/K would cost **~120 K** ⇒ **~90 % uncounted**.
  ⚠ "Free" overstates it: the physical annulus budget is of order 4 kW/quarter, so
  a budgeted march **reduces the widening, not removes it**. The global
  `ledger_err` cannot see this — energy is conserved domain-wide; the **jet** is
  over-drawn. **Adopt the wellbore formalism** (Ramey 1962; Hasan & Kabir) for the
  axial annulus balance rather than inventing one.
- **REPORT EXCAVATED VOLUME ON BOTH PRO-RATING BASES.** By pit depth `M1` is
  0.82x Meier; **by burner descent `M1` is ~0.97x and `M0` ~1.31x**. The basis
  changes the verdict, so print both or neither.
- **REFERENCE ERRATUM:** the Meier burner is **methane/air at λ = 1.2**
  (2.47 kg/h CH₄, 52 kg/h air), **not** methane–oxygen. Correct any write-up that
  says otherwise.

- **⚠ `jet_free_surface = 1` DOES NOT ABORT UNDER PRESCRIBED DESCENT — it silently
  misbehaves.** The free-surface test is `r > surface_patch.foot_r_outer`
  (`:2223`), and `foot_r_outer` **defaults to 0.0** (`:5962`) and is parsed **only**
  in the feet branch (`:6317`) ⇒ **every column is treated as free surface.**
  `foot_body_clearance` *does* abort (`:6287`), and `nozzle_z0`/`nozzle_feed` are
  **forbidden with feet** (`:6312`) hence required without. Any prescribed-descent
  leg built from a feet key set must **strip `jet_free_surface` and every
  `foot_*` key explicitly** — the rest are strict-checked as unused.
- **⚠ THE TAPER'S DIAMETER FLOOR APPLIES AT SHALLOW DEPTH ONLY.** A wall column at
  r = 46 mm still sees ~430–455 W/m²K **below** the plane and recedes at ~1.16 m/h
  against a 1.5 m/h feed, so it stays below the plane and cuts ~160 mm deep. The
  ~96–100 mm floor therefore governs 10–50 mm depth; deeper the hole can reach
  88–96 mm. **`M0` measures 86 mm at 200 mm depth with no blend at all**, so a
  blanket "cannot cut below 92 mm" clause fires on the control. **Gate the floor
  shallow; report deep Ø.**

- **⚠ THE FLOOR/WALL h MAP MUST BRANCH ON `s`, NOT BLEND ACROSS IT.** A single
  radial blend leaves the impingement term (clamped to its 15 mm stand-off value)
  reaching above-plane columns: at r = 40 / 42 / 45 / 47 mm a 40→56 blend gives
  **557 / 467 / 350 / 282 W/m²K**, i.e. **483 / 405 / 303 / 245 kW/m²** against a
  floor `q_pin` of ~520. **That drives the hole out to the hand-over radius
  regardless of physics** — the packet would decide its own answer. Floor branch =
  blend; **wall branch = flat `h_ann`**; mouth taper multiplies both.
- **⚠ THE INNER BLEND MUST START AT THE SKIRT RADIUS `foot_r_outer = 0.040`.** The
  pads sit at 28–40 mm and the burner rides the highest rock under them. Starting
  at 36 mm costs the outer pad rock **10 % at 38 mm and 20 % at 40 mm**, which
  slows its recession and slows the feed — **D2o-0's jam mechanism in mild form.**
- **⚠ AN 8-COARSE-CELL TAPER IMPOSES A DIAMETER FLOOR.** A 16 mm blend from 40 to
  56 mm hands over at a centre of 48 mm ⇒ **Ø cannot go below ~96–100 mm at
  2 mm**, while Meier reads 88 at 50 mm depth and 85 at 200. **The model stays
  wider than Meier at depth BY CONSTRUCTION**; narrowing the taper needs a finer
  mesh (8 mm on 1 mm). **Taper width is a numerical parameter, never a dial.**
- **⚠ `input_feet_d2c` IS NOT THE FROZEN D2i CONFIGURATION.** It differs in two
  keys: `pinned_idle_cycles` **2.0 vs 1.0e9** and `foot_body_clearance` **1 vs 0**.
  Any baseline meant to reproduce the frozen Ø(z) table must import the
  `d2i_gate` harness key set, and the check should be **byte identity of the event
  log**, not a millimetre tolerance.
- **⚠ A FEET-DESCENT RUN NEVER EXECUTES THE ABOVE-PLANE BRANCH.** The nozzle rides
  50 mm above the feet, so at t = 0 **every column has `s = +50 mm > 0`**. Any
  short trace or smoke test of above-plane behaviour needs **prescribed descent
  with the plane already below the original surface**, enthalpy closure on,
  collision guard off.
- **Recession beyond r = 76 mm on the frozen run is mean 0.08 mm but MAX 4.00 mm**
  (a stray two-cell event), and **exactly 0 beyond 90 mm**. A blanket guard must
  therefore be scored on the **mean**; a max-based bar fails on the baseline.

- **⚠⚠ THE DOMINANT MEIER SHAPE ERROR IS A FUNNEL — THE MODEL IS ~45 mm TOO WIDE
  AT THE MOUTH** (measured 2026-09-23 with the interpolated metric on frozen runs
  `d2i_gate/LI0_2mm` at 572 s, `d2m_feet/R40` agrees at 250 s):

  | depth [mm] | 10 | 25 | 50 | 75 | 100 | 150 | 200 |
  |---|---|---|---|---|---|---|---|
  | model | **141** | 128 | 113 | 105 | 99 | 92 | 86 |
  | Meier | 96 | 92 | 88 | 87 | 87 | 86 | 85 |

  Mean Ø over 50–200 mm = 99 mm, never below 86. ⚠ **Do NOT confuse this with
  D2m's "min Ø tracks 2 x foot_r_outer = 80 mm"** — min Ø is the *narrowest*
  section, near the drilling front (`R40` reads 55 mm at 150 mm depth). Any
  mechanism that **widens** the hole points the wrong way at every depth the
  reference covers.
- **⚠ A SCALAR ABOVE-PLANE `h` IS A BLANKET — 89 % OF COLUMNS ARE ABOVE THE
  PLANE.** At 572 s (`nozzle_z` 214 mm, `foot_z` 164 mm): **3210 of 3600 columns
  above the plane, 2294 never touched**, all in-patch. A uniform `h = 142` at
  ~1690 K delivers ~110 kW/m² against ~24 of losses and fires undisturbed block
  top. D2j-0b's "warm tops are inert" does **not** protect it — `C2` ran h = 100
  at 1000 K ~ 20 kW/m², *below* losses. **The above-plane state must be the same
  radially confined expression.**
- **⚠ `side_face_factor_above` MUST BE SET TO 0.0 FOR A CLEAN CONTROL.** Unset,
  the sentinel falls back to `side_face_factor`; with `h = 0` against 293 K the
  Newton solve returns `q_in < 0`, so the leg **cools** its above-plane faces
  instead of leaving them alone.
- **⚠ `input_feet` IS STALE: 0.3 m deep at 4 ms.** The frozen D2i/D2m runs use
  `prob_hi = 0.12 0.12 0.4` with **`timestep = 0.016`** (proven at 2 mm), and
  **`input_feet_d2c` (0.4 m) already exists**. Using `input_feet` makes any
  runtime estimate ~4x too low and runs out of domain.
- **THE MEIER MAGNITUDE BAND IS 85–95 ± 5 mm, and the csv is SHAPE ONLY.** The
  digitised width is a **visible width = a LOWER BOUND** on Ø; the thesis says
  85 mm (Tab. 8.2) while 3.42 L / 0.5 m implies **93 mm**.
- **`side_face` APPEARS ZERO TIMES IN `input_feet`** — `side_face_flux` defaults
  off, so **neither D2k's side-face path nor D2r-1's lever is enabled on Meier**.
- **⚠ INTERPOLATE THE CROSSING RADIUS BEFORE ANY Ø CLAIM.** D2t's `D(z) =
  2*r_centre` snaps to a **4 mm lattice**; 4 mm is **4.7 % of an 85 mm hole** and
  equals the Fig. 8.8 reference's own **±5 mm** precision — the entire
  discrepancy budget, spent before any physics.
- **THE MEIER SHAPE REFERENCE** is `validation/meier/meier_fig8_8_hole_profile.csv`
  (digitised D2c 2026-09-17): 96 mm at 10 mm depth, 92 at 25, 88 at 50, 87 from 75
  to 125, 86 at 150, 85 at 175–200, 84 at 250–300. ⚠ It is a **VISIBLE WIDTH, so a
  LOWER BOUND on Ø**, at **±5 mm**. Thesis D = 85 mm (Tab. 8.2); 3.42 L / 0.5 m
  implies a mean of **93 mm**. A **mild collar** (~10 mm over the top 75 mm), not
  a funnel. **Never quote agreement tighter than ±5 mm.**
- **⚠ THE s-SWITCH IS NOW THE SHARPEST CORNER IN THE MODEL.** On Meier `h` steps
  **0 -> jet_h_above** at the nozzle plane — infinitely sharp in relative terms
  against the arbiter's 4.9x — and D2t measured the latch directly (crossings
  collapse 146 -> 21 and 124 -> 5; columns cross once and stop). This is D2o-0b's
  sharp-corner pathology in `s`. **Measure it before smoothing it**; smoothing and
  scoring shape in one packet is two variables.
- **A SINGLE MESH CANNOT SEE THE LEVER'S EFFECT.** It appears as a *mesh gap*
  (−49 % vs −6 % in D2t's band). Any 2 mm-only Meier run is a **screening** run:
  no convergence claim, no shape verdict.

- **⚠ THE CORE'S MARGIN OVER THE NOZZLE FEED IS ONLY ~5 %, AND `h` STEPS 4.9x AT
  THE PLANE ⇒ COLUMNS LATCH.** Measured (block 150–200 s, `r < 18 mm`): core
  **1.6349 m/h at 1 mm / 1.6667 at 2 mm** against feed **1.5502**. Below the plane
  `h = 700` and a column drills faster than feed, so it **stays** below; above the
  plane `h = 142` (`if(s>0.0,700.0,142.0)`) and it drills far under feed, so it
  **falls further behind**. The two states are separated by the feed rate. ⇒ **any
  reduction of `h` beyond ~5 % creates a permanent above-plane wall**, and an
  above-plane **column count** cannot distinguish "stalled at a cliff" from
  "lagging on a taper". This is D2o-0b's sharp-corner pathology in `s` rather than
  `r`, and it is what the frozen `B_1mm` rim stall (111 columns) already was.
- **⚠⚠ THE RESIDUAL MESH GAP IS THE MASK EDGE, NOT THE INTERIOR (2026-09-23,
  verified against the frozen logs, no new runs). This REFUTES the two notes that
  stood here before it** (a √3 area bound implying ">= x3.4 true-area growth", and
  a "surface roughening" reading of the slope data). Both are withdrawn; the
  record is in `ARCHIVE_DONE.md`. Decisive decomposition, 1 mm / 2 mm − 1:

  | leg | r < 30 (scored) | r < 25 core | 25–30 rim |
  |---|---|---|---|
  | D2j-0b `C1`, sides off | −13.5 / −28.7 / −45.5 | −7.5 / −27.5 / −47.3 | stalled both meshes |
  | D2l `f = 1.0` | +3.4 / +1.5 / −3.2 | +4.0 / +1.2 / −3.8 | +1.9 / +2.2 / −1.7 |
  | D2l `f = 0.2` (= D2r-1 `B`) | −6.4 / −8.0 / −14.0 | **−0.2 / +3.0 / −1.8** | **−21.4 / −33.3 / −44.1** |
  | D2l `f = 0.1` | −10.6 / −20.2 / −27.2 | −1.4 / −3.7 / −6.8 | −39.0 / −73.0 / −97.6 |
  | D2r-1 `T`, sides on | −2.3 / +0.7 / −2.0 | **−0.0 / +3.3 / −1.5** | **−7.8 / −5.2 / −3.7** |

  **D2k/D2l's side-face flux at `f = 0.2` cured the interior cascade.** Everything
  since — D2r-0, D2r-0b, D2r-0c, D2r-1 and the withdrawn D2s — chased the outer
  5 mm annulus.
- **"ABOVE-PLANE" IS THE MASK RIM, NOT THE DRILLING FRONT.** At t = 250 s
  `z_nozzle` = 62.4 mm and the front sits **18–21 mm BELOW** it; nothing in the
  interior is above the plane. In-patch columns with `s <= 0`: **1 at 2 mm**
  (r = 29.83) and **13 at 1 mm** (r = 29.71–29.91), **11 never spall**.
  `side_cols_above` 103 -> 2441 is **mask-edge quantisation** — 1 vs 13 column
  centres in the last 0.3 mm of the r = 30 mm mask, times a cliff 44 cells deep at
  2 mm and 88 at 1 mm (1x44x2 ~ 103, 13x88x2 ~ 2441). In the baseline `B_1mm`,
  **111** in-patch columns are above-plane (r = 27.3–29.9 mm) — a stalled annulus
  that D2r-1's key un-stalled.
- **⚠ D2r-1's CRITERION (e) WAS MISREAD — THE NOTE THAT STOOD HERE IS WITHDRAWN.**
  `T` is **flat at both meshes** (`T_1mm` 1.5634 / 1.5682 / 1.5861 / 1.5846;
  `T_2mm` 1.6001 / 1.5571 / 1.6185 / 1.5969). `B_1mm` **decays** 1.4957 / 1.4314 /
  1.3920 / **1.3309** as its rim stalls. **The "+19.06 % monotone climb" is the
  denominator sinking, not the test drilling faster.** Nobody bought convergence
  by drilling faster; the ceiling was scored against a decaying baseline.
- **⚠ A STAIRCASE FACE NEEDS A >= 1-CELL DROP, so any column owning one has local
  slope >= 1.** Face-weighted `w = (sqrt(1+gx^2+gy^2)-1)/(|gx|+|gy|)` on D2r-1's
  own surfaces is **0.595 (2 mm) / 0.683 (1 mm)**, per-edge 0.73 / 0.91, **min
  exactly 0.414** — not the 0.11–0.18 obtained by averaging over flat columns that
  own no faces. **And central differences vanish on exactly the peak column that
  owns the faces**, so that stencil is mis-specified regardless. Any effective-area
  revival needs a **per-edge** weight from the drop across each face.
- **⚠ `f = 0.2` IS RETIRED AS A MESH-GAP QUESTION BUT NOT AS A PHYSICS ONE.** The
  core convergence is f-dependent (`f = 0.1` -> −1.4 / −3.7 / −6.8; `f = 1.0` ->
  +4.0 / +1.2 / −3.8); `f = 0.2` sits in the good spot and **was never derived**.
  It becomes a live magnitude question the moment a wall is scored. Retire it by
  derivation, never by a scan.
- **EVERY METRIC IN THIS PROJECT IS A DISK-MEAN RATE, AND THAT IS HOW THE ARTEFACT
  SURVIVED FOUR PACKETS.** A single averaged number cannot say *where* two meshes
  disagree. D2t builds the first shape metric: profile `z(r)` on **fixed physical**
  bins (a mesh-dependent binning makes the metric self-referential) plus
  diameter-versus-depth. **Do not choose the scoring radius after seeing which
  radius converges** — fix bins and bar in `CRITERION.md` first.
- **The heated side faces are the STEP FACES OF THE DRILLING FRONT, not a pit
  wall.** The side-face branch sits in the full 3D `(i,j,k)` loop
  (`MMWSpalling.H:1596–1658`), not restricted to `k == k_top_c`, but a face needs
  a **void lateral neighbour**: inside the disk all columns recede together and
  outside it (r > 30 mm) nothing is ever removed, so the disk's outer boundary
  yields **no faces at all**. The 103 faces at 2 mm are the front's roughness —
  exactly where the drain flows. **"Above-plane" columns are the LAGGING ones**
  (`s = z_nozzle − z_face`, `:1759`), i.e. the taller neighbours the drain pours
  into.
- **D2r-1's (e) ratio does reach +19.06 % in the unscored 200–250 s block** (the
  runs go to 250 s; the frozen metric stops at 200 s), so **every follow-on must
  score that block**. ⚠ **But the reading "the test drills progressively faster"
  is WITHDRAWN** — see the corrected note above: `T` is flat at both meshes and it
  is `B_1mm` that decays as its rim stalls.
- **⚠ D2r-1 REDUCED THE DRAIN BY 6.8 % IN WATTS, not by eliminating it.** At 1 mm
  the remedy **adds 197 W** and lowers the drain **75.4 → 70.3 W**; at 2 mm it
  adds 34.5 W for 50.0 → 49.1 W (−1.9 %). The gap closes by **compensating**. The
  per-face mean `wall_step_q` *rises* on the test leg only because the face
  population changes — **compare power, never the per-face mean.**
- **Hash `CRITERION.md` BEFORE THE FIRST LEG, not before the scored ones.** D2r-1
  hashed at 14:23:55 with its three 2 mm legs already finished; disclosed, and
  mitigated by the scored criterion depending only on the 1 mm legs and by the
  implementer pre-registering a failing prediction. The rule exists so the
  question never arises.
- **The AZIMUTH is unused and deliberately deferred.** A 3D surface normal has
  two degrees of freedom; `√(1+gx²+gy²)` absorbs the tilt, but **which way a face
  points is used nowhere** — `:1656` treats +x, −x, +y and −y faces identically,
  so a face looking into a dead pocket is heated like one facing the flow. On a
  rough drilling **front** the flow direction is genuinely ambiguous, so any
  projection is a new modelling assumption and deserves its own packet.
- **The 0.5 mm third mesh is PARKED by the user (2026-09-23).** It remains the
  decisive convergence-vs-cancellation test for D2r-1; two meshes cannot separate
  the two. Expect 2–6 h. Revisit when an expensive slot is available.

- **⚠ STANDING PRE-FLIGHT ITEM 0, bought by THREE consecutive dead packets
  (2026-09-23): name the call site of every new or changed behaviour and
  demonstrate, from each leg's key set, that it is REACHED.** D2r-0 anchored a
  criterion "against C0" without checking C0's keys; D2r-0b specified a cap the
  code ignores in the mandated state; D2r-0c specified a pass inside
  `JetEnthalpyMarch`, which the arbiter never calls. **All three are the same
  error: specifying behaviour without verifying the code path executes under the
  configuration being mandated.**

- **ABOVE-PLANE TOP HEATING IS ALREADY KNOWN INERT.** The arbiter's own
  `CRITERION.md` defines **`C2` = h 100 / T 1000 above the plane, run at BOTH
  meshes**, purpose *"does keeping the exposed rock warm remove the cascade?"*,
  with P4 pre-registering that it should behave like C0. **It did not.** But
  D2j-0b predates D2k, so **`C2` heated only the TOPS** — above-plane **side**
  faces have never been heated by anything. **Do not re-test warm tops.**

- **THE IDENTITY SET DOES NOT EXERCISE `side_face_flux`.** No input in it turns
  the key on (every match in the tree is inside `output/` artefacts), and
  D2q-2b/2c both modified that path. ⇒ **2639-file byte-identity does NOT validate
  the frozen D2l `A_f02` baseline.** A 3.2-min re-run and byte-compare is the
  proof, and it is cheap.

- **`f` AND THE CAP ARE TWO FUDGES ON A RAW FACE COUNT.** The lateral flux is
  `div_kgrad += side_f * (s_in − s_out) * (nfx*inv_dx + nfy*inv_dy)` (`:1633`),
  where `nfx + nfy` **counts exposed cell faces**. The staircase has **too much
  lateral area and too little top area** (opposite signs, neither improving with
  refinement). D2k's note says *"`f` corrects the staircase area, not `h`"*, but
  D2l scanned it as a free dial. **The real fix is an effective surface area from
  the surface normal — `sn_cos_col` already exists (`:4278`, with a cos θ floor)
  and is applied only to the K_I sampling depth. That is packet D2s**, deferred
  because it changes every existing run and would destroy the frozen baselines.

- **CHECK A PREDICTION AGAINST MEASUREMENTS ALREADY ON DISK BEFORE HASHING IT.**
  D2r-0's model predicted a 9× drain reduction at h ≈ 150; **D2l had already
  measured f = 0.2 (≈ 140) → gaps 6/8/14 %, still failing**, and D2q-2c at 80 →
  11/22/30 %. Reconcilable only if **above-plane heating is the decisive term**,
  which nothing has tested — supported by D2q-2c's **no-plane control converging
  to 3.7/4.5/2.0 %** at the same h, differing in nothing else.

- **THE WATER-COOLED SHIELD BIASES BOTH SCORED QUANTITIES, SO IT GOES IN.**
  ~0.042 m² at ~350 K vs rock's ~0.058 m² at ~700 K ⇒ **roughly half the annulus
  heat really goes to the cooling water.** Omitting it inflates wall temperature
  (helping the drain result) **and** absorbed power (hurting the rate ceiling) —
  it is not a neutral simplification.

- **DERIVE A BAND'S ENDPOINTS FROM THE RUN'S OWN STATE, NOT FROM AN INPUT KEY
  (D2r-0, 2026-09-23).** The planner read `jet_mdot = 1e-3` off `input_feet`; the
  key is **3.7826e-3** and the **marched** mass with entrainment is **9.751e-3**
  (mcp = 12.188 W/K, m_ratio 2.58). That put the annulus in the wrong regime
  (laminar instead of **Re = 5568**) and set an indefensible band floor, which
  STOPped a packet on arithmetic rather than physics.

- **AN ANCHOR MUST NAME EVERY KEY THAT DIFFERS (D2r-0).** "Against C0" was not a
  specification: `d2i_gate` does **not** set `side_face_flux` (only `d2l_scan`
  does), so a drain criterion anchored to it would have scored the **D2k
  side-face path**, which already gives 2.23×, rather than the new mechanism.
  **A control must differ from the test leg in exactly one key.**

- **THE ENTHALPY BOUND DOES NOT BIND INSIDE THE HOLE (D2r-0).** Over the whole
  measured wall (0.0310 m² quarter, a **360 mm cone**, mouth r ≈ 78 mm) the
  stream gives up only **8.5–33 %** and reaches the mouth at **1223–1571 K**. Any
  restriction on which rock the gas may heat is therefore **a modelling choice
  that must be argued**, not a cost bound — the argument being that past the
  mouth the flow is no longer annular and Meier's sealed wellhead ducts it away.

- **CALIBRATE A WALL MODEL ON A MEASUREMENT THE PROJECT ALREADY HAS (D2r-0).**
  Two closed forms failed (semi-infinite Robin without losses → 915–1270 K, above
  firing; lumped balance → 790 K where D2o-0b measured 418 K). Anchoring to
  D2q-2c's **630–650 K wall at h = 80, T_gas = 1600** gives a quasi-steady rock
  sink **S = 190 W/m²K** that reproduces the anchor. **`T_fire` = 822.2 K is
  reached at h = 186**, so below that the wall is **warmed, not spalled**.

- **THE REFRAMING (2026-09-23): warm the wall, do not remove it — and widening
  needs no undercut.** 418 K → 700 K cuts the drain from 302 to ~90 kW/m²
  (58 % → 17 % of q_pin) with no removal mechanism. Meier's 80 mm skirt in an
  85–93 mm hole never covers rock outside r = 40 mm, so that rock comes off by
  ordinary **vertical** recession; D2m already measured recession at r = 45–50 mm
  (merely too slow). ⇒ **D2q's per-cell removal stays built, gated and OFF; its
  justification was probably wrong, its by-products stand.**

- **`jet_P_exhaust` IS NOT LOST, IT IS UNCONTACTED.** Under
  `jet_T_ent_mode = exhaust` the leftover stream recirculates into the next
  step's jet, so the ledger closes. The defect is that it never touches the wall
  on the way out. **Do not describe it as an energy leak.**

- **(d)'s `k·ΔT/dx` diagnostic is STILL OWED** from D2q-2c, against D2o-0b's
  302 kW/m² / 58 %. It was requested regardless of verdict and is recoverable
  from existing plotfiles.

- **THE WATER-COOLED SHIELD IS A MISSING COMPETING SINK.** Meier's burner body
  runs 160 L/h of cooling water in the same annulus: ~0.042 m² at ~350 K against
  rock's ~0.058 m² at ~700 K ⇒ **roughly half the annulus heat would go to the
  cooling water**. Any wall-heating result is therefore an **upper bound** until
  this is modelled.

- **`test_feet`, `test_feet_d2c` and `unit/amr_microstructure_regrid` are
  CHECK-ONLY** — they re-read stored plotfiles and run no simulation, so they do
  not exercise the current binary. Do not count them as coverage.

- **CONTIGUITY TESTS ARE THIS LINE'S RECURRING BUG CLASS (D2q-2b, 2026-09-22).**
  "The cell above is void" was used as a proxy for "this is the column top" in
  **two** separate places — beam invariant (a) and the `P_inc` tally — and they
  break in **opposite** directions once an undercut exists: one over-reports
  violations, the other **double-counts incident beam power by 14.29 %**. The
  second **shipped byte-identical in D2q-1 and was wrong**, which is the sharpest
  available demonstration that identity proves only "nothing changed while void
  stayed contiguous". **Grep for any remaining `rem(i,j,k+1)` used as a
  top-of-column test before trusting any per-cell result.**

- **A probe needs a guard that the state under test actually exists (D2q-2b).**
  `queryarr("seed_void", …)` lacked the `spall.` prefix every key in
  `ParseSpallSettings` carries; **AMReX only WARNS on an unconsumed key**, so the
  feature was silently off and two sharp assertions passed against a state that
  did not exist. Only the "is the seeded box actually void" check caught it.

- **The "don't rebuild during a sweep" rule needs harness enforcement, not
  memory (D2q-2b takeaway 7).** It was violated at the checkpoint and invalidated
  a sweep; the run was killed, cleaned and repeated, and nothing was reported
  from the corrupted one.

- **`surface.follow_mask = 1` is ALREADY SET everywhere this line scores**
  (`d2j0_hotdisk/input_hotdisk:135`, `sp_meier_pilot/input_feet:181`,
  `input_feet_d2c:187`), so D2q-1's phi-window `surface_mf` defect is live **only
  in `dev2d`**. Verify by grep; do not spend characterisation runs on it. Related:
  **`beam.P0 = 0.0`** in both the arbiter and Meier inputs, so
  `ColumnSolidAboveCount`'s `ncols·nz` beam-on gather is **not** on this line's
  critical path — time it in whatever packet next scores a beam run.

- **The 2639-file identity set is SINGLE-LEVEL, so the regrid-repair site is
  exercised by no run in it** — and `unit/amr_microstructure_regrid` is **not**
  the direct check, because it has no removal. A D2q-1 review statement to the
  contrary is withdrawn. AMR is parked, so this is a **recorded gap**, not a task.

- **A pointwise corner flux must not be a convergence gate once removal is
  per-cell.** The reflex-corner solution is *physically* singular (`r^(-1/3)` for
  an isotherm boundary), so `k·ΔT/dx` at one cell will never converge however
  correct the model is. **Gate on integrated quantities** — removed volume per
  unit perimeter, disk-mean rate — and report the pointwise flux as a diagnostic.

- **A default-off key keeps byte-identity available as the plumbing oracle.** A
  planner claim that a scaffolding step was "the last packet that can use
  identity" is **withdrawn**: only the mechanism *with the key on* is beyond it,
  and that is what the arbiter is for. Do not spend a separate build-and-sweep
  cycle to preserve an oracle that is not at risk.

- **THE BINARY sha256 IS NOT AN IDENTITY ORACLE ON THIS MACHINE (D2q-1,
  2026-09-22).** Rebuilding *unchanged, git-clean* source moved it
  (`62dea451…` → `83603563…`), and three passing 2639-file sweeps produced three
  different binaries. **Every packet that says "record sha before/after" means
  provenance only — the 2639 output files are the oracle.** Corollary from the
  same step: **batch the identity sweep** into bisectable runs (three at ~12 min
  localised one failure to a single site immediately), and **never rebuild while
  a sweep is running** — it swaps the binary mid-sweep and silently corrupts the
  result.

- **A pre-check that detects one direction of a bug does not clear the other
  (D2q-1 takeaway 4).** `surface_missing_cols` = 0 looked like clearance for the
  mask conversion and was wrong — that counter is blind to a surface cell that is
  not the top solid cell, which is exactly what a `phi_top ≈ dz` tie produces.
  **If a pre-check and the gate disagree, the gate wins.**

- **THERMAL LENGTH CORRECTED: κ/v ≈ 1.6 mm, not 3.2 (2026-09-22).** The 3.2 mm
  that propagated through the roadmap since D2j-0b was textbook granite
  (k = 2.5, ρc = 2.2e6, v = 1.3 m/h → 3.147 mm). The run's own keys
  (`kappa = 1.5`, `rho = 2750`, `Cp = 790`) give **α = 6.9045e-7** and, on the
  prescribed feed, **1.604 mm**. Both knock-ons make things **worse**: more steps
  exceed the thermal length, and **2 mm spans 0.8 cells of the ablation layer,
  1 mm spans 1.6** — a single thermal length is unresolved at 2 mm. **Always
  derive α and κ/v from the run's own input keys.**

- **Anchor a ceiling to an artefact-free reference (D2p-1, planner error).** A
  gate measured against a quantity the bug corrupts cannot be passed by fixing
  the bug: capping the 1 mm rate at +15 % over the *cascade-suppressed* 1 mm run
  made criteria (a) and (e) jointly unsatisfiable, with empty windows in two of
  three blocks — arithmetic available at plan time. Use the **no-plane control C0
  at the same mesh**; removing heating cannot make drilling faster, so it is a
  physical bound, and it still retro-rejects D2k at +52 %.

- **A criterion a fix satisfies by construction is not a criterion.** D2p-0's
  slab rule levelled every step as it formed; D2p-1's clamped 1 mm and 0.5 mm
  legs agreed because they shared a length. Both would have passed with nothing
  left to measure. **Check joint satisfiability, and check what the fix makes
  true for free, before hashing.**

- **Corner regularisation points the WRONG WAY and must not be revisited.**
  Reducing the lateral drain hands heat back to the floor: floor faster, wall
  colder, hole narrower — while Meier's wall is a few mm wider than his skirt
  *because it received that heat*. **The drain is heat delivered where the wall
  ought to spall; a fix must SPEND it there.** That is the acceptance test for
  every successor.

- **The height function cannot express lateral recession — and the continuous
  front does NOT fix that (D2p-0, 2026-09-22, design-time finding).** Void is
  derived from a per-column `phi` (`MMWSpalling.H:2818–2821`), so **solid above
  void in one column — an undercut — cannot exist**; and there is **no per-cell
  firing criterion** (`Sp_field_mf` written only at `k == ksurf`, `:4002–4004`).
  At a floor/wall corner the wall column is heated **at its base**, by the very
  conduction the drain delivers, so any faithful lateral rule needs exactly the
  undercut the state forbids: **it fails at removal, not at triggering.**
  **Refinement is not the escape** — the continuum solution is singular at the
  corner. **And the D2j ablation front is per-column** (`2026-09-21b.md` §1) and
  disclaims this in §0: *"The front does not fix under-resolved lateral
  conduction."* It fixes the **event** artefacts and is worth doing on its own
  merits — it is not a fix for the step face. **Do not name it as the fallback
  for a lateral problem again.**

- **Before a packet claims a mechanism "fits the existing representation", name
  the state that holds it and the line that writes it.** Two packets (D2o-0,
  D2o-0b) and one design fork were spent on assumptions that a grep would have
  refuted in an hour. Both facts above were one grep away.

- **G1-style field gates must be ONE-SIDED (D2o-0b review, 2026-09-22).** A
  two-sided "total within ±20 % of the control" cannot coexist with a packet that
  mandates a 20× wall cut: D2o-0b measured **0.49–0.72×**, a *deficit*, entirely
  outside r = 40 mm with the floor bit-identical. The gate exists to catch an
  **excess** (D2o-0's 2.9× blanket). **Write it as: no excess, no far-field
  blanket, far-field share not above the control's** — and use **one estimator on
  both fields** so bias cancels in the ratio (reading a plotfile's top solid cell
  is *not* the pinned surface temperature and runs **4.2–4.6× high**).

- **The mesh-divergence trigger must watch `k·ΔT/dx` at the transition against
  q_pin, not the wall's distance from the firing temperature (D2o-0b).** A
  "within 50 K of 821 K" test never fires when the wall sits at 417 K, yet the
  divergence risk was real and came from the lateral gradient at the corner.

- **Depth-mean Ø is not a single number for a cone (D2o-0b).** The control reads
  112.3 mm at a 97 mm averaging depth and 100.7 at 236 mm; the packet's quoted
  110.5 corresponds to ≈ 105 mm. **Score it at a matched depth with the control
  recomputed there**, or the target drifts with how deep the test drills. Fourth
  metric trap in this line, with D2m's three and D2l's "profiles cross".

- **The 85–95 mm width carries ~a quarter of the volume and is not clearance
  (D2o-0b).** A flat 80 mm hole at Meier's 1.30 m/h gives **1.82 cm³/s against
  his 2.47 — 26 % low**; 95 mm gives 2.56. So lateral recession is load-bearing
  for the **volume** score, not only the shape.

- **AMR is not a cheap parallel win, and will not buy back the 6 h fine legs.**
  Measured wall ratio is **1.00** — it buys no speed. The energy ledger aborts
  under AMR and beam+AMR silently loses deposition below the band. Parked by the
  user. Running the Meier AMR gates is a correctness exercise, not a cost fix.

- **Fig. 8.8 colour original: request it as a CHECK, not an anchor.** D2n closed
  for two reasons and the second is the binding one — the altered band is an
  **output** a correct model predicts, not an input to calibrate against.
  Treating it as an anchor was planner error. The closure must have **no free
  wall parameter**; the band is then a validation check on it.

- **Ø(z) may only be scored on rock the nozzle plane has already passed (D2o-0,
  2026-09-22).** `s ≤ 0` shuts the flux off above the plane, so rock below it is
  still being worked and its diameter will keep changing. **Any packet that
  scores a diameter must first report how much finished wall the run produced**
  (= the nozzle-plane depth at end of run), and **a run that finishes less than
  the scoring range has no shape result at all** — at any depth, by any
  estimator, matched time or matched geometry. D2o-0's jammed runs finished
  **0–5 mm** against the control's **186 mm**, and a planner re-read that missed
  this briefly claimed the shape had improved. Report the **pit-vs-feet gap**
  too (`patch_min_standoff`, and max hole depth against feet depth): D2o-0 ran to
  250 mm of hole with the feet at 55 mm, which is a runaway pit, not a wellbore.
  Sits alongside D2l's "profiles cross, so never score one depth" and D2m's
  moving-normaliser traps.

- **A prescribed field is not the control's field unless it is bounded in
  radius (D2o-0, 2026-09-22).** `jet_closure = none` removes the enthalpy march,
  and **the march is the only thing that both cooled the gas as it spread and
  bounded it in radius**. Hand-written fields must supply both. D2o-0's did
  neither: a quadratic fitted on r ≤ 68 mm, with a `max(1400, …)` added purely so
  the parser could not evaluate a negative, became a **1400 K gas blanket over
  every radius beyond 85 mm** — 73 % of the columns, **1.7–1.8 kW**, more than
  the control's entire 1639 W, and a total 2.9× the control's. Two rules follow:
  **never apply a fit outside the radius it was fitted on**, and **never let a
  parser-safety clamp set a physical value** (a guard against a non-positive
  branch is the *ambient* temperature, not a working gas temperature). Any
  prescribed-field packet must gate its launch on a pre-flight that prints total
  power and its radial split against the control's own march
  (**2711 W shallow / 1391 W deep**; shares 41.9 / 29.9 / 22.9 / 5.3 % over
  r 0–40 / 40–70 / 70–120 / 120–200 mm shallow, 64.1 / 35.9 / 0 / 0 deep).

- **The floor/wall corner must be radial, not on stand-off (D2o-0, 2026-09-22).**
  `FeetDescent` rests the pads on the **0.9 nearest-rank quantile** of annulus
  heights, so ~10 % of annulus rock always sits *above* the pad plane. **Any
  transition defined on stand-off therefore cuts through the carrying
  population**, and a carrier that lands on the cold side freezes (444–482 K
  against 821 K) and holds the burner up for the rest of the run — ROP
  1.491 → 0.108 m/h. A corner at `r = foot_r_outer` puts every carrier in the
  floor zone by construction. **This is a constraint on the D2o closure itself,
  not only on prescribed-field tests:** the transition must sit clear of the
  annulus height spread, or be smooth enough that no carrying column can freeze.

- **A shape criterion must bind on the whole profile, not its mean (D2o-0).**
  The control already satisfies a shaft-mean of **90.0 mm**, inside Meier's
  85–95, while flaring to 145 mm at the mouth. Likewise the control's **steady
  band gain is +5.0 mm/side, already inside** the 2.5–7.5 target — so band gain
  does not discriminate either. The metrics that do: **whole-profile Ø(z)**,
  **depth-mean Ø** (control 110.5 vs Meier ≈ 95), **mouth Ø** (control 145) and
  **volume rate** (control 4.58 vs Meier 2.47).

- **Wall metrics must be defined on geometry, not on stand-off (D2o-0).**
  `band = (s > 0) & (s ≤ 45 mm)` with no radius limit put **3298 columns of
  untouched far-field rock** into one case's "wall temperature", which is the
  whole of why it read 333 K. Band counts came out 239 / 3298 / 400 across three
  otherwise-identical variants purely according to whether the nozzle sat a few
  mm above or below the original surface at the snapshot. Require: outside the
  skirt radius, **at least one cell below the original rock surface**, and within
  the band height.

- **The "do not rebuild" block from D2k/D2l is LIFTED (2026-09-22).** The
  2639-file reference sweep passed on the current working tree
  (`studies/tree_sweep_0922/`, **2639/2639 byte-identical**, binary swapped
  rather than rebuilt and verified restored to `62dea451…`). The AMR delta is
  key-off exact. Do not re-queue the sweep; re-run it only after the next source
  change.


- **Idle clock (`pinned_idle_cycles`) is mesh-dependent (D2h, 2026-09-20).** Its
  timeout is measured in cell-heating times (t_cell ∝ dz), and the unpinned face
  form is self-limiting, so at finer meshes edge columns park sub-critical and a
  silent front walks inward. Set the key per input (1e9 = off); **changing the
  source default rebases every pinned-closure result since C1** and needs its own
  packet with the 2639-file reference set re-checked.

- **D2e retroactive review (2026-09-19, accepted) — carry to the commit.** The
  2639-file witness set is single-level MMWSpalling only; the Perf edits touch
  shared code (`BC::Constant`, `Numeric/Material/Table.H`) and an AMR regrid
  path, so the deferred full sweep must include the AMR microstructure test and
  other `BC::Constant` users. Nits: the `Table.H` seed can in principle pick a
  bracket one to the right at an ulp tie (step `lo` back while `m_H[lo] > H`);
  `Constant.cpp` `Domain().contains(box)` mixes index types (debug-build
  `sameType` assert on node-centred FABs).

- **How Meier ran the pilot test vs the model (2026-09-18, full Ch. 8
  reread; `Claude_markdowns/2026-09-18c.md`).** Decisions for D2g, before any
  scoring:
  1. **The mouth was inside a sealed, water-cooled wellhead** (Fig. 8.1,
     p. 213), glued to the rock, with exhaust to a cyclone and a bore of
     roughly twice the burner Ø (visual estimate; no dimension given). In
     phase 1 the gas over the surface was confined, cooled exhaust, not
     293 K room air. D2c's free-surface ambient dilution, which D2d credits
     with about half the mouth fix, rests on an open surface that did not
     exist. The wellhead also caps any funnel. **Confirmed by the user against
     independent thesis notes, 2026-09-20: this is a gas-side finding, not
     only a mouth-scoring one — it is the same physics as the tube-wall
     confinement, so Stage 2 must carry at minimum a free-surface-off variant
     as a sensitivity rather than treat dilution as settled.**
  2. **Two drilling phases** (647 s + 736 s, p. 219). Phase 2 restarted in a
     cooled, part-drilled hole (about 0.25 m, estimate) with the hole in
     rough vacuum and room air drawn in through the burner–wellhead gap. The
     model is one continuous run at 1 atm.
  3. **The thesis's ROP figures disagree with each other:** 0.5 m / 1383 s =
     **1.30 m/h**, against the text's 1.5. TSE 15.4 J/mm³ confirms 1383 s.
     Quote 1.30 as the time average, and state the discrepancy.
  4. **The measured ROP is an operator choice at or below the contact
     limit** (crane advance, hook load, clearance set by advance rate,
     pp. 214–223). The model's contact rule is an upper bound. Read the ROP
     one-sided: a model ROP below the measurement is a miss, while a
     somewhat higher ROP with a somewhat narrower hole is consistent.
  5. **The cooling-water ΔT was instrumented** (thermocouples up- and
     downstream, p. 208) but not reported, so the request to Meier is for
     recorded data. Also ask for the wellhead bore and height, the depth at
     the end of phase 1, and what the 1.5 m/h is.
  - **D2f wording corrected 2026-09-18** (the gate is unaffected):
    - "jolts clear isolated high rock" and "gas reached under the rim" are
      our inferences;
    - "wider than the tube at every depth" is overstated (Tab. 8.2: 85 mm;
      Fig. 8.8 visible width 80 → 78 mm near the exit, a lower bound).
- **D2f planning decisions (2026-09-18, revised same day).** A shielded
  "ledge" rule (rock under the rim gets no gas; broken off when undercut) was
  planned and **withdrawn before implementation**: Meier's hole is wider than
  the tube at every depth and he calls the tool contactless, so the rock under
  the feet was removed thermally. The support is the contact rule with the rock
  heated and no clip; the per-pad quantile is read as the operator-jolt
  allowance and must be shown not to control the rate. The model never lifts
  the burner (Meier's holds and jolts shape the hole locally).
- **D2e planning decisions (2026-09-18).** 09-18b's P-b as written
  (`jet_T_ent = 1550` under `fixed`) would also move the march floor and the
  free-surface ambient, so D2e adds `rec_fixed` (T_mix only). P-c
  (`foot_rule = mean`) must run with clearance off (parse abort), so it is read
  against P-a. 8 ms at 1 mm is extrapolated (K/step ∝ dt/dz); B0 checks it
  against D2d's 4 ms R7. The Perf edits are unreviewed until D2e Goal 0.
  Regression sweep 09-17: `heterogeneous_kappa`, `amr_microstructure_regrid`
  fail under `mpirun -np 4` (k_eff diagnostic not reduced per rank; one-line
  `ReduceLongSum` fix when that file is next open). `pkill -f <run name>` kills
  shell waiters too. Advisor's Hutchinson–Suo buckling release (09-18) is the
  one import worth considering later, not now.

- **D2d planning decisions (2026-09-17, planner, before any run).** The
  steady-window mesh criterion is replaced by a matched transient window on a
  0.40 m full-domain pair (a steady 1 mm run needs ≈ 24 h). Matrix domain
  0.70 m (hand needs 0.65 m, runs 13 % fast). Min Ø scored to the feet depth
  (D2c's nozzle-plane range missed the 50 mm feet zone). T sweep lower end
  (1.3 m/h) is below 1600 K and reported open-ended. Run 3 (n = 0.5) is hand
  only. "D2b closures + n = 1" keeps exponential/count/clearance.

- **`surface_mf` is φ-windowed** (`UpdateSurfaceMaskFromPhi`, φ ∈ [0, dz),
  refreshed each step at Removal.H ~1295). At a φ_top = dz tie a column has NO
  surface cell, so radiation/convection losses drop out for that step.
  **Measured by S1 as negligible**: on Meier 2D, 1 of 280 columns for ~1.1 s,
  7e-5 of column-steps, ~0.03 J of 77 kJ. Use `surface.follow_mask = 1` whenever
  a source must track the receding surface.
- **The S1 criterion finding (2026-09-15).**
  - `sp_weibull` + `spall.sample = top_cell` + `local_thermoelastic` fires at
    a fixed **cell-average** T: 656 K at a = 20 µm, p = 1 MPa; ≈ 812 K for the
    Meier Weibull draw (to be confirmed by A1).
  - **h_col ≈ dz/2 − a_f is a removal increment, not a flake.** Two firings
    empty one cell. Flake/PSD claims are retired for self-heated drilling.
    Rossi (prescribed_T) keeps its isotherm-depth reading.
  - At fixed flux, ROP = q/(ρ·Cp·ΔT_fire) at any dz.
  - Robin at 2 mm is not converged (cell +42 %, face −82 %). The converged
    1-D answer is the pinned surface q = h(T_gas − T_fire) − εσ(T_fire⁴ −
    T_a⁴).
  - **Quote Robin as a bracket.** The finest-dz data bracket is 14.0–18.6 m/h
    at (1e4, 1000 K); "15–16 m/h" is an extrapolation.
  - Do not reconstruct a sub-cell T from the 2 mm cell average for the Robin
    flux (it gives ~3000 K).
- **Overshoot time-step rule for removal runs** (binding over the conduction
  limit): ΔT_step = q·dt/(ρ·Cp·dz) ≤ 5 K. That means dt ≤ 11 ms at 2 mm and
  2 MW/m²; ≤ 1.4 ms at 2 mm and 15 MW/m²; ≤ 0.7 ms at 0.125 mm and 2 MW/m².
- **`spall.h_col_events_csv` is edge-triggered.** It logs only the first firing
  per new top cell, so never sum it for recession (≈ 2× undercount). Rossi's
  test depends on these semantics. Use `spall.removal_events_csv` (A1) for
  recession and T_fire: one row per (step, column) removal. Key on `regime`,
  because `h_scan` can be nonzero on regime 2/3 rows.
- **`weibull.enabled = 0` + `spall.enabled` segfault — FIXED by A1.** The
  scalar path's per_face gate is `Sp_field(k_top) >= 1`. `connected_cluster`
  still needs Weibull. The S1 degenerate-draw workaround is now optional.
- **Weibull `V0` (user decision 2026-09-15): use `weibull.V0 = V_cell` in all
  C/D inputs.**
  - With the default 1e-9, T_fire drifts ~+10 K per halving and never
    converges (S1b).
  - The default is left unchanged, so Kant/Rossi/Meier existing inputs are
    unchanged.
  - The IG vol_factor applies only to IG cells (GB = per-facet).
  - The Meier 2D 812 K / 6.3 m/h numbers are V0 = 1e-9; with V_cell, 2 mm
    fires at ≈ 820 K.
- **Pinned Robin closure (C1) caveats.**
  - The closed-form match is mostly a consistency check; the independent
    evidence is that the result lies inside S1's fine-mesh bracket.
  - A sharp uniform patch edge leaves a never-firing rim column pair, so score
    non-rim columns or use a smooth profile.
  - Time to first firing is face-form and mesh-dependent: score by windows and
    `removal_events.csv` `h_applied` fits, never by end depth or plotfile
    staircases.
  - Pin state is not checkpointed.
  - `pinned_idle_cycles` (default 2) is not a knob: 1-D refire gap
    1.004·t_cell, 2-D up to 1.32·t_cell.
- **Burner drawings (planner, 2026-09-16) — the D2b blocker is RESOLVED.**
  Source: `thierrymeier.ch/documents/drillBabydrill_VT5_R4.pdf`
  ("Zusammenstellung Drill baby drill VT5", CATIA V5, 19.11.2015, 17 A2
  sheets, 797 mm long, ca. 14 kg). **Not in the repo** — re-download to
  re-check.
  - Sheet 15/17 **"Abstand Halter F VT5"** = the feet: **OD Ø 80 mm**, bore
    Ø 45 (+0.2/+0.05), plus Ø 56 and a Ø 68 bolt circle with Ø 6.2 holes;
    overall length **58 mm** with an **8 mm** collar, so it protrudes
    **50 mm**; the isometric shows **3 slots**, so the bearing surface is an
    interrupted annulus (a tripod), not a continuous ring.
  - Sheet 14/17 **"Flame jet nozzle D VT5"**: flange **OD Ø 80 mm**, Ø 68
    bolt circle, Ø 6.2 holes at 120° (it bolts to the holder); central bore
    **Ø 7.5 mm** at the outlet; converging contour R12/R20/R10 from a Ø 46
    bore, Ø 56 f7 spigot.
  - Sheet 1/17 assembly, section A-A: confirms the Ø 7.5 bore and a **50 mm**
    nozzle-exit-to-foot-tip dimension. The mantle above the nozzle is a
    67 × 2.5 tube (Ø 67), narrower than the Ø 80 drilling end.
  - **Derived:** burner Ø at the drilling end **80 mm** (r = 0.040 m); foot
    annulus **[0.028, 0.040] m** — section N-N, re-read at higher resolution,
    shows the feet are a **tube OD Ø 80 / ID Ø 56** running 50 mm below an
    8 mm flange (bore Ø 45), with 3 slots at the bolt positions (the earlier
    [0.034, 0.040] "nominal" misread the Ø 68 bolt circle). The jet expands
    inside this collar for its full 50 mm and exhausts through the slots;
    `foot_standoff` **0.050 m**; `D = 7.5e-3 m`, so `SOD/D = 6.67` (not 7).
    Meier's Ø 85–93 mm hole is 5–13 mm larger than the burner, consistent
    with thesis p. 223.
  - **Two discrepancies, do not silently resolve.** (1) The drawing's
    **D = 7.5 mm** vs the reference `.md`'s **7.1 mm** (thesis text): the
    drawing wins; consequences are Re ≈ 5.66e4, h_Martin ~8 % lower
    (h ∝ D^−1.5), potential core 5D = 37.5 mm. D2a's `jet.py` used 7.1 mm.
    (2) The drawing's flows (50 + 5 m³ₙ/h ≈ 62–64 + 3.4–3.6 kg/h) are
    **design** values, above the thesis's **measured** 52 + 2.47 kg/h. Use
    the measured `ṁ = 0.01513 kg/s`.
- **D2a2 lessons that constrain D2b and later jet work (2026-09-16).**
  - **Keep the stagnation decay.** Without it the centre runs away from any
    feed. Baseline treatment for D2b: `jet_stagnation = decay` +
    `jet_entrained_mass = 1` (energy-conserving; nozzle mass discards 12–17 kW
    of the 30.4 kW nozzle budget as `jet_P_decay`).
  - **The centre leads the nozzle** (s_c 60–90 mm, T_stag 1000–1270 K at
    1900 K). Hand predictions must use the centre's equilibrium stand-off
    (`walljet.py` `s_equilibrium`), never s_c = SOD.
  - **Planner hand-table errors, recorded.** D2a2's first table used a
    non-computable criterion (T_gas decays exponentially toward T_s and never
    reaches T_fire) and a linear in-core depletion; the revised table's decay
    rows were computed at **D = 7.1 mm** (`jet.py`), not 7.5. Hand numbers go
    through `walljet.py` (D = 7.5 mm) and are written by the implementer
    before runs.
  - **Cold entrainment lowers flux at every radius** (reviewer, integrated):
    neglecting it is an upper bound on flux and on removal reach.
  - `jet_mdot` is the modelled sector's share (quarter domain = mdot/4).
  - **Wall treatment is open:** `input_drilling`'s
    `spall.flake_coherence_length = 0.004` caps walls at 63–65° and capped
    cells overheat under the pinned flux; `spall.surface_normal = 1` recedes by
    h/cos θ with no matching flux on horizontal area. D2b decides by a
    pre-stated energy-consistency rule.
  - **Most of the march is face-form uptake by unfired rock** (pinned share
    0.01–0.12), so T_gas(r) and the ring flux carry the 2 mm face-form mesh
    sensitivity.
  - Under a prescribed feed the flank freezes once the nozzle plane passes
    it; read prescribed-feed runs before that time. The feet rule removes it.
  - `jet_dr = dx` is first-order in the march (16 % excess-T error at 2 mm in
    a strip) but moves the quarter-disk ROP only 0.4–1.4 % per halving.
  - `patch_P_robin − jet_P_face` is not a clean overdraw measure; a direct
    counter is deferred.
- **Process (planner decision 2026-09-16, from D2a2 takeaway 12).** Per step:
  targeted bit-identity witnesses for the touched paths, with the full
  26k-file sweep once at commit; no 2-D slab parts for jet closures (the slab
  needs an invented jet share); a few 3-D runs rather than full grids; sparse
  `amr.plot_int` (D2a2 wrote 11 GB). Implementers re-read `ACTIVE_STEP.md`
  before long runs (D2a2 was revised mid-phase).
- **`weibull.V0 = V_cell` switches the size effect off, not just its mesh
  dependence** (reviewer, 2026-09-16). `vol_factor ≡ 1`, so `a0_ig` becomes
  "the characteristic flaw of a cell". The right choice, but **D3 must not
  describe the model as size-scaling.**
- **Jet h magnitude conflict (planner, 2026-09-15; D2a2 narrowed it, D2b settles).** Under the D2a2 closure every anchor drilled the r = 40 mm ring at only 0.5–1.3 m/h, so the Martin-vs-bracket difference is now a factor ~2 in ring ROP, not in face power. Martin
  (1977) with Meier's ṁ = 0.0151 kg/s, D = 7.1 mm, Re ≈ 5×10⁴ gives
  h ≈ 1.5 kW/m²K at 2.5 D, SOD 7 D. The reference bracket is 5–30 kW/m²K, and
  Ch. 7's 10 kW/m²K comes from Potter Drilling, not measured. Correlation
  constants were quoted from memory; D2a must verify them against the source.
- **AMReX rank-count identity:** Cell_D is written one file per rank. Compare
  the rank-ordered concatenation or the field arrays, not a file-by-file `cmp`.
- **Energy ledger (`energy_ledger.enabled`)** is exact only with Neumann-0
  outer BCs, and its state is not checkpointed (it re-baselines on restart).
- **Kant p = 0 onset now prints 446.1 °C** (Eq. 15 445.6 °C) in
  `sp_kant_pressure_sweep` (reported by E1), both pre and post E1 on the
  current tree, vs the 461.3 °C quoted in older notes. The cause has not been
  identified (it predates E1). Treat the live test output as authoritative
  until reconciled.
- **P4 cuttings/debris shielding — WITHDRAWN 2026-09-15 (design kept for
  reference; premise invalid, Σ₀ uncalibratable on Meier).** Default-off `shield.*`.
  Per-column debris areal mass Σ(i,j) [kg/m²]: `Σ ← Σ·exp(−dt/τ) + ρ·h_spall`
  each step (SPALL only, not vapor — vaporized rock leaves as gas); the beam is
  attenuated `Pi → Pi·exp(−Σ/Σ₀)` in `AdvanceMicrostructure` only (the sole
  surface-following/removal path). Removal runs BEFORE the beam in `Advance()`
  (747 vs 1229), so the beam reads fresh debris (stable — `h_spall` derives from
  the prior-step temperature, no implicit loop). Σ stored broadcast down each
  column so the per-cell beam read is z-decomposition-safe; `h_s_cand_col` is
  already MPI-gathered. Single-level only (max_level=0 dev harness); AMR debris
  repair deferred. Steady η is set by τ/Σ₀ (degenerate at steady state; τ also
  sets the smoothing timescale — pick τ physical ~1 s, calibrate Σ₀).
- **Refactor track (R-series) is decoupling `src/Integrator/MMWSpalling.H`
  (~4377 lines) into per-concern pieces WITHOUT changing its virtual
  multiple-inheritance graph** (CLAUDE.md Lessons #1 — base-class edits caused
  AMR vtable crashes before). Planned decomposition by field-ownership cluster:
  Thermal/enthalpy, Microstructure, ContinuousDamage (`damage_law`), Cohesive,
  LefmSp (`sp_weibull`), Removal. **R1 (active)** extracts the ~957-line
  `UpdateRemovalAfterCohesive` into `src/Integrator/MMWSpalling/Removal.H`,
  `#include`d *inside the class body* — the method stays an inline member, pure
  code-motion, byte-identical regressions. The partial-header trick is safe
  because the build (Makefile) compiles objects only from `*.cpp` (mindepth ≥2)
  and top-level `*.cc`; `*.H` files are include-only and never standalone-
  compiled, and incremental rebuilds still trigger via the generated `.d` deps.
  Such a fragment must contain NO `#include`/`namespace`/`class` (it is pasted
  into `class MMWSpalling { ... }`). Data-member / helper relocation is
  deferred to R2+. New pieces go under `src/Integrator/MMWSpalling/` (mirrors the
  existing `src/Integrator/Base/` subdir convention).
- The long plan still contains an old Zhang/Oglesby note with `omega_0 = 0.20 m`.
  The corrected validation reference says `omega_0 = 0.02 m`.
- Step numbering contains inserted/moved validation steps such as 4b, 6d, 7b,
  7c, 7d, 8b, and 11b. Keep those labels because tests and notes refer to them.
- Homogeneous Hu reproductions should stay homogeneous. Voronoi heterogeneity is
  a model-extension study for Hu, not the direct reproduction.
- Step 8 changed AMR microstructure repair: passive `D` is still re-evaluated
  from `damage.ic.*` when `damage.enabled = 0`, but evolved `D` is no longer
  reset from ICs when `damage.enabled = 1`; repair refreshes derived fields
  such as `kappa_eff`.
- Step 8b is complete and must remain an indicator validation with
  `damage.enabled = 0`: it validates `f_b = sigma_v / sigma_s` from
  homogeneous thermoelastic stress, not the DP damage law.
- Step 8b fixed a mechanics coupling issue: `surface_patch` prescribed
  temperature must feed top-patch nodal eigenstrain, not only the cell-centred
  heat flux. Preserve that coupling for future surface-T / moving-surface work.
- Step 8b passes Hu surface `f_b` targets but relaxes L6 deep `f_b` to `<0.25`
  because the bottom clamp inflates deeper stresses; strict `<0.1` remains
  blocked on a stable pure-roller elastic solve.
- Quartz alpha-beta transformation strain is implemented only on the
  microstructure mechanics path, default-off per phase via
  `microstructure.phaseN.alpha_beta.*`, and verified by
  `tests/MMWSpalling/alpha_beta_transition/`.
- Step 9 CZM is default-off under `cohesive.enabled`. Its fields
  `gb_delta`, `gb_delta_max`, `gb_traction`, and `gb_cohesive_damage` are
  registered only when enabled and key on `is_grain_boundary`, not `is_gb`.
  Mechanics traction feedback is not active yet.
- Step 11 removal is default-off under `spall.enabled || vapor.enabled`.
  Shared fields `phi`, `removed`, `regime`, and `RoP` register whenever either
  mode is enabled; per-mode diagnostics register only for their active mode.
  Per-column removal is MPI-safe for arbitrary BoxArray partitioning, including
  z-decomposed columns. AMR moving-surface validation is still deferred.
- Step 11 void handling exists only on the microstructure heat path; the
  constant-alpha, global-material explicit, and global-material implicit paths
  do not yet honor `removed`.
- For Hu onset/removal work that needs damage or spall fields, use homogeneous
  one-phase `microstructure.mode = expression` inputs rather than Voronoi. The
  direct Hu reproduction stays homogeneous; heterogeneity is a later extension.
- Step 11b is a calibrated v1 onset validation, not a fully predictive damage
  closure. It uses single-level, spall removal off, `damage.stiffness_floor
  = 1.0`, and material-specific `A_D` values to match Hu onset/LRST. Sandstone
  still develops a deep DP damage zone; the test reports this as a warning.
  Re-enable AMR and strict shallow-depth failure only after the mechanics/regrid
  and DP localization issues are resolved. As of Step 12, `hu_spall_onset` runs
  at `amr.n_cell = 32 32 32` (was `40 40 40`) — the 32³ partition is a clean
  16+16 box split under `max_grid_size = 16` and dodges the 40³ intermittent
  SIGBUS that BUS_ERROR_ANALYSIS.md documents. The sandstone tolerance window
  was retargeted from the Hu-average onset (76.7 s) to the Hu-primary onset
  (89.0 s), keeping ±20%. The live 32³ result is granite 37.6 s / sandstone
  93.0 s, both inside the tolerance windows; ARCHIVE_DONE.md keeps the
  historical 40³ numbers (37.0 s / 87.0 s) with an explicit post-Step-12
  update note.
- AMR multi-level mechanics for `hu_spall_onset` is a known limitation, not a
  bug to chase. Five hypotheses (Jacobi smoother bandwidth, bottom-solver type,
  CF interface placement via `n_error_buf`, RBM null-space deflation, F0
  cross-level inconsistency) were refuted with quantitative numbers. See
  `tests/MMWSpalling/hu_spall_onset/AMR_FAILURE_ANALYSIS.md` §9–§10 and
  `AMR_SMOOTHER_LIMITATION.md` before re-opening; surviving candidates are
  per-level coefficient consistency and vector-elastic multigrid machinery vs
  `MLNodeLaplacian`. The in-tree diagnostic toolkit (`DiagnoseAMRResidual`,
  `ApplyAMRF0Consistency`, `CheckAMRMechanicsFieldsFinite` and the
  `Check*FieldFinite` family, plus the `mmw.debug_amr_mechanics` /
  `ALAMO_AMR_DIAG` / `ALAMO_LEV1_F0_FROM_COARSE` gates) was removed during
  post-Step-18 cleanup — see [git history of `src/Integrator/MMWSpalling.H`]
  for the prior implementation if reopening this work.
- `revised-model-approach.md` (Stages 1–5 of the LEFM Sp + Weibull
  re-architecture) is in the repository root and is the authoritative spec
  for the `sp_weibull` branch. The long plan §12–§18 is the implementation
  staging; cross-check both when planning sp_weibull work.
- Step 12 ships the LEFM Sp + K_Ic infrastructure (toggle, Tada-weight
  K_I integral via trig substitution, Nasseri K_Ic table, optional
  prescribed σ_xx / T_surface parser overrides for verification) under
  `spallation.model = sp_weibull`. Default `damage_law` preserves
  Step 8–11b behavior byte-for-byte. EOS work is dropped indefinitely;
  no `eos.type`, no `src/Numeric/EOS/*`, no `tests/MMWSpalling/eos/`.
  DP "logged-only under sp_weibull" (revised-model-approach §10 #4) is
  NOT implemented yet — `damage.enabled = 1` is a hard parser error
  alongside `sp_weibull`.
- **Kant onset closed form is Kant Eq. 15, NOT Eq. 14 `Sp_red`.** Eq. 14
  `Sp_red = Sp_rock·Sp_fluid + Sp_conf` is an `h_fl`-dependent
  rock-classification number, not an onset temperature. The onset
  equation is Kant Eq. 15
  `Θ_onset = K_Ic·(1−ν)/(2·1.763·√(a/π)·E·α)` at `p=0` —
  `h_fl`-independent and `λ`-independent, a material threshold. `h_fl`
  sets only the onset *time*. Step 18 (confining-pressure sweep) must
  use Eq. 15 with the `−p·ν/(1−ν)` confining term, NOT Eq. 14.
  `kant-2017-central-aare-validation-reference.md` §5 was corrected
  during Step 15b (Kant onset) accordingly.
- Kant's Table-2 390–560 °C band was computed with **frozen room-T
  `K_Ic⁰ = 1.5`** (Kant Table 1 / §4.3), not `K_Ic(T)`. The
  `sp_kant_onset/input_verification` runs use frozen `K_Ic⁰` and are
  the ones scored against the Kant band; `input_model` uses the
  Nasseri-`×1.05` `K_Ic(T)` table and is reported as the production
  model's honest prediction, not banded. Honest agreement at Kant
  onset is ~3.2–3.5% Eq.15-vs-FEM; the FEM `σ_xx` residual is
  ~0.96–0.97× the thin-film ideal (O(dz/thermal-layer) discretization,
  approaches 1 as the thermal layer thickens with slower `h_fl`) —
  hence the 8% target-3 tolerance.
- The **lateral roller box** (`el.bc.type = constant`, all 26 regions:
  `u_x = 0` x-faces, `u_y = 0` y-faces, `u_z = 0` zlo, +z free; all 12
  edges + 8 corners specified consistently) + **full-face**
  `convective_flame` heating is the canonical Kant-style 1-D-confinement
  setup — `ε_xx = ε_yy = 0` exactly, `σ_xx = -EαΔT/(1-ν)` per cell,
  MLMG well-conditioned, coarse cubic 1 mm mesh sufficient. Reuse this
  layout for Step 18 (confining-pressure sweep). The patch + cold-rim
  layout was abandoned during Step 15b Kant validation because the
  sub-cell thermal layer under-resolves `σ_xx` (FEM only ~0.45× ideal)
  and fixing it requires anisotropic fine-`dz` cells that the elastic
  MLMG solver does not converge on (see next note).
- **The elastic MLMG solver fails on anisotropic cells** — including a
  modest 2:1 aspect ratio. Keep mechanics-active meshes cubic. Sub-key
  trap with `el.bc.type = constant`: BC-type goes in the prefix, so the
  sub-keys are `el.bc.constant.type.<region>` (one per all 26 regions).
- **Weibull first-fire is governed by order statistics, not the scale
  parameter.** Under `weibull.enabled = 1` with `N` surface cells, the
  cell that first satisfies `Sp ≥ 1` holds the realized maximum flaw
  length `a_max ≈ a0·(ln N)^(1/m)` (Galambos extreme-value asymptotic),
  NOT `a0` itself. Sp scales as `√a`, so first-fire onset ΔT is
  systematically *below* the deterministic Kant Eq. 15 prediction at
  `a = a0` (~7% below for `N = 576, m = 15`; ~80–110 °C below the Kant
  `[390, 560] °C` band that was derived from a single deterministic
  crack). This is why the Step 15c v1 sweep PASS gate was rejected —
  it scored an order-statistic quantity against a deterministic band.
  Resolved in §15c v2: score Weibull-on FEM first-fire against
  `Eq.15(a_max_realized)` computed from `max(flaw_a_mf)` on the t=0
  plotfile, with the deterministic Kant band kept as a heatmap
  colorbar reference for the Weibull-off verification run
  (`sp_kant_onset/input_verification`, 461.3 °C ∈ [390, 560],
  unchanged). **Always score Weibull-on first-fire against
  `Eq.15(a_max_realized)`, never `Eq.15(a0)` or the Kant band.**
- The Sp criterion stays sign-agnostic in `src/Numeric/SpCriterion.H`:
  `K_I` consumes whatever the `sigma_at_depth` callable returns; the
  `-sigma_xx` convention lives in the `UpdateSpAfterMechanics` caller, not
  the header. Do not move it into the header. Step 16's column h_spall
  scan likewise builds its `sigma_at_depth` callable in the integrator
  with the compression-positive convention; the header stays generic.
- **`h_spall_field_mf` and Rossi's measured observable are different
  physical quantities.** Step 16's diagnostic is a per-column running
  max of the deepest K_I-vs-K_Ic crossing across all steps — under a
  monotonic ramp this captures the LATE-ramp value, not the
  integrated event population. Rossi 2018 measures the spatial
  distribution of distinct cracks across all events during heating.
  Step 16's Rossi test reports the bands as informational only; the
  quantitative validation is Step 16b's scope using a new
  time-integrated 3-D `crack_event_count_mf` field. Edge-triggered
  counting (first crossing per cell only) and stress-relief
  option (ii) (freeze re-firing, no σ touch) are user-approved and
  documented in `rossi-validation-diagnostic-design.md` (repo root).
- **Step 16 geometry was forced off-packet by MLMG conditioning + free-
  lateral physics.** Working config for sp_weibull + spall.enabled +
  prescribed-T heating: `32^3 cubic cells in 50×50×50 mm + Kant-style
  lateral roller box + full-face heating` (heat patch radius >=
  face half-diagonal). Avoid: thin slabs (e.g. 30 mm deep at 1 mm
  cells — intermittent MLMG failures), non-power-of-2 meshes (can't
  coarsen multigrid past 1-2 levels — MLMG stalls at residual ~1e-5),
  free-lateral BCs with finite heated patch (central-column σ_xx
  drops to ~0.3× of 1-D ideal — Sp never fires within the ramp). The
  lateral-roller-box config is a "central column of an infinite
  heated zone" reading of Rossi, not a literal patch-plus-rim
  reproduction.
- **Sub-cell h_col events are scan-quantization-suppressed.** Step 16
  walks the depth scan in integer-cell strides. A column where K_I ≥
  K_Ic only at `s = 0` (h_col < dz) registers as h_col = 0 and is
  NOT counted as firing. At 1.5625 mm cells, this biases the GB/IG
  firing ratio to 100% GB (IG cells with `a0_ig = 4 µm` have Sp ≥ 1
  at the surface but fail the scan at `s = 1`). Step 16b added an
  8-iter sub-cell binary search bracketing `s_last_pass` and
  `s_last_pass + 1`; refines `h_col` to dz/256 (~6 µm). This partly
  lifts the quantization — Step 16b's GB fraction dropped from 100 %
  (Step 16) to 96 % (21 IG firings of 528). The remaining 96:4 ≈
  24:1 GB/IG ratio is **not** Rossi's 5:1 and **will not be fixed by
  Step 16c's depth-binning**; the path to 5:1 is higher peak σ_xx
  (longer/hotter ramp), a re-examined `a0_ig`, or eventually finer
  mesh (option (a)). Step 16c calls this out as a scoped-out
  limitation.
- **σ_xx(z) is local-per-z under 1-D confinement** (Kant-style
  lateral roller box + full-face heating, as used in Step 16's
  working geometry). `σ_xx = -E·α·(T(z) - T_ref)/(1-ν)` depends only
  on the local cell's temperature — zeroing σ at fired cells (design
  note §3 option (i)) does NOT change σ at deeper cells on this
  geometry, so option (i) is ineffective as an R1 fix here. Option
  (i) might help on patch+free-lateral geometries (stress
  redistribution through the cold rim) but that reopens Step 16's
  MLMG conditioning failures.
- **The `crack_event_count_mf` field is cell-centred and is
  insufficient for the Rossi depth profile at dz = 1.5625 mm.** Step
  16b's `cec` flips 0→1 on the first K_I ≥ K_Ic crossing per cell
  (edge-triggered, option (ii) freeze, no σ touch). But at this
  mesh, K_I crashes from above-K_Ic at z = 0 to below-K_Ic at z = dz
  in one cell, so only the top cell ever flags. The Rossi-comparable
  observable is Step 16c's per-step **h_col CSV log + post-process
  histogram** instead of a cell-centred event field. Sub-cell
  crack-tip depths from the bisection (`h_col ∈ (0, dz)`) are
  recorded in `h_spall_field_mf` as a running max per column, and
  Step 16c logs them per step for the histogram. Sub-cell σ_xx
  reconstruction (linear cell-to-cell vs T-derived vs node-based
  interp) is the key modelling knob Step 16c locks down explicitly.
- **Step 16c's option (L) σ_xx reconstruction has a sampling
  artefact at `z = dz/2 − a_f`.** Under linear cell-to-cell σ_xx,
  the top cell has no adjacent cell ABOVE the surface to
  interpolate from, so σ_xx is flat in `z ∈ [0, dz/2]`. K_I
  integration over `(z_crack, z_crack + a_f)` is then constant
  for `z_crack ∈ [0, dz/2 − a_f]` and starts dropping only past
  that. The bisection brackets `(0, dz)` and converges to the
  START of the K_I drop — h_col ≈ dz/2 − a_f for every firing
  column, regardless of a_f scale or firing time. At dz = 1.5625
  mm and a_f ∈ [4, 22] µm: h_col ∈ [759, 779] µm — verified
  arithmetically against 16c's observed [757, 787] µm range
  (matches to within bisection tolerance dz/256 ≈ 6 µm).
  Physically, option (L) discards the 20 % surface-T jump: at
  first-firing time t ≈ 80 s, `T_face = 693 K` (from
  prescribed-T `T_expr`) but `T_top_cell_centre ≈ 625 K`
  (thermal-diffusion gap across dz/2). **Step 16d picks option
  (a) — T-derived σ_xx with the prescribed-T face BC** — to
  address the root cause.
- **`rossi-validation-diagnostic-design.md` is the canonical
  reference for the Rossi validation chain.** §1-§5 cover the
  16b/16c design decisions (edge-triggered counting, stress-
  relief options i/ii/iii, sub-cell bisection, open questions).
  §6 is the Step 16 mechanics status. §7-§10 (added during 16c)
  cover the (L) artefact mechanism, Step 16d scope (option (a)
  only, with the K_I-not-K_Ic-driven argument against (c)), the
  mesh-resolution limit reasoning (two reconstruction points
  between face and top-cell centre means NO sub-cell strategy
  can manufacture features at <100 µm depths — only anisotropic
  z-refinement can), and anti-patterns (no sub-cell tweaks
  beyond (a), no K_Ic-only smoothing, no stop_time tuning past
  Rossi conditions). **§11 (added post-16d) marks the chain as
  CLOSED with R1 DEFERRED**: empirical convergence of 16b/16c/
  16d on three different reconstructions, all peaks outside
  [100, 200] µm; revisit conditions enumerated (anisotropic
  z-refinement, different reference experiment, or accepting
  the shape-without-band-location reading). The Rossi mechanics
  test stays as an ongoing P1-P4 regression; R1/R2/R3 reported
  but not project-gated. Step 17 (per-cell v_n + Zhang/Oglesby)
  and Step 18 (Kant confining-pressure) advance without the
  Rossi R1 gate.
- **Rossi R1 deferral checklist (do not re-litigate without
  reopening one of the §11 revisit paths)**: do NOT chase the
  Rossi peak on this mesh with further σ_xx reconstruction
  tweaks; do NOT smooth K_Ic-vs-z to chase the peak (K_I is
  σ_xx-driven, K_Ic continuity is second-order — verified by
  the 16c/16d empirical chain); do NOT extend stop_time past
  130 s; do NOT modify Weibull parameters (`a0_gb`, `a0_ig`, m)
  to shift R3 toward 5:1 — those are calibrated against Kant.
  The 24:1 R3 ratio reflects late-ramp σ_xx skew on this
  geometry, not a model bug.
- **Step 16d's σ_xx reconstruction is gated to `surface_patch.mode =
  prescribed_T`** because the T(z) lerp uses `surface_patch.T_f(time)`
  to pin T at z=0. Step 19 (active) forks this lambda: under
  prescribed_T it keeps the byte-identical T-derived form; under
  any other heat source (MMW beam, convective_flame, surface_patch
  disabled) it reads σ_xx from `stress_mf` via NodeToCellAverage —
  same pattern as the diagnostic Sp_field path in
  UpdateSpAfterMechanics. The mechanics solve must be active
  (`el.time_evolving = 1`) for the non-prescribed_T branch to work;
  Step 19 adds a fail-fast abort if it isn't. The prescribed_T
  branch byte-identity is the canonical witness for Step 19
  correctness (Step 18 `sp_kant_pressure_sweep` p=0 = 461.3 °C
  must stay unchanged).
- **Step 17 per-cell v_n + regime_field + local-Q infrastructure.**
  Under `spallation.model = sp_weibull && (spall.enabled ||
  vapor.enabled)`, the vapor-RoP path is a three-branch ladder:
  (1) `vapor_prescribed` (synthetic-test priority, unchanged),
  (2) `sp_weibull` local-Gaussian `Q(x,y,t) = (2P/(π·w0²))·
  exp(-2·r²/w0²)`, `rop_v = oneR·Q/(rho·L_tot)` (Step 17 NEW),
  (3) damage_law bulk-averaged `oneR·P/(rho·L_tot·A_beam)` (Step
  11 byte-identical). Column-winner precedence under sp_weibull:
  spall-priority (if Sp fires, regime=1; else if H≥H_vap,
  regime=2); damage_law's winner-on-h preserved. New
  `regime_field_mf` cell-centred plotfile diagnostic carries
  the canonical per-cell regime under sp_weibull; the scalar
  `regime_mf` stays for damage_law back-compat (same per-cell
  value, different field name). Beam class NOT refactored —
  Gaussian formula inlined at the vapor-RoP call site. No new
  parser keys. Local-Q at r=0 exactly equals the bulk-RoP value
  by construction (`A_beam = 0.5·π·w0²`); the local-Q signature
  is the radial decay (Pearson r = +0.988 verified in
  sp_v_n_regime). `hu_end_to_end` is the canonical damage_law
  simultaneous-spall+vap test — must stay byte-identical to
  Step 16d baseline (granite 40.0 s, sandstone 90.5 s) for any
  future change in this area.
- §15c shipped `DetachMode { PerFace, ConnectedCluster }` +
  `weibull.detach_mode` / `weibull.A_crit` parser, `Sp_cluster_id_mf`
  (top-z plane only, 0 = not firing or below-A_crit cluster, positive
  int = labelled), `UpdateSpClusters` 4-connected union-find with
  ReduceIntMax across MPI ranks, and per-step
  `<plot_file>_clusters.csv` (time, cluster_count, max_cluster_area,
  total_firing_area). `Sp_field_mf` is never modified by the labeller;
  cluster gating moves *which* cells detach but not *when* first-fire
  happens. Cross-level AMR connectivity is deferred (current `sp_*`
  runs are `max_level = 0`).
- Step 14b patched `PrincipalStressRatio` (Step 10 production spall
  path) to return a `1e30` sentinel on mixed/tensile principal-stress
  states instead of 0 — physically, a free surface with tensile sigma_zz
  is more spall-prone, so the Hoek-Brown anisotropy gate should pass.
  Production Hu / Kant / Rossi validations no longer need the
  `spall.prescribed_ratio` workaround; that hook stays only for
  synthetic unit-style spall tests (`spall_event`,
  `regime_low_high_power`).
- The Hu validation chain (6d → 7b/7d → 8b → 11b → 14) is closed on
  the `damage_law` branch. Future Hu-style validations should use
  Step 14's `tests/MMWSpalling/hu_end_to_end/` as the canonical setup
  rather than spinning up new test directories.
- ARCHIVE_DONE.md and the Step 11b prose still mention
  `A_D = 0.000667` for sandstone, but the live
  `hu_spall_onset/input_sandstone2` and `hu_end_to_end/input_sandstone2`
  ship `0.005`. Treat the live test inputs as authoritative; the
  archive prose is historical.
- `src/BC/Operator/Elastic/ZloRoller321.H` is untracked in git but is
  load-bearing for `hu_spall_onset` (`el.bc.type = zlo_roller_321`)
  shipped by Step 11b. Treat as "intentionally untracked, awaiting a
  proper checkin"; do not delete from the working tree.
- The kept working-tree changes to `src/Solver/Nonlocal/Newton.H` and
  `src/Integrator/Base/Mechanics.H` (after Step 12's hunk-by-hunk
  cleanup) are: the `ZloRoller321` BC registration; the simplified
  `solver.solve(disp, rhs, model)` API with default-tol fallback in
  Newton; the post-solve `Util::RealFillBoundary(*disp_mf[lev],
  geom[lev])`. These are net improvements and stay. The reverted
  pieces were the `bx → domain` stencil swap, the `domain` arg to
  `GetBC()`, the rhs_mf defensive zeroing, and the AMR-investigation
  escape hatches in `Operator/`.
- Run test and validation simulations with 4 MPI ranks by default on this
  machine, e.g. `mpirun --oversubscribe --bind-to none -np 4 ...`, unless the
  test metadata or user request explicitly requires a different rank count.
- Step 6b fixed the old microstructure path that evolved `Temp` directly and
  skipped `InvertHtoT`; the microstructure path now advances conserved `H`.
- Step 6c split mineral topology from grain topology: `phase` is a mineral ID,
  `grain_id` is a grain label, `is_grain_boundary` is the future CZM/damage
  flag, and `is_gb` remains only a phase-boundary compatibility alias.
- `microstructure.grain_expr` is optional in expression mode; avoid
  comma-bearing parser expressions in inputs because ParmParse splits on
  commas.
- Old Step 4c is now Step 6d so Hu conduction validates the thermal driver
  before mechanics.
- The current Hu validation text has constants, BC trendlines, and qualitative
  field-comparison numbers, but not full digitized hole-1 time-series arrays.
  Do not invent Hu reference curves; add real digitized CSVs if curve checks
  are needed.
- Step 6d added `surface_patch.*` on the explicit global-material path;
  Step 15a extended the convective-flame branch to `AdvanceMicrostructure`
  as well. `AdvanceMaterialImplicit` still does not consume `surface_patch`
  at all (out of scope, ROADMAP Known Note).
- `surface_patch.T_expr` must be one ParmParse token with no whitespace, e.g.
  `638.22+3.74*t`.
- Step 15a wired `surface_patch.mode = prescribed_T | convective_flame`
  (default `prescribed_T`). Under `convective_flame` the Robin patch flux
  `h_fl*(T_flame-Tc)*inv_dz` is applied on the in-patch top row;
  `surface_patch.h_conv` is reused as `h_fl` and `surface_patch.T_flame` is
  a required new key. `surface_patch.T_ext` remains parsed-but-unused
  (reserved for a future non-patch convective exterior). The Step 8b
  top-node eigenstrain override is gated to `prescribed_T`; under
  `convective_flame` the top node reads solved `temp_mf`.
- Step 7b is a homogeneous Hu reproduction and now passes single-level. It ships
  with a bottom clamp instead of Hu's ideal bottom roller because the pure
  roller elastic solve lost convergence by `t = 6.2 s`.
- Do not tighten the Step 7b L3 stress ordering without revisiting the
  mechanical BC; the bottom-clamp case reverses Hu's L3 ordering. Keep stress
  sampling on node plotfiles for L3/L4.
- Step 7d is complete with shipped `20^3 + max_level=1` AMR Hu reruns. Keep the
  bottom-clamp and L3 reporting-only caveats unless the elastic
  BC/preconditioner changes.
- Exploratory `16^3 + max_level=2` Hu AMR runs showed static mechanics
  displacement reuse/prolongation can blow up after regrid; setting
  `el.zero_out_displacement = 1` stabilized Granite. Sandstone still needed
  `amr.abort_on_nan = 0` because dirty `Temp` ghost data triggered plotfile
  warnings, so do not promote that two-level setup until the ghost-fill warning
  is understood or fixed.
- For yt/numpy regressions on this machine, use
  `/Users/tzetze20/Desktop/code/.venv/bin/python`; `/usr/bin/python3` lacks
  `numpy`.

- **THE STEP-FACE MECHANISM — the cause of the D2d–D2i mesh failures
  (D2j-0b, 2026-09-21; supersedes the D2d and D2h diagnoses).** Flux and losses
  act on a column's **top face only**; an exposed vertical step face is
  adiabatic. A firing top cell beside a taller, colder column therefore loses
  `k(T_top − T_wall)/dx` into that column's *interior at the step depth*
  (≈ 240 kW/m² at 2 mm, ≈ 480 at 1 mm, against q_pin ≈ 520). It stalls, goes
  cold, and becomes the next wall — an inward cascade. **Reproduced with no
  feet, no jet, no gas closure and without the `s ≤ 0` exclusion** (C0 cascades
  too, three times slower).
  - **Why the mesh dependence.** A **one-cell** step has ΔT ∝ dz, the two cancel,
    and the loss converges — which is why the static D2j-0 disk looked bounded.
    Once the step exceeds the thermal length (**κ/v ≈ 1.6 mm — corrected
    2026-09-22 from a textbook-granite 3.2**, so MORE steps qualify) the neighbour sits at
    bulk temperature at any mesh, ΔT is pinned, and the flux density goes as
    1/dx. This is the continuum **corner singularity**: a pointwise firing
    criterion evaluated at a geometric corner is mesh-divergent by construction,
    and **refinement makes it worse**.
  - **Measured, not argued:** per-bin `v_1mm/v_2mm` is 0.996–1.007 wherever the
    local slope < 0.7 and falls to 0.95 … 0.53 in the six bands where the slope
    runs 0.7 → 8.2. The whole ≈ 10 % disk-mean offset lives where
    slope × dx ≳ 1.5 mm. **The half-cell removal increment is not a second
    mechanism** (the flat interior is where it would show, and S1 already had the
    flat 1-D case converging).
  - **Dead remedies, each on measurement, not argument:** wall heating above the
    plane (C1 ≡ C2 ≡ C3 band for band — the sink is 10–20 mm below the exposed
    surface); `spall.surface_normal = 1` (inert at 2 and 1 mm, the third
    independent measurement of its inertness); true-area flux × 1/cos θ on the
    top face (deposits the energy in the column's own top cell, never reaching
    the sink); a longer idle clock; more mesh pairs on the same configuration.
  - **The fix must act on the step**, i.e. heat the sink cell directly. Hence
    D2k. Design constraints found in review and binding on it: the side power
    must come **out of the enthalpy march's budget** (the D2a failure D2a2 exists
    to prevent, invisible in a fixed-h hot disk); a side-heated cell must obey
    `T_s = min(T_f, T_pin)` or an unremovable mid-wall cell becomes a heat
    reservoir that flatters the fix; the `s ≤ 0` exclusion must apply identically
    to side faces or the C1/C2/C3 result is silently invalidated; and `f` is
    pre-registered, never fitted (a staircase of one-cell steps has area
    `dx + Δz` against a true `√(dx² + Δz²)`, over-heating by ≤ 41 % at 45° and
    ≈ 16 % at the observed 80° cone).
  - **The arbiter already exists and is cheap:** `d2j0b_plane` C1 at 2 and 1 mm,
    3 min + 26 min. Any candidate remedy is tested there first. An acceptance of
    "the gap stands" is not enough — it must be **within 5 % in every block**, on
    the same estimator as the Meier gate, or a pass will not transfer.

- **D2k's outcome and the two live explanations (2026-09-22).** Side-face flux
  fixes the step-face cascade (arbiter −45 % → −3 %) but drives the Meier rate to
  4.07 m/h. **The decomposition is the thing to carry:** absorbed power 1639 →
  3204 W (2.0×) while the rate went 2.7×, because warming the wall also removes
  the lateral conduction sink (240–480 of ~520 kW/m²), so far more of the absorbed
  power converts to recession.
  - **The second factor is correct physics and must survive any closure.** Only
    the first — a stagnation-impingement `h` on a near-vertical face — is suspect.
  - **Explanation 2 is live and has not been tested:** over a *smooth* steep cone
    the wall jet is attached and the area really is ~6× the projection, so large
    absorbed power is not obviously wrong. What may be wrong is that `h(r, s)` is
    too peaked, so the model digs a pit with an 80° skirt instead of Meier's "half
    elongated spheroid" at Ø 85–93 mm. **Its discriminating signature is rate down
    *and* Ø up together**, which no cut to side flux can produce.
  - **The coupling that may kill explanation 1:** the wall stays warm because it
    is heated, so cascade suppression and the power excess ride the *same*
    parameter in opposite directions. If the f that restores 1.3–1.6 m/h is below
    the f at which the cascade returns, the wall closure cannot satisfy both gates
    and no closure should be written. D2l Gate 0 measures both bounds in ~1 h,
    key-only.
  - **`f` is an area correction and a diagnostic dial only.** It scales `q_side`
    linearly where an `h` scale would, which is why it maps the window — but no
    value found by scanning it may be kept as a setting. A closure value must be
    derived from a named heat-transfer source and pre-registered before the Meier
    run, then checked against Gate 0's window.
  - **Score three axes, never rate alone:** rate one-sided against 1.30 m/h, **Ø
    at matched feet height** against 85–93 mm, and **J/mm³ against Meier's 15.4**.
    Cutting side flux to recover ~1.4 m/h would give the right rate for the wrong
    reason — a narrower hole drilled more slowly, which is D2b's failure again.
  - **Hole shape is compared at matched feet height, never at matched time.**
    D2k's hole looks 30 mm wider at matched time purely from drilling 3× further.

- **Shape, not rate, is the honest score (D2l, 2026-09-22).** With the side-face
  term off — the D2i/D2d configuration — the model removes **1.9× too much rock**
  (4.58 vs 2.47 cm³/s) while the burner rate reads 1.49 m/h, in band. **The rate
  agreement was partly a coincidence of two errors.** Cut Ø is 78 mm against
  Meier's 85–93; the mouth flares to 145 mm. Freshly cut rock is ~78 mm at every
  setting — everything wider is removed *afterwards*, and scales with dwell.
  - **The excess volume and the too-cheap removal are one thing seen twice.** The
    excess 2.11 cm³/s costs ≈ 3.0 kW, ~46 % of absorbed power; remove it and the
    cost lands near 2.6 J/mm³, about what Meier's numbers imply at a comparable
    delivery. In his hole that heat soaked into the wall (the **thermally altered
    band** in the sawn section); in ours it removes rock.
  - **Scoring rules now standing:** score Ø(z) and volume, never rate alone;
    compare at matched feet height, never matched time; quote the volume (water
    filled) and the cut Ø, and keep the **mouth reported rather than scored** (it
    formed under the wellhead across two phases).
  - **Anchors named for the delivery-vs-cost question** (one equation, two
    unknowns, not separable by rate data): the altered-band width, and the
    mouth-flank behaviour with and without dilution.
  - **Reopened:** cuttings absorbing heat was withdrawn as degenerate with
    delivery efficiency *under a fixed-flux BC*; the enthalpy march made the
    budget finite, so that premise has expired.
  - **Largest unverified input:** the ~1900 K driving gas temperature. Meier's
    burner cooling water (160 L/h, thermocouples both ends, ΔT never reported) is
    186 W per kelvin — a 10–30 K rise is 1.9–5.6 kW, i.e. **100–300 K off the
    nozzle gas temperature**, and heat flow scales with (T_gas − 821 K). Cheapest
    question to ask Meier, largest leverage.

- **How the feet reach the heat distribution (code path, for D2m).**
  `FeetDescent` reads only live columns with r in [`foot_r_inner`,
  `foot_r_outer`] (28–40 mm), takes the 0.9 quantile per pad, and sets
  `z_target = z_foot + foot_standoff` with monotone descent.
  `BuildFlameColumns` then derives every column's `s = z_nozzle − z_face` and
  `h = h_expr(r, s)`, with `h = 0` where `s ≤ 0`. **The feet touch heating through
  exactly one scalar — the nozzle height.** Consequence: the 28–40 mm ring is the
  only rock held at the design stand-off for the whole run, so pad recession and
  ROP are locked together; the centre runs away (s_c 157–178 mm against a design
  50 mm) and rock outside 40 mm is left behind with a growing stand-off and
  falling heat. That predicts an erosion-rate transition at r ≈ `foot_r_outer`,
  which is what D2m tests.

- **The burner is a body in the hole, not a point source (second reading of
  Ch. 8 and the A.2 drawings, 2026-09-22).** This is the geometry the gas closure
  has never had:
  - the last **58 mm is an 80 mm skirt with a 56 mm bore and three exhaust
    slots**, resting on the rock. Inside the band the hole wall faces ~12 mm of
    steel across a **2.5–6.5 mm** gap, and **rock beyond 40 mm radius is reached
    only through the slots**;
  - above the skirt, a **67 mm water-cooled shield** in an 85–95 mm hole: a ~1 cm
    annulus with **cold steel at 300–350 K** on one side, over half a metre;
  - the 50 mm stand-off is **7 nozzle diameters, chosen as the impingement
    optimum**, through a **7.1 mm Laval nozzle**;
  - the wellhead is water-cooled and **held near room temperature** by design
    (p. 208, infrared-monitored).
  - **Consequence:** impingement transfers ~700–1500 W/m²K; flow sliding up the
    annulus transfers ~**40 W/m²K** — a factor of ~20 the model does not apply.
    And **Meier's hole is his skirt plus 2.5–7.5 mm per side while ours is
    exactly the skirt** (D2m: 79.3 vs 80), so the missing mechanism is whatever
    reaches *past* the skirt — on this reading, the jet leaving through the slots.
  - **Where the funnel forms:** in the ~50 mm band **below** the nozzle plane, not
    above it. `s ≤ 0` already zeroes flux above the plane, so any test that only
    changes heating up there is inert. **The floor→wall switch must happen at the
    corner, inside the band.**

- **Nozzle temperature — what is and is not sourced (2026-09-22).**
  **1900 K is about 8 % below equilibrium adiabatic**, not "an adiabatic
  chamber": for CH₄/air at λ = 1.2, T_ad ≈ 2040 K and dissociation costs ≈ 40 K,
  so 1900 K implies ≈ 2–3 kW of chamber loss.
  - **The water jacket brackets it from the hardware:** 160 L/h is 186 W/K, so a
    return staying well below boiling caps the loss near 11–13 kW (nozzle above
    ≈ 1450 K), and a designer's 20–40 K rise gives 4–7.5 kW (≈ 1650–1830 K).
    **So 1450–1850 K is defensible and no central value is sourced.**
  - **Do not cite the 29 % figure from §4.2.2 as a chamber-loss prior.** That is
    CW 2, the **pressure-vessel** mantle around a 260 bar bath of products already
    mixed with 350 kg/h of injected CW 1 — a different heat path, which Meier
    himself discards from his jet balance. Planner error, withdrawn.
  - **1436 K is not the drilling-phase reading.** Fig. 8.7 shows the igniter
    thermocouple at its 1163 °C peak only near 300 s with the flame attached;
    through both drilling phases it reads 0.4–0.5 of that (≈ 740–850 K), upstream
    of a lifted flame among cold reactants. It constrains the exit temperature in
    neither direction.
  - **Two routes are closed for good:** Meier no longer holds the cooling-water
    data, and chamber pressure is unreported (only the 6 bar(g) air-line limit),
    so the choked-throat route to stagnation temperature is unavailable.
  - **Keep the current input and sweep 1500–2000 K** in the queued
    nozzle-temperature study rather than adopting a point.
  - **AMENDMENT to "nozzle temperature is a rate lever, not a shape lever":**
    once wall physics exists, duct heating at ~40 W/m²K sustained for ~1000 s
    puts the wall in the **600–870 K** range depending on annulus gas
    temperature, against a **821 K** firing threshold. Near a threshold the
    nozzle temperature becomes a **shape** lever too. The deferral still stands
    for now, but it expires when the wall closure lands.

- **⚠ A WITNESS MUST BE SCORED ON THE SAME POPULATION AS THE THING IT MEASURES
  (2026-09-28, D2z — three instances in one study).** (1) `delta_mean` (`:3315`)
  accumulates `StepCellWeight` (the per-cell share, bounded near 1) while the
  kernel's own `g` (`:1687`) and the offline estimator use `PairDelta` (the pair's
  full height, unbounded): identical at `n = 1`, **−15.4 %** at 8.3 % multi-cell,
  **−87.5 %** at 26.6 %. Its normaliser is wrong too — `delta_n += (nfx + nfy)`
  counts every exposed face, including unstepped ones that contribute nothing.
  (2) `WallStepDiagnostic` (`:2109`) has **no radius gate**, so the patch edge
  (outer column never fired ⇒ `phi = 1`, `col_dt = 0` ⇒ `a = 0` ⇒ **`g = 1`
  exactly**) outnumbers the real faces: mean `g` read **0.843** where the in-patch
  value is **0.335**, so criterion (b) scored **3.10x** against an 8x bar with the
  per-face term contributing **1.03x** — the population term (3.00x) was doing all
  the work and the per-face term was simply unmeasured. (3) Scoring the estimator
  at r < 25 mm while the live population is the r < 30 mm patch produced a spurious
  **2.7x failure on correct code**. **Pin the disc, and say which one, every time.**

- **⚠ A CRITERION BAR MUST NOT BE TIGHTER THAN THE INSTRUMENT (2026-09-28, D2z).**
  The power-identity bar was set at ±0.01 W; the thermo log prints **6 significant
  figures**, so three ~10³ W columns quantise at ±0.005 W each ⇒ a **0.015 W
  floor**. `S0` measured **0.0141 W** — a rounding artefact that reads as a failure.
  Use **≥0.02 W**. ⚠ And `jet_P_face` exists only under `jet_closure`, which the
  hot-disk arbiter does not run, so **the only check that catches a kernel↔march
  split cannot be evaluated on any arbiter leg** — it is Meier-only.

- **⚠ `side_face_area_law` DEFAULTS TO `step` (`:6819`), and `run.py`'s `SHARED`
  block does not set it (2026-09-28).** A Meier leg that omits the key silently
  runs the law the D2z tie rule did **not** pick. Cheap witness:
  `riser_area_ratio < 1` ⟺ `wetted` (measured 0.325 / 0.433 on wetted legs against
  5.765 / 2.340 on step legs). Same class as `side_face_factor`'s 1.0 default
  (`:6748`), which makes a missed call site numerically invisible.

- **✅ `f` IS GONE, AND THE RIM ARTEFACT WENT WITH IT (2026-09-28, D2z).** Both
  area laws cleared the 5 % gate and the 2 % "as good as the fudge" bar at
  r < 25 mm (`step` +1.39 %, `wetted` +1.06 %, reference +0.38 %). At r < 30 mm the
  25–30 mm rim annulus — the artefact that survived D2r-0, D2r-0b, D2r-0c, D2r-1
  and D2s — fell from **−9.11 %** under `f = 0.2` to **−0.71 %**. And the side
  heating **self-regulates**: `wetted` charges 5.8x less area than `step` yet
  delivers **more** side power (25.3 vs 21.4 W), because less area ⇒ taller steps
  ⇒ more area. Both deliver less than `f = 0.2`'s 33.0 W while healing at least as
  well. ⚠ **Every one of these numbers was measured on risers one cell tall or
  less.** No claim about multi-cell risers is supported until `SW` runs.

- **`delta_fallback_cols` — BOTH READINGS ARE RIGHT, ON DIFFERENT WINDOWS
  (2026-09-28).** **0 on every leg in the scored 50–200 s window.** Over the full
  run, which includes the 0–50 s startup where columns have fired fewer than twice,
  means are 0.78 / 6.77 / 2.56 / 32.03 and maxima 39 / 73 / 156 / 331 on
  `S2` / `S2w` / `S1` / `S1w`. `PREFLIGHT.md:35`'s "0.00 %" is a full-run claim and
  is refuted; the scored-window claim holds. **State the window.**

- **⚠⚠ EVERY MEIER THERMO POWER IS A QUARTER-HOLE SUM (2026-09-28).**
  `surface_patch.x0 = 0.0`, `y0 = 0.0` on a 0.12 x 0.12 m domain — the burner sits
  in the **corner**, and `d2x_areasplit/analyze.py:104` already multiplies volume by
  **4.0**. A planner energy check compared 2.1 kW against a 1–2 kW ceiling and was
  **off by 4x**. Full-hole: `A` **4.18 kW**, `SW` **8.26 kW**, `SW`'s side path
  **3.54 kW**. ⚠ **AUDIT D2l's "17.3 % of HHV absorbed, needs 9.1 %"** for the same
  factor before quoting it again.

- **⚠⚠ `A`'s IN-BAND RATE WAS THE STEP-DRAIN THROTTLE — do not build on it
  (2026-09-28).** D2x's own analysis, line 108: **`C` 2.09, `E` 2.11 m/h** is the
  un-throttled pad-ring top-face rate. `A`'s **1.251** was the drain. D2z removed
  that drain, so **any leg that removes the wall heat un-throttles the pad ring.**
  Two arguments have now been built on `A`'s in-band rate and both were wrong.

- **⚠⚠ `side_faces` COUNTS FACES THAT CARRY NO AREA — and it cost the planner a
  premise (2026-09-28).** `side_n` increments for **every** exposed lateral face,
  including above-plane ones whose `sf_c = side_face_factor_above = 0` gives
  `Aside = 0`. On `SW`: **6031 exposed, 5081 (84.3 %) above the plane and zeroed,
  950 below** — and `delta_faces` is **also 950**, so the step branch *is* the
  heated set. ⚠ **AND the area law is applied to `StepCellWeight` (≈1, factor
  0.414), NOT to `delta_mean`** (2.174, factor 1.393) — using the latter is **3.4x
  too large**. Together those two errors produced a spurious "h ≈ 35 W/m²K, so the
  lever is area". **Corrected: the charged area is ~16–18 cm² and the effective
  lateral coefficient is ~660–750 W/m²K — the map's own impingement value (685 at
  r = 40 mm).** ⇒ **THE WALL LEVER IS `h`, AND THERE IS NO AREA EXCESS.** Von Rohr's
  grazing band with the march's exhaust temperature is a **7–9x** cut, 3.54 kW →
  **0.4–0.5 kW**, inside the derived target. ⚠ To size a flux from thermo
  aggregates, use `side_faces − side_cols_above` and the per-cell weight, never the
  raw count and never the pair delta.
  ⚠ `surface_patch.side_face_h_max` (default **80**) is separately **INERT**:
  `:3257`/`:1562` gate it on `spall_per_cell_removal`, which defaults to **0**
  (`:4787`).
  **Measured split of the charged side area, per quarter: bore 1, pads 10, wall
  7 cm²** — more than half is **under the pads**, not on the wall.

- **The model's specific energy is 1.35 J/mm³, and that is a SEPARATE finding if the
  literature is right (2026-09-28).** `SW`: 8.26 kW over 6.1 cm³/s — the sensible
  heat to its own pinned surface temperature, as it must be for energy-set
  recession. With that physics the model **cannot** reach Meier's 2.47 cm³/s at
  8 kW. ⚠ If rock really takes **4–6 J/mm³** the spall energetics are wrong by
  **~3x** — testable against Wilkinson and Tester's measured surface temperatures.
  **It is not a licence to keep the power.** ⚠ The 4–6 J/mm³ band remains an
  **assumption, not a measurement** (D2l).

- **⚠ A FALSE PASS IS AVAILABLE ON ANY SHIELDED-PAD LEG (2026-09-28).** Under the
  pads a third of the multi-cell side area sits on faces that get the **full**
  conductance because they are real steps, and after shielding would get **little
  heat** — the **`f = 0` configuration on that ring**, which can re-cascade. **An
  in-band rate with the 30–40 mm rings falling the way `A`'s did (125 → 90 → 69 mm)
  and the pad-ring multi-cell fraction RISING is the throttle returning, not a
  pass.** Score it, do not report it.

- **✅ SPECIFIC ENERGY IS CONSTANT AT ~1.32 J/mm³, SO HEAT CANNOT BE MOVED TO BUY
  RATE (2026-09-28, D2aa).** Identical on `A`, `SW`, `WL`, `WL81` — it is the
  sensible heat to the pinned surface temperature, which energy-set recession
  requires. ⇒ **rate follows the TOTAL power (`rate ∝ P^1.59`, four points, ±4 %);
  only shape follows the distribution.** Meier's 1.3–1.6 m/h band needs
  **4.17–4.75 kW**. ⚠ The published **4–6 J/mm³** band that would make this a ~3x
  error is recorded in D2l as **an assumption, not a measurement** — settle it with
  the Rauenzahn/Tester flat-face check, not by citation.

- **⚠⚠ THE SIDE TERM DEFENDS ITSELF AGAINST A LOWER `h` (2026-09-28, D2aa).**
  Cutting the lateral coefficient **8.5x** (685 → 81 W/m²K) cut the side power only
  **4.3x**, because the charged riser area **doubled**: 68 → 133 → 143 cm² as `h`
  fell 685 → 142 → 81. Less heating ⇒ steps heal less ⇒ more exposed riser area ⇒
  more power. Same feedback D2z saw at 2 mm (`wetted` delivered **more** side power
  than `step` on 5.8x less area). **Never assume side power scales with `h`.**

- **⚠⚠ A FROZEN GAS TEMPERATURE IS SELF-INCONSISTENT — THE MARCH CONSERVES ENTHALPY
  (2026-09-28, D2aa).** `WL`/`WL81` fed the wall **1463.5 K** measured on `SW`,
  while their **own** exhaust ran at **1601.3 / 1631.0 K** — and `A`, which absorbs
  nothing on the wall, at **1678.8**. **Gas that gives up less heat to the wall
  stays hotter.** Driving ΔT was **22–27 % too small**, and correcting it moves
  `WL81` from 1.78 to **1.91 m/h** — the wrong way. **Any leg that changes how much
  heat a surface absorbs must read the gas temperature LIVE, never frozen from
  another leg.**

- **⚠ THE RATE IS NOT STEADY — report it as a profile against depth (2026-09-28).**
  Disk-mean removal inside r < 25 mm on `A` falls **2.35 → 1.92 → 1.70 m/h** across
  100–200 / 200–300 / 300–400 s. The scored metric is the **depth-based** 100–150 mm
  descent stretch, which is window-independent by design — but a single quoted rate
  hides this and must not stand alone.

- **⚠ `PD` JAMMED: shielding the pad ring is BLOCKED BY LATERAL RECESSION
  (2026-09-28, D2aa, pre-registered and confirmed).** With the ring shielded,
  rings 30–38 mm read **0.0 at every one**, the ring froze, and the burner sat on it
  and stopped descending — **the model's lateral-recession limit (undercut is
  impossible by construction), NOT physics.** ⚠ `pad_shield_open = 0.0` was the
  **bounding** fully-closed case because **the slot width is on no record** (needs
  drawing sheet 15). ⚠ `side_P_pads` goes **negative**: shielded faces receive
  nothing and still radiate. **Do not re-attempt the pad ring until lateral
  recession exists; report it as a limitation.**

- **⚠⚠ `side_face_h` REACHES EVERY LATERAL FACE, AND TWO-THIRDS OF THEM ARE NOT ON
  THE WALL (2026-09-28, D2ab pre-flight).** Measured from D2aa's region columns on
  `WL`: side **power** bore 13.9 % / pads **50.7 %** / wall 35.5 %, side **area**
  13.9 / 49.5 / 36.7 %. The bore-floor staircase risers and the pad-ring step faces
  face into the bore, in **impinging** flow — an annulus coefficient has no business
  on them. ⇒ **D2aa's `WL`/`WL81` are the "all-lateral-faces" lever, NOT the wall
  lever**, and two-thirds of the reduction they measured was misapplied. Fix is a
  radial gate (`side_face_h_r_min`), and **the correction makes the rate WORSE** —
  side power to ~2.5 kW, derived legs near 2.3–2.9 m/h. ⇒ **The binding term is the
  bore and pad-ring lateral faces, which need LATERAL RECESSION — now the critical
  path.**

- **✅ THE 4–6 J/mm³ THREAD IS CLOSED — it was a units mismatch (2026-09-28).** The
  model's **1.32 J/mm³ is delivered-to-rock**; every literature figure in the repo is
  **fuel-basis**: Meier's TSE **15.4 J/mm³ on HHV** (Tab. 8.2,
  `meier_thesis_notes.md:83`), Browning/Fleming **8.7–60**. On a common HHV basis the
  model runs **38 kW of fuel for 6.1 cm³/s = 6.2 J/mm³ against Meier's 15.4** — a
  factor of 2.5 that is **the rate error restated, not a separate energetics error.**
  ⚠ **And a prescribed-flux flat face cannot test it**: at prescribed flux the
  recession rate is `q/(rho c_p dT_spall)` by energy conservation, which the constant
  1.32 J/mm³ already proves the model obeys, so such a run returns the rate law **by
  construction** and tests only the **spall temperature**. ⚠ **No flux-versus-rate
  dataset exists anywhere in the repo.** What does exist is Rauenzahn's compilation of
  measured **surface-temperature rise** (Fig. 2.12, `:123`), and the model's **528 K
  already sits inside it** (`:222`). Do not reopen this as an energetics packet.

- **⚠⚠ MEIER'S RATE IS ONE-SIDED AND THE MODEL COMPUTES A DIFFERENT QUANTITY
  (2026-09-29 — recorded 09-18, then scored against for eight packets).** The burner
  "hung from the hall crane on coil tubing, **keeping the coil essentially in
  tension**" (p. 209); "the **hook load** … giving evidence whether the burner
  **hangs or lies**"; the hole is "only slightly underreamed allowing the necessary
  **clearance for the contactless burner to advance under the effect of gravity**"
  and "the **clearance to the hole is adjusted by changing the advancing rate**"
  (p. 223); **the operator set the advance rate** (notes `:214`). Every Meier leg
  runs `surface_patch.nozzle_descent=feet` — **the CONTACT LIMIT**. ⇒ **Meier's
  1.5 m/h is a LOWER BOUND on rock capability; 3.9 ≥ 1.5 is CONSISTENT, not
  refuted.** **Never score burner descent against 1.3–1.6 as if it were two-sided.**
  Under a prescribed feed the **rate is an input** and the discriminator is the
  **diameter**.

- **⚠ THE DECAY-CLIP WARNING WAS BACKWARDS AND IS RETIRED (2026-09-29).** D2aa item
  0(c) and the planner's 09-28 advice said the stagnation decay is "clipped at twelve
  nozzle diameters" so no equilibrium standoff could form. Read from `A`'s own argv
  the term is **`pow(0.09/max(0.09,s),1.0)`**: the clip is at the **NEAR** end (no
  decay inside 90 mm) and beyond it the flux falls as **1/s without limit**. **An
  equilibrium standoff exists** — ≈222 mm at a 1.5 m/h feed, scaling from `SWc`'s
  3.694 m/h at 90 mm.

- **⚠ `side_faces` vs the WALL BAND: removing wall heat LATCHES it off (2026-09-29,
  D2ab).** `GL`/`GH`/`WLr` all report `side_P_wall = 0.0 W` and zero recession beyond
  r = 38 mm. Less wall heat ⇒ narrower hole (Ø79.6 vs 87.6) ⇒ the wall moves **inside**
  the 40 mm gate ⇒ no faces for the coefficient to act on ⇒ narrower still. Bore and
  pad power agree to **0.35 %/1.17 %** across an `h` span of **4.4x**, so at the
  contact limit **the annulus coefficient is inert by geometry, not by physics.**

- **✅ A RATE RELATION THAT WORKS, AND ITS FAILURE MODE (2026-09-29, D2ab).**
  `rate = (P / J_specific) / (pi D^2/4)` predicts burner descent **within 3 %** on
  every leg whose pit lead is small (`GL` 3.78 vs 3.897; `SWc` 3.60 vs 3.694) and
  fails progressively as the lead grows (`WL81` +34 %, `A` +64 % at a 128 mm lead).
  **Volume that goes into pit depth does not go into burner descent.** Use it to
  predict, and read a large residual as a deep pit rather than as a rate error.

- **⚠⚠ FIG. 8.8 BOUNDS THE PIT LEAD, AND IT KILLS THE "CATEGORY ERROR" EXCUSE
  (2026-09-29, verified by the planner from the PDF).** Rendered pdf p. 251 at 200
  and 500 dpi and contrast-enhanced: the hole is a **straight, essentially
  constant-width channel** terminating in a **rounded tip of roughly its own
  radius**, just above the bottom of the 0.5 m block. **There is NO narrow pit ahead
  of the cylinder**, and the saw plane passes through the axis, so one would have
  been sectioned; a 17–21 cm lead would also have broken out of the block bottom
  while the burner was at ≈30 cm. ⇒ **the floor led the burner by ≤ ≈5 cm after
  1383 s.** Meier, p. 222, verified verbatim: *"it is likely that the penetration
  rate was limited by the power of the combustion and not by poor heat transfer."*
  ⇒ **The rate gap is REAL.** Under `pow(0.09/max(0.09,s),1)` the flux at s ≤ 100 mm
  is ≥90 % of its s = 90 value where `SW` recedes 3.694 m/h, so the model's floor
  runs **3.32 m/h** against Meier's **1.28–1.5** — a **2.2–2.6x floor-flux excess**.
  ⚠ **Still true: burner descent is one-sided as a KINEMATIC quantity** (the
  operator set the feed) — but the FLOOR position is not free, so the flux
  comparison is two-sided.
  ⚠ **The WIDTH is still not measurable from this figure** (edge-finding gave 92–98
  mm on some rows, noise on others — the 09-22 note stands). **85 mm comes from
  Tab. 8.2, never from the photo.**

- **⚠ MEIER'S 3.42 L IS NOT A DIAMETER (2026-09-29).** The water-filled volume
  includes the ≈10 cm-wide **mouth**, the **pause underream** at ≈8–10 cm (Meier:
  "induced by holding the burner steady for several seconds", p. 223) and the
  **rounded tip**. Converting it to a 93.3 mm equivalent bore and scoring diameter
  against it is wrong — a uniform 93 mm bore **contradicts** Fig. 8.8. Score
  **Ø 85 ± 5 mm on the calibrated 15–45 cm section** and **volume ≤ 3.42 L**, with
  the shortfall attributed to the mouth and the pause. ⚠ Do **not** score Ø at
  100 mm depth: that is the pause underream.

- **⚠⚠ P0b (h) FAILS AT THE PERIPHERY ONLY, AND NO BOUND CAN FIX IT (2026-09-30).** Cap grip
  by radius, 60–160 s: **0.0 %** inside 20 mm, **4.3 %** at 20–30, **37.9 %** at 30–40,
  **99.8 %** at 40–50, **100 %** beyond 50 mm. The floor under the burner is fully resolved
  and physical; the cap is the whole thickness law only on the **side-heated** wall, whose
  profile genuinely is flat. Those columns carry **36.2 % of all rock removed** ⇒ **the (h)
  failure and the over-wide hole are ONE defect.**
  ⚠⚠ **The unbounded quantity DIVERGES, it does not merely overshoot:** `ell_eff` p99 is
  **2.9 m**, max **1.6 km**, and **3.0 %** of firings exceed the whole 600 mm block, because
  `ell = dz/ln(ratio)` and the ratio → 1. (`RESULTS_REV7.md`'s "p99 15.0 mm" is `h_f` **after**
  the search-depth cap.) So ~4 % of firings are flat to round-off and **even a 50 mm bound
  grips 4.3 %** — the 20 % bar is unreachable by any bound inside Wilkinson's 0.1–1 mm band.
  ⚠ The grip keeps growing past the reported span: **7.3 % → 47.7 %** overall across the run,
  and **4.7 % → 30.3 %** inside 40 mm, so "interior passes" is a scored-window statement.

- **⚠⚠ THE AXIS RATE IS NOT AN INDEPENDENT TEST UNDER PRESCRIBED FEED (2026-09-30).** The
  burner is forced down at 1.5 m/h. The gap changes at (floor rate − feed) and a wider gap
  cuts floor power — one stable resting point, at floor rate **= feed rate by construction**.
  Measured over rev 7's 500 s: rate **5.84 → 2.01 m/h**, gap **72 → 223 mm**, floor power
  **1690 → 631 W**, side share **0.22 → 0.58**, all still moving at the end. ⇒ a windowed
  rate is a point on a decaying transient; the comparison to revision 2 is valid **because
  both were measured the same way**, but it is not "the model's rate". **The free quantities,
  and therefore the only real tests, are the GAP (223 mm vs Meier's ≤100) and the WIDTH.**

- **⚠ REMOVAL IS ENERGY-LIMITED — VOLUME IS THE WRONG INVARIANT (2026-09-30).** P0b rev 7
  reported "P6 refuted" because volume fell 34 % while power fell 15 %. That is a misreading.
  Flake enthalpy rate **1192 → 1590 W** while delivered power fell **2457 → 2080 W**; energy
  per cubic metre removed went **4.40e8 → 8.87e8 J/m³**, i.e. **exactly doubled**. Revision 5
  shed thick cold slabs (184 K mean superheat), revision 7 sheds thin hot ones (372 K). Less
  rock, hotter rock, same energy. **Hold enthalpy fixed, never volume.**

- **⚠ QUARTER DOMAIN: THE ×4 IS REQUIRED FOR ANY MEIER VOLUME COMPARISON (2026-09-30).** The
  burner sits at the domain corner with mirror planes on two sides, and the reduction scripts
  report the quadrant. `RESULTS_REV7.md` quotes 1.793e-6 m³/s and never compares it to Meier.
  Full hole: **7.18e-6 m³/s in the 60–160 s window (2.90×)**, **5.79e-6 over the whole run
  (2.34×)** against Meier's 2.473e-6. Across revisions 2/5/7 the volume is **2.80× → 3.57× →
  2.34×** — rev 7 is the **best of the three**, which the write-up does not say.

- **⚠ A DIAGNOSTIC THAT HAS NEVER BEEN NON-ZERO HAS NOT BEEN TESTED (2026-09-30).** P0b
  revision 2's "the `h_f` cap binds on 85 % of columns but **0 % of firings**" was reading an
  `rc_flag` bit that **was never written** — the flag encoded flat/`ell`/`T_surf` only. The
  conclusion survived two revisions. ⚠ Related: **grep for a key name before adding one** —
  `surface_patch.side_face_T` already existed (D2ab's side gas temperature) and only a smoke
  run caught the collision.

- **⚠ SIDE-POWER BINS DECOMPOSE *GROSS*, NOT NET (2026-09-30).** `surf_side_P_rin + _rout =
  _top + _buried = 982.4344 W` (agreeing to 1.0e-5 W) `= side_P_face + side_P_loss =
  906.647 + 75.787`. No double count, but the live-gate "38 % inside / 62 % outside" is a
  share of **982 W**, not of the 906.6 W entering the rock. Also: **"min standoff" is the
  window MEAN** of `patch_min_standoff` (itself a per-step minimum over columns); the true
  window minimum is 106.5 mm. The mean reproduces all four revisions consistently, so
  comparisons are sound — but relabel it.

- **✅ DEEPEN A BLOCK DOWNWARD (2026-09-30).** Extend `geometry.prob_lo z`, not `prob_hi`, so
  the surface, the standoff and `dz` do not move. Verified as exactly 100 mm of added rock,
  first firing 0.288 vs 0.304 s.

- **⚠⚠ THE KANT ANCHOR SCORES A VACUOUS VARIABLE (measured 2026-10-01).** It scores the firing
  **superheat**, which is the criterion's own threshold, so mesh-independence there is almost
  tautological. The physical observable — **when** the rock spalls — is 150× worse:

  | leg | firing TIME, 2 / 1 / 0.5 mm | spread | scored SUPERHEAT spread |
  |---|---|---|---|
  | `R150` | 19.000 / 21.775 / 23.100 s | **21.6 %** | **0.16 %** |
  | `R500` | 1.437 / 1.450 / 1.812 s | **26.1 %** | **0.12 %** |
  | `C150` (control) | 33.187 / 28.637 / 26.412 s | 25.7 % | 14.99 % |
  | `C500` (control) | 5.113 / 3.563 / 2.850 s | 79.4 % | 43.86 % |

  Same class as the standing note "a threshold criterion's own variable is a VACUOUS mesh
  witness". ⚠ Note it is vacuous **only in the scored leg** — the controls spread in both — so
  it is a real blind spot, not a useless metric. The 21–26 % is the **cold-start reconstruction
  error** (the second cell has not warmed ⇒ ratio too large ⇒ layer read too thin ⇒ surface
  over-read ⇒ premature firing), which is exactly what the flux-anchored estimator would fix.
  ⇒ **Score the anchor on firing time**, and fix the scoring **before** judging the estimator.

- **⚠ ONLY 29 % OF REMOVED ROCK LIES BEYOND MEIER'S WALL (2026-10-01).** Meier's bore is
  Ø 85 mm ⇒ r = 42.5; the model's is Ø 101 ⇒ r = 50.5. Cumulative share of removed rock by
  radius: 15.9 % inside 20 mm, 36.0 inside 30, 56.6 inside 38, **70.7 inside 42**, 86.5 inside
  48, 94.2 inside 54. ⇒ the bore error is worth **29 %** of the volume, not the 36 % the
  40 mm gate suggested, and the capped population the criterion suppresses **straddles** the
  wall (38–42 mm is inside the bore, 14 % of the rock, already 83 % capped at fraction 0.234).
  ⇒ **the criterion cannot cut cleanly at 42.5 mm**; the map's radial tail must do the rest.

- **⚠⚠ THE 40–50 mm BAND IS FLOOR, NOT WALL — AND P1's GAIN THERE IS ARTEFACT-ASSISTED
  (2026-10-01).** Measured delivered top flux: **1.19–1.45 MW/m² at 38–42 mm** (the highest in
  the hole, because the map peaks where the local standoff is shortest, at the shallow rim of
  the bowl), **1.02 → 0.65 at 42–50 mm**, 0.45 → 0.24 at 50–56, ≤0.2 beyond 56. ⇒ those columns
  are **floor under floor-intensity heat**, and their flat two-cell reading is the fit's
  blindness to a **thin hot skin on a substrate pre-warmed from the side** — not a side-heated
  wall. Three consequences: (i) holding the flux-anchored estimator out of P1 is right and the
  measured cancellation is correct physics — anchoring to 1 MW/m² reads the skin and fires those
  columns, **as rock under that flux must**; (ii) **the flake reference's suppression at
  40–50 mm is the right answer for the wrong reason** and its contribution to (e), (e′) and (f)
  must be recorded as **artefact-assisted, not validated shape physics** — the genuine content
  is the band **beyond 50 mm** plus the slab and Kant behaviour; (iii) the width at 40–50 mm is
  **the map's radial tail**, which must run **before** the estimator fix or the bore reopens and
  the estimator is blamed for it.

- **✅ THE ESTIMATOR'S ROOT CAUSE: BOTH FITS ASSUME DECAY TO AMBIENT (2026-10-01).** The
  two-cell fit and the flux-anchored fit share that assumption, and it is why both misread the
  wall band. The right generalisation is a **three-parameter fit with the substrate temperature
  as the base** — top cell, cell below, and the flux: three data, three unknowns (surface rise,
  decay length, base). It resolves a thin skin on pre-warmed rock **and fires it**, and reads
  near-zero drive on a genuinely side-heated column with little top flux. The adopted residual
  (**0.78 on the axis → 0.10 at 40–50 mm**) is exactly the signature of a non-ambient base.

- **⚠⚠ NO MEIER RATE, VOLUME, WIDTH OR STANDOFF IS QUOTABLE FROM P1 (2026-10-01).** The only
  converged shape number is the **collar: Ø 93.8 mm at 10 mm depth and 88.2 at 25 mm** against
  Meier's 96 and 92 — inside Fig. 8.8's ±5 mm, and **the first correct mouth** (rev 7: 124.8 /
  112.9). ⚠⚠ It must never be quoted without the cone (**72 mm at 100 mm, 56 at 150**) and the
  stall beside it. Volume per mm of axis advance **3.1 mL/mm**; `ledger_P_loss` 135 → **1331 W**;
  side power inside r = 40 mm 37.8 → **93.7 %**.

- **⚠ MEIER'S VOLUME PER mm OF AXIS IS 5.93 mL/mm, NOT 7.1 (2026-10-01).** 3.42 L over the
  **576 mm** the burner descended (1.5 m/h × 1383 s), which back-solves to an **86.9 mm** bore,
  matching the reported 85–87. The 7.1 used in P1's write-up implies a 482 mm hole. ⇒ **rev 7
  was +13 % OVER Meier on this measure, not 6 % under**, and P1 is **−47 %**, not −56 %. The
  overshoot runs from +13 % to −47 %.

- **⚠ A CRITERION WHOSE DENOMINATOR THE TREATMENT REMOVES CANNOT BE SCORED (2026-10-01).**
  P1's (e′) asked for cap grip < 20 % beyond r = 40 mm; the criterion **deletes 10,771 of
  10,809** firings there and the cap still binds on **94.7 % of the 38 survivors**. Score the
  substantive claim — **grip over ALL radii, 30.6 % → 8.6 %** — and report the survivors beside
  it. ⚠ Deliverable 5's wall-band witness has the same shape (asks for ≪ 0.1; lowest band median
  **0.191**, nothing beyond 50 mm fires). **Became unobservable rather than satisfied.**

- **⚠ (c)'s SHIPPED-ONSET CLAIM HOLDS ON ONE LEG OF THE PAIR (2026-10-01).** Measured at 2 mm on
  the P1 binary: `S150` **461.76 K** superheat (reproduces rev 7's 463.3, −1.5 K) but `S500`
  **478.33** and `R500` **482.28** — the quoted "465.8" is **rev 7's `R500`**, not a P1 leg. ⇒
  the 500 leg moves **+12.5 to +16.5 K**; "where it fires on the anchor does not change" is true
  for the 150 leg only, which weakens the "chosen on a witness, not to pass" defence by that
  much. (c) still passes comfortably and the refit is still in band.

- **✅ A DIFFERENCE IS ONLY WELL DEFINED ON ONE PROFILE (2026-10-01).** Leaving the crack side on
  the piecewise-linear `T_at_zv` while the reference read the exponential cost **2.4 % mesh
  spread at `R500`**, falling to **0.97 %** once both ends read one profile. ⚠ `K_Ic`'s tip
  temperature **still** reads `T_at_zv` — the named cause of (d)'s 2.85 % — and fixing it changes
  **every** reconstructing leg, so it is a successor, not a patch.

- **⚠⚠ THE REFIT AND THE CRITERION VERY NEARLY CANCEL — AND THAT WAS NEVER PRE-REGISTERED
  (2026-10-01).** `A80`'s rate **level** is **6.0 m/h against `A08`'s 2.5**: the 80 µm flaw alone
  **multiplies the removal rate by 2.4**. The differential criterion pulls it back to 2.476,
  i.e. `A08`'s level to within **−1.2 %**. ⇒ (d′)'s "the refit restored the level" was **not
  luck and not neutrality** — enlarging the flaw makes firing much easier and applying only the
  0.330 stress jump makes it correspondingly harder. **Neither change is separately defensible
  at these settings; they are defensible only together.** The successor must treat `a0` and the
  drive fraction as **one** calibration, not two.

- **⚠ THE EXPOSED-FACE CENSUS REPRODUCES EXCEPT ONE ROW (2026-10-01).** `r < 40 mm` with
  `t < 120 s` is the **unique** cut returning both headline rows (152/316 frozen, 78 % with ≥ 2
  exposed sides); at that same cut **live columns with 0 exposed sides is 51 %, not 69 %**. Six
  exposure definitions and a radius/threshold sweep were tried; none returns 69 % while also
  returning the headline rows. ⚠ The correction **strengthens** the roughness reading — even
  *live* columns stand proud of a neighbour about half the time. Reproduce with
  `studies/p1_gradient_criterion/census.py` against event-log sha `26d8996f26ae2c90…`.

- **⚠⚠ `a0` AND THE DRIVE FRACTION ARE ONE CALIBRATION KNOB, NOT TWO (measured 2026-10-01).**
  P1's control legs: the flaw enlargement **alone** multiplies the column rate by **2.46×**
  (`A08` absolute/8 µm **2.506** → `A80` absolute/80 µm **6.158 m/h**) and the differential
  criterion divides it back by **2.49×** (→ `S` **2.476**, **−1.2 %** of `A08`). The cancellation
  is **structural, not luck**: holding Kant's onset fixed requires
  `√(a_new/a_old) · 0.3297 = 1`, i.e. a ratio of **9.20**, and 80/8 = 10.0. On a self-similar
  front the drive fraction is the same 0.3297 on Kant and on the column, so the fit transfers
  exactly. ⇒ **(d′) is a check of the refit arithmetic, not independent evidence that the change
  is benign**, and **(c) carries two things of which only one is refit-independent**: the form's
  exactness (1 − e^{−c_f} to 16 digits on two physically different legs) is a property of the
  flake reference alone; the onset *location* in temperature is not.
  ⇒ **ANY future packet that changes how sharp the model reads the near-surface profile must
  refit `a0` against Kant IN THE SAME PACKET and re-measure the column level** — otherwise the
  rate moves by a factor for no physical reason. The two queued successors (the `K_Ic` tip
  temperature on the exponential; the three-parameter estimator) both do this. It is also the
  right frame for the outstanding `a0_gb` = 400 µm justification.
  ⚠ Not an alarm about the rate **level** — that is P2's item (blowing `h_eff = h·B/(e^B − 1)`,
  `cp(T)`), with literature magnitudes, and was never P1's to fix.

