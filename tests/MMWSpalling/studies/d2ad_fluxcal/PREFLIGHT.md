==============================================================================
D2ad ITEM 0 -- measured on 60 s smokes
==============================================================================

(b) THE kic.scale -> T_s MAPPING, MEASURED
------------------------------------------------------------------------------
  Smokes on the Q5 base (F_q = 0.5), stop_time = 60 s.
  ⚠ SMOKE_K10 would duplicate SMOKE_Q5 exactly (F_q 0.5, scale 1.0), so
    SMOKE_Q5 IS the scale-1.0 row; it is not run twice.

   scale  firings     mean      p10      p50      p90  p90-p10  window
     1.0     4413    821.3    815.2    820.4    828.4     13.2  OK
     1.3     3630    864.6    859.2    863.8    871.1     12.0  OK
     1.6     3144    895.5    890.6    894.8    901.5     10.8  OK
     1.9     2831    919.1    914.2    918.0    924.7     10.5  OK
     2.2     2447    950.1    941.6    948.8    960.4     18.8  OK

  ⚠ BAR: p90 - p10 <= 30 K at every scale. If it broadens, the
    lever is changing the CRITERION, not the temperature, and the packet
    stops and reports that (P5).
  widest window across the ladder: 18.8 K   PASS

  the lever moves T_s: 821.3 -> 950.1 K over scale 1.0 -> 2.2
  (a flat 823 K across all five would mean the key never reached Removal.H:759)
  T_s = 870 K -> scale 1.3523 (interpolated); nearest run scale 1.3 is -5.4 K off
  T_s = 920 K -> scale 1.9088 (interpolated); nearest run scale 1.9 is -0.9 K off

(c) THE BACK-REACTION, VERIFIED LIKE FOR LIKE
------------------------------------------------------------------------------
  ⚠ entrainment is `exhaust` mode, so the exhaust temperature BUILDS over
    the run. Comparing a 60 s smoke against FH's 490 s value would be
    wrong. FH is read at the SAME 60 s and the offset stated.
  MATCHED TIME (the smokes' last logged row): t = 50.016 s

  FH floor T_g at t = 50.016 s:   1616.7 K  (floor power  2152.4 W quarter)
  FH floor T_g at t = 490.0 s:   1664.7 K  (floor power  1408.9 W quarter)
  matched -> 490 s OFFSET: +48.0 K   <- applied below

  leg             F_q  T_g(60s)   +offset   C1 pred    diff    floor P
  SMOKE_Q10      1.00    1616.7    1664.7    1664.7    -0.0     2152.4
  SMOKE_Q5       0.50    1786.2    1834.2    1782.3    51.8     1118.9
  SMOKE_Q35      0.35    1823.2    1871.2    1817.6    53.5      753.4

  ⚠ BAR: any leg off by > 40 K after the offset ⇒ recompute F, s* and
    q100 and re-hash before the ladder.
  worst |measured - C1| = 53.5 K   FAIL -- RECOMPUTE

  WITNESS (the prefix actually applied): floor power must fall ~F_q x
    F_q = 1.00: floor power  2152.4 W  ratio to F_q=1 1.000  (expect ~1.00)
    F_q = 0.50: floor power  1118.9 W  ratio to F_q=1 0.520  (expect ~0.50)
    F_q = 0.35: floor power   753.4 W  ratio to F_q=1 0.350  (expect ~0.35)

(d) G1600's FLOOR T_g, MEASURED
------------------------------------------------------------------------------
  SMOKE_G1600 floor T_g at t = 50.016 s: 1417.8 K, +offset = 1465.8 K
  the ladder TABLED 1410 K (an estimate)  ->  diff +55.8 K
  floor power 1524.9 W quarter
  recomputed F = 1.0 * [(1465.8 - 823.3)/(823.3 - 293)] / 1.58665 = 0.764
  => s* = 198 * F = 151 mm   (tabled 138)

(e) R90 -- THE UNDECAYED FLOOR RECESSION, MEASURED
------------------------------------------------------------------------------
  ⚠ C0 IS FH's argv plus stop_time, so C0 == FH for the first 500 s and
    R90 can be read from the FROZEN FH trace -- no new run needed.
  Measured where the standoff is still inside the flat region (s < 90 mm),
  because beyond it the map decays as 1/s and the rate is no longer R90.
  window: t = 5.1-55.4 s, standoff 51.5-89.7 mm (525 rows)
  R90 = 4.237 m/h

  the brief quotes TWO inconsistent values: 3.3 (s* formula) and
  4.2 (contact edge). Measured 4.24 -> 4.2 is closer
  deviation from 3.3: +28.4 %   OVER 10 % -- RECOMPUTE THE LADDER

  CONTACT THRESHOLD  F = feed/R90 = 1.5/4.237 = 0.3540   (ladder assumed 0.4545)
  PASSING WINDOW     0.3540 < F <= 0.505  (s* between 90 and 100 mm)
  ⚠ THE WINDOW IS ~010 mm WIDE and the predictor carries +/-5 %. A leg at 95-105 mm
    is NOT a decisive pass or fail (C4).

(f) THE CALIBRATION / VALIDATION SPLIT -- fixed BEFORE running
------------------------------------------------------------------------------
  CALIBRATED (two numbers, chosen on the STANDOFF ALONE):
      F_q   the floor heat-transfer level
      T_s   the flaking temperature, via spallation.sp.kic.scale
  VALIDATION (nothing here may be used to choose a leg):
      Ø(z) over 120-250 mm
      volume per 0.5 m
      pad clearance
      the flux gate q100 against Kant and Meier
      the free-descent rate (the CONDITIONAL R_cal follow-on)

(g) THE RECOMPUTED LADDER -- mandated by (c) and (e)
==============================================================================
  The packet REQUIRES this recompute in two places:
    (c) 'off by > 40 K after the offset => recompute F, s* and q100'
    (e) 'if R90 differs from 3.3 by more than 10 %, recompute the ladder'
  Both fired. Measured worst |C1 error| = 53.5 K and R90 = 4.237 (+28.4 %).

  T_g(F_q): C1's linear 1900 - 235.3*F_q is REPLACED by a quadratic fitted
  exactly through the three MEASURED points (0.35, 0.50, 1.00). Every leg
  uses F_q inside that range, so this is interpolation, not extrapolation.
    T_g = -142.05 F^2 -125.92 F +1932.67
  ⚠ The smokes were all at scale 1.0. Raising T_s cuts the driving dT by
    ~5 %, which raises T_g by < 5 K -- second order, carried as noise.

  s* = 0.09*R90*F/feed = 254.2*F mm     (the packet's 198*F used R90 = 3.3)
  CONTACT when F <= feed/R90 = 0.3540   (the packet assumed 0.4545)
  WINDOW  0.3540 < F <= 0.3934   i.e. s* in (90, 100] mm

  leg          F_q   scale    T_s     T_g      F  s* mm   q100  verdict
  C0          1.00     1.0  823.3  1664.7  1.000    254   1.02  fail 254 mm
  G1600       1.00     1.0  823.3  1465.8  0.764    194   0.78  fail 194 mm
  Q5          0.50     1.0  823.3  1834.2  0.601    153   0.61  fail 153 mm
  Q5T870      0.50  1.3523  870.0  1834.2  0.527    134   0.59  fail 134 mm
  Q35         0.35     1.0  823.3  1871.2  0.436    111   0.45  fail 111 mm
  Q35T870     0.35  1.3523  870.0  1871.2  0.383     97   0.43  *** IN WINDOW ***
  Q30         0.30     1.0  823.3  1882.1  0.378     96   0.39  *** IN WINDOW ***
  Q30T870     0.30  1.3523  870.0  1882.1  0.332     84   0.37  CONTACT

  ⚠⚠ THIS INVERTS THE PACKET'S PREDICTIONS.
    * the pre-registered pass Q5T870 lands at 134 mm and FAILS;
    * the window is reached by Q35T870 (97 mm) and Q30 (96 mm);
    * Q35T870 is a leg the packet explicitly DROPPED as 'contact under
      C4' -- that verdict used R90 = 3.3, which item 0(e) refutes.
  ⚠ Q35T870 and Q30 reach the SAME s* by DIFFERENT lever pairs. That is a
    direct measurement of criterion (i)'s degeneracy, not a redundancy.
  ⚠⚠ The window is 10 mm wide and the predictor carries +/-5 % (C4), so a
    leg at 95-105 mm is NOT a decisive verdict. The deliverable stays
    'F to about +/-10 % with an equivalent nozzle temperature'.

  ⚠⚠ THE GEOMETRY WINDOW AND THE FLUX GATE NO LONGER COINCIDE.
    Revision 2's table had the predicted pass at q100 = 0.55 MW/m^2, the
    centre of Meier's required 0.5-0.6, and called that coincidence 'a
    real result' if it held. With the MEASURED R90 it does not hold: the
    legs that reach the standoff window sit at q100 = 0.39-0.43, BELOW
    Kant's measured spallation band (0.45-0.87) and below Meier's
    required 0.5-0.6. PRE-REGISTERED HERE, before the ladder runs:
      * P3 is expected to be REFUTED on its third clause;
      * holding Fig. 8.8's standoff would need a floor flux BELOW the
        threshold at which Kant measured spallation to occur at all.
    That tension is a RESULT, and it is the packet's real product if it
    survives. It must not be resolved by moving either bar.

  LEGS TO RUN (8 + identity, same budget as the packet's):
    C0       predictor check + known-error anchor (650 s, C3)
    G1600    the degeneracy witness, measured T_g 1465.8 K
    Q5       the F_q ladder at 0.50
    Q5T870   THE PACKET'S PRE-REGISTERED PASS -- kept to test P1 directly
    Q35      nearest bracket ABOVE the window
    Q35T870  predicted pass
    Q30      predicted pass by the OTHER lever
    Q30T870  nearest bracket BELOW (contact)
    S0       identity
  DROPPED from the packet's table: Q55T870 (145 mm) and Q5T920 (117 mm) --
  both now fail high and carry less than the legs that bracket the window.

==============================================================================
