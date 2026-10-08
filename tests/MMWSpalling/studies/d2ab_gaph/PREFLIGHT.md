==============================================================================
D2ab ITEM 0 -- the gap-derived annulus coefficient
==============================================================================

0. WHICH DIAMETER SETS THE GAP
   Meier's MEASURED 85-93 mm, NOT the leg's own diameter. The leg's
   diameter is a model OUTPUT; using it would make h depend on the
   result it produces. Burner Ø80 mm = 2*foot_r_outer = the skirt that
   rests on the rock.
   gap = (D_hole - 80)/2 = 2.5 - 6.5 mm
   ⚠ The packet says 'a gap of 2.5-6 mm'; (93-80)/2 = 6.5, so the top of
     the packet's range is 6.5 not 6. Reported, not silently replaced.

1. HYDRAULIC DIAMETER, REYNOLDS NUMBER, REGIME
   D_h = 4A/P = 2(R_o - R_i) = D_o - D_i for a concentric annulus.

   mdot basis              D_hole     D_h         A   T_gas      Re  regime
   nozzle feed (15.13 g/s)     85     5.0      6.48    1464    2134  transitional
   nozzle feed (15.13 g/s)     93    13.0     17.66    1464    2035  transitional
   march, SW (20.87 g/s)      85     5.0      6.48    1464    2944  transitional
   march, SW (20.87 g/s)      93    13.0     17.66    1464    2808  transitional
   march, A (34.48 g/s)       85     5.0      6.48    1679    4396  turbulent
   march, A (34.48 g/s)       93    13.0     17.66    1679    4192  turbulent

   ⚠ Re = 2 mdot / (pi (R_o+R_i) mu) -- the hole diameter cancels almost
     exactly, so the 85-93 mm uncertainty moves Re by < 5 %. What moves
     Re is the MASS FLOW BASIS, by a factor of 2.3.

   WHICH MASS FLOW. At steady state only the nozzle feed crosses the
   exit plane: jet_T_ent_mode = exhaust means the march's entrained gas
   IS the hole's own recirculated exhaust, so counting it in the annulus
   double-counts mass that never leaves. PRIMARY = the nozzle feed.
   The marched value is carried as the upper sensitivity.

2. NUSSELT CORRELATIONS AND THEIR SOURCES
   laminar   Nu = 5.385   Incropera Table 8.2, annulus D_i/D_o -> 1,
                          outer wall heated, inner insulated. These run
                          D_i/D_o = 0.86-0.94, i.e. parallel plates.
   turbulent Nu = 0.023 Re^0.8 Pr^0.4   Incropera Eq. 8.60 (Dittus-
                          Boelter, heating), on the annulus D_h.
   ⚠ Re = 2035-2134 on the primary basis is TRANSITIONAL. The choice
     between these two IS a factor of two and is NOT derivable. The
     bracket is their UNION, exactly as D2aa ran von Rohr at both ends.
   ⚠ Entry length: L_e ~ 0.05 Re D_h = 0.5 m laminar, far longer than
     the 50-200 mm annulus, so the real flow is in the ENTRANCE region
     where Nu exceeds the fully-developed value. This biases the true h
     toward the UPPER end of the bracket. Stated, not corrected for:
     an entry-length correction would be a third free choice.

3. h_low AND h_high

   basis                   D_hole    h_lam   h_turb
   nozzle feed (15.13 g/s)     85     99.7    169.3
   nozzle feed (15.13 g/s)     93     38.4     62.7
   march, SW (20.87 g/s)      85     99.7    219.0
   march, SW (20.87 g/s)      93     38.4     81.1
   march, A (34.48 g/s)       85    111.9    337.6
   march, A (34.48 g/s)       93     43.1    125.0

   PRIMARY BRACKET (nozzle feed, Ø85-93, laminar U turbulent):
     h_low  =   38.4 W/m^2K
     h_high =  169.3 W/m^2K
   >>> GL uses side_face_h = 38, GH uses side_face_h = 169 <<<

   FULL union including the marched mass flow: 38 - 338
   ⚠ NO point value inside the bracket is taken. Both ends are run.

   SENSITIVITIES, reported and NOT folded into the bracket:
     Ø67 water-cooled shield instead of the Ø80 skirt, Ø85: D_h 18 mm, Re 2316, h 28-50
     Ø67 water-cooled shield instead of the Ø80 skirt, Ø93: D_h 26 mm, Re 2201, h 19-33
     film properties at (T_gas+T_wall)/2 = 1157 K: k 0.0763 vs 0.0926 -> h lower by 18 %

4. WHAT EACH END PREDICTS
==============================================================================
   TOP-FACE power is invariant to the side treatment:
     WL    P_top = 4.16 kW
     WL81  P_top = 4.11 kW
     SW    P_top = 3.96 kW
     mean 4.08 kW  <- used below

   BORE + PADS revert to the INHERITED map, so their power is SW's.
   ⚠ SW was frozen WITHOUT the region columns (the D2aa item-0(a) debt),
     so SW's own split is not on disk. Bracketed by WL's wall share:
       WL   wall = 468.6/1321.6 = 35.5 % of side power
       WL81 wall = 275.6/ 825.5 = 33.4 %
     so bore+pads = 60-70 % of SW's 4.09 kW side total.
       bore+pads = 2.45 - 2.86 kW, central 2.66

   WALL power scales with h, measured on the WL pair:
     h=142 -> 468.6 W  (3.30 W per W/m^2K)
     h= 81 -> 275.6 W  (3.40 W per W/m^2K)
     slope 3.35 W per (W/m^2K), used linearly

   RATE via the packet's rate ~ P^1.59, anchored on SW (8.05 kW,
   3.694 m/h). ⚠ On the frozen legs this anchor UNDER-predicts by
   3-4 % (D2aa: WL -3.3 %, WL81 -4.4 %), which is inside the 15 % bar.

   leg        h   P_wall   P_rock kW   rate m/h
   GL      38.4     129 W     6.86          2.87   [2.73 - 3.00]
   GH     169.3     567 W     7.30          3.16   [3.02 - 3.31]

   PRE-REGISTERED, criterion (a) scores RATE against these:
     GL: h = 38 W/m^2K -> rate 2.87 m/h  (15 % bar: 2.44 - 3.30)
     GH: h = 169 W/m^2K -> rate 3.16 m/h  (15 % bar: 2.69 - 3.64)

   ⚠ DISAGREEMENT WITH THE PACKET, reported not replaced: P2 pre-registers
     'GL/GH land 2.3-2.9 m/h'. This derivation gives 2.87-3.16, above that band. The
     difference is the bore+pads estimate: the packet's 'roughly 2.5 kW'
     side total implies P_rock 6.6 kW and 2.72 m/h with the same P^1.59,
     whereas SW's measured 4.09 kW side total at 60-70 % gives 2.45-2.86
     for bore+pads ALONE, before the wall is added. Criterion (a) is
     scored against ITEM 0 (this file), as the packet specifies; P2 is
     reported alongside.

   SHAPE, REPORTED not scored -- by interpolation in rate on the frozen
   legs (A 1.251, WL81 1.780, WL 2.073, SW 3.694):
     D@100 [mm]     GL   87.1   GH   87.3
     D@125 [mm]     GL   84.7   GH   85.1
     pit lead [mm]  GL   63.5   GH   54.3
     volume [cm3]   GL    1.2   GH    1.1
   ⚠ Both predicted rates are ABOVE WL's 2.073 and below SW's 3.694, so
     these are interpolations, not extrapolations.

==============================================================================
D2l FACTOR-OF-FOUR AUDIT -- the debt owed by three packets. DISCHARGED.
==============================================================================
   d2l_scan/analyze.py:
     vdot()      returns 4.0 * sum(h_applied) * dz^2 / dt   -> FULL hole
     tse_rock    = 4.0 * p_face / vmm3    -> 4 on BOTH sides, cancels
     tse_hhv     = P_HHV(38 kW, full burner) / vmm3(full hole)  consistent
     plot :397   4 * p_face, labelled 'absorbed (full hole)'
   RESULTS.md:190 'absorbed, full hole' f=1 -> 12.82 kW.
   Reproduced from d2l_scan/output/M_f10/thermo.dat over 150-250 s:
     jet_P_face (quarter) 3204.2 W  x4 = 12.82 kW   EXACT MATCH
   The f = 0 row (6.56 kW) is the linear-in-f extrapolation of
     f=0.1 7.35, f=0.2 7.92, f=0.5 10.10, f=1.0 12.82 kW (full hole).
   VERDICT: D2l DID apply the factor of four. '17.3 % of 38 kW HHV' and
   'needs 9.1 %' both stand as written. Nothing to restate. Debt closed.

