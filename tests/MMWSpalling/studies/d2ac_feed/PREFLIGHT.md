==============================================================================
D2ac ITEM 0 -- prescribed feed at Meier's kinematics
==============================================================================

(a) THE STANDOFF CEILING, THE DOMAIN, AND nozzle_z0
------------------------------------------------------------------------------
  Fig. 8.8 caps the pit lead at ~50 mm, so the standoff ceiling is
  ~100 mm (foot standoff 50 mm + lead 50 mm).
  Block height 0.5 m  ->  geometry.prob_hi z = 0.5
  dz must stay 2 mm AND cubic (side_face_flux enforces dx=dy=dz):
      amr.n_cell z = 0.5/0.002 = 250   (was 200 at lz = 0.4)
      dx = 0.12/60 = 2.0 mm = dz  OK
  cost: 250/200 = 1.25x the cells of the lz = 0.4 legs.

  nozzle_z0: at t = 0 the block is flat at z = 0.5, the feet rest on it,
  and the nozzle sits foot_standoff above the tips (VT5 drawing).
      nozzle_z0 = 0.5 + 0.05 = 0.550 m
  Initial standoff 50 mm is inside the flat region
  (decay = 1.000), so the legs start at full flux.
  ⚠ The dominant convention on disk (nozzle_z0 = lz + 0.020) is a
    DIFFERENT study's 20 mm standoff; Meier's tool is 50 mm.

  EQUILIBRIUM STANDOFF -- where floor recession equals the feed:
      R(s) = R90 * decay(s) = R90 * 0.09/s   for s > 90 mm
      R(s) = feed   =>   s = 0.09 * R90 / feed

  feed [m/h]     s (fit 3.694)   s (endpt 3.578)
  1.50                   222 mm             215 mm
  1.30                   256 mm             248 mm
  1.28                   260 mm             252 mm

  ⇒ the model settles at 215-222 mm at 1.5 m/h (FM/FL/FH) and 248-256 mm at 1.3 (FM13).
    Meier's ceiling is <= 100 mm. P1 (180-260 mm) is reproduced.

  THE DELIVERABLE -- floor-flux excess at a 100 mm standoff:
      decay(100 mm) = 0.90
      R90 = 3.694 (fit     ) -> floor 3.32 m/h  =>  2.22x (vs 1.5) .. 2.60x (vs 1.28)
      R90 = 3.578 (endpoint) -> floor 3.22 m/h  =>  2.15x (vs 1.5) .. 2.52x (vs 1.28)

  >>> PRE-REGISTERED FLOOR-FLUX EXCESS: 2.15 - 2.60 x <<<
  The packet states 2.22-2.60x from the fit anchor alone; carrying the
  endpoint anchor widens the low end to 2.15. Reported, not replaced.
  ⚠ Both are INSIDE the packet's P2 band of 2.0-3.0x.

  ⚠ Domain-floor check (the packet's own inertness witness): at the feed
    the burner descends 1.5*500/3600 = 208 mm in 500 s; with a 222 mm
    standoff the floor sits 70 mm
    above the domain bottom. The calibrated 150-300 mm section is
    reached (430 mm of floor depth), so stop_time = 500 stands.

(b) THE PAD-CONTACT TEST, OFFLINE
------------------------------------------------------------------------------
  ⚠ nozzle_collision_radius CANNOT be the test. In A's argv it is not
    set at all, so it defaults to 0.0 and the guard can never fire; even
    at the packet's 0.002 it would need a live column with r <= 2 mm to
    reach the NOZZLE PLANE, 50 mm above the foot tips, and the axis is
    the fastest-receding column in every leg. Vacuous at any feed.

  THE REAL TEST -- the feet resting on the pad ring (hook load, p. 221):
      foot tips     z_tip(t) = nozzle_z0 - nozzle_feed*t - 0.05
      pad ring      28 <= r <= 40 mm
      contact       max top-face height over that annulus >= z_tip
  Computed offline from the removal log via the SAME zface() the shape
  metric uses, per output time. Reported as: first contact time, total
  contact fraction, and the worst-case interference in mm.
  ⚠ It is a genuine test: D2x put the pad-ring top-face rate at ~2.1 m/h,
    so a 1.5 m/h feed is marginal, not obviously clear.

(c) THE D2l FACTOR-OF-FOUR AUDIT
------------------------------------------------------------------------------
  ⚠ ALREADY DISCHARGED -- in D2ab item 0, one packet ago. It is listed
    here as 'owed by five packets', but D2ab reproduced it and closed it.
  Evidence, from d2ab_gaph/PREFLIGHT.md (sha 9840fa2b...):
      d2l_scan/analyze.py vdot()   returns 4.0 * sum(...)  -> FULL hole
      tse_rock = 4.0*p_face/vmm3   -> the 4 is on BOTH sides, cancels
      tse_hhv  = P_HHV(38 kW full) / vmm3(full hole)        consistent
      plot :397 uses 4*p_face, labelled 'absorbed (full hole)'
      reproduced: M_f10 quarter jet_P_face 3204.2 W x 4 = 12.82 kW,
      matching d2l_scan/RESULTS.md:190 'absorbed, full hole' EXACTLY.
  VERDICT: D2l DID apply the factor. '17.3 % of HHV' and 'needs 9.1 %'
  stand as written. Nothing to restate. NOT re-derived here -- the rule
  is to report a frozen pre-flight, not to recompute it each packet.

(d) FORCED KEY DROPS -- a gap in the packet's leg spec
------------------------------------------------------------------------------
  The packet says the shared block 'replaces nozzle_descent=feet'. That
  is not sufficient: A's argv carries ten more feet-only keys, and the
  code REFUSES them under prescribed. Measured, not reasoned:
      swap nozzle_descent only  -> rc=6, ABORT at MMWSpalling.H:7593
         'surface_patch.jet_free_surface = 1 requires
          surface_patch.nozzle_descent = feet'
      + drop jet_fs_*           -> rc=6 (foot_body_clearance, :7596)
      + drop the whole block    -> rc=0, unused-key list EMPTY

  FORCED DROP LIST (no choice -- the code will not run otherwise):
      surface_patch.jet_free_surface
      surface_patch.jet_fs_aspect
      surface_patch.jet_fs_exponent
      surface_patch.foot_r_inner
      surface_patch.foot_r_outer
      surface_patch.foot_standoff
      surface_patch.foot_rule
      surface_patch.foot_npads
      surface_patch.foot_pad_quantile
      surface_patch.foot_body_clearance

  ⚠ CONFOUND, stated up front: FM therefore differs from A by MORE than
    the feed. Two consequences, both bounded:
    1. free-surface dilution is gone. Measured on the frozen legs, it
       barely fired: A jet_fs_cols mean 0.019 (max 5), jet_m_ratio max
       1.040; SWc exactly 0 and 1.000. So the drop is near-inert, and it
       is FORCED -- the model does not support it under a prescribed feed.
    2. foot_r_inner/outer go to their 0.0 defaults, so the D2aa region
       split (side_P_bore/pads/wall) bins EVERY charged face as 'wall'.
       The region columns are therefore MEANINGLESS on these legs and
       are not used; criterion (b) is computed offline from geometry,
       which is why the packet put it there.
    ⚠ foot_r_inner also silently sets nozzle_collision_radius, so with
      it dropped the collision guard is doubly inert -- see (b).

(e) THE LZ TRAP IN THE IMPORTED MACHINERY
------------------------------------------------------------------------------
  d2v_hmap/analyze.py:30 has LZ = 0.4 and d2x_areasplit/analyze.py:47
  has LZ = 0.4, both MODULE-LEVEL. zface() returns LZ - nv*dz, so
  importing them 'unchanged' against a 0.5 m block would misplace every
  depth by 100 mm and silently corrupt Ø(z), the pit and the volume.
  ⇒ analyze.py sets V.LZ and X.LZ per leg from that leg's own .done
    metadata before each call, and restores them after. The functions
    are untouched; only the constant they read is bound per leg.

==============================================================================
PRE-REGISTERED, for CRITERION.md
==============================================================================
  nozzle_z0            0.550 m
  lz / n_cell z        0.5 m / 250
  equilibrium standoff 215-222 mm at 1.5 m/h, 248-256 at 1.3
  floor-flux excess    2.15 - 2.60 x   <- THE DELIVERABLE
  Meier rate band      1.28 - 1.5 m/h
  Ø target             85.0 +/- 5 mm (Tab. 8.2 TEXT, not photo)
  volume ceiling       3.42 L per 0.5 m
