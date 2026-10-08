#!/usr/bin/env python3
"""D2ac item 0 -- on paper, hashed BEFORE the first leg.

Nothing here reads a D2ac output. Inputs are the packet's quoted sources, the
frozen SWc run, and the source tree. Run it, hash it, then launch.

  (a) the standoff ceiling, the domain, and nozzle_z0
  (b) the pad-contact test, specified offline
  (c) the D2l factor-of-four audit
  (d) the FORCED key drops -- a gap in the packet's leg spec, resolved by the
      code rather than by choice
"""
import numpy as np

# ---- quoted sources (ACTIVE_STEP.md "Sources", one row each) --------------
D_HOLE_TAB82 = 85.0        # mm, Tab. 8.2 TEXT -- never photo-derived
VOL_WATER = 3.42           # L per 0.5 m, water filling
MEIER_MEAN = 1.28          # m/h, 0.49 m / 1383 s
MEIER_OP = 1.5             # m/h, operator's stated advance between pauses
LZ = 0.5                   # m, block height (Fig. 8.8 caption)
D_BURNER = 80.0            # mm, burner outer diameter
FOOT_STANDOFF = 0.050      # m, VT5 drawing: nozzle sits 50 mm above foot tips
FOOT_RI, FOOT_RO = 0.028, 0.040   # m, the pad ring

# ---- the map's decay term, LIFTED VERBATIM from A.log's h_expr ------------
#   pow(0.09/max(0.09,s),1.0)   ->  flat for s <= 90 mm, then 1/s
S_KNEE = 0.090             # m
decay = lambda s: S_KNEE / max(S_KNEE, s)

# ---- the anchor, measured on frozen SWc (d2ab_gaph), NOT assumed ----------
# Between 100 and 150 mm burner descent the pit lead is CONSTANT at 38.0 mm, so
# the floor recession rate EQUALS the burner rate in steady state and the
# packet's burner-rate anchor is also the floor-rate anchor. Axial standoff is
# foot_standoff 50 + lead 38 = 88 mm, which is INSIDE the flat region
# (decay = 1.000), so no de-scaling is needed.
R90_FIT = 3.694            # m/h, SWc ROP fitted over the 100-150 mm stretch
R90_ENDPOINT = 3.578       # m/h, same stretch by endpoint differencing
S_AXIAL_SWC = 0.088        # m


def main():
    P = print
    P("=" * 78)
    P("D2ac ITEM 0 -- prescribed feed at Meier's kinematics")
    P("=" * 78)

    P("\n(a) THE STANDOFF CEILING, THE DOMAIN, AND nozzle_z0")
    P("-" * 78)
    P(f"  Fig. 8.8 caps the pit lead at ~50 mm, so the standoff ceiling is")
    P(f"  ~100 mm (foot standoff {FOOT_STANDOFF*1e3:.0f} mm + lead 50 mm).")
    P(f"  Block height {LZ} m  ->  geometry.prob_hi z = {LZ}")
    P(f"  dz must stay 2 mm AND cubic (side_face_flux enforces dx=dy=dz):")
    P(f"      amr.n_cell z = {LZ}/0.002 = {int(LZ/0.002)}   (was 200 at lz = 0.4)")
    P(f"      dx = 0.12/60 = {0.12/60*1e3:.1f} mm = dz  OK")
    P(f"  cost: {int(LZ/0.002)}/200 = {LZ/0.002/200:.2f}x the cells of the lz = 0.4 legs.")
    P("")
    P(f"  nozzle_z0: at t = 0 the block is flat at z = {LZ}, the feet rest on it,")
    P(f"  and the nozzle sits foot_standoff above the tips (VT5 drawing).")
    P(f"      nozzle_z0 = {LZ} + {FOOT_STANDOFF} = {LZ+FOOT_STANDOFF:.3f} m")
    P(f"  Initial standoff {FOOT_STANDOFF*1e3:.0f} mm is inside the flat region")
    P(f"  (decay = {decay(FOOT_STANDOFF):.3f}), so the legs start at full flux.")
    P(f"  ⚠ The dominant convention on disk (nozzle_z0 = lz + 0.020) is a")
    P(f"    DIFFERENT study's 20 mm standoff; Meier's tool is 50 mm.")
    P("")
    P("  EQUILIBRIUM STANDOFF -- where floor recession equals the feed:")
    P("      R(s) = R90 * decay(s) = R90 * 0.09/s   for s > 90 mm")
    P("      R(s) = feed   =>   s = 0.09 * R90 / feed")
    P(f"\n  {'feed [m/h]':<12}{'s (fit 3.694)':>16}{'s (endpt 3.578)':>18}")
    for feed in (MEIER_OP, 1.3, MEIER_MEAN):
        P(f"  {feed:<12.2f}{1e3*S_KNEE*R90_FIT/feed:>14.0f} mm"
          f"{1e3*S_KNEE*R90_ENDPOINT/feed:>16.0f} mm")
    P(f"\n  ⇒ the model settles at {1e3*S_KNEE*R90_ENDPOINT/MEIER_OP:.0f}-"
      f"{1e3*S_KNEE*R90_FIT/MEIER_OP:.0f} mm at {MEIER_OP} m/h (FM/FL/FH) and "
      f"{1e3*S_KNEE*R90_ENDPOINT/1.3:.0f}-{1e3*S_KNEE*R90_FIT/1.3:.0f} mm at 1.3 (FM13).")
    P(f"    Meier's ceiling is <= 100 mm. P1 (180-260 mm) is reproduced.")
    P("")
    P("  THE DELIVERABLE -- floor-flux excess at a 100 mm standoff:")
    P(f"      decay(100 mm) = {decay(0.100):.2f}")
    for R, lab in ((R90_FIT, "fit"), (R90_ENDPOINT, "endpoint")):
        m = R * decay(0.100)
        P(f"      R90 = {R:.3f} ({lab:8s}) -> floor {m:.2f} m/h"
          f"  =>  {m/MEIER_OP:.2f}x (vs {MEIER_OP}) .. {m/MEIER_MEAN:.2f}x (vs {MEIER_MEAN})")
    lo = R90_ENDPOINT * decay(0.100) / MEIER_OP
    hi = R90_FIT * decay(0.100) / MEIER_MEAN
    P(f"\n  >>> PRE-REGISTERED FLOOR-FLUX EXCESS: {lo:.2f} - {hi:.2f} x <<<")
    P(f"  The packet states 2.22-2.60x from the fit anchor alone; carrying the")
    P(f"  endpoint anchor widens the low end to {lo:.2f}. Reported, not replaced.")
    P(f"  ⚠ Both are INSIDE the packet's P2 band of 2.0-3.0x.")
    P("")
    P("  ⚠ Domain-floor check (the packet's own inertness witness): at the feed")
    P(f"    the burner descends {MEIER_OP}*500/3600 = "
      f"{MEIER_OP*500/3600*1e3:.0f} mm in 500 s; with a {1e3*S_KNEE*R90_FIT/MEIER_OP:.0f} mm")
    P(f"    standoff the floor sits {LZ*1e3 - MEIER_OP*500/3600*1e3 - 1e3*S_KNEE*R90_FIT/MEIER_OP:.0f} mm")
    P(f"    above the domain bottom. The calibrated 150-300 mm section is")
    P(f"    reached ({MEIER_OP*500/3600*1e3 + 1e3*S_KNEE*R90_FIT/MEIER_OP:.0f} mm of floor depth), so stop_time = 500 stands.")

    P("\n(b) THE PAD-CONTACT TEST, OFFLINE")
    P("-" * 78)
    P("  ⚠ nozzle_collision_radius CANNOT be the test. In A's argv it is not")
    P("    set at all, so it defaults to 0.0 and the guard can never fire; even")
    P("    at the packet's 0.002 it would need a live column with r <= 2 mm to")
    P("    reach the NOZZLE PLANE, 50 mm above the foot tips, and the axis is")
    P("    the fastest-receding column in every leg. Vacuous at any feed.")
    P("")
    P("  THE REAL TEST -- the feet resting on the pad ring (hook load, p. 221):")
    P(f"      foot tips     z_tip(t) = nozzle_z0 - nozzle_feed*t - {FOOT_STANDOFF}")
    P(f"      pad ring      {FOOT_RI*1e3:.0f} <= r <= {FOOT_RO*1e3:.0f} mm")
    P("      contact       max top-face height over that annulus >= z_tip")
    P("  Computed offline from the removal log via the SAME zface() the shape")
    P("  metric uses, per output time. Reported as: first contact time, total")
    P("  contact fraction, and the worst-case interference in mm.")
    P("  ⚠ It is a genuine test: D2x put the pad-ring top-face rate at ~2.1 m/h,")
    P("    so a 1.5 m/h feed is marginal, not obviously clear.")

    P("\n(c) THE D2l FACTOR-OF-FOUR AUDIT")
    P("-" * 78)
    P("  ⚠ ALREADY DISCHARGED -- in D2ab item 0, one packet ago. It is listed")
    P("    here as 'owed by five packets', but D2ab reproduced it and closed it.")
    P("  Evidence, from d2ab_gaph/PREFLIGHT.md (sha 9840fa2b...):")
    P("      d2l_scan/analyze.py vdot()   returns 4.0 * sum(...)  -> FULL hole")
    P("      tse_rock = 4.0*p_face/vmm3   -> the 4 is on BOTH sides, cancels")
    P("      tse_hhv  = P_HHV(38 kW full) / vmm3(full hole)        consistent")
    P("      plot :397 uses 4*p_face, labelled 'absorbed (full hole)'")
    P("      reproduced: M_f10 quarter jet_P_face 3204.2 W x 4 = 12.82 kW,")
    P("      matching d2l_scan/RESULTS.md:190 'absorbed, full hole' EXACTLY.")
    P("  VERDICT: D2l DID apply the factor. '17.3 % of HHV' and 'needs 9.1 %'")
    P("  stand as written. Nothing to restate. NOT re-derived here -- the rule")
    P("  is to report a frozen pre-flight, not to recompute it each packet.")

    P("\n(d) FORCED KEY DROPS -- a gap in the packet's leg spec")
    P("-" * 78)
    P("  The packet says the shared block 'replaces nozzle_descent=feet'. That")
    P("  is not sufficient: A's argv carries ten more feet-only keys, and the")
    P("  code REFUSES them under prescribed. Measured, not reasoned:")
    P("      swap nozzle_descent only  -> rc=6, ABORT at MMWSpalling.H:7593")
    P("         'surface_patch.jet_free_surface = 1 requires")
    P("          surface_patch.nozzle_descent = feet'")
    P("      + drop jet_fs_*           -> rc=6 (foot_body_clearance, :7596)")
    P("      + drop the whole block    -> rc=0, unused-key list EMPTY")
    P("")
    P("  FORCED DROP LIST (no choice -- the code will not run otherwise):")
    for k in ("jet_free_surface", "jet_fs_aspect", "jet_fs_exponent",
              "foot_r_inner", "foot_r_outer", "foot_standoff", "foot_rule",
              "foot_npads", "foot_pad_quantile", "foot_body_clearance"):
        P(f"      surface_patch.{k}")
    P("")
    P("  ⚠ CONFOUND, stated up front: FM therefore differs from A by MORE than")
    P("    the feed. Two consequences, both bounded:")
    P("    1. free-surface dilution is gone. Measured on the frozen legs, it")
    P("       barely fired: A jet_fs_cols mean 0.019 (max 5), jet_m_ratio max")
    P("       1.040; SWc exactly 0 and 1.000. So the drop is near-inert, and it")
    P("       is FORCED -- the model does not support it under a prescribed feed.")
    P("    2. foot_r_inner/outer go to their 0.0 defaults, so the D2aa region")
    P("       split (side_P_bore/pads/wall) bins EVERY charged face as 'wall'.")
    P("       The region columns are therefore MEANINGLESS on these legs and")
    P("       are not used; criterion (b) is computed offline from geometry,")
    P("       which is why the packet put it there.")
    P("    ⚠ foot_r_inner also silently sets nozzle_collision_radius, so with")
    P("      it dropped the collision guard is doubly inert -- see (b).")

    P("\n(e) THE LZ TRAP IN THE IMPORTED MACHINERY")
    P("-" * 78)
    P("  d2v_hmap/analyze.py:30 has LZ = 0.4 and d2x_areasplit/analyze.py:47")
    P("  has LZ = 0.4, both MODULE-LEVEL. zface() returns LZ - nv*dz, so")
    P("  importing them 'unchanged' against a 0.5 m block would misplace every")
    P("  depth by 100 mm and silently corrupt Ø(z), the pit and the volume.")
    P("  ⇒ analyze.py sets V.LZ and X.LZ per leg from that leg's own .done")
    P("    metadata before each call, and restores them after. The functions")
    P("    are untouched; only the constant they read is bound per leg.")
    P("")
    P("=" * 78)
    P("PRE-REGISTERED, for CRITERION.md")
    P("=" * 78)
    P(f"  nozzle_z0            {LZ+FOOT_STANDOFF:.3f} m")
    P(f"  lz / n_cell z        {LZ} m / {int(LZ/0.002)}")
    P(f"  equilibrium standoff {1e3*S_KNEE*R90_ENDPOINT/MEIER_OP:.0f}-"
      f"{1e3*S_KNEE*R90_FIT/MEIER_OP:.0f} mm at 1.5 m/h, "
      f"{1e3*S_KNEE*R90_ENDPOINT/1.3:.0f}-{1e3*S_KNEE*R90_FIT/1.3:.0f} at 1.3")
    P(f"  floor-flux excess    {lo:.2f} - {hi:.2f} x   <- THE DELIVERABLE")
    P(f"  Meier rate band      {MEIER_MEAN} - {MEIER_OP} m/h")
    P(f"  Ø target             {D_HOLE_TAB82} +/- 5 mm (Tab. 8.2 TEXT, not photo)")
    P(f"  volume ceiling       {VOL_WATER} L per 0.5 m")


if __name__ == "__main__":
    main()
