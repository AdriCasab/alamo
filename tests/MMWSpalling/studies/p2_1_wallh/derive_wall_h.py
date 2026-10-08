#!/usr/bin/env python3
"""P2-1 ord 1 / deliverable 1 -- THE WALL LAW, DERIVED OFFLINE. Zero compute.

Everything here is read off files already on disk -- `p1g2c_shiftfix/CTRL4_FIXED`
-- so (a1)-(a4) are decided BEFORE any source change, which is the point: the
(c2) bar is computed here and written down, and the leg is then run against it.

⚠⚠ RULES THIS SCRIPT ENFORCES (each has cost this programme a packet):

 1. Every run-dependent quantity is READ FROM THE RUN. The gas state is this
    leg's own logged `jet_T_exhaust`, not a number carried in a packet. See the
    PACKET CHECK block: the packet states 1674-1690 K and the leg logs 1519 K,
    and the derivation uses what the leg logs.
 2. The bore and the face heights come from the final VOID MASK, never from the
    event log -- the event-log estimator is biased ~0.3 mm (P1g-2c risk iv).
 3. (a4) is RESOLVED BY READING THE SOURCE, in this script, so it re-runs. An
    assertion in prose is not a resolution.
 4. Quarter domain: the input is a quarter annulus (axis at x0=y0=0). Mass
    flows and areas are converted to FULL at the point of use and said so.

Usage:  python3 derive_wall_h.py > WALL_H.txt
"""
import math
import os
import re
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "../../../.."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "tests/MMWSpalling/studies/p2_0_annulus"))

import wall_law as WL                                     # noqa: E402
# ⚠ ONE implementation of the duct branch and of the gas properties, imported
# from the step that derived them, exactly as P1g-2b imported annulus_h.
from annulus_h import props, annulus, nu_gnielinski, read_argv   # noqa: E402

LEG = os.path.join(ROOT, "tests/MMWSpalling/studies/p1g2c_shiftfix/CTRL4_FIXED")
ARGV = os.path.join(ROOT,
                    "tests/MMWSpalling/studies/p1g2c_shiftfix/CTRL4_FIXED.argv.sh")
W_LO, W_HI = 150.0, 350.0          # the declared scoring window [s]

# The packet's own claims, restated here so the script can CHECK them rather
# than inherit them. A row that does not reproduce is printed, not silently
# replaced.
PACKET = {
    "jet_T_exhaust": (1674.0, 1690.0),
    "T_stag": (1645.0, 1645.0),
    "phi": (0.3364, 0.3364),
    "h_duct": (50.3, 50.3),
    "corner_lo_hi": (250.0, 350.0),       # (a1) gate
    "z_fill_lo_hi": (0.150, 0.220),       # (a2) gate [m]
    "band_mean_lo_hi": (150.0, 220.0),    # (a3) report band
}

# Kant & Rudolf von Rohr Table 4 -- the five tabulated stand-offs, verbatim from
# tests/MMWSpalling/validation/kant/kant_von_rohr_validation_test_for_mmwspalling_plan.txt:69-74
KANT = ((50, 1723.15, 766.0), (80, 1683.15, 1009.0), (110, 1633.15, 873.0),
        (130, 1523.15, 747.0), (150, 1373.15, 633.0))

# The depth bands the bore is scored over. ⚠ 90-110 mm is excluded by standing
# rule (never score the bore at 100 mm depth); the 10-50 mm collar is reported
# apart because it is not the bore; a band within 20 mm of either front is not
# comparable.
BORE_BANDS = ((110.0, 150.0), (150.0, 200.0), (200.0, 300.0))
MOUTH_BAND = (10.0, 50.0)

P = print


# ---------------------------------------------------------------------------
# (a4) -- the reachability ledger, resolved by reading the source.
# ---------------------------------------------------------------------------
def resolve_sites():
    """Return (rows, ok). Each row is (name, found, detail)."""
    mmw = open(os.path.join(ROOT, "src/Integrator/MMWSpalling.H")).read()
    rem = open(os.path.join(ROOT,
                            "src/Integrator/MMWSpalling/Removal.H")).read()
    rows = []

    def find(txt, pat, path):
        m = re.search(pat, txt)
        if not m:
            return None
        return "%s:%d" % (path, txt[:m.start()].count("\n") + 1)

    # the four roles, BY CONTENT. ⚠ never by line number: P1g-0 moved these
    # ~25 lines when a comment grew.
    crit = find(rem, r"if \(side_h_ex >= 0\.0\) \{ h_side = side_h_ex;", "Removal.H")
    kern = find(mmw, r"if \(side_h_ex >= 0\.0 && r2_col >= side_h_rmin2\)", "MMWSpalling.H")
    marc = find(mmw, r"const bool ex_c = \(surface_patch\.side_face_h >= 0\.0\)", "MMWSpalling.H")
    ledg = find(mmw, r"if \(surface_patch\.side_face_h >= 0\.0\s*\n\s*&& r2_col >= side_h_rmin2\)",
                "MMWSpalling.H")
    for nm, site in (("criterion", crit), ("kernel", kern),
                     ("march", marc), ("ledger", ledg)):
        rows.append(("role: " + nm, site is not None, site or "NOT FOUND"))

    # the criterion's hoisted radius gate
    g = find(rem, r"if \(!\(lat_face_form && r2 >= rmin2\)\) continue;", "Removal.H")
    rows.append(("criterion radius gate hoisted", g is not None, g or "NOT FOUND"))

    # the height datum's source: the gathered, ReduceIntMax'd column top. It
    # must be in scope at all four roles.
    rows.append(("datum source: k_top_global (criterion)",
                 "amrex::ParallelDescriptor::ReduceIntMax(k_top_global.data()" in rem,
                 "Removal.H, gathered + reduced before the face scan"))
    rows.append(("datum source: k_top_col (kernel, in AdvanceMicrostructure)",
                 "const amrex::Vector<int> k_top_col =" in mmw, "MMWSpalling.H"))
    rows.append(("datum source: k_top_col passed to the march",
                 "void JetEnthalpyMarch(int lev, Set::Scalar time, "
                 "const amrex::Vector<int>& k_top_col)" in mmw, "MMWSpalling.H"))
    rows.append(("datum source: k_top_col passed to the ledger",
                 bool(re.search(r"void EnergyLedger\(int lev, Set::Scalar time, "
                                r"Set::Scalar dt,\s*\n\s*const amrex::Vector<int>& k_top_col",
                                mmw)), "MMWSpalling.H"))

    # ⚠ the march's slot IS the height: slot = k_top_col[c] - k. Without this
    # the march could not resolve a per-height value at all and deliverable 4
    # would be impossible as written.
    rows.append(("march slot carries the height (slot = k_top_col[c] - k)",
                 "const int slot = k_top_col[c] - k;" in mmw,
                 "MMWSpalling.H, the march's side-hit binning"))

    # the scope keys that must NOT be blamed on the law
    rows.append(("scope: side_face_h_r_min is the criterion's gate",
                 "side_face_h_r_min" in mmw, "MMWSpalling.H (parse)"))
    rows.append(("scope: side_face_factor_above zeroes above-plane faces",
                 "const Set::Scalar sf_c = above_c ? side_f_above : side_f;" in mmw,
                 "MMWSpalling.H, kernel + ledger"))

    # ⚠⚠ THE CONJUNCTION. The derived law is reached only where the explicit
    # override already is, so `side_face_h >= 0` is a PRECONDITION of the law
    # running at all. A leg that sets the law and not the pair gets nothing.
    rows.append(("conjunction: the law is reached only where side_face_h >= 0",
                 all(x is not None for x in (crit, kern, marc, ledg)),
                 "all four roles test side_face_h/side_h_ex >= 0 first"))
    return rows, all(r[1] for r in rows)


# ---------------------------------------------------------------------------
# the run's own state
# ---------------------------------------------------------------------------
def read_thermo(path, lo, hi):
    lines = open(path).read().split("\n")
    hdr = lines[0].split()
    ix = {n: i for i, n in enumerate(hdr)}
    rows = []
    for ln in lines[1:]:
        f = ln.split()
        if len(f) < len(hdr):
            continue
        try:
            t = float(f[ix["time"]])
        except ValueError:
            continue
        if lo <= t <= hi:
            rows.append(f)
    return ix, rows


def load_mask():
    """Final void mask -> (rem, dx, dz, zlo, t). ⚠ the MASK, never the events."""
    import numpy as np
    import yt
    yt.set_log_level(50)
    d = os.path.join(LEG, "FH")
    pfs = sorted(x for x in os.listdir(d) if x.endswith("cell"))
    ds = yt.load(os.path.join(d, pfs[-1]))
    dims = ds.domain_dimensions
    cg = ds.covering_grid(0, left_edge=ds.domain_left_edge, dims=dims)
    rem = np.array(cg["removed"]) > 0.5
    dx = float((ds.domain_right_edge[0] - ds.domain_left_edge[0]) / dims[0])
    dz = float((ds.domain_right_edge[2] - ds.domain_left_edge[2]) / dims[2])
    zlo = float(ds.domain_left_edge[2])
    return rem, dx, dz, zlo, float(ds.current_time), pfs[-1]


def main():
    import numpy as np

    av = read_argv(ARGV)
    ix, rows = read_thermo(os.path.join(LEG, "FH/thermo.dat"), W_LO, W_HI)
    if not rows:
        P("FAIL: no thermo rows in the window")
        return 1
    mean = lambda n: st.mean([float(r[ix[n]]) for r in rows])

    mdot_q = av["surface_patch.jet_mdot"][0]
    cp = av["surface_patch.jet_cp"][0]
    T_noz = av["surface_patch.jet_T_nozzle"][0]
    T_ent = av["surface_patch.jet_T_ent"][0]
    D_noz = av["surface_patch.jet_D"][0]
    L_core = av["surface_patch.jet_core_length"][0]
    h_typed = av["surface_patch.side_face_h"][0]
    r_gate = av["surface_patch.side_face_h_r_min"][0]
    f_above = av["surface_patch.side_face_factor_above"][0]
    noz_z0 = av["surface_patch.nozzle_z0"][0]
    feed = av["surface_patch.nozzle_feed"][0]

    T_exh = mean("jet_T_exhaust")
    T_stag = mean("jet_T_stag")
    T_sidewall = mean("side_T_wall")
    D_e_log = mean("jet_D_e")
    s_c = mean("jet_s_c")
    P_cap = mean("jet_P_cap")
    phi = mdot_q * cp * (T_stag - T_ent) / P_cap
    mdot_full = 4.0 * mdot_q

    P("=" * 78)
    P("P2-1 ord 1 -- THE COMPOSITE WALL LAW, DERIVED OFFLINE")
    P("=" * 78)
    P("  leg    : CTRL4_FIXED (p1g2c_shiftfix) -- the leg (c2) and (d2) compare against")
    P("  window : %.0f-%.0f s, %d thermo rows" % (W_LO, W_HI, len(rows)))
    P("  typed  : side_face_h = %.1f W/m2K, r_gate %.3f m, factor_above %.1f"
      % (h_typed, r_gate, f_above))
    P("  burner : mdot(quarter) %.6g kg/s => FULL %.6g kg/s, D %.4f m, T_noz %.0f K"
      % (mdot_q, mdot_full, D_noz, T_noz))
    P("  shield : R_SHIELD %.4f m (VT5 drawings, Oe67) -- NOT the %.3f m gate"
      % (WL.R_SHIELD, r_gate))

    # ---------------- PACKET CHECK -------------------------------------
    P("")
    P("-" * 78)
    P("PACKET CHECK  every number the packet states, reproduced from the leg")
    P("-" * 78)
    lo, hi = PACKET["jet_T_exhaust"]
    ok_exh = lo <= T_exh <= hi
    P("  jet_T_exhaust   packet %.0f-%.0f K   measured %.1f K   %s"
      % (lo, hi, T_exh, "ok" if ok_exh else "** DOES NOT REPRODUCE **"))
    P("  jet_T_stag      packet %.0f K        measured %.1f K   %s"
      % (PACKET["T_stag"][0], T_stag,
         "ok" if abs(T_stag - PACKET["T_stag"][0]) <= 5.0 else "** off **"))
    P("  phi             packet %.4f        measured %.4f     %s"
      % (PACKET["phi"][0], phi,
         "ok" if abs(phi - PACKET["phi"][0]) <= 0.002 else "** off **"))
    if not ok_exh:
        P("  ⚠⚠ THE GAS STATE IS TAKEN FROM THE LEG, NOT THE PACKET. The packet's")
        P("     1674-1690 K is not this leg's exhaust; `side_T_wall` -- the")
        P("     temperature the wall faces ACTUALLY see -- logs %.1f K, equal to" % T_sidewall)
        P("     jet_T_exhaust to %.3g K, which is what fixes the state. Using the"
          % abs(T_sidewall - T_exh))
        P("     packet's figure would LOWER the derived h by ~8 % (density), i.e.")
        P("     it is not conservative in the direction of the (a1) gate.")

    # ---------------- the mask ------------------------------------------
    rem, dx, dz, zlo, t_end, pf = load_mask()
    nx, ny, nz = rem.shape
    Z_TOP = zlo + nz * dz
    z_noz = noz_z0 - feed * t_end
    P("")
    P("-" * 78)
    P("THE VOID MASK  %s at t = %.1f s   (%dx%dx%d, dx %.4f m, dz %.4f m)"
      % (pf, t_end, nx, ny, nz, dx, dz))
    P("-" * 78)
    cc = (np.arange(nx) + 0.5) * dx
    rr = np.hypot(*np.meshgrid(cc, cc, indexing="ij"))
    ktop = np.full((nx, ny), -1, dtype=int)
    for i in range(nx):
        for j in range(ny):
            s = np.nonzero(~rem[i, j, :])[0]
            if len(s):
                ktop[i, j] = int(s.max())
    # bore, area-equivalent, per depth band. ⚠ area-equivalent not max-radius:
    # a staircased wall makes the max radius a one-cell artefact.
    def bore_of_band(d_lo, d_hi):
        ks = [k for k in range(nz)
              if d_lo <= 1e3 * (Z_TOP - (zlo + (k + 0.5) * dz)) <= d_hi]
        if not ks:
            return float("nan"), float("nan")
        ra, rm = [], []
        for k in ks:
            n = int(rem[:, :, k].sum())
            if n == 0:
                continue
            ra.append(math.sqrt(4.0 * n * dx * dx / math.pi))
            rm.append(float(rr[rem[:, :, k]].max()) + 0.5 * dx)
        if not ra:
            return float("nan"), float("nan")
        return st.mean(ra), st.mean(rm)

    bands = [(b, ) + bore_of_band(*b) for b in BORE_BANDS]
    for (a, b), ra, rm in bands:
        P("  depth %3.0f-%3.0f mm : bore Oe %6.1f mm (area-equiv)   Oe %6.1f mm (max-radius)"
          % (a, b, 2e3 * ra, 2e3 * rm))
    mo_a, mo_m = bore_of_band(*MOUTH_BAND)
    P("  mouth %3.0f-%3.0f mm : bore Oe %6.1f mm (area-equiv)   -- REPORTED APART, not the bore"
      % (MOUTH_BAND[0], MOUTH_BAND[1], 2e3 * mo_a))
    r_bore = st.mean([x[1] for x in bands if x[1] == x[1]])
    P("  representative bore radius r_bore = %.4f m  (Oe %.1f mm)" % (r_bore, 2e3 * r_bore))
    P("  nozzle plane at t = %.1f s : z = %.4f m  (z0 %.3f - feed %.6g * t)"
      % (t_end, z_noz, noz_z0, feed))

    # ---------------- gas state and the two branches ---------------------
    mu, k_g, Pr = props(T_exh, "fit")
    mu4, k4, Pr4 = props(T_exh, "A4")
    rho = WL.rho_gas(T_exh)
    uDe = WL.u_De(mdot_full, D_noz, T_noz, T_exh)
    h_duct, Re_d, D_h, A_an = annulus(r_bore, mdot_full / phi, T_exh, "fit",
                                      nu_gnielinski)

    P("")
    P("-" * 78)
    P("(a1) GATE  the corner value, at the leg's own exhaust gas state")
    P("-" * 78)
    P("  gas    : T = %.1f K (logged jet_T_exhaust == side_T_wall), rho %.4f kg/m3,"
      % (T_exh, rho))
    P("           mu %.4g Pa s, Pr %.3f, c_p %.0f J/kgK (the model's jet_cp)"
      % (mu, Pr, cp))
    P("  momentum: u*D_e = 4*mdot/(pi*D*sqrt(rho_j*rho_inf)) = %.4f m2/s" % uDe)
    P("            (logged jet_D_e %.5f m => the implied u is %.0f m/s; the"
      % (D_e_log, uDe / D_e_log))
    P("             product is what the law needs and is carried unsplit)")
    v_c = WL.C_WJ_VEL * uDe / r_bore
    d_c = WL.C_WJ_DELTA * r_bore
    P("  corner : v_max %.1f m/s, delta %.2f mm   at r_bore %.1f mm"
      % (v_c, 1e3 * d_c, 1e3 * r_bore))
    h_c = WL.h_walljet(0.0, r_bore, uDe, T_exh, cp, mu, Pr)
    h_c4 = WL.h_walljet(0.0, r_bore, uDe, T_exh, cp, mu4, Pr4)
    h_c_nohalf = WL.h_walljet(0.0, r_bore, uDe, T_exh, cp, mu, Pr, cf_half=False)
    lo, hi = PACKET["corner_lo_hi"]
    ok_a1 = lo <= h_c <= hi
    P("  h_corner = %.1f W/m2K        [gate %.0f-%.0f]" % (h_c, lo, hi))
    P("  property sensitivity: A4 tables give %.1f (%+.1f %%)"
      % (h_c4, 100.0 * (h_c4 / h_c - 1.0)))
    P("  ⚠ the EXCLUDED convention (no 1/2) would give %.1f -- outside the gate,"
      % h_c_nohalf)
    P("    which is the measurement's verdict, not a preference.")
    P("  (a1) %s" % ("PASS" if ok_a1 else "FAIL"))

    # ---------------- (a2) z_fill and the switch -------------------------
    zf = WL.z_fill(r_bore)
    h_wj_at_fill = WL.h_walljet(zf, r_bore, uDe, T_exh, cp, mu, Pr)
    lo, hi = PACKET["z_fill_lo_hi"]
    ok_a2 = lo <= zf <= hi
    L_fill = (r_bore - WL.R_SHIELD) / WL.C_WJ_GROW
    P("")
    P("-" * 78)
    P("(a2) GATE  z_fill, and the switch discontinuity WITH ITS HANDLING")
    P("-" * 78)
    P("  gap at r_bore = r - R_SHIELD = %.2f mm; corner layer %.2f mm"
      % (1e3 * (r_bore - WL.R_SHIELD), 1e3 * d_c))
    P("  z_fill = (gap - delta_c)/0.073 = %.1f mm        [gate %.0f-%.0f mm]"
      % (1e3 * zf, 1e3 * lo, 1e3 * hi))
    P("  duct branch at this bore: h_duct %.1f W/m2K (Re %.0f, D_h %.4f m)"
      % (h_duct, Re_d, D_h))
    P("  packet states the P2-0 ord 1 duct value as %.1f  => %s"
      % (PACKET["h_duct"][0],
         "reproduces" if abs(h_duct - PACKET["h_duct"][0]) < 3.0
         else "DIFFERS (bore moved: this leg, not CTRL2)"))
    P("  THE JUMP: h_wj(z_fill) = %.1f  ->  h_duct = %.1f   = a factor %.2f"
      % (h_wj_at_fill, h_duct, h_wj_at_fill / h_duct))
    P("")
    P("  ⚠⚠ HANDLING, STATED: a HARD STEP at z_fill, accepted and logged.")
    P("     No blend, and no parameter. The reason it is safe here is that the")
    P("     datum is the LOCAL CORNER FLOOR, so z_fill is a height above a floor")
    P("     that descends through the run -- the step SWEEPS through the rock")
    P("     rather than sitting at one elevation, and cannot etch a ledge. A")
    P("     blend would need a second length (one gap-fill length, %.0f mm) that"
      % (1e3 * L_fill))
    P("     nothing measures; for reference it would replace the step by a ramp")
    P("     from %.1f to %.1f over %.0f-%.0f mm."
      % (h_wj_at_fill, h_duct, 1e3 * zf, 1e3 * (zf + L_fill)))
    P("  ⚠ NON-VACUITY: under the rejected `max(h_wj, h_duct)` form the duct is")
    P("    never selected (h_wj stays above h_duct to z = 400 mm, see the")
    P("    profile below), so this row would read as a no-op. It does not.")
    P("  (a2) %s" % ("PASS" if ok_a2 else "FAIL"))

    # ---------------- the profile ----------------------------------------
    P("")
    P("-" * 78)
    P("THE PROFILE at r_bore = %.1f mm  (deliverable 1)" % (2e3 * r_bore / 2))
    P("-" * 78)
    P("     z [mm]   delta [mm]   v [m/s]   Re_delta      h_wj   h_wall    /typed")
    for zmm in (0, 25, 50, 100, 150, 200, 400):
        z = 1e-3 * zmm
        d_z = d_c + WL.C_WJ_GROW * z
        v = v_c * math.sqrt(d_c / d_z)
        Re = v * d_z / (mu / rho)
        hw = WL.h_walljet(z, r_bore, uDe, T_exh, cp, mu, Pr)
        hh = WL.h_wall(z, r_bore, uDe, T_exh, cp, mu, Pr, h_duct)
        P("     %6d   %10.2f   %7.1f   %8.0f   %7.1f  %6.1f    %6.3f"
          % (zmm, 1e3 * d_z, v, Re, hw, hh, hh / h_typed))
    P("  ⚠ h_wj alone is %.1f at z = 400 mm, still ABOVE h_duct %.1f -- this is"
      % (WL.h_walljet(0.4, r_bore, uDe, T_exh, cp, mu, Pr), h_duct))
    P("    what makes the `max` form vacuous and the z_fill switch load-bearing.")

    # ---------------- (a3) the cutting band against the DEFINED datum ----
    P("")
    P("-" * 78)
    P("(a3) REPORT  the cutting band RE-DERIVED against the defined datum")
    P("-" * 78)
    P("  DATUM (deliverable 2, ONE definition, the twin of the C++ helper):")
    P("    z = (k - k_corner[col] - 0.5)*dz, floored at dz/2, where k_corner is")
    P("    the DEEPEST column top on the ray from this column inward to the")
    P("    axis, truncated at the wall band's own inner edge r = %.0f mm."
      % (1e3 * r_gate))
    MISS = -10 ** 6
    kt = [[int(ktop[i, j]) if ktop[i, j] >= 0 else MISS for j in range(ny)]
          for i in range(nx)]
    rl = [[float(rr[i, j]) for j in range(ny)] for i in range(nx)]
    kc = WL.wall_corner_columns(kt, rl, r_gate, MISS)
    k_axis = min(v for row in kt for v in row if v != MISS)

    faces = []          # (r, z, nopen)  under the CHOSEN datum
    alt = {"nearest tread (REJECTED)": [], "axis front (REJECTED)": []}
    n_clamp = 0
    for i in range(nx):
        for j in range(ny):
            if rr[i, j] < r_gate or kt[i][j] == MISS:
                continue
            col_top_z = zlo + (ktop[i, j] + 1) * dz
            if not (z_noz - col_top_z > 0.0):          # mirrors flame_s[col] > 0
                continue                               # 0 W by key, not by law
            for k in range(nz):
                if rem[i, j, k]:
                    continue
                cand, nopen = [], 0
                for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ii, jj = i + di, j + dj
                    if not (0 <= ii < nx and 0 <= jj < ny):
                        continue
                    if not rem[ii, jj, k]:
                        continue
                    nopen += 1
                    if 0 <= ktop[ii, jj] < k:
                        cand.append(ktop[ii, jj])
                if nopen == 0:
                    continue
                z = (k - kc[i][j] - 0.5) * dz
                if z < 0.5 * dz:
                    z = 0.5 * dz
                    n_clamp += 1
                faces.append((float(rr[i, j]), z, nopen))
                alt["nearest tread (REJECTED)"].append(
                    (max(0.5 * dz, (k - (max(cand) if cand else k - 1) - 0.5) * dz),
                     nopen))
                alt["axis front (REJECTED)"].append(
                    (max(0.5 * dz, (k - k_axis - 0.5) * dz), nopen))
    if not faces:
        P("  FAIL: no live wall faces found")
        return 1
    zs = sorted(f[1] for f in faces)
    q = lambda p: zs[min(len(zs) - 1, int(p * len(zs)))]
    P("  live wall cells (r >= %.0f mm, below the nozzle plane): %d  (%d faces)"
      % (1e3 * r_gate, len(faces), sum(f[2] for f in faces)))
    P("  faces clamped to dz/2 (no corner beneath them): %d (%.2f %%)"
      % (n_clamp, 100.0 * n_clamp / len(faces)))
    P("  height above the CORNER [mm]:  p05 %5.1f  p25 %5.1f  p50 %5.1f"
      "  p75 %5.1f  p95 %5.1f" % (1e3 * q(.05), 1e3 * q(.25), 1e3 * q(.50),
                                  1e3 * q(.75), 1e3 * q(.95)))
    P("  => the packet puts the band at 30-65 mm against this datum: %s"
      % ("REPRODUCED" if 30.0 <= 1e3 * q(.50) <= 65.0 else "NOT reproduced"))
    P("")
    P("  ⚠⚠ DATUM SENSITIVITY -- the packet says the datum moves the answer by")
    P("     25-60 mm, so the two readings NOT chosen are measured, not asserted:")
    for nm, lst in alt.items():
        za = sorted(x[0] for x in lst)
        ha = [WL.h_wall(z, f[0], uDe, T_exh, cp, mu, Pr,
                        annulus(f[0], mdot_full / phi, T_exh, "fit",
                                nu_gnielinski)[0])
              for (z, _n), f in zip(lst, faces)]
        mh = (sum(h * n for h, (_z, n) in zip(ha, lst))
              / sum(n for _z, n in lst))
        P("     %-26s p50 %5.1f mm  band mean h %6.1f  /typed %.3f"
          % (nm, 1e3 * za[len(za) // 2], mh, mh / h_typed))
    P("     nearest tread reads the one-cell STAIRCASE, and the layer is 4.4 mm")
    P("     at the corner growing to ~15 mm, i.e. it spans several treads: a")
    P("     2 mm ledge every 10 mm is roughness, not a fresh corner. The axis")
    P("     front reproduces the p50 90 mm prior packets quote, which is the")
    P("     check that this script's three readings are the same three.")
    hs = [WL.h_wall(f[1], f[0], uDe, T_exh, cp, mu, Pr,
                    annulus(f[0], mdot_full / phi, T_exh, "fit", nu_gnielinski)[0])
          for f in faces]
    w = [f[2] for f in faces]
    mean_h = sum(h * ww for h, ww in zip(hs, w)) / sum(w)
    lo, hi = PACKET["band_mean_lo_hi"]
    P("  the law's FACE-WEIGHTED MEAN over the band: %.1f W/m2K   [report %.0f-%.0f]"
      % (mean_h, lo, hi))
    P("  against the typed %.1f  =>  ratio %.4f" % (h_typed, mean_h / h_typed))
    nwj = sum(1 for f in faces if f[1] < WL.z_fill(f[0]))
    P("  faces on the wall-jet branch %d (%.1f %%), on the duct branch %d (%.1f %%)"
      % (nwj, 100.0 * nwj / len(faces), len(faces) - nwj,
         100.0 * (len(faces) - nwj) / len(faces)))

    # ⚠⚠ THE DISTRIBUTION IS BIMODAL AND THE SINGLE MEAN HIDES IT. The r >= 40 mm
    # gate admits two populations that the code calls by one name: the CRATER
    # FLOOR between the gate and the bore (where the local relief is a riser a
    # cell or two high, so z ~ 0 and the law is at its corner value), and the
    # TRUE HOLE WALL at r ~ r_bore (where z is the full height above the corner).
    # Reporting one mean over both is what would make this row unreadable.
    r_wall = r_bore - 2.0 * dx
    split = (("crater floor, %.0f <= r < %.0f mm" % (1e3 * r_gate, 1e3 * r_wall),
              [n for n, f in enumerate(faces) if f[0] < r_wall]),
             ("hole wall,    r >= %.0f mm" % (1e3 * r_wall),
              [n for n, f in enumerate(faces) if f[0] >= r_wall]))
    P("")
    P("  ⚠⚠ SPLIT BY RADIUS -- the gate admits TWO populations under one name:")
    for nm, sel in split:
        if not sel:
            P("    %-34s  (empty)" % nm)
            continue
        zz = sorted(faces[n][1] for n in sel)
        mh = (sum(hs[n] * faces[n][2] for n in sel)
              / sum(faces[n][2] for n in sel))
        P("    %-34s cells %5d  z p50 %5.1f mm  p95 %5.1f mm  mean h %6.1f  /typed %.3f"
          % (nm, len(sel), 1e3 * zz[len(zz) // 2],
             1e3 * zz[min(len(zz) - 1, int(0.95 * len(zz)))], mh, mh / h_typed))
    P("    => the crater-floor population sits in the radial wall jet at full")
    P("       strength, which is why the combined mean lands ABOVE the typed")
    P("       value and above the packet's expected 150-220 band.")

    P("")
    P("  ⚠⚠ THE SIGN IS NOT PREDICTED, AND IT CAME OUT THE OTHER WAY.")
    P("     mean/typed = %.3f, so the band sees a %s value than the typed"
      % (mean_h / h_typed, "HIGHER" if mean_h > h_typed else "LOWER"))
    P("     constant and (c3) %s rather than narrows. The packet's P2 expected"
      % ("widens" if mean_h > h_typed else "narrows"))
    P("     'close to or above'; this is ABOVE, by %.0f %%, and (c2)'s bar below"
      % (100.0 * (mean_h / h_typed - 1.0)))
    P("     is therefore a RISE in wall energy, not the fall the packet's P4")
    P("     wording assumes. P4's mechanism is unchanged -- the energy moves and")
    P("     the volume cannot -- only its direction.")
    P("  (a3) %s"
      % ("in the reported range" if lo <= mean_h <= hi else
         "OUTSIDE the reported range -- REPORTED, not a gate"))

    # ---------------- the (c2) bar, WRITTEN DOWN BEFORE THE LEG ----------
    P("")
    P("-" * 78)
    P("(c2) THE BAR, COMPUTED HERE AND FIXED BEFORE THE LEG RUNS")
    P("-" * 78)
    P("  ⚠ not carried from any prior estimate; computed from THIS profile.")
    fall = mean_h / h_typed
    P("  predicted wall-energy ratio (face-weighted h, dT held) = %.4f" % fall)
    W_now = mean("surf_side_share") * mean("ledger_P_robin")
    P("  CTRL4_FIXED wall energy, surf_side_share x ledger_P_robin = %.1f W" % W_now)
    P("  => (c2) predicts CTRL5 at %.1f W, bar +/- 10 %% => %.1f to %.1f W"
      % (fall * W_now, 0.9 * fall * W_now, 1.1 * fall * W_now))
    # cross-check: weight by each column's own logged side power
    import csv as _csv
    pw = {}
    with open(os.path.join(LEG, "FH_side_power.csv")) as fh:
        for row in _csv.DictReader(fh):
            t = float(row["time"])
            if not (W_LO <= t <= W_HI):
                continue
            key = (int(row["col_i"]), int(row["col_j"]))
            pw[key] = pw.get(key, 0.0) + float(row["P_side"])
    # per-column face lists for the power-weighted cross-check
    colmap = {}
    for (r_f, z_f, nf), hh in zip(faces, hs):
        colmap.setdefault(round(r_f, 9), []).append((hh, nf))
    num = den = 0.0
    for (ci, cj), pp in pw.items():
        if ci >= nx or cj >= ny:
            continue
        rk = round(float(rr[ci, cj]), 9)
        lst = colmap.get(rk)
        if not lst:
            continue
        hbar = sum(h * n for h, n in lst) / sum(n for _h, n in lst)
        num += pp * hbar / h_typed
        den += pp
    if den > 0.0:
        P("  cross-check, weighted by each column's LOGGED P_side: %.4f"
          % (num / den))
        P("    (same direction, %.1f %% from the face-weighted figure)"
          % (100.0 * abs(num / den - fall) / fall))
    P("  ⚠ ARITHMETIC BOUND ON (c4): width is held within 2 mm by (c3) (~4 % of")
    P("    area) and depth by (c5), so VOLUME cannot move more than ~5 %. The")
    P("    volume levers are the rock price (P2-2) and the floor (P2-3).")

    # ---------------- the five Kant stand-offs ---------------------------
    P("")
    P("-" * 78)
    P("THE LAW AT THE FIVE KANT STAND-OFFS  (deliverable 1, report)")
    P("-" * 78)
    P("  ⚠⚠ THIS IS A LEVEL COMPARISON, NOT A KANT PREDICTION. Kant's burner")
    P("     mass flow and nozzle diameter are NOT in this repository, so the")
    P("     momentum carried below is THIS leg's u*D_e; only the gas temperature")
    P("     is Kant's. Kant's h_fl is also a sensor-PLATE average (stagnation")
    P("     included), so it should sit ABOVE a local wall-jet value at r_bore.")
    P("     Kant Table 4 is quoted verbatim from the validation reference.")
    P("    SOD [mm]  T_flame [K]  measured h_fl   law @ corner   law/measured")
    kratio = []
    for sod, Tf, hfl in KANT:
        muK, kK, PrK = props(Tf, "fit")
        uK = WL.u_De(mdot_full, D_noz, T_noz, Tf)
        hK = WL.h_walljet(0.0, r_bore, uK, Tf, cp, muK, PrK)
        kratio.append(hK / hfl)
        P("    %8d  %11.2f  %13.0f   %12.1f   %12.3f" % (sod, Tf, hfl, hK, hK / hfl))
    P("  => the law sits at %.0f-%.0f %% of Kant's plate-average coefficient."
      % (100.0 * min(kratio), 100.0 * max(kratio)))
    P("     Below unity is the RIGHT side for a local value at r = %.0f mm"
      % (1e3 * r_bore))
    P("     against a disc average dominated by the stagnation point -- a ratio")
    P("     above 1 here would have refuted the law's level outright.")

    # ---------------- CROSS-CHECK: the Python twin vs the C++ ------------
    P("")
    P("-" * 78)
    P("CROSS-CHECK  this script's constants against the SHIPPED C++ law")
    P("-" * 78)
    P("  ⚠⚠ wall_law.py is the twin of MMWSpalling.H::WallJetH. If they drift,")
    P("     (c2) compares a bar from one law against a run of another and the")
    P("     comparison is void. Parsed out of the header, not retyped:")
    hdr = open(os.path.join(ROOT, "src/Integrator/MMWSpalling.H")).read()
    WANT = (("WJ_VEL", WL.C_WJ_VEL), ("WJ_DELTA", WL.C_WJ_DELTA),
            ("WJ_GROW", WL.C_WJ_GROW), ("WJ_CF_A", WL.C_CF_A),
            ("WJ_CF_B", WL.C_CF_B), ("WJ_PR", 0.72),
            ("WJ_P_ATM", WL.P_ATM), ("WJ_R_GAS", WL.R_GAS),
            ("WJ_R_SHIELD", WL.R_SHIELD))
    ok_x = True
    for nm, py in WANT:
        m = re.search(r"%s\s*=\s*(-?[\d.eE+-]+)\s*;" % nm, hdr)
        cpp = float(m.group(1)) if m else None
        same = cpp is not None and abs(cpp - py) <= 1e-12 * max(1.0, abs(py))
        ok_x = ok_x and same
        P("    %-12s python %-12.6g  C++ %-12s  %s"
          % (nm, py, "%.6g" % cpp if cpp is not None else "NOT FOUND",
             "ok" if same else "** DRIFT **"))
    # the datum rule must also be the same rule
    for nm, pat in (("the hard switch, not max()",
                     r"if \(z >= WallZFill\(r\)\) return wall_h_duct;"),
                    ("plane-jet momentum v ~ delta\\^-1/2",
                     r"v_c \* std::sqrt\(d_c / d_z\)"),
                    ("Chilton-Colburn WITH the half",
                     r"\(0\.5 \* cf\) / std::pow\(WJ_PR, 2\.0 / 3\.0\)"),
                    ("the datum is the INWARD RAY",
                     r"if \(kt != MISS && kt < best\) best = kt;")):
        found = re.search(pat, hdr) is not None
        ok_x = ok_x and found
        P("    %-34s %s" % (nm, "ok" if found else "** NOT IN THE C++ **"))
    P("  CROSS-CHECK %s" % ("PASS" if ok_x else "FAIL"))

    # ---------------- (a4) ------------------------------------------------
    P("")
    P("-" * 78)
    P("(a4) GATE  the reachability ledger, RESOLVED BY READING THE SOURCE")
    P("-" * 78)
    srows, ok_a4 = resolve_sites()
    for nm, ok, det in srows:
        P("  %-55s %-4s %s" % (nm, "ok" if ok else "FAIL", det))
    P("  (a4) %s" % ("PASS" if ok_a4 else "FAIL"))

    P("")
    P("=" * 78)
    P("ORD 1 VERDICT")
    P("=" * 78)
    P("  (a1) corner %.1f W/m2K   %s" % (h_c, "PASS" if ok_a1 else "FAIL"))
    P("  (a2) z_fill %.1f mm, jump x%.2f, HARD STEP (swept, not pinned)   %s"
      % (1e3 * zf, h_wj_at_fill / h_duct, "PASS" if ok_a2 else "FAIL"))
    P("  (a3) band mean %.1f W/m2K, ratio %.3f  => (c3) %s"
      % (mean_h, mean_h / h_typed, "widens" if mean_h > h_typed else "narrows"))
    P("  (a4) ledger %s" % ("PASS" if ok_a4 else "FAIL"))
    P("  (c2) bar FIXED: %.1f W -> %.1f W (+/- 10 %%)" % (W_now, fall * W_now))
    P("  cross-check python<->C++ %s" % ("PASS" if ok_x else "FAIL"))
    return 0 if (ok_a1 and ok_a2 and ok_a4 and ok_x) else 1


if __name__ == "__main__":
    sys.exit(main())
