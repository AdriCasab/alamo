#!/usr/bin/env python3
"""P0b rev5 -- the Meier FH reduction.

Covers criterion (b) (the power identity on a march-exercising leg),
deliverable 9 (the side-power bins, OVER THE RATE WINDOW) and deliverable 10
(the ledger residual).  Takes a run directory containing <tag>/thermo.dat and
<tag>_removal_events.csv.

usage: analyze_fh.py <dir> <tag>          e.g. analyze_fh.py output_fh_rev5 FH
"""
import csv, os, sys
import numpy as np

W_LO, W_HI = 60.0, 160.0          # THE RATE WINDOW -- everything is quoted here
IDENT_BAR  = 0.02                 # W, criterion (b)

def thermo(p):
    L = open(p).read().split("\n")
    hd = L[0].split()
    rs = [x.split() for x in L[1:] if x.strip()]
    d = {}
    for i, k in enumerate(hd):
        try:    d[k] = np.array([float(r[i]) for r in rs])
        except (IndexError, ValueError): pass
    return d

def main():
    D, TAG = sys.argv[1], sys.argv[2]
    th = thermo(os.path.join(D, TAG, "thermo.dat"))
    t  = th["time"]
    win = (t >= W_LO) & (t <= W_HI)
    late = t >= 400.0
    print("=== %s : %d thermo rows, t = %.3f .. %.3f s ===" % (D, len(t), t[0], t[-1]))

    # ---------- (b) THE POWER IDENTITY ----------
    print("\n--- (b) SCORED: power identity, jet_P_face == patch_P_robin + "
          "side_P_face + side_P_loss ---")
    if all(k in th for k in ("jet_P_face","patch_P_robin","side_P_face","side_P_loss")):
        res = np.abs(th["jet_P_face"] - (th["patch_P_robin"] + th["side_P_face"]
                                         + th["side_P_loss"]))
        for lab, m in (("rate window 60-160 s", win), ("t >= 150 s", t >= 150.0),
                       ("whole run", np.ones_like(t, bool))):
            if m.any():
                e = float(res[m].max())
                print("  %-22s max residual %9.3e W   [bar %.2f]  %s"
                      % (lab, e, IDENT_BAR, "PASS" if e <= IDENT_BAR else "FAIL"))
        # the march must actually be carrying side power, or (b) is vacuous
        sp = th["side_P_face"][win]
        print("  march-exercising: side_P_face over the window = %.1f .. %.1f W, "
              "jet_P_face = %.1f .. %.1f W"
              % (sp.min(), sp.max(), th["jet_P_face"][win].min(), th["jet_P_face"][win].max()))
    else:
        print("  (columns absent)")

    # ---------- deliverable 9: THE SIDE-POWER BINS, OVER THE WINDOW ----------
    print("\n--- deliverable 9: the side share, binned OVER THE RATE WINDOW ---")
    if "surf_side_share" in th:
        s = th["surf_side_share"]
        for lab, m in (("60-160 s (the rate window)", win), ("t >= 400 s", late),
                       ("end of run", np.array([len(t)-1]))):
            mm = m if m.dtype == bool else np.isin(np.arange(len(t)), m)
            print("  surf_side_share  %-26s %.4f" % (lab, float(s[mm].mean())))
    def binshare(num, den, lab):
        if num in th and den in th:
            a, b = th[num][win].mean(), th[den][win].mean()
            print("  %-34s %8.2f / %8.2f W = %.4f" % (lab, a, b, a / b if b else 0.0))
    print("  bin 1 -- by RADIUS at the 40 mm gate (side_P_wall is r >= foot_r_outer):")
    if all(k in th for k in ("side_P_bore","side_P_pads","side_P_wall")):
        b_, p_, w_ = (th[k][win].mean() for k in ("side_P_bore","side_P_pads","side_P_wall"))
        tot = b_ + p_ + w_
        print("    bore %8.2f  pads %8.2f  wall %8.2f W   wall share %.4f"
              % (b_, p_, w_, w_ / tot if tot else 0.0))
    print("  bin 2 -- does the RECEIVING CELL fire next? (top cell vs buried):")
    if "surf_side_P_top" in th and "surf_side_P_buried" in th:
        a, b = th["surf_side_P_top"][win].mean(), th["surf_side_P_buried"][win].mean()
        print("    top-cell %8.2f  buried %8.2f W   top share %.4f"
              % (a, b, a / (a + b) if (a + b) else 0.0))
    else:
        print("    (surf_side_P_top absent -- rev2 run)")
    print("  bin 3 -- MULTI-CELL RISER share (7b is exact only for a one-cell step):")
    for k in ("wall_step_f_multi", "wall_step_f_multi_pads"):
        if k in th: print("    %-24s %.4f" % (k, float(th[k][win].mean())))

    # ---------- deliverable 10: THE LEDGER RESIDUAL ----------
    print("\n--- deliverable 10: the ledger residual ---")
    for k in ("ledger_err", "ledger_err_solid"):
        if k in th: print("  %-18s window %+.4f   end %+.4f"
                          % (k, float(th[k][win].mean()), float(th[k][-1])))
    if all(k in th for k in ("ledger_E_H_solid","ledger_E_src","ledger_E_flake")):
        EH, ES, EF = (th[k] for k in ("ledger_E_H_solid","ledger_E_src","ledger_E_flake"))
        dEH = EH - EH[0]
        # what the balance SHOULD be:  dE_H_solid = E_src - E_flake
        miss = dEH - (ES - EF)
        print("  dE_H_solid %+.4e J   E_src %+.4e J   E_flake %+.4e J"
              % (dEH[-1], ES[-1], EF[-1]))
        print("  unaccounted = dE_H_solid - (E_src - E_flake) = %+.4e J "
              "(%.1f %% of E_src)" % (miss[-1], 100 * miss[-1] / ES[-1] if ES[-1] else 0))
        # does the gap track the firings, the side power, or the beam?
        dt = np.diff(t)
        rate = np.diff(miss) / np.where(dt > 0, dt, np.nan)      # W of unaccounted
        for k in ("side_P_face", "patch_P_robin", "ledger_P_loss", "jet_P_face"):
            if k in th:
                x = th[k][1:]
                g = np.isfinite(rate) & np.isfinite(x)
                if g.sum() > 10:
                    c = float(np.corrcoef(rate[g], x[g])[0, 1])
                    print("    corr(unaccounted rate, %-14s) = %+.3f  "
                          "(mean %8.1f W vs %8.1f W)"
                          % (k, c, np.nanmean(rate[g]), x[g].mean()))
    if "ledger_surface_missing_cols" in th:
        print("  ledger_surface_missing_cols max %g"
              % float(th["ledger_surface_missing_cols"].max()))

    # ---------- (f) the Meier numbers ----------
    ev = os.path.join(D, TAG + "_removal_events.csv")
    if not os.path.exists(ev):
        print("\n(no event log)"); return
    print("\n--- (f) REPORTED: the Meier numbers ---")
    ti, ha, tt, ts, fl, ci, cj = [], [], [], [], [], [], []
    with open(ev) as f:
        for r in csv.DictReader(f):
            if r.get("regime") != "1": continue
            ti.append(float(r["time"])); ha.append(float(r["h_applied"]))
            tt.append(float(r["T_top"])); ci.append(int(r["col_i"])); cj.append(int(r["col_j"]))
            ts.append(float(r.get("T_surf") or 0.0)); fl.append(float(r.get("rc_flag") or 0.0))
    ti, ha, tt, ts, fl = map(np.array, (ti, ha, tt, ts, fl))
    ci, cj = np.array(ci), np.array(cj)
    ax = (ci == int(np.median(ci))) & (cj == int(np.median(cj)))   # the axis column
    m = ax & (ti >= W_LO) & (ti <= W_HI)
    print("  axis column (%d,%d): %d firings in the window"
          % (int(np.median(ci)), int(np.median(cj)), int(m.sum())))
    print("  axis ROP 60-160 s        %.3f m/h" % (ha[m].sum() / (W_HI - W_LO) * 3600))
    print("  axis depth at t = %.1f s  %.1f mm" % (ti[ax].max(), 1e3 * ha[ax].sum()))
    w = (ti >= W_LO) & (ti <= W_HI)
    print("  firing window p10-p90     %.1f K  (T_top)"
          % (np.percentile(tt[w], 90) - np.percentile(tt[w], 10)))
    print("  flat fallback at firings  %.1f %%"
          % (100 * np.mean(np.array([int(x) & 1 for x in fl])[w])))
    print("  T_surf - T_top at firing  %+.1f K" % float(np.mean(ts[w] - tt[w])))
    print("  firings total %d, window %d" % (len(ti), int(w.sum())))

main()
