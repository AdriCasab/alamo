#!/usr/bin/env python3
"""P0b revision 7 -- the Meier reduction.

Scores (h) and (i)'s Meier half, reports (f), deliverable 2's side split at the
LIVE gate, deliverable 6's closure comparison and deliverable 5's ledger.

usage: analyze_fh7.py <dir> <tag> [<control dir> <control tag>]
"""
import csv, os, sys
import numpy as np

W_LO, W_HI = 60.0, 160.0
CELL = 0.002 ** 3

def thermo(p):
    L = open(p).read().split("\n"); hd = L[0].split()
    rs = [x.split() for x in L[1:] if x.strip() and len(x.split()) == len(hd)]
    return {k: np.array([float(r[i]) for r in rs]) for i, k in enumerate(hd)}

def events(p):
    o = {k: [] for k in ("t", "h", "i", "j", "T", "Ts", "f", "hf", "ell", "nv")}
    for r in csv.DictReader(open(p)):
        if r.get("regime") != "1": continue
        o["t"].append(float(r["time"]));  o["h"].append(float(r["h_applied"]))
        o["i"].append(int(r["col_i"]));   o["j"].append(int(r["col_j"]))
        o["T"].append(float(r["T_top"])); o["Ts"].append(float(r.get("T_surf") or 0))
        o["f"].append(int(float(r.get("rc_flag") or 0)))
        o["hf"].append(float(r.get("h_f") or 0)); o["nv"].append(int(float(r["n_voided"])))
        o["ell"].append(float(r.get("ell_eff") or 0))
    return {k: np.array(v) for k, v in o.items()}

def summarise(D, TAG, label):
    th = thermo(os.path.join(D, TAG, "thermo.dat"))
    ev = events(os.path.join(D, TAG + "_removal_events.csv"))
    t = th["time"]; win = (t >= W_LO) & (t <= min(W_HI, t[-1]))
    w = (ev["t"] >= W_LO) & (ev["t"] <= min(W_HI, ev["t"].max()))
    ax = (ev["i"] == 0) & (ev["j"] == 0)
    span = min(W_HI, ev["t"].max()) - W_LO
    out = {
        "label": label, "t_end": t[-1],
        "rop_axis": ev["h"][ax & w].sum() / span * 3600,
        "vol_rate": ev["nv"][w].sum() * CELL / span,
        "hf_p50": 1e3 * np.percentile(ev["hf"][w], 50),
        "hf_p90": 1e3 * np.percentile(ev["hf"][w], 90),
        "ell_p50": 1e3 * np.percentile(ev["ell"][w], 50),
        "bound": np.mean((ev["f"][w] & 16) > 0),
        "search": np.mean((ev["f"][w] & 8) > 0),
        "flat": np.mean((ev["f"][w] & 1) > 0),
        "Ttop_p50": np.percentile(ev["T"][w], 50),
        "n_fire": int(w.sum()),
        "depth": 1e3 * ev["h"][ax].sum(),
    }
    for k, col in (("P_floor", "patch_P_robin"), ("P_side", "side_P_face"),
                   ("P_rin", "surf_side_P_rin"), ("P_rout", "surf_side_P_rout"),
                   ("P_top", "surf_side_P_top"), ("P_bur", "surf_side_P_buried"),
                   ("standoff", "patch_min_standoff"), ("err", "ledger_err_solid")):
        out[k] = float(th[col][win].mean()) if col in th else float("nan")
    if all(c in th for c in ("jet_P_face", "patch_P_robin", "side_P_face", "side_P_loss")):
        out["ident"] = float(np.abs(th["jet_P_face"] - (th["patch_P_robin"]
                        + th["side_P_face"] + th["side_P_loss"]))[win].max())
    return out

def main():
    runs = [summarise(sys.argv[1], sys.argv[2], "closure FOLLOW (rev7)")]
    if len(sys.argv) > 4:
        runs.append(summarise(sys.argv[3], sys.argv[4], "closure CELL (control)"))
    print("=== P0b rev7, rate window %.0f-%.0f s ===" % (W_LO, W_HI))
    keys = [("t_end", "run reached [s]", "%10.1f"),
            ("rop_axis", "axis ROP [m/h]", "%10.3f"),
            ("vol_rate", "removed volume rate [m3/s]", "%10.3e"),
            ("depth", "axis depth at end [mm]", "%10.1f"),
            ("hf_p50", "h_f p50 [mm]", "%10.3f"),
            ("hf_p90", "h_f p90 [mm]", "%10.3f"),
            ("ell_p50", "ell_eff p50 [mm]", "%10.3f"),
            ("bound", "bound by flake_h_max", "%10.4f"),
            ("search", "bound by h_spall_max", "%10.4f"),
            ("flat", "flat fallback", "%10.4f"),
            ("Ttop_p50", "T_top p50 at firing [K]", "%10.1f"),
            ("P_floor", "floor power [W]", "%10.1f"),
            ("P_side", "side power [W]", "%10.1f"),
            ("P_rin", "  side, r <  40 mm [W]", "%10.1f"),
            ("P_rout", "  side, r >= 40 mm [W]", "%10.1f"),
            ("P_top", "  side onto top cells [W]", "%10.1f"),
            ("P_bur", "  side onto buried [W]", "%10.1f"),
            ("standoff", "min standoff [m]", "%10.4f"),
            ("err", "ledger_err_solid", "%10.2e"),
            ("ident", "power identity max [W]", "%10.3e")]
    print("%-28s" % "" + "".join("%24s" % r["label"] for r in runs))
    for k, lab, fmt in keys:
        row = "".join(("%24s" % (fmt % r[k])) if k in r else "%24s" % "-" for r in runs)
        print("%-28s%s" % (lab, row))
    r = runs[0]
    print()
    print("(h) capped fraction %.1f %% [< 20 %%]  %s ; h_f p50 %.3f mm [0.3-1.0]  %s"
          % (100 * r["bound"], "PASS" if r["bound"] < 0.20 else "FAIL",
             r["hf_p50"], "PASS" if 0.3 <= r["hf_p50"] <= 1.0 else "FAIL"))
    print("(i) Meier half: binds on %.1f %% [> 10 %%]  %s"
          % (100 * r["bound"], "PASS" if r["bound"] > 0.10 else "FAIL"))
    print("witness -- the re-bin is NOT at the dead gate: side r<40mm = %.1f W (must be > 0)"
          % r["P_rin"])
    if len(runs) > 1:
        a, b = runs[0], runs[1]
        print("\ndeliverable 6 -- the wall closure's real effect, same flake bound:")
        print("  removed volume  %.3e -> %.3e m3/s   %+.1f %%"
              % (b["vol_rate"], a["vol_rate"], 100 * (a["vol_rate"] / b["vol_rate"] - 1)))
        print("  axis ROP        %.3f -> %.3f m/h     %+.1f %%"
              % (b["rop_axis"], a["rop_axis"], 100 * (a["rop_axis"] / b["rop_axis"] - 1)))
        print("  delivered power %.1f -> %.1f W       %+.1f %%"
              % (b["P_floor"] + b["P_side"], a["P_floor"] + a["P_side"],
                 100 * ((a["P_floor"] + a["P_side"]) / (b["P_floor"] + b["P_side"]) - 1)))

main()
