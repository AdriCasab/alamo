#!/usr/bin/env python3
"""D2c deliverable 10: the corrected Meier scoring (ACTIVE_STEP §Meier targets),
applied in D2d. Reuses d2b_feet_rop/analyze.py (Run, analyze: steady-window
finder, ROP, ΔT_fire, R, ledger, pinned shares) and adds the D2c metrics.

  score.py RUN_DIR [RUN_DIR ...]      print the scored table of each run
  score.py --dry-run RUN_DIR          same, into RESULTS.md between <!-- DRYRUN --> markers
  score.py --mesh RUN_2MM RUN_1MM     full-domain mesh verdict (both on steady windows)
  score.py --sweep RUN_1600 RUN_1750 RUN_1900   ROP(T_nozzle) and the T range giving 1.3-1.6 m/h
  score.py --p6 RUN [RUN ...]         evaluate every PREDICTIONS.md P6 line (held / refuted /
                                      not tested), runs identified by their .done name
  score.py --window-a RUN [RUN ...]   the window (a) table (D2d amendment iii)

A RUN_DIR is a plot_file prefix with <prefix>.done (run metadata, as written by
the D2b/D2d run.py), <prefix>/thermo.dat and <prefix>_removal_events.csv.

Definitions (quarter domain; volumes and powers x4 = full hole):
  steady window   analyze.Run.steady(150): the earliest t0 with [t0, t_end]
                  >= 150 s, 50 s descent sub-fits within ±10 % of their mean
                  and |d jet_s_c/dt| < 0.02 mm/s. None -> "no steady window":
                  the metrics are shown over analyze's fallback window (the
                  last max(150 s, t_end/2)) and every verdict is FAIL/"not
                  scored" for lack of a window.
  ROP             nozzle_z descent rate over the window.
  whole volume    4 * sum h_applied dx dy over every removal row in the window
                  (thermal and mechanical, regime 4) / window length.
                  Inside-Ø 93 is reported, not scored.
  mechanical      4 * (foot_mech_vol(t1) - foot_mech_vol(t0)) / (whole volume
                  in the window); > 5 % is a finding.
  wall profile    Ø(z) = 2 x azimuth-mean r_wall(z) at the window end: the
                  quarter is split into N_SEC = 9 equal azimuth sectors, r_wall
                  of a sector = max r of its columns deeper than z (0 if none),
                  mean over the sectors. z = depth below the original top.
  drilled depth   0 to the nozzle-plane depth at the window end (lz - nozzle_z).
  feet depth      0 to the feet depth at the window end (lz - foot_z), i.e. 50 mm
                  deeper than the drilled depth. D2d amendment (i): the minimum Ø
                  is scored over this range, because that band is where the burner
                  body passes; the nozzle-plane minimum is also reported.
  window (a)      D2d amendment (iii): [150 s, t_a], t_a = the first time the centre
                  depth lz - (nozzle_z - jet_s_c) reaches 0.36 m, i.e. where a 0.40 m
                  domain would have stopped; t_a = t_end (flagged) if never reached.
                  The P6 bands of PREDICTIONS.md are stated on this window.
  mouth           Fig. 8.8 (validation/meier/meier_fig8_8_hole_profile.csv,
                  visible width = lower bound): PASS iff at z = 25 and 50 mm
                  0 <= Ø - width <= 10 mm. Every CSV row inside the drilled depth
                  enters the RMS; any such row with Ø < width is listed (and
                  fails the profile, since the width is a lower bound).
  depth-mean Ø    mean of Ø(z) over the drilled depth: 85-93 mm.
  minimum Ø       min of Ø(z) over the FEET depth: >= 80 PASS; 76-80
                  "unresolved at 2 mm" (one cell = 4 mm in Ø); < 76 FAIL.
  removal power   rho Cp ΔT_fire x whole volume rate (x4 included) vs 2.83 kW;
                  R = that / (4 mean jet_P_face) <= 1.05.
The script also writes output/<run name>_fig88_profile.csv in this study
(depth, model Ø azimuth-mean, model Ø max-r, Meier visible width); nothing is
written next to the scored run.
"""
import argparse
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUDIES = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(STUDIES, "d2b_feet_rop"))
import analyze as an  # noqa: E402

wj = an.wj
FIG = os.path.join(os.path.dirname(STUDIES), "validation", "meier", "meier_fig8_8_hole_profile.csv")
RESULTS = os.path.join(HERE, "RESULTS.md")
OUT = os.path.join(HERE, "output")
N_SEC = 9
BAND_ROP = (1.04, 1.92)
BAND_V = (1.98, 2.96)
BAND_D = (85.0, 93.0)
BAND_DT = (500.0, 560.0)
MOUTH_Z = (25.0, 50.0)
Z_A = 0.36                    # D2d amendment (iii): the 0.40 m domain's bottom-stop depth
MOUTH_TOL = 10.0
MIN_D, MIN_D_UNRES = 80.0, 76.0
MECH_MAX = 0.05
R_MAX = 1.05
P_MEIER = wj.P_MEIER_ROCK
MESH_TOL = 0.05
LEDGER_TOL = 1e-10


def load_fig():
    rows = []
    for line in open(FIG):
        if line.startswith("#") or not line.strip():
            continue
        parts = [p.strip() for p in line.split(",")]
        rows.append((float(parts[0]), float(parts[1]), parts[2] if len(parts) > 2 else ""))
    return rows


def load_run(pf):
    m = an.meta_of(pf)
    if m is None:
        raise SystemExit(f"{pf}.done missing (run metadata)")
    return an.Run(pf, m)


def az_rwall(run, depth, z, nsec=N_SEC):
    xc = (np.arange(run.nx) + 0.5) * run.dz
    X, Y = np.meshgrid(xc, xc, indexing="xy")
    a = np.arctan2(Y, X).ravel()
    sec = np.minimum(nsec - 1, (a / (0.5 * math.pi) * nsec).astype(int))
    out = np.zeros(nsec)
    deep = depth > z
    for s in range(nsec):
        sel = deep & (sec == s)
        out[s] = float(run.r[sel].max()) if sel.any() else 0.0
    return float(out.mean())


def verdict(ok, steady=True):
    if not steady:
        return "FAIL (no window)"
    return "PASS" if ok else "FAIL"


def score(pf, write_profile=True):
    run = load_run(pf)
    r = an.analyze(run, "B")
    th, m = run.th, run.m
    W = r["W"]
    steady = r["steady"] is not None
    dur = W[1] - W[0]
    t_end_w = W[1]
    nz_end = float(np.interp(t_end_w, th["time"], th["nozzle_z"]))
    drilled = run.lz - nz_end
    d_w = run.depth(t_end_w)
    feet = run.lz - float(np.interp(t_end_w, th["time"], th["foot_z"]))
    zg = np.arange(0.0, max(feet, run.dz), run.dz / 2.0)
    D_az = np.array([2.0 * az_rwall(run, d_w, z) for z in zg]) * 1e3
    D_mx = np.array([2.0 * run.r_wall(d_w, z) for z in zg]) * 1e3
    in_drilled = zg <= drilled
    fig = load_fig()
    fig_in = [(z, w, n) for z, w, n in fig if z * 1e-3 <= drilled]
    model_at = lambda z_mm: float(np.interp(z_mm * 1e-3, zg, D_az))
    mouth = {z: (model_at(z), next(w for zz, w, _ in fig if zz == z)) for z in MOUTH_Z}
    mouth_ok = all(0.0 <= mv - w <= MOUTH_TOL for mv, w in mouth.values())
    below = [(z, model_at(z), w) for z, w, _ in fig_in if model_at(z) < w]
    rms = float(np.sqrt(np.mean([(model_at(z) - w) ** 2 for z, w, _ in fig_in]))) if fig_in else float("nan")
    # depth-mean and the mouth keep the nozzle-plane range; the minimum is
    # taken to the feet depth (D2d amendment (i)), with both reported
    dmean = float(np.mean(D_az[in_drilled]))
    dmin_drill = float(np.min(D_az[in_drilled]))
    dmin = float(np.min(D_az))
    # volumes
    vdot = r["vdot"]
    mech = float("nan")
    if "foot_mech_vol" in th:
        mv = np.interp(W, th["time"], th["foot_mech_vol"])
        mech = 4.0 * float(mv[1] - mv[0]) / max(1e-30, vdot * dur)
    P_rm = wj.RHOCP * r["dT"] * vdot
    rop = r["rop"]
    rows = [
        ("steady window", "≥ 150 s, |ds_c/dt| < 0.02 mm/s, 50 s sub-fits ±10 %",
         (f"{W[0]:.0f}–{W[1]:.0f} s" if steady else
          f"none (fallback {W[0]:.0f}–{W[1]:.0f} s; s_c drift {r['sc_drift'] * 1e3:+.2f} mm/s; sub-fits "
          + ", ".join(f"{x:.2f}" for x in r["subs"]) + " m/h)"),
         "exists" if steady else "FAIL: no steady window"),
        ("ROP", "[1.04, 1.92] m/h", f"{rop:.2f} m/h (ring columns {r['ring_rop']:.2f})",
         verdict(BAND_ROP[0] <= rop <= BAND_ROP[1], steady)),
        ("whole volume rate", "[1.98, 2.96] cm³/s",
         f"{vdot * 1e6:.2f} cm³/s (inside Ø 93, reported: {r['vdot_hole'] * 1e6:.2f})",
         verdict(BAND_V[0] <= vdot * 1e6 <= BAND_V[1], steady)),
        ("mechanical share", "≤ 5 % (above = finding)",
         "no clearance in this run" if not np.isfinite(mech) else f"{mech * 100:.2f} %",
         "—" if not np.isfinite(mech) else ("ok" if mech <= MECH_MAX else "FINDING")),
        ("wall profile (mouth)", "0 ≤ Ø − width ≤ 10 mm at 25 and 50 mm (Fig. 8.8)",
         "; ".join(f"z {z:.0f}: Ø {mv:.1f} vs {w:.0f} ({mv - w:+.1f})" for z, (mv, w) in mouth.items())
         + f"; RMS over {len(fig_in)} rows inside {drilled * 1e3:.0f} mm: {rms:.1f} mm; rows below the width: "
         + (", ".join(f"{z:.0f} mm ({mv:.0f} < {w:.0f})" for z, mv, w in below) if below else "none"),
         verdict(mouth_ok and not below, steady)),
        ("depth-mean Ø (drilled depth)", "85–93 mm",
         f"{dmean:.1f} mm over 0–{drilled * 1e3:.0f} mm (max-r Ø {float(np.mean(D_mx[in_drilled])):.1f})",
         verdict(BAND_D[0] <= dmean <= BAND_D[1], steady)),
        ("minimum Ø (feet depth)", "≥ 80 mm; 76–80 unresolved at 2 mm",
         f"{dmin:.1f} mm at {zg[int(np.argmin(D_az))] * 1e3:.0f} mm over 0–{feet * 1e3:.0f} mm "
         f"(max-r Ø {float(np.min(D_mx)):.1f}); over the drilled depth only: {dmin_drill:.1f} mm",
         ("PASS" if dmin >= MIN_D else "unresolved at 2 mm" if dmin >= MIN_D_UNRES else "FAIL")
         if steady else "FAIL (no window)"),
        ("ΔT_fire", "500–560 K", f"{r['dT']:.0f} K (T_fire p10–p90 {r['Tf'][1]:.0f}–{r['Tf'][2]:.0f} K)",
         verdict(BAND_DT[0] <= r["dT"] <= BAND_DT[1], steady)),
        ("removal power / R", "vs 2.83 kW; R ≤ 1.05",
         f"{P_rm / 1e3:.2f} kW ({P_rm / P_MEIER:.2f}× Meier); R = {r['R']:.3f} "
         f"(face power ×4 {r['P_face'] / 1e3:.2f} kW)",
         verdict(r["R"] <= R_MAX, steady)),
        ("ledger", "round-off", f"{r['ledger']:.1e}", "PASS" if r["ledger"] < LEDGER_TOL else "FAIL"),
        ("mesh", "full domain 1 mm vs 2 mm, ROP within 5 %", "score.py --mesh (D2d run 7)", "—"),
    ]
    if write_profile:
        os.makedirs(OUT, exist_ok=True)
        with open(os.path.join(OUT, os.path.basename(pf) + "_fig88_profile.csv"), "w") as f:
            f.write("# D2c score.py: model hole profile at the window end "
                    f"({W[0]:.1f}-{W[1]:.1f} s{'' if steady else ', NOT steady'}) vs Meier Fig. 8.8\n")
            f.write("depth_mm,model_D_azimuth_mean_mm,model_D_max_r_mm,meier_visible_width_mm\n")
            for z, w, _ in fig:
                if z * 1e-3 <= drilled:
                    f.write(f"{z:.0f},{model_at(z):.2f},{float(np.interp(z * 1e-3, zg, D_mx)):.2f},{w:.0f}\n")
    info = dict(run=run, r=r, W=W, steady=steady, drilled=drilled, feet=feet, zg=zg, D_az=D_az,
                D_mx=D_mx, mouth=mouth, rms=rms, below=below, dmean=dmean, dmin=dmin,
                dmin_drill=dmin_drill, vdot=vdot, mech=mech, status=m.get("status", "?"),
                t_end=run.t_end, name=m.get("name", os.path.basename(pf)))
    return rows, info


def table(name, rows, info):
    L = [f"### {name} ({info['status']} at {info['t_end']:.1f} s)", "",
         "| criterion | target | simulation | verdict |", "|---|---|---|---|"]
    L += [f"| {a} | {b} | {c} | {d} |" for a, b, c, d in rows]
    return L


def centre_depth(run):
    th = run.th
    return run.lz - (th["nozzle_z"] - th["jet_s_c"]), th["time"]


def t_a_of(run):
    """First time the centre depth reaches Z_A (the 0.40 m domain's bottom
    stop). Returns (t_a, reached)."""
    zc, t = centre_depth(run)
    ok = (t > 0.0) & (zc >= Z_A)
    if np.any(ok):
        return float(t[np.argmax(ok)]), True
    return float(t[-1]), False


def window_a(pf):
    """D2d amendment (iii): metrics over [150 s, t_a]."""
    run = load_run(pf)
    r = an.analyze(run, "B")
    th = run.th
    t_a, reached = t_a_of(run)
    W = (150.0, t_a)
    if t_a <= W[0] + 10.0:
        return dict(run=run, name=run.m.get("name"), W=W, reached=reached, ok=False)
    sel = (th["time"] >= W[0]) & (th["time"] <= W[1])
    rop = -float(np.polyfit(th["time"][sel], th["nozzle_z"][sel], 1)[0]) * 3600.0
    dur = W[1] - W[0]
    ev = run.ev
    inw = (ev["time"] > W[0]) & (ev["time"] <= W[1])
    vdot = 4.0 * float(np.sum(ev["h_applied"][inw])) * run.dz * run.dz / dur
    mech = float("nan")
    if "foot_mech_vol" in th:
        mv = np.interp(W, th["time"], th["foot_mech_vol"])
        mech = 4.0 * float(mv[1] - mv[0]) / max(1e-30, vdot * dur)
    nz_end = float(np.interp(W[1], th["time"], th["nozzle_z"]))
    drilled = run.lz - nz_end
    feet = run.lz - float(np.interp(W[1], th["time"], th["foot_z"]))
    d_w = run.depth(W[1])
    zg = np.arange(0.0, max(feet, run.dz), run.dz / 2.0)
    D_az = np.array([2.0 * az_rwall(run, d_w, z) for z in zg]) * 1e3
    ind = zg <= drilled
    at = lambda z_mm: float(np.interp(z_mm * 1e-3, zg, D_az))
    fig = {z: w for z, w, _ in load_fig()}
    dT = r["Tf"][0] - 293.15
    sel_t = (th["time"] >= W[0]) & (th["time"] <= W[1])
    return dict(run=run, name=run.m.get("name"), W=W, reached=reached, ok=True, rop=rop,
                vdot=vdot, mech=mech, drilled=drilled, feet=feet, zg=zg, D_az=D_az,
                d25=at(25.0), d50=at(50.0), fig=fig, dmean=float(np.mean(D_az[ind])),
                dmin=float(np.min(D_az)), dmin_drill=float(np.min(D_az[ind])), dT=dT,
                sc=float(np.mean(th["jet_s_c"][sel_t])),
                sc_drift=float(np.polyfit(th["time"][sel_t], th["jet_s_c"][sel_t], 1)[0]),
                T_rec=float(np.mean(th["jet_T_rec"][sel_t])) if "jet_T_rec" in th else float("nan"),
                steady=an.analyze(run, "B")["steady"] is not None)


def window_a_table(pfs):
    L = ["| run | window (a) | centre 0.36 m reached | ROP | whole volume | mech share | Ø 25 | Ø 50 | "
         "depth-mean Ø | min Ø (feet) | mean s_c | s_c drift | mean T_rec |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    out = {}
    for pf in pfs:
        w = window_a(pf)
        out[w["name"]] = w
        if not w["ok"]:
            L.append(f"| {w['name']} | too short ({w['W'][0]:.0f}–{w['W'][1]:.0f} s) | "
                     f"{w['reached']} | — | — | — | — | — | — | — | — | — | — |")
            continue
        mech_txt = "—" if not np.isfinite(w["mech"]) else f"{w['mech'] * 100:.2f} %"
        L.append(f"| {w['name']} | {w['W'][0]:.0f}–{w['W'][1]:.0f} s | {'yes' if w['reached'] else 'NO (t_end)'} | "
                 f"{w['rop']:.2f} m/h | {w['vdot'] * 1e6:.2f} cm³/s | {mech_txt} | "
                 f"{w['d25']:.0f} | {w['d50']:.0f} | {w['dmean']:.0f} | {w['dmin']:.0f} | "
                 f"{w['sc'] * 1e3:.0f} mm | {w['sc_drift'] * 1e3:+.3f} mm/s | {w['T_rec']:.0f} K |")
    return L, out


# PREDICTIONS.md P6 (hashed a28f1987...): the bands are stated on window (a).
P6_BANDS = dict(rop_scored=(1.31, 2.27), rop_clamp=(1.21, 2.10), rop_1600=(0.91, 1.57),
                rop_1750=(1.11, 1.91), vol_scored=(3.17, 5.48), dmean_scored=(91.0, 131.0),
                sc_scored=(160.0, 200.0), d25_scored=134.0, d50_scored=118.0,
                d25_run6=174.0, rop_run6=1.71, mech_max=5.0)


def p6(pfs):
    """Evaluate every P6 line mechanically: held / refuted / not tested."""
    _, W = window_a_table(pfs)
    L = ["| P6 item | prediction | measurement | verdict |", "|---|---|---|---|"]

    def band_row(item, pred, w, key, band, unit="m/h", scale=1.0):
        if w is None or not w.get("ok"):
            L.append(f"| {item} | {pred} | run not available | not tested |")
            return
        v = w[key] * scale
        L.append(f"| {item} | {pred} | {v:.2f} {unit} (window (a) {w['W'][0]:.0f}–{w['W'][1]:.0f} s) | "
                 f"{'held' if band[0] <= v <= band[1] else 'REFUTED'} |")

    r1, r2 = W.get("R1_scored"), W.get("R2_clamp")
    r4, r5, r6 = W.get("R4_1600"), W.get("R5_1750"), W.get("R6_d2bjet_n1")
    # P1 scored
    if r1 and r1.get("ok"):
        sc = r1["sc"] * 1e3
        settles = abs(r1["sc_drift"]) * 1e3 < 0.02
        parts = [f"steady window: {'yes' if r1['steady'] else 'no'}",
                 f"mean s_c {sc:.0f} mm, drift {r1['sc_drift'] * 1e3:+.3f} mm/s"]
        ref = (r1["steady"] and not r1["reached"]) or (settles and not
                                                       P6_BANDS["sc_scored"][0] <= sc <= P6_BANDS["sc_scored"][1])
        L.append(f"| P1 scored | no steady window before the 0.40 m bottom depth; s_c settles in "
                 f"{P6_BANDS['sc_scored'][0]:.0f}–{P6_BANDS['sc_scored'][1]:.0f} mm | {'; '.join(parts)} | "
                 f"{'REFUTED' if ref else 'held'} |")
    else:
        L.append("| P1 scored | as above | run not available | not tested |")
    # P1 clamp
    if r2 and r2.get("ok"):
        eq = abs(r2["sc_drift"]) * 1e3 < 0.02
        L.append(f"| P1 clamp | no centre equilibrium (|ds_c/dt| stays ≥ 0.02 mm/s) | "
                 f"drift {r2['sc_drift'] * 1e3:+.3f} mm/s, steady window "
                 f"{'yes' if r2['steady'] else 'no'} | {'REFUTED' if eq else 'held'} |")
    else:
        L.append("| P1 clamp | as above | run not available | not tested |")
    L.append("| P1 n = 0.5 | needs a ≥ 1.55 m domain | R3_n05 not run (decision) | not tested |")
    band_row("P2 scored ROP", f"{P6_BANDS['rop_scored'][0]}–{P6_BANDS['rop_scored'][1]} m/h",
             r1, "rop", P6_BANDS["rop_scored"])
    band_row("P2 clamp ROP", f"{P6_BANDS['rop_clamp'][0]}–{P6_BANDS['rop_clamp'][1]} m/h",
             r2, "rop", P6_BANDS["rop_clamp"])
    band_row("P3 1600 K ROP", f"{P6_BANDS['rop_1600'][0]}–{P6_BANDS['rop_1600'][1]} m/h",
             r4, "rop", P6_BANDS["rop_1600"])
    band_row("P3 1750 K ROP", f"{P6_BANDS['rop_1750'][0]}–{P6_BANDS['rop_1750'][1]} m/h",
             r5, "rop", P6_BANDS["rop_1750"])
    if all(x and x.get("ok") for x in (r4, r5, r1)):
        v = [r4["rop"], r5["rop"], r1["rop"]]
        mono = v[0] < v[1] < v[2]
        L.append(f"| P3 ordering | ROP increases with T_nozzle | {v[0]:.2f} → {v[1]:.2f} → {v[2]:.2f} m/h "
                 f"(1600 → 1750 → 1900 K) | {'held' if mono else 'REFUTED'} |")
    else:
        L.append("| P3 ordering | ROP increases with T_nozzle | runs not available | not tested |")
    # P4
    if r1 and r1.get("ok"):
        w25, w50 = r1["fig"][25.0], r1["fig"][50.0]
        mouth_pass = (0.0 <= r1["d25"] - w25 <= 10.0) and (0.0 <= r1["d50"] - w50 <= 10.0)
        L.append(f"| P4 mouth | stays too wide: {P6_BANDS['d25_scored']:.0f} / {P6_BANDS['d50_scored']:.0f} mm "
                 f"±20 at 25 / 50 mm | {r1['d25']:.0f} / {r1['d50']:.0f} mm vs Fig. 8.8 {w25:.0f} / {w50:.0f} | "
                 f"{'REFUTED (the mouth passes)' if mouth_pass else 'held'} |")
        dm = r1["dmean"]
        L.append(f"| P4 depth-mean Ø | {P6_BANDS['dmean_scored'][0]:.0f}–{P6_BANDS['dmean_scored'][1]:.0f} mm | "
                 f"{dm:.0f} mm | {'held' if P6_BANDS['dmean_scored'][0] <= dm <= P6_BANDS['dmean_scored'][1] else 'REFUTED'} |")
        band_row("P4 whole volume", f"{P6_BANDS['vol_scored'][0]}–{P6_BANDS['vol_scored'][1]} cm³/s",
                 r1, "vdot", P6_BANDS["vol_scored"], unit="cm³/s", scale=1e6)
        mech = r1["mech"] * 100.0 if np.isfinite(r1["mech"]) else float("nan")
        L.append(f"| P4 mechanical share | ≤ 5 % | {mech:.2f} % | "
                 f"{'REFUTED' if mech > P6_BANDS['mech_max'] else 'held'} |")
    else:
        for item in ("P4 mouth", "P4 depth-mean Ø", "P4 whole volume", "P4 mechanical share"):
            L.append(f"| {item} | see PREDICTIONS P6 | R1 not available | not tested |")
    # P5 run 6
    if r6 and r6.get("ok") and r1 and r1.get("ok"):
        wider = r6["d25"] > r1["d25"]
        rop_ok = abs(r6["rop"] / P6_BANDS["rop_run6"] - 1.0) <= 0.27
        L.append(f"| P5 run 6 mouth | wider than run 1 (hand {P6_BANDS['d25_run6']:.0f} vs "
                 f"{P6_BANDS['d25_scored']:.0f} mm at 25 mm) | {r6['d25']:.0f} vs {r1['d25']:.0f} mm | "
                 f"{'held' if wider else 'REFUTED'} |")
        L.append(f"| P5 run 6 ROP | {P6_BANDS['rop_run6']:.2f} m/h ±27 % | {r6['rop']:.2f} m/h | "
                 f"{'held' if rop_ok else 'REFUTED'} |")
    else:
        L.append("| P5 run 6 | see PREDICTIONS P6 | runs not available | not tested |")
    oa, ob = W.get("Oa_core5"), W.get("Ob_ricou")
    if oa and ob and r1 and all(x.get("ok") for x in (oa, ob, r1)):
        order = ob["rop"] < oa["rop"] < r1["rop"]
        L.append(f"| P5 core 5 / Ricou ordering | Ricou < core 5 < scored (hand 1.72 < 1.75 < 1.79) | "
                 f"{ob['rop']:.2f} < {oa['rop']:.2f} < {r1['rop']:.2f} m/h | "
                 f"{'held' if order else 'REFUTED'} |")
    else:
        L.append("| P5 core 5 / Ricou ordering | Ricou < core 5 < scored | optional runs not available | "
                 "not tested |")
    L.append("| P5 aspect / exponent | hand model only | no D2d run (by design) | not tested |")
    return L


def mesh(pf2, pf1):
    """D2d amendment (ii): a MATCHED transient window, not steady windows.
    The window is [150 s, min(t_end of the pair)]; the criterion is the burner
    ROP within 5 %. s_c, T_rec, the pinned share and the centre depth are
    reported at the same times."""
    run2, run1 = load_run(pf2), load_run(pf1)
    t1 = min(run2.t_end, run1.t_end)
    W = (150.0, t1)
    L = []
    if t1 <= W[0] + 10.0:
        return f"matched window too short (t_end {run2.t_end:.0f} / {run1.t_end:.0f} s)"

    def over(run):
        th = run.th
        sel = (th["time"] >= W[0]) & (th["time"] <= W[1])
        rop = -float(np.polyfit(th["time"][sel], th["nozzle_z"][sel], 1)[0]) * 3600.0
        zc, t = centre_depth(run)
        return dict(rop=rop, sc=float(np.interp(W[1], th["time"], th["jet_s_c"])),
                    T_rec=float(np.interp(W[1], th["time"], th["jet_T_rec"])),
                    zc=float(np.interp(W[1], t, zc)),
                    pin=float(np.mean(th["pinned_cols_pinned"][sel] /
                                      np.maximum(th["pinned_cols_pinned"][sel] + th["pinned_cols_face"][sel], 1e-30))))

    a2, a1 = over(run2), over(run1)
    d = a1["rop"] / a2["rop"] - 1.0
    ok = abs(d) <= MESH_TOL
    L.append(f"Matched window {W[0]:.0f}–{W[1]:.0f} s (t_end {run1.t_end:.0f} s at 1 mm, "
             f"{run2.t_end:.0f} s at 2 mm; steady windows are NOT required, D2d decision 2).")
    L.append("")
    L.append("| quantity | 2 mm | 1 mm | 1 mm vs 2 mm |")
    L.append("|---|---|---|---|")
    L.append(f"| burner ROP over the window | {a2['rop']:.3f} m/h | {a1['rop']:.3f} m/h | {d * 100:+.1f} % |")
    L.append(f"| s_c at the window end | {a2['sc'] * 1e3:.1f} mm | {a1['sc'] * 1e3:.1f} mm | "
             f"{(a1['sc'] / a2['sc'] - 1) * 100:+.1f} % |")
    L.append(f"| T_rec at the window end | {a2['T_rec']:.0f} K | {a1['T_rec']:.0f} K | "
             f"{(a1['T_rec'] / a2['T_rec'] - 1) * 100:+.1f} % |")
    L.append(f"| centre depth at the window end | {a2['zc'] * 1e3:.0f} mm | {a1['zc'] * 1e3:.0f} mm | "
             f"{(a1['zc'] / a2['zc'] - 1) * 100:+.1f} % |")
    L.append(f"| pinned share | {a2['pin']:.2f} | {a1['pin']:.2f} | {(a1['pin'] - a2['pin']):+.2f} |")
    L.append("")
    L.append(f"**Mesh verdict: {'PASS' if ok else 'FAIL'}** (burner ROP within {MESH_TOL * 100:.0f} %).")
    return "\n".join(L)


def sweep(pfs):
    """ROP(T_nozzle). Steady-window ROPs are used only if EVERY run in the sweep
    has a steady window; otherwise every point is the window (a) ROP, labelled
    transient. The two bases are never mixed in one interpolation (ACTIVE_STEP
    D2d §Scoring)."""
    pts = []
    for pf in pfs:
        rows, info = score(pf, False)
        w = window_a(pf)
        pts.append((info["run"].m["T_nozzle"], info["r"]["rop"], info["steady"],
                    w["rop"] if w.get("ok") else float("nan"), w.get("W"), info["W"]))
    pts.sort()
    all_steady = all(p[2] for p in pts)
    basis = "steady window" if all_steady else "window (a), TRANSIENT"
    idx = 1 if all_steady else 3
    T = np.array([p[0] for p in pts])
    v = np.array([p[idx] for p in pts])
    out = [f"Basis: {basis}" + ("" if all_steady else
                                " (not every run has a steady window, so no steady ROP enters the fit)"),
           "",
           "| T_nozzle | steady-window ROP | window (a) ROP | window used for the fit |",
           "|---|---|---|---|"]
    for t, rs, st, ra, wa, ws in pts:
        out.append(f"| {t:.0f} K | {rs:.2f} m/h {'(steady)' if st else '(no window: fallback)'} | "
                   f"{ra:.2f} m/h | {(ws if all_steady else wa)[0]:.0f}–{(ws if all_steady else wa)[1]:.0f} s |")
    out.append("")
    if np.all(np.diff(v) > 0):
        lo = float(np.interp(1.3, v, T)) if v[0] <= 1.3 <= v[-1] else None
        hi = float(np.interp(1.6, v, T)) if v[0] <= 1.6 <= v[-1] else None
        lo_s = f"{round(lo)} K" if lo is not None else ("at or below %.0f K" % T[0] if 1.3 < v[0]
                                                        else "above %.0f K" % T[-1])
        hi_s = f"{round(hi)} K" if hi is not None else ("at or below %.0f K" % T[0] if 1.6 < v[0]
                                                        else "above %.0f K" % T[-1])
        out.append(f"Meier's 1.3–1.6 m/h corresponds to T_nozzle {lo_s} to {hi_s} on the {basis} basis "
                   "(linear interpolation between the swept points). This is an inference conditional on "
                   "J-M, n = 1, free-surface dilution and core 8; it is never reused as an input.")
    else:
        out.append("ROP is not monotone in T_nozzle: no T range reported")
    return out


def dry_run(pf):
    rows, info = score(pf)
    name = os.path.basename(pf)
    L = ["## 2. Dry run: D2c scoring on D2b's B_JM_A2 (deliverable 10)", "",
         "The D2b scored run, re-scored with the D2c rules. No D2b table is rewritten; the D2b",
         "configuration has no clearance and no free-surface dilution.", ""]
    L += table(name, rows, info)
    zz = [10, 25, 50, 75, 100, 150, 200, 250]
    L += ["", "Profile at the window end (Ø in mm; model azimuth-mean / max-r / Fig. 8.8 width):", "",
          "| depth | " + " | ".join(f"{z} mm" for z in zz) + " |", "|---|" + "---|" * len(zz)]
    fig = {z: w for z, w, _ in load_fig()}
    L.append("| model (az-mean) | " + " | ".join(
        f"{float(np.interp(z * 1e-3, info['zg'], info['D_az'])):.0f}" if z * 1e-3 <= info["drilled"] else "—"
        for z in zz) + " |")
    L.append("| model (max r) | " + " | ".join(
        f"{float(np.interp(z * 1e-3, info['zg'], info['D_mx'])):.0f}" if z * 1e-3 <= info["drilled"] else "—"
        for z in zz) + " |")
    L.append("| Fig. 8.8 | " + " | ".join(f"{fig[z]:.0f}" if z in fig else "—" for z in zz) + " |")
    L += ["", "Reading:",
          f"- **Steady window:** none, as expected. The D2b run hit the bottom with s_c still drifting "
          f"({info['r']['sc_drift'] * 1e3:+.2f} mm/s), so every scored row reads FAIL (no window). "
          "The numbers are shown over the fallback window.",
          f"- **Whole-excavation volume:** {info['vdot'] * 1e6:.2f} cm³/s, "
          f"{info['vdot'] * 1e6 / 2.96:.1f}× the top of the band. D2b's in-Ø 93 score was "
          f"{info['r']['vdot_hole'] * 1e6:.2f} cm³/s and is now only reported.",
          f"- **Profile metric:** the mouth reads {info['mouth'][25.0][0]:.0f} / {info['mouth'][50.0][0]:.0f} mm "
          "at 25 / 50 mm against 92 / 88 mm (Fig. 8.8), and the RMS over the drilled depth is "
          f"{info['rms']:.1f} mm. This is the D2b funnel, the target of the D2c free-surface closure.",
          f"- **Drilled-depth Ø:** depth-mean {info['dmean']:.1f} mm, minimum {info['dmin']:.1f} mm, over "
          f"0–{info['drilled'] * 1e3:.0f} mm (nozzle-plane depth at the window end). The minimum applies the "
          "2 mm rule.",
          f"- **Profile CSV:** `output/{name}_fig88_profile.csv` in this study; no D2b output was changed.", ""]
    txt = open(RESULTS).read()
    a, b = "<!-- DRYRUN -->", "<!-- PROBES -->"
    head, rest = txt.split(a, 1)
    tail = rest.split(b, 1)[1]
    open(RESULTS, "w").write(head + a + "\n\n" + "\n".join(L) + "\n" + b + tail)
    print("\n".join(L))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="*")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--mesh", nargs=2, metavar=("RUN_2MM", "RUN_1MM"))
    ap.add_argument("--sweep", nargs="+")
    ap.add_argument("--p6", nargs="+")
    ap.add_argument("--window-a", nargs="+", dest="window_a")
    args = ap.parse_args()
    if args.p6:
        print("\n".join(p6([os.path.abspath(x) for x in args.p6])))
        return
    if args.window_a:
        L, _ = window_a_table([os.path.abspath(x) for x in args.window_a])
        print("\n".join(L))
        return
    if args.mesh:
        print(mesh(*args.mesh))
        return
    if args.sweep:
        print("\n".join(sweep(args.sweep)))
        return
    if args.dry_run:
        if len(args.runs) != 1:
            sys.exit("--dry-run takes one run")
        dry_run(os.path.abspath(args.runs[0]))
        return
    for pf in args.runs:
        rows, info = score(os.path.abspath(pf))
        print("\n".join(table(os.path.basename(pf), rows, info)))
        print()


if __name__ == "__main__":
    main()
