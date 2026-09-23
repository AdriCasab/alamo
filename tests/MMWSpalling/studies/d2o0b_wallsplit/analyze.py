#!/usr/bin/env python3
"""D2o-0b analysis -> analysis.md + d2o0b_profiles.png + d2o0b_vprofile.png.

Copied and fixed from studies/d2o0_wallsplit/analyze.py (frozen, not edited).

Estimators reused by path, so every number matches the study it is compared to:
  d2l_scan/analyze.py   mrun / az_rwall / vdot / common_feet_depth
  d2m_feet/analyze.py   vmetrics / binned / cross / CEN (D2j-0's v(r), ABSOLUTE)

THE FIX that matters here: **the wall band is defined on GEOMETRY**, not on `s`.
D2o-0 used `band = (s > 0) & (s <= 45 mm)` with no radius limit and no
requirement that the column be below the original rock surface; band counts came
out 239 / 3298 / 400 across three otherwise-identical variants purely according
to whether the nozzle sat a few mm above or below the original surface at the
snapshot, and one case's "wall" was 3298 columns of untouched rock. Here a wall
column must satisfy ALL of:
    r > R_CORNER, face at least one cell BELOW the original rock surface,
    and within BAND_H below the nozzle plane.
"""
import argparse
import importlib.util
import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


M = load("d2m_analyze", os.path.join(STUD, "d2m_feet", "analyze.py"))
AZ = M.AZ
R = load("d2o0b_run", os.path.join(HERE, "run.py"))
OUT = R.OUT
CTRL_DIR = R.CTRL_DIR

WIN = (150.0, 250.0)
Z_SHAFT = 0.100                 # P1/P2 scored below this depth
FEET_MIN = 0.150                # a case shallower than this is unscoreable
BAND_H = 0.045                  # band height below the nozzle plane
P1_WHOLE = (80.0, 100.0)
P1_MEAN = (85.0, 95.0)
P2_MEAN = (85.0, 100.0)
P2_MOUTH = 110.0
P3_VDOT = (2.0, 3.2)
P4_T = (600.0, 800.0)
P4_P = (1500.0, 3500.0)
T_FIRE, P5_NEAR = 821.0, 50.0
Q1_STALL, Q1_ROP_TOL = 10.0, 0.30
CTRL_ROP = 1.491
L_TRUNC = [None]                # t_end of a deliberately truncated 1 mm leg, if any
ORDER = ["W_2mm", "Wz_2mm", "Wsens_2mm", "W_1mm"]


def base_row(run, z_m):
    """No jet_* columns exist under jet_closure = none; absorbed power comes
    from patch_P_robin, computed identically for control and test."""
    th = run.th
    sel = (th["time"] >= WIN[0]) & (th["time"] <= WIN[1])
    d = dict(rop=-float(np.polyfit(th["time"][sel], th["nozzle_z"][sel], 1)[0]) * 3600.0,
             p_face=float(np.mean(th["patch_P_robin"][sel])),
             ledger=float(np.max(np.abs(th["ledger_err"][sel]))),
             t_end=float(th["time"][-1]),
             stall_end=float(th["foot_stall_time"][-1]))
    v = AZ.vdot(run, WIN)
    d["vdot_cm3"] = v * 1e6
    fz = run.lz - th["foot_z"]
    d["feet_end"] = float(fz[-1])
    i = int(np.searchsorted(fz, z_m))
    d["t_match"] = float(th["time"][min(i, len(fz) - 1)])
    d["reached"] = bool(fz[-1] >= z_m - 1e-9)
    dep = run.depth(d["t_match"])
    zs = np.arange(0.01, z_m + 1e-9, 0.01)
    d["prof_z"], d["prof"] = zs, np.array(
        [2.0 * AZ.az_rwall(run, dep, z) * 1e3 for z in zs])
    d["dia_mean"] = float(d["prof"].mean())
    d["dia_min"] = float(d["prof"].min())
    d["mouth"] = float(d["prof"][0])
    deep = zs >= Z_SHAFT
    d["has_shaft"] = bool(deep.any())
    d["shaft_min"] = float(d["prof"][deep].min()) if deep.any() else float("nan")
    d["shaft_mean"] = float(d["prof"][deep].mean()) if deep.any() else float("nan")
    d["shaft_max"] = float(d["prof"][deep].max()) if deep.any() else float("nan")
    return d


def mrun_any(out, name):
    """AZ.mrun, but also able to read a leg that was stopped early. A killed run
    writes no .done, so its metadata is rebuilt from run.py's own job() rather
    than fabricating a completion record."""
    run = AZ.mrun(out, name)
    if run is not None:
        return run, False
    pf = os.path.join(out, name)
    if not os.path.exists(os.path.join(pf, "thermo.dat")):
        return None, False
    _, _, meta = R.job(name)
    return M.AZ.AN.Run(pf, meta), True


def case_row(name, out, z_m):
    run, trunc = mrun_any(out, name)
    if run is None:
        return None
    d = base_row(run, z_m)
    d.update(run=run, case=name, out=out, dirname=name, truncated=trunc)
    d["v"] = M.vmetrics(run, (max(0.0, d["t_match"] - 100.0), d["t_match"]))
    return d


def wall_state(d):
    """P4 on a GEOMETRIC band (see the module docstring)."""
    run = d["run"]
    case = d["dirname"].rsplit("_", 1)[0]
    st = R.surface_state(os.path.join(d["out"], d["dirname"]), 1e9)
    if st is None:
        return None
    t, Ts, zf = st
    th = run.th
    j = int(np.argmin(np.abs(th["time"] - t)))
    zn = th["nozzle_z"][j]
    s = zn - zf
    ok = np.isfinite(s) & np.isfinite(Ts)
    band = (ok & (run.r > R.R_CORNER)                       # outside the skirt
            & (zf <= run.lz - run.dz + 1e-12)               # at least one cell below the top
            & (s > 0.0) & (s <= BAND_H))                    # within the band height
    if not band.any():
        return dict(t=t, n=0)
    h = R.h_py(case, run.r[band], s[band])
    Tg = R.T_py(run.r[band])
    q = h * (Tg - Ts[band])
    P = 4.0 * float(np.sum(q)) * run.dz * run.dz
    depth = run.lz - zf[band]
    return dict(t=t, n=int(band.sum()), P=P, depth=depth, T=Ts[band],
                Tmin=float(Ts[band].min()), Tmax=float(Ts[band].max()),
                Tmean=float(Ts[band].mean()), near_fire=float(T_FIRE - Ts[band].max()),
                prof=[(lo, float(Ts[band][(depth >= lo) & (depth < lo + 0.05)].mean())
                       if ((depth >= lo) & (depth < lo + 0.05)).any() else float("nan"))
                      for lo in np.arange(0.0, 0.40, 0.05)])


def main():
    argparse.ArgumentParser().parse_args()
    names = [n for n in ORDER if os.path.exists(os.path.join(OUT, n + ".done"))]
    if not names:
        raise SystemExit("no completed runs in " + OUT)
    two = [n for n in names if n.endswith("_2mm")]
    ctl = AZ.mrun(CTRL_DIR, "LI0_2mm")
    z_m = AZ.common_feet_depth([AZ.mrun(OUT, n) for n in two] + [ctl])
    # the 1 mm leg may have been stopped early on purpose; it is read anyway
    if not os.path.exists(os.path.join(OUT, "W_1mm.done")):
        names = names + (["W_1mm"] if os.path.exists(os.path.join(OUT, "W_1mm", "thermo.dat"))
                         else [])
    rows = [r for r in (case_row(n, OUT, z_m) for n in two) if r]
    cr = case_row("LI0_2mm", CTRL_DIR, z_m)
    cr["case"] = "control"

    L = ["# D2o-0b analysis — prescribed floor/wall heating, field bounded in radius", "",
         f"Binary `{R.J.binary_sha()[:8]}…`; predictions hashed before launch (RESULTS §0). "
         f"Compared at **matched feet depth {z_m*1e3:.0f} mm**, never matched time. "
         "**The control's depth-mean is recomputed at that depth, not taken as the packet's "
         "110.5** (PREDICTIONS §0: it is 112.3 at 97 mm and 100.7 over the whole hole).", "",
         "## 1. Q1 — the gate: did the radial corner stop the jam?", "",
         "| case | h_wall | ROP 150–250 s [m/h] | vs control | feet depth at 600 s | "
         "`foot_stall_time` end [s] | Q1 |", "|---|---|---|---|---|---|---|"]
    for d in [cr] + rows:
        nm = d["case"]
        hw = "—" if nm == "control" else f"{0.0 if nm.startswith('Wz') else (R.H_WALL_SENS if 'sens' in nm else R.H_WALL):g}"
        q1 = (d["stall_end"] < Q1_STALL and d["feet_end"] >= FEET_MIN
              and abs(d["rop"] / CTRL_ROP - 1.0) <= Q1_ROP_TOL)
        L.append(f"| {nm} | {hw} | **{d['rop']:.3f}** | "
                 f"{'—' if nm == 'control' else f'{d[chr(114)+chr(111)+chr(112)]/cr[chr(114)+chr(111)+chr(112)]:.2f}×'} | "
                 f"{d['feet_end']*1e3:.0f} mm | {d['stall_end']:.0f} | "
                 f"{'—' if nm == 'control' else ('**PASS**' if q1 else '**FAIL**')} |")
    L += ["", f"Q1 requires `foot_stall_time` < {Q1_STALL:g} s, feet depth ≥ {FEET_MIN*1e3:.0f} mm "
          f"by 600 s, and ROP within {Q1_ROP_TOL*100:.0f} % of {CTRL_ROP:g} m/h. "
          "D2o-0 (stand-off corner) gave ROP 0.108 m/h and stall 219 s.", ""]

    L += ["## 2. Shape — P1 and P2", "",
          "| case | shaft Ø min/mean/max (z>100mm) | depth-mean Ø | mouth Ø | V̇ [cm³/s] | "
          "P1 whole | P1 mean | P2 mean | P2 mouth | P3 |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for d in [cr] + rows:
        nm = d["case"]
        if not d["has_shaft"] or not d["reached"]:
            sh = "— unscoreable"
            p1w = p1m = "—"
        else:
            sh = f"{d['shaft_min']:.1f} / **{d['shaft_mean']:.1f}** / {d['shaft_max']:.1f}"
            p1w = "PASS" if (d["shaft_min"] >= P1_WHOLE[0] and d["shaft_max"] <= P1_WHOLE[1]) else "FAIL"
            p1m = "PASS" if P1_MEAN[0] <= d["shaft_mean"] <= P1_MEAN[1] else "FAIL"
        p2m = "PASS" if P2_MEAN[0] <= d["dia_mean"] <= P2_MEAN[1] else "FAIL"
        p2o = "PASS" if d["mouth"] <= P2_MOUTH else "FAIL"
        p3 = "PASS" if P3_VDOT[0] <= d["vdot_cm3"] <= P3_VDOT[1] else "FAIL"
        L.append(f"| {nm} | {sh} | **{d['dia_mean']:.1f}** | **{d['mouth']:.1f}** | "
                 f"**{d['vdot_cm3']:.2f}** | {p1w} | {p1m} | {p2m} | {p2o} | {p3} |")
    L += ["", f"P1 whole profile {P1_WHOLE[0]:.0f}–{P1_WHOLE[1]:.0f} mm and mean "
          f"{P1_MEAN[0]:.0f}–{P1_MEAN[1]:.0f}; P2 depth-mean {P2_MEAN[0]:.0f}–{P2_MEAN[1]:.0f} and "
          f"mouth ≤ {P2_MOUTH:.0f}; P3 V̇ {P3_VDOT[0]}–{P3_VDOT[1]} cm³/s (Meier 2.47, "
          "control 4.58).", "",
          "## 3. Ø(z) at matched geometry [mm], 20 mm intervals", "",
          "| case | " + " | ".join(f"{z*1e3:.0f}" for z in rows[0]["prof_z"][1::2]) + " |",
          "|---|" + "---|" * len(rows[0]["prof_z"][1::2])]
    for d in [cr] + rows:
        L.append(f"| {d['case']} | " + " | ".join(f"{v:.0f}" for v in d["prof"][1::2]) + " |")

    L += ["", "## 4. v(r) — absolute, 2 mm annuli (never normalised)", "",
          "| case | v_c [m/h] | " + " | ".join(f"{r} mm" for r in (25, 35, 45, 50, 55, 60)) + " |",
          "|---|---|---|---|---|---|---|---|"]
    for d in [cr] + rows:
        p = d["v"]["prof"]
        L.append(f"| {d['case']} | {d['v']['v_c']:.3f} | " + " | ".join(
            f"{p[int(np.argmin(np.abs(M.CEN*1e3-r)))]:.3f}" for r in (25, 35, 45, 50, 55, 60))
            + " |")

    L += ["", "## 5. Wall — P4, P5, on a GEOMETRIC band", "",
          f"A wall column has r > {R.R_CORNER*1e3:.0f} mm, its face **at least one cell below the "
          f"original rock surface**, and 0 < s ≤ {BAND_H*1e3:.0f} mm. **Not `s` alone** — that is "
          "what put 3298 columns of untouched rock into D2o-0's P3.", "",
          "| case | t [s] | band cols | T_wall min/mean/max [K] | 821 − T_max | wall heat [W] | "
          "P4 T | P4 P | P5 |", "|---|---|---|---|---|---|---|---|---|"]
    ws = {}
    for d in rows + [cr]:
        w = wall_state(d)
        ws[d["case"]] = w
        if w is None or w.get("n", 0) == 0:
            L.append(f"| {d['case']} | — | 0 | — | — | — | — | — | — |")
            continue
        ctlrow = d["case"] == "control"
        pw = "— (impinging field)" if ctlrow else f"**{w['P']:.0f}**"
        L.append(f"| {d['case']} | {w['t']:.0f} | {w['n']} | "
                 f"{w['Tmin']:.0f} / **{w['Tmean']:.0f}** / {w['Tmax']:.0f} | "
                 f"{w['near_fire']:+.0f} K | {pw} | "
                 f"{'PASS' if P4_T[0] <= w['Tmean'] <= P4_T[1] else 'FAIL'} | "
                 f"{'—' if ctlrow else ('PASS' if P4_P[0] <= w['P'] <= P4_P[1] else 'FAIL')} | "
                 f"{'**FIRES**' if w['near_fire'] <= P5_NEAR else 'no'} |")
    L += ["", "Wall temperature by depth [K], 50 mm bins:", "",
          "| case | " + " | ".join(f"{z*1e3:.0f}–{z*1e3+50:.0f}" for z in np.arange(0, 0.40, 0.05))
          + " |", "|---|" + "---|" * 8]
    for d in rows + [cr]:
        w = ws.get(d["case"])
        if w and w.get("n", 0):
            L.append(f"| {d['case']} | " + " | ".join(
                "—" if not np.isfinite(v) else f"{v:.0f}" for _, v in w["prof"]) + " |")

    # The 1 mm leg was stopped early on purpose and does not reach the 150-250 s
    # scoring window, so it yields no row. Its single usable datum (the corner
    # step at the common time t = 50 s) is reported in RESULTS.md, not here.
    d1 = None
    w1 = os.path.join(OUT, "W_1mm", "thermo.dat")
    if os.path.exists(w1):
        t1 = float(np.loadtxt(w1, skiprows=1, ndmin=2)[-1, 0])
        if t1 >= WIN[1]:
            d1 = case_row("W_1mm", OUT, z_m)
        else:
            L_TRUNC[0] = t1
    d2 = next((d for d in rows if d["case"] == "W_2mm"), None)
    L += ["", "## 6. Mesh — 2 mm vs 1 mm at matched feet depth", ""]
    if d1 and d2:
        if d1.get("truncated"):
            L += [f"**The 1 mm leg was stopped at t = {d1['t_end']:.0f} s by design** (the scoring "
                  "window is 150–250 s and shape is unscoreable because Q1 failed); it is a "
                  "mechanism test, not a shape run. See RESULTS §0.", ""]
        L += ["| quantity | 2 mm | 1 mm | gap |", "|---|---|---|---|",
              f"| ROP 150–250 s | {d2['rop']:.3f} | {d1['rop']:.3f} | "
              f"**{(d1['rop']/d2['rop']-1)*100:+.1f} %** |",
              f"| depth-mean Ø | {d2['dia_mean']:.1f} | {d1['dia_mean']:.1f} | "
              f"{d1['dia_mean']-d2['dia_mean']:+.1f} mm |",
              f"| V̇ | {d2['vdot_cm3']:.2f} | {d1['vdot_cm3']:.2f} | "
              f"{(d1['vdot_cm3']/d2['vdot_cm3']-1)*100:+.1f} % |"]
    elif L_TRUNC[0]:
        L += [f"**The 1 mm leg was stopped at t = {L_TRUNC[0]:.1f} s of 600** and never reached "
              "the 150–250 s scoring window, so it produces no row here. It was ended "
              "deliberately: Q1 already fails at 2 mm and refinement can only make the step-face "
              "sink worse (it scales as 1/dx), the packet requires both meshes only for a "
              "**positive** claim, and shape is unscoreable regardless. Its one usable datum — the "
              "corner step at the common time t = 50 s — is in RESULTS.md §2. "
              "**This study is single-mesh and nothing in it is converged.**"]
    else:
        L += ["_1 mm leg absent — **single mesh is not quotable** and a positive result "
              "may not be claimed without it._"]

    L += ["", "## 7. Ledger and absorbed power", "",
          "| case | max abs ledger_err | patch_P_robin [W] |", "|---|---|---|"]
    for d in [cr] + rows + ([d1] if d1 else []):
        L.append(f"| {d['case']} | {d['ledger']:.1e} | {d['p_face']:.0f} |")

    txt = "\n".join(L)
    open(os.path.join(HERE, "analysis.md"), "w").write(txt + "\n")
    print(txt)
    figs(cr, rows, d1, ws, z_m)


def figs(cr, rows, d1, ws, z_m):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    allr = [cr] + rows + ([d1] if d1 else [])
    f, ax = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for d in allr:
        ax[0].plot(d["prof"], d["prof_z"] * 1e3, "-o", ms=3, label=d["case"])
    ax[0].axvspan(85, 95, color="g", alpha=0.15, label="Meier shaft 85–95 mm")
    ax[0].axvline(80, color="k", ls="--", lw=1, label="skirt Ø 80 mm")
    ax[0].axhline(Z_SHAFT * 1e3, color="gray", ls=":", lw=1)
    ax[0].invert_yaxis()
    ax[0].set(xlabel="hole Ø [mm]", ylabel="depth below original top [mm]",
              title=f"Ø(z) at matched feet depth {z_m*1e3:.0f} mm")
    ax[0].grid(alpha=0.3)
    ax[0].legend(fontsize=8)
    for d in allr:
        w = ws.get(d["case"])
        if w and w.get("n", 0):
            zz = [z * 1e3 + 25 for z, v in w["prof"] if np.isfinite(v)]
            tt = [v for _, v in w["prof"] if np.isfinite(v)]
            ax[1].plot(tt, zz, "-o", ms=4, label=d["case"])
    ax[1].axvspan(*P4_T, color="g", alpha=0.15, label="P4 600–800 K")
    ax[1].axvline(T_FIRE, color="r", ls="--", lw=1.2, label="spall 821 K")
    ax[1].axvline(570, color="C1", ls=":", lw=1.2, label="alteration ≈570 K")
    ax[1].invert_yaxis()
    ax[1].set(xlabel="wall surface temperature [K]", ylabel="depth [mm]",
              title="Wall temperature by depth (geometric band)")
    ax[1].grid(alpha=0.3)
    ax[1].legend(fontsize=8)
    f.suptitle("D2o-0b: radial floor/wall corner, field bounded in radius", y=0.99)
    f.tight_layout()
    f.savefig(os.path.join(HERE, "d2o0b_profiles.png"), dpi=110)
    plt.close(f)

    f2, a2 = plt.subplots(figsize=(7.2, 4.8))
    a2.axvspan(42.5, 46.5, color="g", alpha=0.15, label="Meier implied floor half-width")
    for d in allr:
        a2.plot(M.CEN * 1e3, d["v"]["prof"], "-", lw=1.8, label=d["case"])
    a2.axvline(R.R_CORNER * 1e3, color="k", ls="--", lw=1, label="radial corner r = 40 mm")
    a2.set(xlabel="radius [mm]", ylabel="recession rate [m/h]", xlim=(0, 70),
           title="v(r), absolute — does the wall erode at all?")
    a2.grid(alpha=0.3)
    a2.legend(fontsize=8)
    f2.tight_layout()
    f2.savefig(os.path.join(HERE, "d2o0b_vprofile.png"), dpi=110)
    plt.close(f2)
    print("wrote d2o0b_profiles.png, d2o0b_vprofile.png")


if __name__ == "__main__":
    main()
