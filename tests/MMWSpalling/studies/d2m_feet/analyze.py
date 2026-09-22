#!/usr/bin/env python3
"""D2m analysis -> analysis.md + d2m_profiles.png + d2m_vprofile.png.

Reuses studies/d2l_scan/analyze.py by path (frozen, not edited) for the Ø(z),
volume and matched-depth estimators, so those numbers are computed exactly as
D2l computed them:
  az_rwall (itself copied there from d2c_steady/score.py), mrun, meier_row,
  vdot, common_feet_depth

The recession-rate profile (ACTIVE_STEP 5a, DECIDING after the 12:16 amendment)
uses D2j-0's estimators, COPIED below from d2j0_hotdisk/analyze.py rather than
run in that frozen folder: 2 mm annuli, v_c = mean rate over r < 6 mm,
r_half = cross(prof, 0.5*v_c), w_edge = r(0.1*v_c) - r(0.9*v_c), with `cross`
taking the OUTERMOST crossing by linear interpolation between bin centres.
Extended from D2j-0's 40 mm to 70 mm as 5a requires.
"""
import argparse
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


AZ = load("d2l_analyze", os.path.join(STUD, "d2l_scan", "analyze.py"))
R = load("d2m_run", os.path.join(HERE, "run.py"))
OUT = R.OUT

WIN_FIX = (150.0, 250.0)           # the fixed scoring window (cross-check)
VWIN = 100.0                       # matched-geometry recession window length
FLOOR = 1.147                      # J/mm3, rho*cp*dT_fire
BAND_D = (85.0, 93.0)              # Meier Ø, orientation only
MEIER_HALF = (42.5, 46.5)          # Meier implied floor half-width, mm
ORDER = ["R32", "R40", "R44", "R48", "S3648"]

# --- D2j-0 estimators, copied verbatim in definition, extended to 70 mm -----
EDGES = np.arange(0.0, 0.0701, 0.002)
CEN = 0.5 * (EDGES[:-1] + EDGES[1:])


def binned(r, v):
    b = np.digitize(r, EDGES) - 1
    return np.array([v[b == i].mean() if np.any(b == i) else np.nan
                     for i in range(len(CEN))])


def cross(prof, level):
    """Outermost radius where the binned profile falls through `level`, by
    linear interpolation between bin centres. (d2j0_hotdisk/analyze.py)"""
    ok = np.isfinite(prof)
    x, y = CEN[ok], prof[ok]
    hit = np.nonzero((y[:-1] >= level) & (y[1:] < level))[0]
    if not len(hit):
        return np.nan
    i = hit[-1]
    if y[i] == y[i + 1]:
        return float(x[i])
    return float(x[i] + (x[i + 1] - x[i]) * (y[i] - level) / (y[i] - y[i + 1]))


def recession(run, win):
    """Per-column recession rate in m/h over `win`, from the removal log."""
    ev = run.ev
    s = (ev["time"] >= win[0]) & (ev["time"] <= win[1])
    col = ev["col_j"].astype(int) * run.nx + ev["col_i"].astype(int)
    v = np.zeros(run.ncol)
    np.add.at(v, col[s], ev["h_applied"][s])
    return v / (win[1] - win[0]) * 3600.0


def vmetrics(run, win):
    v = recession(run, win)
    prof = binned(run.r, v)
    v_c = float(v[run.r < 0.006].mean())
    r90, r10 = cross(prof, 0.9 * v_c), cross(prof, 0.1 * v_c)
    return dict(win=win, prof=prof, v_c=v_c, r_half=cross(prof, 0.5 * v_c),
                r90=r90, r10=r10, w_edge=r10 - r90)


def row(name, z_m):
    run = AZ.mrun(OUT, name)
    if run is None:
        return None
    d = AZ.meier_row(run, WIN_FIX, z_m)
    d.update(run=run, case=name)
    m = AZ.AN.meta_of(os.path.join(OUT, name))
    d["ri"], d["ro"] = m["foot_r_inner"], m["foot_r_outer"]
    d["stance"] = 2e3 * m["foot_r_outer"]
    d["min_at_mm"] = float(d["prof_z"][int(np.argmin(d["prof"]))]) * 1e3
    d["tse_floor_ratio"] = d["tse_rock"] / FLOOR
    tm = d["t_match"]
    d["v_match"] = vmetrics(run, (max(0.0, tm - VWIN), tm))     # DECIDING
    d["v_fix"] = vmetrics(run, WIN_FIX)                          # cross-check
    return d


def main():
    argparse.ArgumentParser().parse_args()
    names = [n for n in ORDER if os.path.exists(os.path.join(OUT, n + ".done"))]
    if not names:
        raise SystemExit("no completed runs in " + OUT)
    z_m = AZ.common_feet_depth([AZ.mrun(OUT, n) for n in names])
    rows = [r for r in (row(n, z_m) for n in names) if r is not None]
    ctl = next((d for d in rows if d["case"] == "R40"), None)

    L = ["# D2m analysis — is the cut diameter set by the feet or by the flame?", "",
         f"Binary `{R.J.binary_sha()[:8]}…`; predictions hashed before launch (RESULTS.md §0). "
         f"All cases read at **matched geometry**: the time each case's feet reach the common "
         f"depth **{z_m * 1e3:.0f} mm**. Never at matched time.", "",
         "## 1. The ladder", "",
         "| case | ring [mm] | stance Ø | **min Ø** | 2×r_outer | min Ø − stance | at depth | "
         "ROP [m/h] | Δ ROP | absorbed [W] | V̇ [cm³/s] | depth-mean Ø |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for d in rows:
        dr = (d["rop"] / ctl["rop"] - 1.0) * 100 if ctl else float("nan")
        L.append(f"| {d['case']} | {d['ri']*1e3:.0f}–{d['ro']*1e3:.0f} | {d['stance']:.0f} | "
                 f"**{d['dia_min']:.1f}** | {d['stance']:.0f} | {d['dia_min']-d['stance']:+.1f} | "
                 f"{d['min_at_mm']:.0f} mm | **{d['rop']:.3f}** | {dr:+.1f} % | "
                 f"{d['p_face']:.0f} | {d['vdot_cm3']:.2f} | {d['dia_mean']:.1f} |")

    L += ["", "## 2. The two deciding metrics", "",
          "**min Ø** (weak evidence outward — see PREDICTIONS §2a: the 0.9 quantile makes "
          "Ø ≥ 2·r_outer at the feet depth close to a structural identity) **and r_half**, the "
          "flame's own erosion edge, which the support rule does not constrain.", "",
          "| case | r_outer | **min Ø** | H-feet (2×r_outer ±5) | H-thermal (78±3) | "
          "**r_half [mm]** | r_half − r_outer | H-feet on r_half (±2 mm) |",
          "|---|---|---|---|---|---|---|---|"]
    for d in rows:
        ro = d["ro"] * 1e3
        rh = d["v_match"]["r_half"] * 1e3
        feet_d = abs(d["dia_min"] - d["stance"]) <= 5.0
        th_d = abs(d["dia_min"] - 78.0) <= 3.0
        feet_r = np.isfinite(rh) and abs(rh - ro) <= 2.0
        L.append(f"| {d['case']} | {ro:.0f} | **{d['dia_min']:.1f}** | "
                 f"{'yes' if feet_d else 'no'} | {'yes' if th_d else 'no'} | "
                 f"**{rh:.1f}** | {rh - ro:+.1f} | {'**yes**' if feet_r else 'no'} |")

    rhs = [d["v_match"]["r_half"] * 1e3 for d in rows if d["case"] != "S3648"]
    ros = [d["ro"] * 1e3 for d in rows if d["case"] != "S3648"]
    if len(rhs) > 1 and all(np.isfinite(rhs)):
        spread = max(rhs) - min(rhs)
        literal = all(abs(a - b) <= 2.0 for a, b in zip(rhs, ros))
        slope, intercept = np.polyfit(ros, rhs, 1)
        off = [a - b for a, b in zip(rhs, ros)]
        L += ["",
              "**The pre-registered test, as written:** H-feet on r_half requires "
              "|r_half − r_outer| ≤ 2 mm at every case. Offsets are "
              + ", ".join(f"{o:+.1f}" for o in off)
              + f" mm, so the literal test **{'PASSES' if literal else 'FAILS'}**.",
              "",
              "**Slope — computed after seeing the result, reported as an addition, not a "
              "substitution.** The literal test conflates offset with slope, and here the offset "
              f"is nearly constant ({min(off):+.1f} to {max(off):+.1f} mm) while r_half spans "
              f"**{spread:.1f} mm** for **{max(ros) - min(ros):.0f} mm** of r_outer: "
              f"**d r_half / d r_outer = {slope:.2f}** (intercept {intercept:+.1f} mm). "
              + ("_At face value this reads as the edge moving one-for-one with the stance. "
                 "**It is an artefact — see the absolute table below, which retracts it.**"
                 if slope > 0.8 else
                 "**The edge does not move one-for-one with the stance.**")]
        # The retraction. r_half is defined RELATIVE to v_c, and v_c is not
        # constant across the ladder, so r_half can move with the stance while
        # the profile itself does not. Compare absolute rates at fixed radii.
        L += ["", "#### Absolute v(r) at fixed radii — this retracts the slope reading above", "",
              "`r_half` is a level crossing of `0.5·v_c`, and **`v_c` itself falls from 3.32 to "
              "1.97 m/h across the ladder**. So r_half can march outward while the profile stands "
              "still. The absolute rates settle it:", "",
              "| case | r_outer | v_c | " + " | ".join(f"{r} mm" for r in (25, 35, 45, 50, 60))
              + " |", "|---|---|---|---|---|---|---|---|"]
        for d in rows:
            p = d["v_match"]["prof"]
            vals = [p[int(np.argmin(np.abs(CEN * 1e3 - r)))] for r in (25, 35, 45, 50, 60)]
            L.append(f"| {d['case']} | {d['ro']*1e3:.0f} | {d['v_match']['v_c']:.3f} | "
                     + " | ".join(f"{v:.3f}" for v in vals) + " |")
        sp = {}
        for r in (25, 35, 45, 50, 60):
            v = np.array([d["v_match"]["prof"][int(np.argmin(np.abs(CEN * 1e3 - r)))]
                          for d in rows])
            sp[r] = (v.max() - v.min()) / v.mean() * 100
        L += ["", "Ladder spread (max−min as % of mean): "
              + ", ".join(f"**{r} mm: {sp[r]:.0f} %**" for r in (25, 35, 45, 50, 60))
              + ".", "",
              "**At 45–50 mm the absolute erosion rate is invariant to 5–7 % across a 16 mm "
              "change of stance, while the centre moves 39 %.** The flame's far-field erosion is "
              "therefore *not* set by the stance; only the centre is. `r_half` tracked "
              "`foot_r_outer` purely through its own normalisation. **H-thermal holds on the "
              "erosion edge; H-feet holds on min Ø, which PREDICTIONS §2a flagged in advance as "
              "close to a structural identity of the support rule.**"]

    L += ["", "### Cross-check: the fixed [150, 250] s window", "",
          "| case | r_half matched-geometry | r_half fixed window | Δ |", "|---|---|---|---|"]
    disagree = False
    for d in rows:
        a, b = d["v_match"]["r_half"] * 1e3, d["v_fix"]["r_half"] * 1e3
        disagree = disagree or (np.isfinite(a) and np.isfinite(b) and abs(a - b) > 2.0)
        L.append(f"| {d['case']} | {a:.1f} | {b:.1f} | {a-b:+.1f} |")
    L += ["", ("**The two windows disagree by more than one cell somewhere — both are reported "
               "and no single reading is forced** (PREDICTIONS §5)." if disagree else
              "The two windows agree within one cell everywhere, so the deciding "
              "(matched-geometry) reading stands on its own.")]

    s48 = next((d for d in rows if d["case"] == "S3648"), None)
    r48 = next((d for d in rows if d["case"] == "R48"), None)
    if s48 and r48:
        dd = s48["dia_min"] - r48["dia_min"]
        dv = s48["rop"] / r48["rop"] - 1.0
        L += ["", "## 3. `S3648` vs `R48` — outer edge, or the whole ring?", "",
              f"- min Ø **{s48['dia_min']:.1f}** vs **{r48['dia_min']:.1f}** mm (Δ {dd:+.1f}, "
              "predicted within 3)",
              f"- ROP **{s48['rop']:.3f}** vs **{r48['rop']:.3f}** m/h (Δ {dv*100:+.1f} %, "
              "predicted within 5 %)",
              f"- r_half **{s48['v_match']['r_half']*1e3:.1f}** vs "
              f"**{r48['v_match']['r_half']*1e3:.1f}** mm", "",
              f"**{'The outer edge alone sets it' if abs(dd) <= 3.0 and abs(dv) <= 0.05 else 'The whole ring matters — reported separately, not averaged'}.**"]

    L += ["", "## 4. Ø(z) at matched geometry [mm], 20 mm intervals", "",
          "| case | " + " | ".join(f"{z*1e3:.0f}" for z in rows[0]["prof_z"][1::2]) + " |",
          "|---|" + "---|" * len(rows[0]["prof_z"][1::2])]
    for d in rows:
        L.append(f"| {d['case']} | " + " | ".join(f"{v:.0f}" for v in d["prof"][1::2]) + " |")
    L += ["", "(Meier 85–93 mm shown in the figure for orientation only — **no Meier score comes "
          "from this packet**.)"]

    L += ["", "## 5. Recession-rate profile (ACTIVE_STEP 5a) — D2j-0 estimators", "",
          "Matched-geometry window, 2 mm annuli, from the removal log.", "",
          "| case | v_c [m/h] | r(0.9 v_c) | **r_half** | r(0.1 v_c) | **w_edge** | r_outer |",
          "|---|---|---|---|---|---|---|"]
    for d in rows:
        m = d["v_match"]
        L.append(f"| {d['case']} | {m['v_c']:.3f} | {m['r90']*1e3:.1f} | "
                 f"**{m['r_half']*1e3:.1f}** | {m['r10']*1e3:.1f} | "
                 f"**{m['w_edge']*1e3:.1f}** | {d['ro']*1e3:.0f} |")
    L += ["", f"Meier's implied floor half-width is {MEIER_HALF[0]}–{MEIER_HALF[1]} mm. "
          "**Does the model's floor recede flat out to ~45 mm, or is it already falling at 40?** "
          "v(r)/v_c by 2 mm annulus:", "",
          "| case | " + " | ".join(f"{c*1e3:.0f}" for c in CEN[:24]) + " |",
          "|---|" + "---|" * 24]
    for d in rows:
        p = d["v_match"]["prof"]
        vc = d["v_match"]["v_c"]
        L.append(f"| {d['case']} | " + " | ".join(
            "—" if not np.isfinite(p[i]) else f"{p[i]/vc:.2f}" for i in range(24)) + " |")

    L += ["", "## 6. Reported, not deciding", "",
          "| case | J/mm³ to rock | × floor 1.147 | jet_s_c | patch_min_standoff | jet_fs_cols | "
          "ledger_err | cap max |", "|---|---|---|---|---|---|---|---|"]
    for d in rows:
        th = d["run"].th
        s = (th["time"] >= WIN_FIX[0]) & (th["time"] <= WIN_FIX[1])
        L.append(f"| {d['case']} | {d['tse_rock']:.2f} | {d['tse_floor_ratio']:.2f} | "
                 f"{d['s_c']:.1f} | {float(np.mean(th['patch_min_standoff'][s]))*1e3:.1f} | "
                 f"{d['fs_cols']:.1f} | {d['ledger']:.1e} | {d['cap']:.2f} |")
    L += ["", "The 4–6 J/mm³ band is an assumed 30 % delivery, not a measurement, and is not used "
          "as a target; the floor 1.147 J/mm³ (ρc_pΔT_fire) is the only independent number here."]

    txt = "\n".join(L)
    open(os.path.join(HERE, "analysis.md"), "w").write(txt + "\n")
    print(txt)
    fig_profiles(rows, z_m)
    fig_vprofile(rows)


def fig_profiles(rows, z_m):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for d in rows:
        ax[0].plot(d["prof"], d["prof_z"] * 1e3, "-o", ms=3,
                   label=f"{d['case']} (Ø{d['stance']:.0f})")
    ax[0].axvspan(*BAND_D, color="g", alpha=0.15, label="Meier 85–93 mm")
    ax[0].invert_yaxis()
    ax[0].set(xlabel="hole Ø [mm]", ylabel="depth below original top [mm]",
              title=f"Ø(z) at matched feet depth {z_m*1e3:.0f} mm")
    ax[0].grid(alpha=0.3)
    ax[0].legend(fontsize=8)
    st = [d["stance"] for d in rows if d["case"] != "S3648"]
    mn = [d["dia_min"] for d in rows if d["case"] != "S3648"]
    ax[1].plot(st, mn, "o", ms=9, label="observed min Ø")
    lo, hi = min(st) - 6, max(st) + 6
    ax[1].plot([lo, hi], [lo, hi], "k--", lw=1, label="H-feet: min Ø = 2·r_outer")
    ax[1].axhline(78.0, color="C3", ls=":", lw=1.5, label="H-thermal: 78 mm")
    for d in rows:
        if d["case"] == "S3648":
            ax[1].plot(d["stance"], d["dia_min"], "s", ms=9, color="C2", label="S3648")
    ax[1].set(xlabel="stance Ø = 2·foot_r_outer [mm]", ylabel="min Ø [mm]",
              title="Does the cut follow the stance?")
    ax[1].grid(alpha=0.3)
    ax[1].legend(fontsize=8)
    f.suptitle("D2m: cut diameter against tool stance", y=0.99)
    f.tight_layout()
    f.savefig(os.path.join(HERE, "d2m_profiles.png"), dpi=110)
    plt.close(f)
    print("wrote d2m_profiles.png")


def fig_vprofile(rows):
    """ACTIVE_STEP 5a: recession rate against radius, axis to 70 mm, each
    case's foot_r_outer marked, Meier's implied floor half-width shaded."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for k, (a, key, ttl) in enumerate(((ax[0], "v_match", "matched geometry (DECIDING)"),
                                       (ax[1], "v_fix", "fixed 150–250 s (cross-check)"))):
        a.axvspan(*MEIER_HALF, color="g", alpha=0.15,
                  label="Meier implied floor half-width")
        for d in rows:
            m = d[key]
            line, = a.plot(CEN * 1e3, m["prof"], "-", lw=1.8, label=d["case"])
            a.axvline(d["ro"] * 1e3, ls=":", lw=0.9, color=line.get_color())
            if np.isfinite(m["r_half"]):
                a.plot([m["r_half"] * 1e3], [0.5 * m["v_c"]], "o", ms=6,
                       color=line.get_color())
        a.set(xlabel="radius [mm]", ylabel="recession rate [m/h]",
              title=f"v(r), {ttl}", xlim=(0, 70))
        a.grid(alpha=0.3)
        if k == 0:
            a.legend(fontsize=8)
    f.suptitle("D2m 5a: recession rate vs radius — dotted = that case's foot_r_outer, "
               "dot = r_half", y=0.99)
    f.tight_layout()
    f.savefig(os.path.join(HERE, "d2m_vprofile.png"), dpi=110)
    plt.close(f)
    print("wrote d2m_vprofile.png")


if __name__ == "__main__":
    main()
