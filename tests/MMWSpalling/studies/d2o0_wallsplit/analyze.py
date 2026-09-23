#!/usr/bin/env python3
"""D2o-0 analysis -> analysis.md + d2o0_profiles.png + d2o0_vprofile.png.

Estimators are reused by path from the frozen studies, so every number is
computed the way the study it is compared against computed it:
  d2l_scan/analyze.py   mrun / meier_row / az_rwall / common_feet_depth (Ø(z),
                        V̇, matched-depth) -- itself carrying d2c_steady's az_rwall
  d2m_feet/analyze.py   recession / binned / cross / vmetrics / CEN (D2j-0's
                        v(r) estimators, 2 mm annuli to 70 mm), ABSOLUTE

New here, because no earlier study needed it: per-column surface temperature
from the plotfiles (yt), which gives the wall temperature by depth and the wall
heat. `jet_closure = none` means there is no enthalpy march and therefore no
in-code wall-power diagnostic, so P4 is computed from the prescribed field and
the wall area -- as ACTIVE_STEP item 5 requires it to be, and says so.
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
R = load("d2o0_run", os.path.join(HERE, "run.py"))
OUT = R.OUT
CTRL_DIR = os.path.dirname(R.CONTROL)

WIN = (150.0, 250.0)
BLOCKS = [(150.0, 250.0), (250.0, 350.0), (350.0, 450.0)]
MEIER_SHAFT = (85.0, 95.0)          # P1
SKIRT_MM = 80.0                     # the stance
GAIN_BAND = (2.5, 7.5)              # P2, mm per side
T_WALL_BAND = (600.0, 800.0)        # P3
P_WALL_BAND = (1500.0, 3500.0)      # P4, W
T_FIRE = 821.0
P5_NEAR = 50.0
Z_SHAFT = 0.10                      # P1/P2 scored BELOW this depth (flat-start transient)
ORDER = ["W_2mm", "Wz_2mm", "Wsens_2mm", "W_1mm"]


# ---------------------------------------------------------------- plotfiles
def plotfiles(pf):
    return sorted(d for d in os.listdir(pf) if d.endswith("cell") and d[:-4].isdigit())


def surface_state(pf, want_t=None):
    """(t, T_surface[col], z_face[col]) from the last plotfile at or before
    want_t. Column order matches AN.Run.r: index = j*nx + i."""
    import yt
    ds_names = plotfiles(pf)
    if not ds_names:
        return None
    best = None
    for n in ds_names:
        ds = yt.load(os.path.join(pf, n))
        t = float(ds.current_time)
        if want_t is None or t <= want_t + 1e-6:
            best = (t, ds)
    if best is None:
        return None
    t, ds = best
    nx, ny, nz = (int(v) for v in ds.domain_dimensions)
    dz = float(ds.domain_right_edge[2]) / nz
    g = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    T = np.asarray(g["Temp"])          # (nx, ny, nz)
    rem = np.asarray(g["removed"])
    solid = rem < 0.5
    any_solid = solid.any(axis=2)
    ktop = np.where(any_solid, nz - 1 - np.argmax(solid[:, :, ::-1], axis=2), -1)
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    Ts = np.where(any_solid, T[ii, jj, np.clip(ktop, 0, nz - 1)], np.nan)
    zf = np.where(any_solid, (ktop + 1) * dz, np.nan)
    return t, Ts.T.ravel(), zf.T.ravel()          # -> index j*nx + i


# ---------------------------------------------------------------- per case
def base_row(run, z_m):
    """D2l's meier_row, minus the jet_* columns. `jet_closure = none` registers
    no enthalpy-march diagnostics (no jet_P_face / jet_P_cap / jet_s_c /
    jet_fs_cols), so absorbed power is taken from `patch_P_robin`, which exists
    in both the control and these runs -- computed the same way for both so the
    comparison stays apples-to-apples. Ø(z) and V̇ still come from d2l_scan's
    az_rwall and vdot by path."""
    th = run.th
    sel = (th["time"] >= WIN[0]) & (th["time"] <= WIN[1])
    rop = -float(np.polyfit(th["time"][sel], th["nozzle_z"][sel], 1)[0]) * 3600.0
    v = AZ.vdot(run, WIN)
    d = dict(rop=rop, p_face=float(np.mean(th["patch_P_robin"][sel])),
             vdot_cm3=v * 1e6, tse_rock=4.0 * float(np.mean(th["patch_P_robin"][sel])) /
             (v * 1e9) if v else float("nan"),
             ledger=float(np.max(np.abs(th["ledger_err"][sel]))),
             t_end=float(th["time"][-1]))
    fz = run.lz - th["foot_z"]
    i = int(np.searchsorted(fz, z_m))
    d["t_match"] = float(th["time"][min(i, len(fz) - 1)])
    d["reached"] = bool(fz[-1] >= z_m - 1e-9)
    dep = run.depth(d["t_match"])
    zs = np.arange(0.01, z_m + 1e-9, 0.01)
    d["prof_z"] = zs
    d["prof"] = np.array([2.0 * AZ.az_rwall(run, dep, z) * 1e3 for z in zs])
    d["dia_min"] = float(d["prof"].min())
    d["dia_mean"] = float(d["prof"].mean())
    return d


def case_row(name, out, z_m):
    run = AZ.mrun(out, name)
    if run is None:
        return None
    d = base_row(run, z_m)
    d.update(run=run, case=name, out=out)
    d["v"] = M.vmetrics(run, (max(0.0, d["t_match"] - 100.0), d["t_match"]))
    # min Ø BELOW the 10 cm shaft depth: the flat-start transient (s = 50 mm > 45
    # everywhere at t = 0) widens the top of the hole and is not wall erosion.
    zs, pr = d["prof_z"], d["prof"]
    deep = zs >= Z_SHAFT
    d["has_shaft"] = bool(deep.any())
    d["shaft_min"] = float(pr[deep].min()) if deep.any() else float("nan")
    d["shaft_mean"] = float(pr[deep].mean()) if deep.any() else float("nan")
    d["shaft_max"] = float(pr[deep].max()) if deep.any() else float("nan")
    # Three gain measures (PREDICTIONS §5, amended). The control's values are
    # +32.4 / +5.0 / +5.0 -- the first reproduces ACTIVE_STEP's "~33 mm/side",
    # which is how the definition was checked.
    d["gain_mouth"] = 0.5 * (float(pr[0]) - SKIRT_MM)          # z = 10 mm
    i_pass = int(np.argmin(np.abs(zs - (z_m - 0.055))))        # fully past the band
    d["gain_steady"] = 0.5 * (float(pr[i_pass]) - float(pr[-1]))
    d["gain_shaft"] = 0.5 * (d["shaft_mean"] - SKIRT_MM)
    d["gain"] = d["gain_mouth"]                                # P2 is scored on the packet's
    d["dirname"] = name
    return d


def wall_state(d):
    """Wall temperature by depth and wall heat, at end of run, from the
    prescribed field + the plotfile surface temperature."""
    run = d["run"]
    case = d["dirname"].rsplit("_", 1)[0]
    st = surface_state(os.path.join(d["out"], d["dirname"]))
    if st is None:
        return None
    t, Ts, zf = st
    th = run.th
    j = int(np.argmin(np.abs(th["time"] - t)))
    zn = th["nozzle_z"][j]
    s = zn - zf
    ok = np.isfinite(s) & np.isfinite(Ts)
    band = ok & (s > 0.0) & (s <= R.S_SWITCH)
    if not band.any():
        return dict(t=t, n=0)
    r, dz = run.r, run.dz
    h = R.h_py(case, r[band], s[band])
    Tg = R.T_py(r[band], s[band])
    q = h * (Tg - Ts[band])                       # W/m2, per prescribed field
    area = dz * dz                                # projected top-face area
    P = 4.0 * float(np.sum(q) * area)             # quarter -> full hole
    depth = run.lz - zf[band]
    return dict(t=t, n=int(band.sum()), T=Ts[band], depth=depth, P=P,
                Tmin=float(Ts[band].min()), Tmax=float(Ts[band].max()),
                Tmean=float(Ts[band].mean()),
                near_fire=float(T_FIRE - Ts[band].max()),
                prof=[(lo, float(Ts[band][(depth >= lo) & (depth < lo + 0.05)].mean())
                       if ((depth >= lo) & (depth < lo + 0.05)).any() else float("nan"))
                      for lo in np.arange(0.0, 0.40, 0.05)])


def main():
    argparse.ArgumentParser().parse_args()
    names = [n for n in ORDER if os.path.exists(os.path.join(OUT, n + ".done"))]
    if not names:
        raise SystemExit("no completed runs in " + OUT)
    runs2 = [AZ.mrun(OUT, n) for n in names if n.endswith("_2mm")]
    ctl = AZ.mrun(CTRL_DIR, "LI0_2mm")
    z_m = AZ.common_feet_depth(runs2 + [ctl])
    rows = [r for r in (case_row(n, OUT, z_m) for n in names if n.endswith("_2mm")) if r]
    cr = case_row("LI0_2mm", CTRL_DIR, z_m)
    cr["case"] = "control"

    L = ["# D2o-0 analysis — prescribed floor/wall heating", "",
         f"Binary `{R.J.binary_sha()[:8]}…`; predictions hashed before launch (RESULTS §0). "
         f"Ø compared at **matched feet depth {z_m*1e3:.0f} mm**, never matched time. "
         "P1/P2 are scored **below 10 cm depth**, where no rock saw the flat-start transient "
         "(PREDICTIONS §4 Q2).", "",
         "## 1. THE RESULT: the prescribed split jams the burner", "",
         "| case | h_wall | skirt | ROP 150–250 s [m/h] | vs control | nozzle descent in 450 s | "
         "foot_stall_time at t_end [s] | patch_P_robin [W] |", "|---|---|---|---|---|---|---|---|"]
    for d in [cr] + rows:
        th = d["run"].th
        nzp = th["nozzle_z"][th["nozzle_z"] > 0]
        desc = (float(nzp.max()) - float(nzp[-1])) * 1e3   # row 0 is pre-init (0)
        nm = d["case"]
        hw = "—" if nm == "control" else f"{R.H_WALL_SENS if 'sens' in nm else R.H_WALL:g}"
        sk = "—" if nm == "control" else ("zero" if nm.startswith("Wz") else "grazing")
        rel = "—" if nm == "control" else f"{d['rop'] / cr['rop']:.2f}×"
        L.append(f"| {nm} | {hw} | {sk} | **{d['rop']:.3f}** | {rel} | "
                 f"{desc:.1f} mm | {th['foot_stall_time'][-1]:.0f} | {d['p_face']:.0f} |")
    L += ["", "**The burner descends ~27 mm and then stops.** `patch_min_standoff` grows from "
          "60 to 181 mm over the run: the centre pit keeps deepening (it is floor-zone and stays "
          "impingement-heated) while the feet are held up. The hole never reaches 10 cm depth, so "
          "**P1 and P2 are unscoreable as defined** — there is no shaft.", "",
          "### Why: a locking ratchet between the 0.9 pad quantile and the s = 45 mm switch", "",
          "The feet rest on the **0.9 nearest-rank quantile** of the annulus face heights, so ~10 % "
          "of annulus rock sits **above** the pad plane, i.e. at s < the 50 mm stand-off. **The "
          "floor/wall switch at s = 45 mm cuts straight through that population** — the switch is "
          "only 5 mm below the plane the feet themselves define.", "",
          "| run | annulus cols in the WALL zone (s ≤ 45 mm) | their surface T | fate |",
          "|---|---|---|---|",
          "| control | 5 of 162 | **680–814 K** | at the 821 K threshold — keep firing, tool descends |",
          "| W_2mm | 15 of 162 | **444–482 K** | 340 K below firing — **frozen permanently, tool jams** |",
          "",
          "`foot_carry_cols` = 18 in W_2mm against 20–26 in the control, and those carriers cannot "
          "spall at h = 40. **Any annulus column that drifts below the switch freezes and holds the "
          "burner up for the rest of the run.** This is a property of where the transition was "
          "placed, not of the grazing coefficient: all three variants (h = 40 grazing, h = 40 with "
          "a zero skirt, h = 20) stall, and their hole shapes are identical to 0.1 mm.", "",
         "## 2. Shape — such as it is", "",
         "| case | h_wall | skirt | shaft Ø min/mean/max (z>10cm) | **mouth gain** | "
         "**steady band gain** | shaft gain | mouth Ø | ROP [m/h] | V̇ [cm³/s] |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for d in [cr] + rows:
        nm = d["case"]
        hw = "—" if nm == "control" else f"{R.H_WALL_SENS if 'sens' in nm else R.H_WALL:g}"
        sk = "—" if nm == "control" else ("zero" if nm.startswith("Wz") else "grazing")
        sh = ("— (no shaft: hole < 10 cm)" if not d["has_shaft"] else
              f"{d['shaft_min']:.1f} / **{d['shaft_mean']:.1f}** / {d['shaft_max']:.1f}")
        L.append(f"| {nm} | {hw} | {sk} | {sh} | "
                 f"**{d['gain_mouth']:+.1f}** | **{d['gain_steady']:+.1f}** | "
                 f"{'—' if not d['has_shaft'] else f'{d[chr(103)+chr(97)+chr(105)+chr(110)+chr(95)+chr(115)+chr(104)+chr(97)+chr(102)+chr(116)]:+.1f}'} | {d['prof'][0]:.1f} | "
                 f"{d['rop']:.3f} | {d['vdot_cm3']:.2f} |")
    L += ["", f"P1 target shaft Ø **{MEIER_SHAFT[0]:.0f}–{MEIER_SHAFT[1]:.0f} mm**; "
          f"P2 target gain **{GAIN_BAND[0]}–{GAIN_BAND[1]} mm/side**; skirt = {SKIRT_MM:.0f} mm.", ""]
    L += ["", "Control values, which is how the definitions were checked: mouth gain "
          f"**{cr['gain_mouth']:+.1f}**, steady band gain **{cr['gain_steady']:+.1f}** mm/side. "
          "**These are not the control's own values** (+32.4 / +5.0, which reproduce ACTIVE_STEP's "
          "~33): the comparison depth collapsed to 48 mm because the test runs stalled, so every "
          "case is read inside its own start-up funnel. That is why the shape table decides "
          "nothing here.", ""]
    for d in rows:
        if not d["has_shaft"]:
            L.append(f"- **{d['case']}**: **P1 and P2 UNSCOREABLE** — the hole never reached "
                     f"10 cm depth (feet at {d['prof_z'][-1]*1e3:.0f} mm at end of run). "
                     "Reporting them from the start-up funnel would be meaningless.")
            continue
        p1_whole = MEIER_SHAFT[0] <= d["shaft_min"] and d["shaft_max"] <= MEIER_SHAFT[1]
        p1_mean = MEIER_SHAFT[0] <= d["shaft_mean"] <= MEIER_SHAFT[1]
        p2 = GAIN_BAND[0] <= d["gain_mouth"] <= GAIN_BAND[1]
        p2s = GAIN_BAND[0] <= d["gain_steady"] <= GAIN_BAND[1]
        L.append(f"- **{d['case']}**: P1 whole-shaft {'PASS' if p1_whole else 'FAIL'} "
                 f"({d['shaft_min']:.1f}-{d['shaft_max']:.1f} mm), P1 shaft-mean "
                 f"{'PASS' if p1_mean else 'FAIL'} ({d['shaft_mean']:.1f}); "
                 f"P2 mouth {'PASS' if p2 else 'FAIL'} ({d['gain_mouth']:+.1f}), "
                 f"P2 steady {'PASS' if p2s else 'FAIL'} ({d['gain_steady']:+.1f} mm/side).")

    L += ["", "## 3. Ø(z) at matched geometry [mm], 20 mm intervals", "",
          "| case | " + " | ".join(f"{z*1e3:.0f}" for z in rows[0]["prof_z"][1::2]) + " |",
          "|---|" + "---|" * len(rows[0]["prof_z"][1::2])]
    for d in [cr] + rows:
        L.append(f"| {d['case']} | " + " | ".join(f"{v:.0f}" for v in d["prof"][1::2]) + " |")

    L += ["", "## 4. v(r) — absolute, 2 mm annuli (D2m estimators, never normalised)", "",
          "| case | v_c [m/h] | " + " | ".join(f"{r} mm" for r in (25, 35, 45, 50, 55, 60)) + " |",
          "|---|---|---|---|---|---|---|---|"]
    for d in [cr] + rows:
        p = d["v"]["prof"]
        vals = [p[int(np.argmin(np.abs(M.CEN * 1e3 - r)))] for r in (25, 35, 45, 50, 55, 60)]
        L.append(f"| {d['case']} | {d['v']['v_c']:.3f} | "
                 + " | ".join(f"{v:.3f}" for v in vals) + " |")
    L += ["", "**Q6 asks whether v(r) becomes flatter (cylinder) or steeper (narrower cone).** "
          "Ratio v(45 mm)/v_c:", "",
          "| case | v(45)/v_c |", "|---|---|"]
    for d in [cr] + rows:
        p = d["v"]["prof"]
        v45 = p[int(np.argmin(np.abs(M.CEN * 1e3 - 45)))]
        L.append(f"| {d['case']} | {v45/d['v']['v_c']:.3f} |")

    L += ["", "## 5. Wall temperature and wall heat (P3, P4, P5)", "",
          "Surface temperature of the band columns (0 < s ≤ 45 mm) at end of run, from the "
          "plotfile. **Wall heat is computed from the prescribed field and the wall area** — with "
          "`jet_closure = none` there is no march and therefore no in-code wall-power "
          "diagnostic.", "",
          "| case | t [s] | band cols | T_wall min/mean/max [K] | 821 − T_max | **wall heat [W]** | "
          "P3 | P4 | P5 |", "|---|---|---|---|---|---|---|---|---|"]
    ws = {}
    for d in rows + [cr]:
        w = wall_state(d)
        ws[d["case"]] = w
        if w is None or w.get("n", 0) == 0:
            L.append(f"| {d['case']} | — | 0 | — | — | — | — | — | — |")
            continue
        ctlrow = d["case"] == "control"
        p3 = T_WALL_BAND[0] <= w["Tmean"] <= T_WALL_BAND[1]
        p4 = P_WALL_BAND[0] <= w["P"] <= P_WALL_BAND[1]
        p5 = w["near_fire"] <= P5_NEAR
        # The control is heated by the IMPINGING field, not this prescribed one,
        # so its wall heat under this field is meaningless and is not shown.
        pw = "— (impinging field)" if ctlrow else f"**{w['P']:.0f}**"
        L.append(f"| {d['case']} | {w['t']:.0f} | {w['n']} | "
                 f"{w['Tmin']:.0f} / **{w['Tmean']:.0f}** / {w['Tmax']:.0f} | "
                 f"{w['near_fire']:+.0f} K | {pw} | "
                 f"{'PASS' if p3 else 'FAIL'} | {'—' if ctlrow else ('PASS' if p4 else 'FAIL')} | "
                 f"{'**FIRES**' if p5 else 'no'} |")
    L += ["", "Wall temperature by depth [K], 50 mm bins (end of run):", "",
          "| case | " + " | ".join(f"{z*1e3:.0f}–{z*1e3+50:.0f}" for z in np.arange(0, 0.40, 0.05))
          + " |", "|---|" + "---|" * 8]
    for d in rows + [cr]:
        w = ws.get(d["case"])
        if not w or w.get("n", 0) == 0:
            continue
        L.append(f"| {d['case']} | " + " | ".join(
            "—" if not np.isfinite(v) else f"{v:.0f}" for _, v in w["prof"]) + " |")

    L += ["", "## 6. Mesh — 2 mm vs 1 mm by 100 s block (a single mesh is not quotable)", ""]
    d1 = case_row("W_1mm", OUT, z_m) if os.path.exists(os.path.join(OUT, "W_1mm.done")) else None
    d2 = next((d for d in rows if d["case"] == "W_2mm"), None)
    if d1 and d2:
        L += ["| block | 2 mm ROP | 1 mm ROP | gap |", "|---|---|---|---|"]
        for b in BLOCKS:
            g = []
            for d in (d2, d1):
                th = d["run"].th
                s = (th["time"] >= b[0]) & (th["time"] <= b[1])
                g.append(-float(np.polyfit(th["time"][s], th["nozzle_z"][s], 1)[0]) * 3600.0
                         if s.sum() > 2 else float("nan"))
            L.append(f"| {b[0]:.0f}–{b[1]:.0f} s | {g[0]:.3f} | {g[1]:.3f} | "
                     f"**{(g[1]/g[0]-1)*100:+.1f} %** |")
        L += ["", f"1 mm shaft mean Ø **{d1['shaft_mean']:.1f} mm** vs 2 mm "
              f"**{d2['shaft_mean']:.1f} mm**; steady band gain {d1['gain_steady']:+.1f} vs "
              f"{d2['gain_steady']:+.1f} mm/side."]
    else:
        L += ["_1 mm leg not complete_"]

    L += ["", "## 7. Ledger and absorbed power", "",
          "| case | max abs ledger_err | patch_P_robin [W] | J/mm³ to rock (floor 1.147) |",
          "|---|---|---|---|"]
    for d in [cr] + rows + ([d1] if d1 else []):
        L.append(f"| {d['case']} | {d['ledger']:.1e} | {d['p_face']:.0f} | "
                 f"{d['tse_rock']:.2f} |")

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
    ax[0].axvspan(*MEIER_SHAFT, color="g", alpha=0.15, label="Meier shaft 85–95 mm")
    ax[0].axvline(SKIRT_MM, color="k", ls="--", lw=1, label="skirt Ø 80 mm")
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
    ax[1].axvspan(*T_WALL_BAND, color="g", alpha=0.15, label="P3 600–800 K")
    ax[1].axvline(T_FIRE, color="r", ls="--", lw=1.2, label="spall threshold 821 K")
    ax[1].axvline(570, color="C1", ls=":", lw=1.2, label="alteration ≈570 K")
    ax[1].invert_yaxis()
    ax[1].set(xlabel="wall surface temperature [K]", ylabel="depth [mm]",
              title="Wall temperature by depth, end of run")
    ax[1].grid(alpha=0.3)
    ax[1].legend(fontsize=8)
    f.suptitle("D2o-0: prescribed floor/wall split — does a cylinder appear?", y=0.99)
    f.tight_layout()
    f.savefig(os.path.join(HERE, "d2o0_profiles.png"), dpi=110)
    plt.close(f)

    f2, a2 = plt.subplots(figsize=(7.2, 4.8))
    a2.axvspan(42.5, 46.5, color="g", alpha=0.15, label="Meier implied floor half-width")
    for d in allr:
        a2.plot(M.CEN * 1e3, d["v"]["prof"], "-", lw=1.8, label=d["case"])
    a2.axvline(R.R_SKIRT * 1e3, color="k", ls="--", lw=1, label="skirt r = 40 mm")
    a2.set(xlabel="radius [mm]", ylabel="recession rate [m/h]", xlim=(0, 70),
           title="v(r), absolute — cone or cylinder?")
    a2.grid(alpha=0.3)
    a2.legend(fontsize=8)
    f2.tight_layout()
    f2.savefig(os.path.join(HERE, "d2o0_vprofile.png"), dpi=110)
    plt.close(f2)
    print("wrote d2o0_profiles.png, d2o0_vprofile.png")


if __name__ == "__main__":
    main()
