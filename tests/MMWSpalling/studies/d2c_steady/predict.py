#!/usr/bin/env python3
"""D2c deliverable 8: hand predictions for the D2d campaign, written to
PREDICTIONS.md and hashed (sha256 + mtime into RESULTS.md) before any D2d run.

  predict.py            write PREDICTIONS.md + RESULTS.md §0 (refuses to overwrite)
  predict.py --print    print only

Model: handmodel.py (steady = d2b_feet_rop/feetmodel with the D2c closures;
transient = burner on feet with free-surface dilution and body clearance).
Scored configuration (ACTIVE_STEP): J-M (Martin h_ref at T_nozzle), 1900 K,
exhaust recirculation, entrained mass, momentum D_e (nozzle.py), core 8, far
law power n = 1, free surface aspect 1 / exponent 1, clearance on.
"""
import argparse
import datetime
import hashlib
import io
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import handmodel as hm  # noqa: E402
import nozzle as nz  # noqa: E402

fm, wj = hm.fm, hm.wj
PRED = os.path.join(HERE, "PREDICTIONS.md")
RESULTS = os.path.join(HERE, "RESULTS.md")
D2B_JM = os.path.join(os.path.dirname(HERE), "d2b_feet_rop", "output", "B_JM_A2")
T_END = 3000.0
HEIGHTS = (0.30, 0.40)
BOTTOM = 0.040
STEADY_STOP = 50.0
BAND = (1.04, 1.92)
MOUTH_TARGET = {25: 92.0, 50: 88.0}          # meier_fig8_8_hole_profile.csv


def case(T=1900.0, far="power", n=1.0, core=8.0, momentum=True, fs=True, aspect=1.0, exponent=1.0,
         clear=True):
    return dict(T=T, cfg=hm.cfg_for(T, far, n, core, momentum),
                fs=hm.FS(fs, aspect, exponent), clear=clear)


CASES = {
    "scored (n = 1)": case(),
    "n = 0.5": case(n=0.5),
    "clamp": case(far="clamp"),
    "T_nozzle 1600 K": case(T=1600.0),
    "T_nozzle 1750 K": case(T=1750.0),
    "core 5": case(core=5.0),
    "Ricou (core 3.125)": case(core=3.125),
    "D2b closures + n = 1": case(momentum=False, core=5.0, fs=False),
    "aspect 0.5": case(aspect=0.5),
    "aspect 2": case(aspect=2.0),
    "exponent 0.8": case(exponent=0.8),
}
D2B_CHECK = dict(T=1900.0, cfg=fm.Cfg(), fs=hm.FS(False), clear=False)


def t_reach(res, depth):
    i = np.nonzero(res["z_c"] >= depth)[0]
    return float(res["t"][i[0]]) if i.size else None


def at(res, key, t):
    return float(np.interp(t, res["t"], res[key]))


def fit(res, t0, t1, key="d"):
    sel = (res["t"] >= t0) & (res["t"] <= t1)
    return float(np.polyfit(res["t"][sel], res[key][sel], 1)[0])


def hole(c, res, w0, w1):
    """Profile at w1 (rerun to that time: z(r) is the transient's end state),
    whole-volume rate, burner ROP and mechanical share over [w0, w1]."""
    r2 = hm.transient(c["T"], c["cfg"], c["fs"], t_end=w1, clear=c["clear"])
    drilled = float(r2["d"][-1]) - fm.STANDOFF
    pm = hm.profile_metrics(r2, drilled)
    dV = at(res, "V", w1) - at(res, "V", w0)
    return dict(win=(w0, w1), drilled=drilled, pm=pm,
                d25=2.0 * float(hm.wall_profile(r2, [0.025])[0]),
                vol=hm.volume_rate(res, w0, w1), rop=fit(res, w0, w1),
                mech=(at(res, "V_mech", w1) - at(res, "V_mech", w0)) / max(1e-30, dV))


def evaluate(c, t_end=T_END):
    st = hm.steady(c["T"], c["cfg"])
    res = hm.transient(c["T"], c["cfg"], c["fs"], t_end=t_end, clear=c["clear"])
    t0, rate = hm.steady_window(res)
    out = dict(st=st, res=res, t0=t0, rate=rate, h_ref=res["h_ref"],
               De_ref=c["cfg"].De_ref, D_e_steady=st["D_e"] if st else None)
    out["t_bottom"] = {H: t_reach(res, H - BOTTOM) for H in HEIGHTS}
    if t0 is not None:
        t_run = t0 + hm.STEADY_LEN + STEADY_STOP
        noz = at(res, "d", t_run) - fm.STANDOFF
        sc = at(res, "s_c", t_run)
        need = noz + sc + BOTTOM
        out.update(t_run=t_run, noz_depth=noz, sc_end=sc, H_need=need,
                   H_rule=0.30 if need <= 0.30 else 0.40, fits=need <= 0.40,
                   steady_hole=hole(c, res, t0, t0 + hm.STEADY_LEN))
    else:
        out.update(t_run=None, noz_depth=None, sc_end=None, H_need=None, H_rule=0.40, fits=False,
                   steady_hole=None)
    # what a 0.40 m run measures: 150 s to the bottom stop
    tb = out["t_bottom"][0.40]
    out["bottom_hole"] = hole(c, res, 150.0, tb) if tb and tb > 200.0 else None
    out["rate_150_bottom"] = out["bottom_hole"]["rop"] if out["bottom_hole"] else None
    return out


def d2b_crosscheck():
    """Hand transient with the D2b closures vs the D2b B_JM_A2 simulation."""
    res = hm.transient(D2B_CHECK["T"], D2B_CHECK["cfg"], D2B_CHECK["fs"], t_end=342.0,
                       clear=False)
    th = None
    p = os.path.join(D2B_JM, "thermo.dat")
    if os.path.exists(p):
        names = open(p).readline().split()
        a = np.loadtxt(p, skiprows=1)
        th = {n: a[:, i] for i, n in enumerate(names)}
    rows = []
    hand_rop = fit(res, 171.0, 341.0) * 3600.0
    hand_zc = at(res, "z_c", 341.0)
    hand_trec = float(np.mean(res["T_rec"][(res["t"] >= 171.0) & (res["t"] <= 341.0)]))
    hand_sc = at(res, "s_c", 341.0)
    if th is not None:
        sel = (th["time"] >= 171.0) & (th["time"] <= 341.6)
        top = 0.30
        sim_rop = -np.polyfit(th["time"][sel], th["nozzle_z"][sel], 1)[0] * 3600.0
        j = np.nonzero(th["time"] <= 341.0)[0][-1]
        sim_zc = top - (th["nozzle_z"][j] - th["jet_s_c"][j])
        sim_trec = float(np.mean(th["jet_T_rec"][sel]))
        sim_sc = float(th["jet_s_c"][j])
        rows = [("burner ROP 171–341 s [m/h]", hand_rop, sim_rop),
                ("centre depth at 341 s [mm]", hand_zc * 1e3, sim_zc * 1e3),
                ("s_c at 341 s [mm]", hand_sc * 1e3, sim_sc * 1e3),
                ("mean T_rec 171–341 s [K]", hand_trec, sim_trec)]
    return rows


def fmt(v, f="{:.2f}", none="—"):
    return none if v is None else f.format(v)


def build():
    lines = []
    w = lines.append
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    w("# D2c hand predictions for the D2d Meier campaign")
    w("")
    w(f"Written {now} by `predict.py`, **before any D2d run**. The sha256 and mtime of this file are")
    w("recorded in `RESULTS.md` §0. Nothing here may be edited after a D2d run starts.")
    w("")
    w("## Model")
    w("")
    w("- **Steady state:** `d2b_feet_rop/feetmodel.py` with the D2c closures (`Cfg`: far law,")
    w("  momentum D_e = jet_De_ref·√(T_mix/T_ent), core length). The D2b geometry is unchanged:")
    w("  the ring at r_q = 38.97 mm at s = 50 mm; s linear from s_c (r ≤ 18.75 mm) to 50 mm;")
    w("  face at T_fire = 821 K; T_rec = T_gas(r_q) as a fixed point.")
    w("  Deep in the hole no free-surface column marches, because the nozzle plane is below")
    w("  the original surface, so dilution does not enter the steady state.")
    w("- **Transient** (`handmodel.transient`, dt 0.25 s):")
    w("  - the ring (feet) and the centre advance at their pinned closed-form rates;")
    w("  - outer columns on the 2 mm mesh grid march while s > 0, with ALAMO's per-bin")
    w("    free-surface dilution;")
    w("  - T_rec = the mouth inlet T (the exhaust T with no free-surface column);")
    w("  - the annulus band r_q–40 mm is cut to the feet depth when proud (clearance).")
    w("  It gives the time to a steady window (the score's rule applied to the smooth")
    w("  series: 150 s with |ds_c/dt| < 0.02 mm/s and 50 s sub-fits within ±10 %), the wall")
    w("  profile, and the whole-excavation and mechanical volume. It covers the full jet;")
    w("  shares and rates per quarter follow by /4.")
    w("- **Anchors per T_nozzle:** Martin h_ref = `walljet.h_anchor('M', T)`; jet_De_ref =")
    w("  `nozzle.de_ref(T)` (γ 1.3, J including the pressure thrust).")
    w("")
    # cross-check
    xc = d2b_crosscheck()
    w("## Model check against D2b (not a prediction)")
    w("")
    w("The transient model with the D2b closures (clamp, nozzle D, core 5, no free surface, no")
    w("clearance), compared with the D2b scored run B_JM_A2:")
    w("")
    w("| quantity | hand | D2b sim | hand / sim − 1 |")
    w("|---|---|---|---|")
    errs = []
    for n, hv, sv in xc:
        e = hv / sv - 1.0
        errs.append(abs(e))
        w(f"| {n} | {hv:.2f} | {sv:.2f} | {e * 100:+.0f} % |")
    rop_err = abs(xc[0][1] / xc[0][2] - 1.0) if xc else 0.25
    tol = max(0.10, 2.0 * rop_err)
    w("")
    w(f"The hand model runs fast and hot against the simulation. The refutation bands below")
    w(f"are **±{tol * 100:.0f} % on rates**, i.e. twice the hand model's D2b ROP error of")
    w(f"{rop_err * 100:.0f} %, rounded; and **±20 mm on Ø and stand-off** (5 mesh cells in Ø).")
    w("")
    ev = {n: evaluate(c) for n, c in CASES.items()}
    sc = ev["scored (n = 1)"]

    def steady_row(n):
        e = ev[n]
        s = e["st"]
        if s is None:
            return f"| {n} | none | — | — | — | — | — |"
        deep = s["s_c"] >= fm.S_C_MAX
        sc_txt = "none (deep-pit limit)" if deep else f"{s['s_c'] * 1e3:.0f} mm"
        return (f"| {n} | {s['v']:.2f} | {sc_txt} | "
                f"{s['T_rec']:.0f} K | {s['T_stag']:.0f} K | {s['phi']:.3f} | {s['D_e'] * 1e3:.2f} mm | "
                f"{s['P_face'] / 4e3:.2f} kW |")

    w("## P1. Centre equilibrium, time to a steady window, domain height")
    w("")
    w("Steady state (full jet; face power per quarter):")
    w("")
    w("| case | ring ROP [m/h] | centre s_c | T_rec | T_stag | φ | D_e | P_face / 4 |")
    w("|---|---|---|---|---|---|---|---|")
    for n in ("scored (n = 1)", "n = 0.5", "clamp"):
        w(steady_row(n))
    w("")
    w("Transient:")
    w("")
    w("| case | steady window from | ROP in window | run end (window + 50 s) | nozzle-plane depth then | s_c then | domain needed | centre 40 mm above the bottom of 0.30 / 0.40 m at |")
    w("|---|---|---|---|---|---|---|---|")
    for n in ("scored (n = 1)", "n = 0.5", "clamp"):
        e = ev[n]
        tb = e["t_bottom"]
        w(f"| {n} | {fmt(e['t0'], '{:.0f} s', 'never (to 3000 s)')} | {fmt(e['rate'] and e['rate'] * 3600, '{:.2f} m/h')} | "
          f"{fmt(e['t_run'], '{:.0f} s')} | {fmt(e['noz_depth'] and e['noz_depth'] * 1e3, '{:.0f} mm')} | "
          f"{fmt(e['sc_end'] and e['sc_end'] * 1e3, '{:.0f} mm')} | {fmt(e['H_need'], '{:.2f} m')} | "
          f"{fmt(tb[0.30], '{:.0f} s')} / {fmt(tb[0.40], '{:.0f} s')} |")
    w("")
    w("**Domain-height rule** (ACTIVE_STEP deliverable 8): use 0.30 m if nozzle-plane depth")
    w("at the end of the steady window + s_c + 40 mm ≤ 0.30 m, else 0.40 m. The run end")
    w("(window + 50 s steady-stop) is used as the end.")
    w("")
    w(f"- **Scored run: the rule gives {sc['H_rule']:.2f} m, but neither height fits.** The")
    w(f"  hand steady window opens at {fmt(sc['t0'], '{:.0f}')} s. By the end of the run the")
    w(f"  nozzle plane is {fmt(sc['noz_depth'] and sc['noz_depth'] * 1e3, '{:.0f}')} mm deep, so the domain")
    w(f"  must be ≥ {fmt(sc['H_need'], '{:.2f}')} m. In a 0.40 m domain the centre reaches the bottom")
    w(f"  watchdog at {fmt(sc['t_bottom'][0.40], '{:.0f}')} s, before any window.")
    w("  **The hand model therefore predicts that D2d run 1 at 0.40 m stops at the bottom")
    w("  with no steady window.** That decision (a taller domain, or scoring without a window)")
    w("  is D2d's; this packet does not trim or change anything.")
    w("- **What delays the window:** the s_c drift criterion. The ROP sub-fits are within")
    w(f"  ±10 % from early on, while |ds_c/dt| < 0.02 mm/s is first met at {fmt(sc['t0'], '{:.0f}')} s.")
    e5 = ev["n = 0.5"]
    w(f"- **n = 0.5:** s_c ≈ {e5['st']['s_c'] * 1e3:.0f} mm; the domain would have to be")
    w(f"  ≥ {fmt(e5['H_need'], '{:.2f}')} m. **n = 0.5 does not fit** either height, so D2d reports")
    w("  n = 0.5 from the hand model only (run 3 is not run).")
    ec = ev["clamp"]
    w(f"- **Clamp:** no centre equilibrium (s_c → the 2 m bisection bound). The")
    w(f"  deep-pit limit is {ec['st']['v']:.2f} m/h, so run 2 ends at the bottom")
    w(f"  ({fmt(ec['t_bottom'][0.40], '{:.0f}')} s at 0.40 m).")
    w("")
    w("## P2. Ring ROP and its sensitivity to n")
    w("")
    w("| case | steady ring ROP | transient ROP, 150 s to the 0.40 m bottom | ROP in the hand steady window |")
    w("|---|---|---|---|")
    for n in ("scored (n = 1)", "n = 0.5", "clamp"):
        e = ev[n]
        w(f"| {n} | {e['st']['v']:.3f} m/h | {fmt(e['rate_150_bottom'] and e['rate_150_bottom'] * 3600, '{:.2f} m/h')} | "
          f"{fmt(e['rate'] and e['rate'] * 3600, '{:.2f} m/h')} |")
    v1, v05, vc = (ev[k]["st"]["v"] for k in ("scored (n = 1)", "n = 0.5", "clamp"))
    w("")
    w(f"Computed sensitivity (steady):")
    w(f"- n = 1 → 0.5 changes the ring ROP by {100 * (v05 / v1 - 1):+.1f} %;")
    w(f"- n = 1 → clamp (n = 0) by {100 * (vc / v1 - 1):+.1f} %.")
    w("")
    w("The ring hardly feels n. The far law moves the centre (s_c), and the ring feels it only")
    w("through T_rec. **The scored steady ring ROP, "
      f"{v1:.2f} m/h, is {'above' if v1 > BAND[1] else 'inside'} Meier's band [1.04, 1.92].**")
    w("")
    w("## P3. T_nozzle sweep (scored closures)")
    w("")
    w("| T_nozzle | Martin h_ref | jet_De_ref | steady D_e (at T_rec) | steady ROP | transient ROP 150 s → 0.40 m bottom | steady T_rec |")
    w("|---|---|---|---|---|---|---|")
    for T, n in ((1600.0, "T_nozzle 1600 K"), (1750.0, "T_nozzle 1750 K"), (1900.0, "scored (n = 1)")):
        e = ev[n]
        w(f"| {T:.0f} K | {e['h_ref']:.0f} | {nz.de_ref(T)!r} m | {e['st']['D_e'] * 1e3:.2f} mm | "
          f"{e['st']['v']:.2f} m/h | {fmt(e['rate_150_bottom'] and e['rate_150_bottom'] * 3600, '{:.2f} m/h')} | "
          f"{e['st']['T_rec']:.0f} K |")
    Ts = np.array([1600.0, 1750.0, 1900.0])
    vs = np.array([ev["T_nozzle 1600 K"]["st"]["v"], ev["T_nozzle 1750 K"]["st"]["v"], v1])
    lo = float(np.interp(1.3, vs, Ts)) if vs.min() <= 1.3 <= vs.max() else None
    hi = float(np.interp(1.6, vs, Ts)) if vs.min() <= 1.6 <= vs.max() else None
    w("")
    w(f"Steady ROP rises {100 * (vs[2] / vs[0] - 1) / 300:.3f} %/K (of the 1600 K value) between 1600 and 1900 K.")
    w(f"Linear interpolation gives 1.3 m/h at T_nozzle ≈ {fmt(lo, '{:.0f} K', 'below 1600 K')} and")
    w(f"1.6 m/h at ≈ {fmt(hi, '{:.0f} K', 'outside 1600–1900 K')}.")
    w("This is conditional on J-M, n = 1, free-surface dilution and core 8, and is never reused")
    w("as an input.")
    w("")
    w("## P4. Hole: mouth, depth-mean, minimum, volume (scored)")
    w("")

    def hole_row(label, e, hh):
        pm = hh["pm"]
        return (f"| {label} | {hh['win'][0]:.0f}–{hh['win'][1]:.0f} s | {hh['drilled'] * 1e3:.0f} mm | "
                f"{pm['d20'] * 1e3:.0f} | {hh['d25'] * 1e3:.0f} | {pm['d50'] * 1e3:.0f} | {pm['mean'] * 1e3:.0f} | "
                f"{pm['min'] * 1e3:.0f} | {hh['vol'] * 1e6:.2f} | {hh['rop'] * 3600:.2f} | {hh['mech'] * 100:.2f} % |")

    w("Two windows:")
    w("- **(a) 150 s to the 0.40 m bottom stop.** This is what D2d run 1 will measure, and the")
    w("  P6 bands use it.")
    w("- **(b) The hand steady window.** It needs a ≥ 0.65 m domain.")
    w("")
    w("The profile is taken at the window end, over the drilled depth (0 to the nozzle-plane")
    w("depth). Ø is in mm. The volume rate is the whole excavation of the full hole, in cm³/s")
    w("(score target 1.98–2.96).")
    w("")
    w("| case | window | drilled depth | Ø at 20 mm | Ø at 25 mm | Ø at 50 mm | depth-mean Ø | min Ø | whole volume rate | burner ROP [m/h] | mechanical share |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    bh = sc["bottom_hole"]
    w(hole_row("scored, (a) to the 0.40 m bottom", sc, bh))
    if sc["steady_hole"]:
        w(hole_row("scored, (b) hand steady window", sc, sc["steady_hole"]))
    w("")
    pm = bh["pm"]
    w("Reading of window (a) against the targets:")
    w(f"- **Mouth:** Fig. 8.8 has 92 mm at 25 mm and 88 mm at 50 mm (±10 mm). The hand gives")
    w(f"  {bh['d25'] * 1e3:.0f} and {pm['d50'] * 1e3:.0f} mm, so the **mouth criterion is predicted to fail**:")
    w("  free-surface dilution narrows the funnel but does not remove it. The D2b closures give")
    w("  the Ø in the P5 row.")
    w(f"- **Minimum Ø:** {pm['min'] * 1e3:.0f} mm. The hand floor is 2·r_q = 78 mm; the wall")
    w("  grid is 2 mm, so Ø moves in 4 mm steps.")
    w(f"- **Depth-mean Ø:** {pm['mean'] * 1e3:.0f} mm against the 85–93 mm target.")
    w(f"- **Volume rate:** {bh['vol'] * 1e6:.2f} cm³/s against 1.98–2.96.")
    w(f"- **Mechanical share:** {bh['mech'] * 100:.2f} % against the 5 % finding threshold.")
    w("  This is a lower bound: the hand model has one smooth clearance band (r_q–40 mm),")
    w("  while the mesh has discrete columns proud of a nearest-rank pad.")
    w("")
    w("## P5. Variants (same rows)")
    w("")
    w("Hole metrics are for window (a), 150 s to the 0.40 m bottom stop.")
    w("")
    w("| case | steady ROP | steady s_c | hand window from | domain needed | bottom stop (0.40 m) | ROP (a) | Ø 20 | Ø 25 | Ø 50 | depth-mean Ø | min Ø | volume rate (a) | mech share (a) |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for n in ("scored (n = 1)", "core 5", "Ricou (core 3.125)", "D2b closures + n = 1", "aspect 0.5",
              "aspect 2", "exponent 0.8"):
        e = ev[n]
        s = e["st"]
        hh = e["bottom_hole"]
        pm = hh["pm"]
        w(f"| {n} | {s['v']:.2f} | {s['s_c'] * 1e3:.0f} mm | {fmt(e['t0'], '{:.0f} s', 'never')} | "
          f"{fmt(e['H_need'], '{:.2f} m')} | {fmt(e['t_bottom'][0.40], '{:.0f} s')} | {hh['rop'] * 3600:.2f} | "
          f"{pm['d20'] * 1e3:.0f} | {hh['d25'] * 1e3:.0f} | {pm['d50'] * 1e3:.0f} | {pm['mean'] * 1e3:.0f} | "
          f"{pm['min'] * 1e3:.0f} | {hh['vol'] * 1e6:.2f} | {hh['mech'] * 100:.2f} % |")
    w("")
    w("## P6. What D2d would refute")
    w("")
    rt = tol * 100
    s1 = sc["st"]
    w(f"Bands: ±{rt:.0f} % on rates, ±20 mm on Ø and stand-off (see the model check).")
    w("")
    w(f"- **P1, scored:** refuted if run 1 finds a steady window before its centre reaches the")
    w("  0.40 m bottom, or if its s_c settles (drift < 0.02 mm/s for 150 s) outside")
    w(f"  {s1['s_c'] * 1e3 - 20:.0f}–{s1['s_c'] * 1e3 + 20:.0f} mm.")
    w("- **P1, clamp (run 2):** refuted if run 2 shows a centre equilibrium, i.e. |ds_c/dt|")
    w("  < 0.02 mm/s for 150 s before the bottom.")
    w("- **P1, n = 0.5:** not tested by D2d (does not fit); refutable only with a ≥ "
      f"{fmt(e5['H_need'], '{:.1f}')} m domain.")
    rb = sc["rate_150_bottom"] * 3600 if sc["rate_150_bottom"] else None
    w(f"- **P2:** refuted if run 1's burner ROP from 150 s to the bottom stop is outside")
    w(f"  {fmt(rb and rb * (1 - tol))}–{fmt(rb and rb * (1 + tol))} m/h.")
    rcl = ev["clamp"]["rate_150_bottom"] * 3600 if ev["clamp"]["rate_150_bottom"] else None
    w(f"- **P2, clamp:** refuted if run 2's burner ROP from 150 s to its bottom stop is outside")
    w(f"  {fmt(rcl and rcl * (1 - tol))}–{fmt(rcl and rcl * (1 + tol))} m/h (hand {fmt(rcl)} m/h).")
    for T, n in ((1600.0, "T_nozzle 1600 K"), (1750.0, "T_nozzle 1750 K")):
        e = ev[n]
        r = e["rate_150_bottom"] * 3600 if e["rate_150_bottom"] else None
        w(f"- **P3, {T:.0f} K:** refuted if the ROP from 150 s to the bottom stop is outside")
        w(f"  {fmt(r and r * (1 - tol))}–{fmt(r and r * (1 + tol))} m/h.")
    w(f"- **P3, ordering:** refuted if ROP is not increasing in T_nozzle across runs 4, 5 and 1.")
    bh = sc["bottom_hole"]
    w(f"- **P4 (window a):** refuted if run 1's Ø at 25 or 50 mm is within ±10 mm of Fig. 8.8, i.e.")
    w(f"  the mouth passes; the hand says {bh['d25'] * 1e3:.0f} / {bh['pm']['d50'] * 1e3:.0f} mm, ±20 mm.")
    w(f"  Also refuted if the depth-mean Ø is outside {bh['pm']['mean'] * 1e3 - 20:.0f}–{bh['pm']['mean'] * 1e3 + 20:.0f} mm,")
    w(f"  the mechanical share exceeds 5 % (hand {bh['mech'] * 100:.2f} %), or the whole volume rate")
    w(f"  is outside {bh['vol'] * 1e6 * (1 - tol):.2f}–{bh['vol'] * 1e6 * (1 + tol):.2f} cm³/s.")
    e6 = ev["D2b closures + n = 1"]
    w(f"- **P5, run 6 (D2b closures + n = 1):** refuted if its Ø at 25 mm is not wider than run 1's")
    w(f"  (hand: {e6['bottom_hole']['d25'] * 1e3:.0f} vs {bh['d25'] * 1e3:.0f} mm; the free-surface dilution is the funnel")
    w("  lever), or if its ROP from 150 s to the bottom differs from the hand value "
      f"{fmt(e6['rate_150_bottom'] and e6['rate_150_bottom'] * 3600)} m/h by more than ±{rt:.0f} %.")
    w("- **P5, core 5 / Ricou (optional runs):** refuted if their ROP ordering against run 1")
    w(f"  differs from the hand's (core 5 {fmt(ev['core 5']['rate_150_bottom'] and ev['core 5']['rate_150_bottom'] * 3600)}, "
      f"Ricou {fmt(ev['Ricou (core 3.125)']['rate_150_bottom'] and ev['Ricou (core 3.125)']['rate_150_bottom'] * 3600)}, "
      f"scored {fmt(rb)} m/h).")
    w("- **P5, aspect / exponent:** hand only (no D2d run). Reported for the size of the")
    w("  mouth lever.")
    w("")
    w("## dt for D2d (overshoot rule, 2 mm)")
    w("")
    T_max = float(np.max(sc["res"]["T_stag"]))
    w(f"The hand transient's maximum T_stag is {T_max:.0f} K (scored). run.py evaluates the")
    w("overshoot rule, ≤ 5 K per step, at q = h(0, s ≥ 50 mm)·(T_stag,max − 821) − loss.")
    w("")
    return "\n".join(lines) + "\n", ev


def results_header(sha, mtime):
    return f"""# D2c results: steady-state and no-funnel closures (pre-registration; no scored runs)

## 0. Pre-registration record

- `PREDICTIONS.md` sha256 `{sha}`
- `PREDICTIONS.md` mtime {mtime}
- Written by `predict.py` before any D2d run and before the D2c cost probes.
- Scored configuration: ACTIVE_STEP §Scored configuration, unchanged.

<!-- NOZZLE -->

<!-- DRYRUN -->

<!-- PROBES -->

<!-- NOTES -->
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--print", action="store_true")
    args = ap.parse_args()
    text, _ = build()
    if args.print:
        print(text)
        return
    if os.path.exists(PRED):
        raise SystemExit(f"{PRED} exists; refusing to overwrite the pre-registration")
    with open(PRED, "w") as f:
        f.write(text)
    sha = hashlib.sha256(open(PRED, "rb").read()).hexdigest()
    mtime = datetime.datetime.fromtimestamp(os.path.getmtime(PRED)).strftime("%Y-%m-%d %H:%M:%S")
    if os.path.exists(RESULTS):
        raise SystemExit(f"{RESULTS} exists; record the hash by hand: {sha} {mtime}")
    with open(RESULTS, "w") as f:
        f.write(results_header(sha, mtime))
    print(text)
    print(f"sha256 {sha}  mtime {mtime}")


if __name__ == "__main__":
    main()
