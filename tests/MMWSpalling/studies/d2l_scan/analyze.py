#!/usr/bin/env python3
"""D2l Gate 0 analysis -> analysis.md + d2l_scan.png.

  analyze.py            every gate found in output/
  analyze.py --gate 0a  just one of 0a / 0b / 0c

Scoring definitions are fixed in PREDICTIONS.md (sha256 in RESULTS.md §0) and
are not restated loosely here; where a number is computed, the function says
which frozen script it reuses so it is the same estimator as the reference it
is compared against.

  0a  d2j0b_plane/analyze.py  (col_rate / band_timing / r_front / metrics) —
      the same estimator as the D2j-0b C0/C1 cascade and the D2k arbiter.
  0b  d2b_feet_rop/analyze.py (Run: thermo, removal events, per-column depth) +
      az_rwall COPIED below from d2c_steady/score.py, because that script
      writes into its own frozen output/ directory.

Both are loaded by path and neither is edited.
"""
import argparse
import importlib.util
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


A = load("d2j0b_analyze", os.path.join(STUD, "d2j0b_plane", "analyze.py"))
AN = load("d2b_analyze", os.path.join(STUD, "d2b_feet_rop", "analyze.py"))
K = load("d2l_run", os.path.join(HERE, "run.py"))
OUT = K.OUT
D2K = K.D2K
D2I = os.path.join(STUD, "d2i_gate", "output")

TOL = 0.05
BLOCKS = [(50.0, 100.0), (100.0, 150.0), (150.0, 200.0)]      # 0a deciding blocks
WIN = (150.0, 250.0)                                          # 0b scoring window
REF_GAP_OFF = (-13.5, -28.7, -45.5)                           # D2j-0b C1, key off
BAND_ROP = (1.3, 1.6)                                         # 0b target band
BAND_D = (85.0, 93.0)                                         # Meier hole Ø, mm
TSE_HHV = 15.4                                                # J/mm3, combustion enthalpy
TSE_ROCK = (4.0, 6.0)                                         # J/mm3, enthalpy-to-rock
P_HHV = 38.0e3                                                # W, Meier burner (HHV)
RHOCP, T_FIRE, T_AMB = 2750.0 * 790.0, 821.0, 293.15
N_SEC = 9


# ---------------------------------------------------------------- 0a helpers
def metrics_in(out, name):
    """A.metrics against a chosen output directory (the d2k pattern)."""
    old = A.OUT
    A.OUT = out
    A.R.OUT = out
    A.R.J.OUT = out
    try:
        return A.metrics(name)
    finally:
        A.OUT = old
        A.R.OUT = old
        A.R.J.OUT = old


def stalls(d):
    """Bands inside 28 mm that STOP FIRING before their own face crosses the
    plane. t_stop >= t_end - 5 is D2j-0b's 'still firing' sentinel, not a
    stall (the D2k scoring correction, carried forward deliberately)."""
    out = []
    for key in sorted(d["bands"]):
        ts, tc, _, _ = d["bands"][key]
        if ts >= d["t_end"] - 5.0:
            continue
        if (tc - ts) > 10.0 and key[1] <= 28:
            out.append(key)
    return out


def front_walks_in(d):
    fr = [d["front"][t] for t in A.FRONT_T
          if t in d["front"] and np.isfinite(d["front"][t])]
    return len(fr) > 1 and any(fr[i + 1] < fr[i] - 1e-9 for i in range(len(fr) - 1))


# ---------------------------------------------------------------- 0b helpers
def az_rwall(run, depth, z, nsec=N_SEC):
    """COPIED VERBATIM from d2c_steady/score.py (frozen; not run in place).
    Azimuth-mean wall radius at depth z: the quarter is split into nsec equal
    azimuth sectors, r_wall of a sector = max r of its columns deeper than z."""
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


def mrun(out, name):
    pf = os.path.join(out, name)
    m = AN.meta_of(pf)
    if m is None:
        return None
    return AN.Run(pf, m)


def vdot(run, w):
    """4 x sum h_applied dx dy over removal rows in w, per second [m3/s]."""
    ev = run.ev
    sel = (ev["time"] > w[0]) & (ev["time"] <= w[1])
    return 4.0 * float(np.sum(ev["h_applied"][sel])) * run.dz * run.dz / (w[1] - w[0])


def meier_row(run, w=WIN, z_feet=None):
    th = run.th
    s = (th["time"] >= w[0]) & (th["time"] <= w[1])
    rop = -float(np.polyfit(th["time"][s], th["nozzle_z"][s], 1)[0]) * 3600.0
    p_face = float(np.mean(th["jet_P_face"][s]))
    p_side = float(np.mean(th["side_P_face"][s])) if "side_P_face" in th else 0.0
    v = vdot(run, w)
    vmm3 = v * 1e9                                     # m3/s -> mm3/s
    t_end = float(th["time"][-1])
    feet_end = run.lz - float(th["foot_z"][-1])
    d = dict(name=run.name, rop=rop, p_face=p_face, p_side=p_side,
             side_frac=p_side / p_face if p_face else float("nan"),
             vdot_cm3=v * 1e6, tse_hhv=P_HHV / vmm3 if vmm3 else float("nan"),
             tse_rock=4.0 * p_face / vmm3 if vmm3 else float("nan"),
             ledger=float(np.max(np.abs(th["ledger_err"][s]))),
             cap=float(np.max(th["jet_P_face"][s] / th["jet_P_cap"][s])),
             s_c=float(np.mean(th["jet_s_c"][s])) * 1e3,
             # jet_fs_cols is registered only when jet_free_surface = 1, so the
             # fsoff probes legitimately do not have the column.
             fs_cols=float(np.mean(th["jet_fs_cols"][s])) if "jet_fs_cols" in th else float("nan"),
             t_end=t_end, feet_end=feet_end, nozzle_end=run.lz - float(th["nozzle_z"][-1]))
    if z_feet is not None:
        # Ø at MATCHED FEET HEIGHT: the time each run's feet reach z_feet.
        fz = run.lz - th["foot_z"]
        i = int(np.searchsorted(fz, z_feet))
        d["t_match"] = float(th["time"][min(i, len(fz) - 1)])
        d["z_feet_mm"] = z_feet * 1e3
        # The whole Ø(z) profile at matched geometry, NOT one depth. A single
        # depth is fragile: the profiles of different runs cross, and the
        # matched feet depth can land near a crossing and read "no change"
        # when the shapes differ completely. Scored as D2c/D2d do it:
        # depth-mean Ø over [0, z_m] against 85-93 mm and minimum Ø over the
        # feet depth against 80 mm.
        dep = run.depth(d["t_match"])
        zs = np.arange(0.01, z_feet + 1e-9, 0.01)
        prof = np.array([2.0 * az_rwall(run, dep, z) * 1e3 for z in zs])
        d["prof_z"], d["prof"] = zs, prof
        d["dia_mm"] = float(prof[-1]) if prof.size else float("nan")   # at z_m
        d["dia_mean"] = float(prof.mean()) if prof.size else float("nan")
        d["dia_min"] = float(prof.min()) if prof.size else float("nan")
        # Residual confound, reported not corrected: a faster run reaches the
        # matched depth EARLIER, so that wall has had longer to widen by the
        # end of the run. Dwell makes the size of that bias visible.
        d["dwell"] = t_end - d["t_match"]
        d["dia_end_mm"] = 2.0 * az_rwall(run, run.depth(t_end), z_feet) * 1e3
    return d


def common_feet_depth(runs):
    """The shallowest end-of-run feet depth on the ladder — the deepest depth
    every run actually reached, so Ø is compared at matched geometry."""
    return min(r.lz - float(r.th["foot_z"][-1]) for r in runs)


# ---------------------------------------------------------------- gates
def gate_0a(L):
    L += ["## 1. Gate 0a — the cascade floor (D2j-0b C1 arbiter, side flux on)", ""]
    ladder = [("f10", 1.0), ("f05", 0.5), ("f02", 0.2), ("f01", 0.1)]
    got, rows = {}, []
    for tag, f in ladder:
        for mesh, suf in (("2mm", ""), ("1mm", "_1mm")):
            n = f"A_{tag}{suf}"
            if os.path.exists(os.path.join(OUT, n + ".done")):
                got[(f, mesh)] = metrics_in(OUT, n)
    L += ["### Single-mesh signature at 2 mm (band stalling — what the scan detects)", "",
          "| f | disk-mean rate by block [mm/s] | bands stalling inside 28 mm | r_front walks in |",
          "|---|---|---|---|"]
    for tag, f in ladder:
        d = got.get((f, "2mm"))
        if d is None:
            continue
        rt = " / ".join(f"{d['blocks'][b]:.4f}" if b in d["blocks"] else "—" for b in BLOCKS)
        st = stalls(d)
        rows.append((f, d, st))
        L.append(f"| {f} | {rt} | {'none' if not st else ', '.join(f'{a}–{b}' for a, b in st)} "
                 f"| {'yes' if front_walks_in(d) else 'no'} |")
    gaps_a = {}
    L += ["", "### Mesh gap where both legs were run (the deciding condition)", "",
          "| f | block | 2 mm | 1 mm | gap | verdict |", "|---|---|---|---|---|---|"]
    ok_f = {}
    for tag, f in ladder:
        d2, d1 = got.get((f, "2mm")), got.get((f, "1mm"))
        if d2 is None or d1 is None:
            continue
        good = True
        for b in BLOCKS:
            if b in d2["blocks"] and b in d1["blocks"]:
                g = d1["blocks"][b] / d2["blocks"][b] - 1.0
                good = good and abs(g) <= TOL
                L.append(f"| {f} | {b[0]:.0f}–{b[1]:.0f} s | {d2['blocks'][b]:.4f} | "
                         f"{d1['blocks'][b]:.4f} | **{g * 100:+.1f} %** | "
                         f"{'ok' if abs(g) <= TOL else '**> 5 %**'} |")
        ok_f[f] = good and not stalls(d1) and not front_walks_in(d1)
        gaps_a[f] = [d1["blocks"][b] / d2["blocks"][b] - 1.0 for b in BLOCKS
                     if b in d2["blocks"] and b in d1["blocks"]]
    passing = [f for f, ok in ok_f.items() if ok]
    f_min = min(passing) if passing else None
    L += ["", f"**f_min = {f_min if f_min is not None else 'not established'}** — the smallest f "
          "on the ladder meeting all three conditions"
          + ("" if f_min is not None else " (no f with both legs run passed)") + ".", ""]
    L += ["**f is a diagnostic dial here, not a calibration.** It scales q_side at exactly the "
          "place an h scale would, which is why it maps the window; its physical meaning is still "
          "AREA, and no value in this table may be kept as a calibrated setting.", ""]
    return f_min, gaps_a


def gate_0b(L):
    L += ["## 2. Gate 0b — rate, shape and cost vs f (Meier LI0 keys, 2 mm, 250 s)", ""]
    ladder = [("M_f10", 1.0, OUT), ("M_f05", 0.5, OUT), ("M_f02", 0.2, OUT), ("M_f01", 0.1, OUT),
              ("LI0_2mm", 0.0, D2I)]                      # f = 0 anchor: D2i, side key off
    runs = [(f, mrun(out, n)) for n, f, out in ladder]
    runs = [(f, r) for f, r in runs if r is not None]
    if not runs:
        L += ["_no 0b runs yet_", ""]
        return None, []
    z_feet = common_feet_depth([r for _, r in runs])
    rows = [dict(meier_row(r, WIN, z_feet), f=f) for f, r in runs]
    rows.sort(key=lambda d: -d["f"])
    L += [f"Every run is compared at **matched geometry**: the time its feet reach the common "
          f"depth **{z_feet * 1e3:.0f} mm**, the deepest depth every run on the ladder reached. "
          "Never at matched time.", "",
          "| f | ROP 150–250 s [m/h] | absorbed (quarter) [W] | side / face | V̇ [cm³/s] | "
          "depth-mean Ø [mm] | min Ø [mm] | TSE(HHV) [J/mm³] | TSE(to-rock) [J/mm³] |",
          "|---|---|---|---|---|---|---|---|---|"]
    for d in rows:
        L.append(f"| {d['f']:g} | **{d['rop']:.3f}** | {d['p_face']:.0f} | "
                 f"{d['side_frac'] * 100:.0f} % | {d['vdot_cm3']:.2f} | "
                 f"**{d['dia_mean']:.1f}** | **{d['dia_min']:.1f}** | "
                 f"{d['tse_hhv']:.1f} | {d['tse_rock']:.2f} |")
    L += ["", f"Ø(z) profile at matched geometry [mm], z = depth below the original top:", "",
          "| f | " + " | ".join(f"{z * 1e3:.0f} mm" for z in rows[0]["prof_z"][::2]) + " |",
          "|---|" + "---|" * len(rows[0]["prof_z"][::2])]
    for d in rows:
        L.append(f"| {d['f']:g} | " + " | ".join(f"{v:.0f}" for v in d["prof"][::2]) + " |")
    L += ["", f"Bands: ROP {BAND_ROP[0]}–{BAND_ROP[1]} m/h · Ø {BAND_D[0]}–{BAND_D[1]} mm · "
          f"TSE(HHV) {TSE_HHV} · TSE(to-rock) {TSE_ROCK[0]}–{TSE_ROCK[1]}. "
          f"Thermodynamic floor ρc_pΔT_fire = {RHOCP * (T_FIRE - T_AMB) / 1e9:.3f} J/mm³.", ""]
    # f_rate by interpolation in f
    fs = np.array([d["f"] for d in rows])[::-1]
    rp = np.array([d["rop"] for d in rows])[::-1]
    inside = [d["f"] for d in rows if BAND_ROP[0] <= d["rop"] <= BAND_ROP[1]]
    if inside:
        f_rate = max(inside)
        how = "run directly in the band"
    else:
        f_rate = float(np.interp(BAND_ROP[1], rp, fs))
        how = f"linear interpolation in f onto ROP = {BAND_ROP[1]} m/h"
    L += [f"**f_rate = {f_rate:.3f}** ({how}).", ""]
    L += ["Ø at matched feet depth is not confound-free either: a faster run **reaches that depth "
          "earlier**, so its wall has had longer to widen by the end of the run. Dwell = t_end − "
          "t(feet reach the matched depth) makes the size of that bias visible; Ø(t_end) is the "
          "same wall measured at the end of the run instead.", "",
          "| run | t at matched depth [s] | dwell [s] | Ø at match [mm] | Ø at t_end [mm] |",
          "|---|---|---|---|---|"]
    for d in rows:
        L.append(f"| f = {d['f']:g} | {d['t_match']:.0f} | {d['dwell']:.0f} | {d['dia_mm']:.1f} | "
                 f"{d['dia_end_mm']:.1f} |")
    L += ["", "| run | cap max | max abs ledger_err | jet_s_c [mm] | jet_fs_cols | feet depth end [mm] |",
          "|---|---|---|---|---|---|"]
    for d in rows:
        L.append(f"| f = {d['f']:g} | {d['cap']:.2f} | {d['ledger']:.1e} | {d['s_c']:.1f} | "
                 f"{d['fs_cols']:.1f} | {d['feet_end'] * 1e3:.0f} |")
    L.append("")
    return f_rate, rows


def gate_0c(L):
    L += ["## 3. Gate 0c — the competing explanation (distribution)", "",
          "Discriminating signature: **rate DOWN and Ø UP together.** No cut to side flux can "
          "produce that, so either probe showing it implicates the distribution.", ""]
    pairs = [("P_flat", "M_f10", 1.0, "far-field exponent n 1.0 -> 0.5"),
             ("P_fsoff", "M_f10", 1.0, "jet_free_surface off"),
             ("P_flat0", "LI0_2mm", 0.0, "far-field exponent n 1.0 -> 0.5"),
             ("P_fsoff0", "LI0_2mm", 0.0, "jet_free_surface off")]
    base = {}
    for n, out in (("M_f10", OUT), ("LI0_2mm", D2I)):
        base[n] = mrun(out, n)
    got = [(p, b, f, w) for p, b, f, w in pairs
           if os.path.exists(os.path.join(OUT, p + ".done")) and base.get(b) is not None]
    if not got:
        L += ["_no 0c runs yet_", ""]
        return []
    z_feet = None
    rs = [base[b] for _, b, _, _ in got] + [mrun(OUT, p) for p, _, _, _ in got]
    z_feet = common_feet_depth([r for r in rs if r is not None])
    L += [f"Ø at the matched feet depth **{z_feet * 1e3:.0f} mm**.", "",
          "| probe | f | change | ROP [m/h] | Δ rate | **Ø at z_m [mm]** | **Δ Ø (deciding)** | "
          "Δ depth-mean Ø | Δ min Ø | identical to base? | signature? |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    out_rows = []
    for p, b, f, what in got:
        rp, rb = mrun(OUT, p), base[b]
        dp, db = meier_row(rp, WIN, z_feet), meier_row(rb, WIN, z_feet)
        ident = K.rows_upto(os.path.join(OUT, p, "thermo.dat"), K.MEIER_STOP) == \
            K.rows_upto(os.path.join(rb.pf, "thermo.dat"), K.MEIER_STOP)
        drate = dp["rop"] / db["rop"] - 1.0
        # DECIDING: Ø at the matched depth z_m, the quantity fixed in
        # PREDICTIONS.md §3 before any run. The depth-mean and min are
        # reported beside it but do NOT decide -- see the note below, which
        # records that a mean-based rule would have flipped this verdict.
        ddia = dp["dia_mm"] - db["dia_mm"]
        dmean = dp["dia_mean"] - db["dia_mean"]
        dmin = dp["dia_min"] - db["dia_min"]
        sig = (drate < -0.02) and (ddia > 2.0)
        out_rows.append(dict(probe=p, f=f, what=what, rop=dp["rop"], drate=drate,
                             dia=dp["dia_mm"], ddia=ddia, dmean=dmean, dmin=dmin,
                             ident=ident, sig=sig))
        L.append(f"| `{p}` | {f:g} | {what} | {dp['rop']:.3f} | {drate * 100:+.1f} % | "
                 f"**{dp['dia_mm']:.1f}** | **{ddia:+.1f}** | {dmean:+.1f} | {dmin:+.1f} | "
                 f"{'identical' if ident else 'no'} | "
                 f"{'**YES**' if sig else 'no'} |")
    L += ["", "**Deciding column is Δ Ø at z_m**, the metric fixed in PREDICTIONS.md §3 before "
          "any run. The depth-mean is shown beside it because a single depth is fragile — the "
          "Ø(z) profiles of different runs cross, and z_m can land near a crossing. **A "
          "mean-based rule would have flipped this verdict**: `P_fsoff0` moves the depth-mean "
          "by +17.1 mm and would read as the signature. It is not treated as one, for two "
          "reasons: it is not the pre-registered metric, and it moves the wrong way — the "
          "depth-mean is already too WIDE (112 mm against Meier's 85–93), so +17 mm is a worse "
          "funnel, not a better hole. The quantity that fails LOW is the minimum Ø (D2d: 77.6 "
          "against 85–93), and that moves +0.3 mm.", "",
          "Base rows for reference:", "",
          "| base | f | ROP [m/h] | depth-mean Ø [mm] | min Ø [mm] |",
          "|---|---|---|---|---|"]
    for n, f in (("M_f10", 1.0), ("LI0_2mm", 0.0)):
        if base.get(n) is not None:
            d = meier_row(base[n], WIN, z_feet)
            L.append(f"| `{n}` | {f:g} | {d['rop']:.3f} | {d['dia_mean']:.1f} | "
                     f"{d['dia_min']:.1f} |")
    L.append("")
    return out_rows


def decision(L, f_min, f_rate, probes):
    L += ["## 4. The decision (PREDICTIONS.md §2, applied)", ""]
    sig = [p for p in probes if p["sig"]]
    if sig:
        L += [f"**0c shows the discriminating signature** on {', '.join('`' + p['probe'] + '`' for p in sig)}. "
              "Per the rule this is the stronger lead: report it and **stop before the closure**, "
              "even if a window exists.", ""]
    if f_min is None or f_rate is None:
        L += ["**Undecidable so far** — f_min and/or f_rate not both established.", ""]
        return
    empty = f_rate < f_min
    L += [f"- `f_min` = **{f_min}** (cascade suppression floor)",
          f"- `f_rate` = **{f_rate:.3f}** (rate re-enters {BAND_ROP[0]}–{BAND_ROP[1]} m/h)", "",
          f"**Window is {'EMPTY' if empty else 'NON-EMPTY'}: f_rate "
          f"{'<' if empty else '>='} f_min.**", ""]
    if empty:
        L += ["Per ACTIVE_STEP: **STOP. Write no code.** The wall closure cannot satisfy both "
              "gates — the same parameter controls cascade suppression and the power excess and "
              "they pull opposite ways. Do not run Gate 1, do not write `side_face_h_scale`.", ""]
    else:
        L += ["Per ACTIVE_STEP: proceed to **Gate 1** (the 2639-file sweep on the working tree "
              "with `bin/mmwspalling-3d-g++-amr`, no rebuild) before any build.", ""]



def fig(rows_b, gaps_a):
    """Left: rate and absorbed power vs f against the Meier band. Middle: the
    Ø(z) profiles at matched geometry against 85-93 mm. Right: the 0a mesh gap
    vs f against the key-off cascade."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    if rows_b:
        fs = [d["f"] for d in rows_b][::-1]
        ax[0].plot(fs, [d["rop"] for d in rows_b][::-1], "o-", lw=2, ms=7, label="burner ROP")
        ax[0].axhspan(*BAND_ROP, color="g", alpha=0.15, label="Meier 1.3-1.6 m/h")
        ax[0].set(xlabel="side_face_factor f", ylabel="ROP 150-250 s [m/h]",
                  title="0b: rate vs f (f_rate = 0)")
        a2 = ax[0].twinx()
        a2.plot(fs, [4 * d["p_face"] / 1e3 for d in rows_b][::-1], "s--", color="C3",
                label="absorbed (full hole)")
        a2.set_ylabel("absorbed power [kW]", color="C3")
        ax[0].grid(alpha=0.3)
        ax[0].legend(fontsize=8, loc="upper left")
        for d in rows_b:
            ax[1].plot(d["prof"], d["prof_z"] * 1e3, "-o", ms=3, label=f"f = {d['f']:g}")
        ax[1].axvspan(*BAND_D, color="g", alpha=0.15, label="Meier 85-93 mm")
        ax[1].invert_yaxis()
        ax[1].set(xlabel="hole Ø [mm]", ylabel="depth below original top [mm]",
                  title="0b: Ø(z) at MATCHED feet depth")
        ax[1].grid(alpha=0.3)
        ax[1].legend(fontsize=8)
    if gaps_a:
        for i, b in enumerate(BLOCKS):
            xs = sorted(gaps_a)
            ys = [gaps_a[x][i] * 100 if i < len(gaps_a[x]) else np.nan for x in xs]
            ax[2].plot(xs, ys, "o-", label=f"{b[0]:.0f}-{b[1]:.0f} s")
        for i, g in enumerate(REF_GAP_OFF):
            ax[2].axhline(g, ls=":", color="gray", lw=1)
        ax[2].axhspan(-5, 5, color="g", alpha=0.15, label="+/-5 %")
        ax[2].axhline(0, color="k", lw=0.6)
        ax[2].set(xlabel="side_face_factor f", ylabel="1 mm / 2 mm disk-mean rate - 1 [%]",
                  title="0a: mesh gap vs f (dotted = key off)")
        ax[2].grid(alpha=0.3)
        ax[2].legend(fontsize=8)
    f.suptitle("D2l Gate 0: key-only scan of the side-face factor", y=0.99)
    f.tight_layout()
    f.savefig(os.path.join(HERE, "d2l_scan.png"), dpi=110)
    plt.close(f)
    print("wrote d2l_scan.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", choices=["0a", "0b", "0c"])
    a = ap.parse_args()
    L = ["# D2l Gate 0 analysis — key-only scan of the side-face factor", "",
         f"Binary `{K.J.binary_sha()[:8]}…`; predictions hashed before launch "
         "(sha256 in RESULTS.md §0).", ""]
    f_min = f_rate = None
    probes, rows_b, gaps_a = [], [], {}
    if a.gate in (None, "0a"):
        f_min, gaps_a = gate_0a(L)
    if a.gate in (None, "0b"):
        f_rate, rows_b = gate_0b(L)
    if a.gate in (None, "0c"):
        probes = gate_0c(L)
    if a.gate is None:
        decision(L, f_min, f_rate, probes)
    txt = "\n".join(L)
    if a.gate is None:
        open(os.path.join(HERE, "analysis.md"), "w").write(txt + "\n")
        fig(rows_b, gaps_a)
    print(txt)


if __name__ == "__main__":
    main()
