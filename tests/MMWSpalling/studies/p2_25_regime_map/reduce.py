#!/usr/bin/env python3
"""P2-25 -- RESULTS.txt: the regime map from thermo.dat, the events CSV and the metadata of each leg (no plotfile but
P400R's last one, for L10's volume check), the gates' logs, and the e0 verdict. Usage: reduce.py > RESULTS.txt

Units: times s; q MW/m2 (the ABSORBED centre flux (1 - R) 2P/(pi w0^2)); P W; alpha 1/m; depths and lengths m printed
mm; T K; volumes m3; energy J; specific energy GJ/m3 (= J/m3 / 1e9).
Timing (SR 7(iii)): a thermo row stamped t is the step that ENDED at t; an events-CSV row stamps the step's START, so
the first regime-1 row of a leg is one step before its first thermo row with removals; the last step has a CSV row
and no thermo row, so the two counts are compared over the steps both record.
OUTCOME (the packet's declared thresholds): S SPALL = a removal, before any melt, and N_rem(stop) >= 20; St STALL =
1 <= N_rem(stop) < 20 with the first removal before any melt; M MELT-FIRST = a melt and no removal before it; B BLOCKED =
fired, no removal, no melt by stop; C COLD = no firing by stop -- UNDECIDED, never "no spall".
Removed VOLUME = sum of h_applied x dx dy over regime-1 rows (the lines run along grid axes, so the along-line
thickness is the volume per face area); on P400R cross-checked against sum (1 - fr_f) V of its last plotfile (fr_f is
the face list's solid fraction; the packet named `phi`, which is the level set z_surface - z_cell: L10, corrected).
Specific energy THREE ways: beam_E_dep(stop) / V (absorbed), ledger_E_src(stop) / V (net: absorbed minus the surface
losses -- the packet's "absorbed" column, which is net) and P stop / V (beam).
"""
import csv
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LEGS = os.environ.get("P2_25_LEGS", os.path.join(HERE, "legs"))   # the env override: a smoke reduce
P400 = os.path.join(HERE, "..", "p2_23_beam_port", "legs", "P400", "run")
DT_ = ["0p5", "1", "2", "4", "8", "16", "71"]
QT_ = ["0p3", "0p5", "1", "2", "4", "8"]
QV = {"0p3": 0.3, "0p5": 0.5, "1": 1.0, "2": 2.0, "4": 4.0, "8": 8.0}
DV = {"0p5": 0.5, "1": 1.0, "2": 2.0, "4": 4.0, "8": 8.0, "16": 16.0, "71": 71.4}
SIGN = {"S": "S ", "St": "St", "M": "M ", "B": "B ", "C": "C "}


def mm(x):
    return "%.2f" % (1e3 * x) if x is not None and math.isfinite(x) else "-"


def fmt(t):
    return "%.2f" % t if t is not None else "-"


def meta(run):
    kv = {}
    p = os.path.join(run, "metadata")
    if os.path.exists(p):
        for ln in open(p):
            if " = " in ln:
                k, v = ln.rstrip("\n").split(" = ", 1)
                kv[k.strip()] = v.strip()
    return kv


class Leg:
    def __init__(self, name):
        self.name = name
        d = os.path.join(LEGS, name)
        self.run = os.path.join(d, "run")
        self.rc = int(open(os.path.join(d, "rc")).read()) if os.path.exists(os.path.join(d, "rc")) else None
        self.ok = self.rc == 0 and os.path.exists(os.path.join(self.run, "thermo.dat"))
        if not self.ok:
            return
        self.m = meta(self.run)
        L = open(os.path.join(self.run, "thermo.dat")).read().splitlines()
        self.h = [c.strip() for c in L[0].split("\t")]
        rows = [dict(zip(self.h, [float(x) for x in r.split("\t")])) for r in L[1:] if r.strip()]
        self.rows = [r for r in rows if r["time"] > 0.0]
        self.dt = float(self.m.get("timestep", "nan"))
        self.stop = float(self.m.get("stop_time", "nan"))
        lo = [float(x) for x in self.m["geometry.prob_lo"].split()]
        hi = [float(x) for x in self.m["geometry.prob_hi"].split()]
        nc = [int(x) for x in self.m["amr.n_cell"].split()]
        self.dx = [(b - a) / n for a, b, n in zip(lo, hi, nc)]
        self.nz = nc[2]
        self.V = self.dx[0] * self.dx[1] * self.dx[2]
        self.rho = float(self.m.get("microstructure.phase0.grimsel.rho", "2650"))
        self.P = float(self.m.get("beam.schedule.P", "nan").split()[0])
        self.w0 = float(self.m.get("beam.omega0", "nan"))
        self.x0, self.y0 = float(self.m.get("beam.x0", "nan")), float(self.m.get("beam.y0", "nan"))
        er, lam, td = (float(self.m[k]) for k in ("beam.eps_r", "beam.lambda", "beam.tan_delta"))
        self.alpha = 2.0 * math.pi * math.sqrt(er) * (td / lam)
        n = math.sqrt(er)
        self.R = ((n - 1) / (n + 1)) ** 2
        self.q = (1 - self.R) * 2 * self.P / (math.pi * self.w0 ** 2) / 1e6
        evp = self.run + "_removal_events.csv"
        ev = list(csv.DictReader(open(evp))) if os.path.exists(evp) else []
        self.e1 = [e for e in ev if e["regime"] == "1"]
        self.e2 = [e for e in ev if e["regime"] == "2"]
        f = lambda r: r["fr_n_floor"] + r["fr_n_wall"] + r["fr_n_roof"]
        self.t_fire = next((r["time"] for r in self.rows if r["bk_n_eval"] > 0), None)
        self.fire_row = next((r for r in self.rows if r["bk_n_eval"] > 0), None)
        self.t_rem = next((r["time"] for r in self.rows if f(r) > 0), None)
        self.t_rem_csv = min((float(e["time"]) for e in self.e1), default=None)
        mr = next((r for r in self.rows if "bc_z_melt" in r and math.isfinite(r["bc_z_melt"])), None)
        self.t_melt = mr["time"] if mr else None
        self.N = len(self.e1)
        t_last = self.rows[-1]["time"]
        self.N_thermo = sum(f(r) for r in self.rows)
        self.N_csv_common = sum(1 for e in self.e1 if float(e["time"]) < t_last - 1e-9)
        self.vol = sum(float(e["h_applied"]) for e in self.e1) * self.dx[0] * self.dx[1]
        self.first = self.rows[0]
        # b2: the HEATING step -- before the melt AND before the first removal: a removal replaces the central column's
        # top, so bc_T_peak drops by the removed cell's excess in one step (~95 K on D0p5_Q2), a jump dt cannot cause
        t_end = min(t for t in (self.t_melt, self.t_rem, float("inf")) if t is not None)
        dT = lambda rs: max((abs(b["bc_T_peak"] - a["bc_T_peak"]) for a, b in zip(rs, rs[1:])
                             if math.isfinite(a["bc_T_peak"]) and math.isfinite(b["bc_T_peak"])), default=float("nan"))
        self.dTmax = dT([r for r in self.rows if r["time"] < t_end])
        self.dTall = dT(self.rows)
        self.Tmax = max((r["bc_T_peak"] for r in self.rows if math.isfinite(r["bc_T_peak"])), default=float("nan"))
        self.last = self.rows[-1]
        # the outcome
        before_melt = self.t_rem is not None and (self.t_melt is None or self.t_rem < self.t_melt)
        if before_melt and self.N >= 20:
            self.out = "S"
        elif before_melt and 1 <= self.N < 20:
            self.out = "St"
        elif self.t_melt is not None:
            self.out = "M"
        elif self.t_fire is not None:
            self.out = "B"
        else:
            self.out = "C"

    def N_at(self, t):
        return sum(1 for e in self.e1 if float(e["time"]) < t - 1e-9)

    def vol_at(self, t):
        return sum(float(e["h_applied"]) for e in self.e1 if float(e["time"]) < t - 1e-9) * self.dx[0] * self.dx[1]

    def floor_depth(self, t=None):
        ks = []
        for e in self.e1:
            if t is not None and float(e["time"]) >= t - 1e-9:
                continue
            k, nv = int(float(e["k"])), int(float(e["n_voided"]))
            ks.append(k - (nv - 1) if e["dir"] == "5" else k)
        return (self.nz - min(ks)) * self.dx[2] if ks else 0.0

    def outcome_at(self, t):
        """the outcome over the prefix [0, t) (for the F twin's comparison at its own stop)"""
        n = self.N_at(t)
        tr = self.t_rem if (self.t_rem is not None and self.t_rem <= t) else None
        tm = self.t_melt if (self.t_melt is not None and self.t_melt <= t) else None
        tf = self.t_fire if (self.t_fire is not None and self.t_fire <= t) else None
        bm = tr is not None and (tm is None or tr < tm)
        return "S" if bm and n >= 20 else "St" if bm else "M" if tm is not None else "B" if tf is not None else "C"


def census(L, r_max):
    """columns whose centre lies within r_max of the beam axis, per 45-degree sector"""
    nx, ny = int(round(0.1 / L.dx[0])), int(round(0.1 / L.dx[1]))
    s = [0] * 8
    for i in range(nx):
        for j in range(ny):
            x = (i + 0.5) * L.dx[0] - L.x0
            y = (j + 0.5) * L.dx[1] - L.y0
            if math.hypot(x, y) <= r_max:
                s[int((math.atan2(y, x) + math.pi) / (math.pi / 4)) % 8] += 1
    return s


def sectors(L, r_max):
    s = [0] * 8
    for e in L.e1:
        x = (int(e["col_i"]) + 0.5) * L.dx[0] - L.x0
        y = (int(e["col_j"]) + 0.5) * L.dx[1] - L.y0
        if math.hypot(x, y) <= r_max:
            s[int((math.atan2(y, x) + math.pi) / (math.pi / 4)) % 8] += 1
    return s


def txt(f):
    p = os.path.join(HERE, f)
    return open(p).read() if os.path.exists(p) else ""


def main():
    out = []
    p = lambda s="": out.append(s)
    names = sorted(n for n in os.listdir(LEGS) if ".old" not in n and os.path.isdir(os.path.join(LEGS, n)))
    legs = {n: Leg(n) for n in names}
    p("P2-25 the regime map (units: s, MW/m2 absorbed at the centre, W, 1/m, mm, K, m3, GJ/m3)")
    p("=" * 118)
    # ---- the gates' logs
    p("GATES (logs written by chain.sh)")
    gl = {
        "z0": re.search(r"^z0 PASS", txt("Z0.txt"), re.M) is not None,
        "z1": re.search(r"^z1 PASS", txt("Z1.txt"), re.M) is not None,
        "z2": "vapour mass 3.28 g" in txt("Z2.txt") and "vapour mass 0 g" in txt("Z2.txt")
              and "vapour mass 14.4 g" in txt("Z2.txt"),
        "a0": "vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref" in txt("GATE_tree.txt"),
        "a1": "VERDICT: IDENTICAL" in txt("GATE_face_cmp.txt"),
        "a2": all(re.search(r"PASS", txt("A2_%s.log" % t)) is not None and "FAIL  " not in txt("A2_%s.log" % t)
                  for t in ("retired_keys", "latent_heat", "beam_tau", "face_beam", "face_melt", "face_vapor")),
    }
    for k, v in gl.items():
        p("  %-3s %s" % (k, "PASS" if v else "FAIL / not on disk"))
    # ---- b0 / b0r
    P4 = legs.get("P400R")
    p("\nb0 / b0r  P400R (D19: alpha 2000, eps_r 1, 400 W, 20 s) against P2-23's P400")
    b0 = False
    if P4 and P4.ok:
        n20 = P4.N_at(20.0 + 1e-9)
        b0 = P4.t_rem_csv is not None and P4.t_rem_csv <= 2.0 and n20 >= 500
        ref = list(csv.DictReader(open(P400 + "_removal_events.csv"))) if os.path.exists(P400 + "_removal_events.csv") else []
        ref1 = [e for e in ref if e["regime"] == "1"]
        vol_ff = None
        try:
            import numpy as np
            import yt
            yt.set_log_level(50)
            pf = sorted(x for x in os.listdir(P4.run) if x.endswith("cell"))[-1]
            ds = yt.load(os.path.join(P4.run, pf))
            ad = ds.all_data()
            vol_ff = float(np.sum((1.0 - ad["fr_f"].to_value()) * ad["cell_volume"].to_value()))
        except Exception as ex:                                       # noqa
            vol_ff = None
        p("  %s b0: the first regime-1 CSV row %s s (bar <= 2.0), N_rem(20 s) %d (bar >= 500)" % (
            "PASS" if b0 else "FAIL", fmt(P4.t_rem_csv), n20))
        p("  b0r: P2-23 P400 %d events, first %s s; P400R %d, first %s s; P400R max bc_T_peak %.1f K (%s 1410.5 K)" % (
            len(ref1), fmt(min((float(e["time"]) for e in ref1), default=None)), P4.N, fmt(P4.t_rem_csv), P4.Tmax,
            "above" if P4.Tmax > 1410.5 else "below"))
        p("  L10 (P400R): sum h_applied dx dy %.4e m3; sum (1 - fr_f) V of the last plotfile %s m3 (%s); N x dz^3 %.4e "
          "(the upper bound)" % (P4.vol, "%.4e" % vol_ff if vol_ff else "-", "rel %+.2f %%" % (
              100 * (P4.vol / vol_ff - 1)) if vol_ff else "-", P4.N * P4.V))
    else:
        p("  P400R not on disk / rc != 0")
    # ---- b1, b2 per leg
    p("\nb1 / b2 PER LEG: rc; beam_alpha_max on the first lit row vs alpha from the argv (%g); beam_P_inc/P (bar 0.8535 "
      "+- 0.002); bc_z_melt in the header; max per-step d(bc_T_peak) before the melt AND the first removal (b2 bar "
      "12 K; the all-rows maximum in brackets: a removal drops the column's peak in one step)")
    bad_b1, bad_b2 = [], []
    for n, L in legs.items():
        if n == "P400R":
            continue
        if not L.ok:
            p("  %-14s rc %s -- EXCLUDED" % (n, L.rc))
            bad_b1.append(n)
            continue
        a_ok = "%g" % L.first["beam_alpha_max"] == "%g" % L.alpha
        pr = L.first["beam_P_inc"] / L.first["beam_P_beam"] if L.first["beam_P_beam"] > 0 else float("nan")
        ok = a_ok and abs(pr - 0.8535) <= 0.002 and "bc_z_melt" in L.h
        if not ok:
            bad_b1.append(n)
        if not (L.dTmax <= 12.0):          # NaN (no heating row) is not a failure
            if math.isfinite(L.dTmax):
                bad_b2.append(n)
        p("  %-14s rc 0  alpha %s (argv %s) %s  P_inc/P %.4f  bc_z_melt %s  dT/step %.2f K (all rows %.1f)%s" % (
            n, "%g" % L.first["beam_alpha_max"], "%g" % L.alpha, "ok" if a_ok else "MISMATCH", pr, "bc_z_melt" in L.h,
            L.dTmax, L.dTall, "" if L.dTmax <= 12.0 else "  > 12 K"))
    p("  b1: %s; b2: %s" % ("PASS (every leg)" if not bad_b1 else "FAIL on " + " ".join(bad_b1),
                             "every leg <= 12 K" if not bad_b2 else "above 12 K on " + " ".join(bad_b2)))
    # ---- per-leg table
    p("\nPER LEG: outcome; t_fire, t_rem (thermo) / first regime-1 CSV row (+dt = t_rem); t_melt; N_rem(stop) (CSV) and "
      "the thermo sum over the steps both record; stop; bc_T_peak, bc_ell_p, bc_ell_pk at t_fire; N_vap")
    for n, L in legs.items():
        if not L.ok:
            continue
        fr = L.fire_row
        p("  %-14s %-2s fire %-7s rem %-7s (csv %-7s) melt %-7s N %5d (thermo %d / csv %d) stop %5.0f  T %s ell_p %s "
          "ell_pk %s  N_vap %d" % (
              n, L.out, fmt(L.t_fire), fmt(L.t_rem), fmt(L.t_rem_csv), fmt(L.t_melt), L.N, L.N_thermo, L.N_csv_common,
              L.stop, "%.0f" % fr["bc_T_peak"] if fr else "-", mm(fr["bc_ell_p"]) if fr else "-",
              mm(fr["bc_ell_pk"]) if fr else "-", len(L.e2)))
    # ---- m1 the map
    p("\nm1 THE MAP (rows 1/alpha, columns q absorbed; each cell: OUTCOME t_fire/t_rem/t_melt N_rem stop)")
    p("  %-8s" % "1/a \\ q" + "".join("%-30s" % ("q %g" % QV[q]) for q in QT_))
    nC = 0
    for d in DT_:
        cells = []
        for q in QT_:
            L = legs.get("D%s_Q%s" % (d, q))
            if not L or not L.ok:
                cells.append("%-30s" % "(not run)")
                continue
            nC += L.out == "C"
            cells.append("%-30s" % ("%s %s/%s/%s N%d s%.0f" % (SIGN[L.out], fmt(L.t_fire), fmt(L.t_rem), fmt(L.t_melt),
                                                               L.N, L.stop)))
        p("  %-8s" % ("%g mm" % DV[d]) + "".join(cells))
    p("  letters only:")
    for d in DT_:
        p("    %-7s %s" % ("%g mm" % DV[d], " ".join("%-2s" % (legs["D%s_Q%s" % (d, q)].out
                                                         if "D%s_Q%s" % (d, q) in legs and legs["D%s_Q%s" % (d, q)].ok
                                                         else "--") for q in QT_)))
    p("  C cells (undecided): %d (P1: <= 6)" % nC)
    # ---- m2 the boundary
    p("\nm2 THE BOUNDARY")
    for q in QT_:
        S = [d for d in DT_ if legs.get("D%s_Q%s" % (d, q)) and legs["D%s_Q%s" % (d, q)].ok
             and legs["D%s_Q%s" % (d, q)].out == "S"]
        nonS = [d for d in DT_ if d not in S and legs.get("D%s_Q%s" % (d, q)) and legs["D%s_Q%s" % (d, q)].ok]
        fn = legs["D%s_Q%s" % (nonS[0], q)] if nonS else None
        p("  q %-4g deepest S: %-8s  bc_ell_pk at t_fire in the S cells %s; the first non-S cell %s: %s, ell_pk %s" % (
            QV[q], "%g mm" % DV[S[-1]] if S else "none",
            " ".join("%g:%s" % (DV[d], mm(legs["D%s_Q%s" % (d, q)].fire_row["bc_ell_pk"])) for d in S) or "-",
            "%g mm" % DV[nonS[0]] if nonS else "-", fn.out if fn else "-",
            mm(fn.fire_row["bc_ell_pk"]) if fn and fn.fire_row else "-"))
    for d in DT_:
        S = [q for q in QT_ if legs.get("D%s_Q%s" % (d, q)) and legs["D%s_Q%s" % (d, q)].ok
             and legs["D%s_Q%s" % (d, q)].out == "S"]
        p("  1/alpha %-5g lowest-q S: %s" % (DV[d], "%g MW/m2" % QV[S[0]] if S else "none"))
    # ---- m3 the K cuts
    p("\nm3 THE K CUTS at q 2 (kic x0.5 / x1 / x2): outcome, N_rem, bc_T_peak at t_fire, t_fire")
    for d in DT_:
        row = []
        for k in ("K05_D%s_Q2", "D%s_Q2", "K2_D%s_Q2"):
            L = legs.get(k % d)
            row.append("%-2s N%-5d T %s f %s" % (L.out, L.N, "%.0f" % L.fire_row["bc_T_peak"] if L.fire_row else "-",
                                                 fmt(L.t_fire)) if L and L.ok else "(not run)")
        p("  %-6s %s" % ("%g mm" % DV[d], "  |  ".join("%-28s" % r for r in row)))
    # ---- m4 the W cuts
    p("\nm4 THE SPOT CUTS at q 2 (P x (w/12.88)^2): outcome, N_rem, P_inc/P, the columns inside r <= w (the code "
      "lights EVERY column under a Gaussian, Pi > 0: the census is the spot's footprint), distinct removal columns")
    for d in ("1", "4"):
        for k in ("W6_D%s_Q2", "D%s_Q2", "W19_D%s_Q2"):
            L = legs.get(k % d)
            if not L or not L.ok:
                p("  %-11s (not run)" % (k % d))
                continue
            cols = len({(e["col_i"], e["col_j"]) for e in L.e1})
            p("  %-11s w %.2f mm: %-2s N %5d  P_inc/P %.4f  columns in r <= w %d  removal columns %d" % (
                k % d, 1e3 * L.w0, L.out, L.N, L.first["beam_P_inc"] / L.first["beam_P_beam"], sum(census(L, L.w0)), cols))
    # ---- m5 the F twin
    p("\nm5 THE 0.78 mm TWIN vs the 1.56 mm leg over the twin's own run (the base read as a PREFIX to the twin's stop)")
    for f, b in (("F_D1_Q2", "D1_Q2"), ("F_D2_Q2", "D2_Q2"), ("F_D4_Q2", "D4_Q2"), ("F_D2_Q1", "D2_Q1"),
                 ("F_D2_Q4", "D2_Q4")):
        F, B = legs.get(f), legs.get(b)
        if not (F and F.ok and B and B.ok):
            p("  %-8s (not run)" % f)
            continue
        T = F.stop
        p("  %-8s to %.0f s: 0.78 mm %-2s rem %s melt %s N %d vol %.3e m3 | 1.56 mm %-2s rem %s melt %s N %d vol %.3e "
          "m3 | vol ratio %s" % (
              f, T, F.out, fmt(F.t_rem), fmt(F.t_melt), F.N, F.vol, B.outcome_at(T),
              fmt(B.t_rem if B.t_rem is not None and B.t_rem <= T else None),
              fmt(B.t_melt if B.t_melt is not None and B.t_melt <= T else None), B.N_at(T), B.vol_at(T),
              "%.3f" % (F.vol / B.vol_at(T)) if B.vol_at(T) > 0 else "-"))
    # ---- m6 the S cells
    p("\nm6 THE S CELLS: rate over [stop - 10, stop] (events/s, m3/s); specific energy (GJ/m3) absorbed beam_E_dep / V, "
      "net ledger_E_src / V, beam P stop / V; V_rem and N x dz^3 (the upper bound); the floor depth at stop")
    for n, L in legs.items():
        if not L.ok or L.out != "S":
            continue
        n10 = sum(1 for e in L.e1 if L.stop - 10.0 <= float(e["time"]) < L.stop)
        v10 = sum(float(e["h_applied"]) for e in L.e1 if L.stop - 10.0 <= float(e["time"]) < L.stop) * L.dx[0] * L.dx[1]
        Edep, Enet = L.last["beam_E_dep"], L.last["ledger_E_src"]
        p("  %-14s %.1f ev/s %.3e m3/s; E/V %.2f (abs) %.2f (net) %.2f (beam); V %.3e m3, N dz^3 %.3e; floor %s mm" % (
            n, n10 / 10.0, v10 / 10.0, Edep / L.vol / 1e9, Enet / L.vol / 1e9, L.P * L.stop / L.vol / 1e9, L.vol,
            L.N * L.V, mm(L.floor_depth())))
    # ---- m7 isotropy
    S_legs = [L for L in legs.values() if L.ok and L.out == "S" and L.name.startswith("D")]
    pick = max(S_legs, key=lambda L: L.N) if S_legs else None
    p("\nm7 ISOTROPY: regime-1 events with r <= w in 8 azimuth sectors, each over its column census inside r <= w")
    for L in [x for x in (pick, legs.get("P400R")) if x is not None and x.ok]:
        ev, ce = sectors(L, L.w0), census(L, L.w0)
        rate = [e / c for e, c in zip(ev, ce) if c > 0]
        p("  %-14s events %s (sum %d of %d) census %s; per-column max/min %.2f (bar <= 2 at >= 100 events)" % (
            L.name, ev, sum(ev), L.N, ce, max(rate) / min(rate) if rate and min(rate) > 0 else float("inf")))
    # ---- m8 vapour
    p("\nm8 VAPOUR on any leg: N_vap, mass, the first vapour time, the outcome")
    nv = [(n, L) for n, L in legs.items() if L.ok and L.e2]
    for n, L in nv:
        p("  %-14s N_vap %d = %.3g g; first %.2f s; outcome %s" % (n, len(L.e2), 1e3 * len(L.e2) * L.V * L.rho,
                                                                min(float(e["time"]) for e in L.e2), L.out))
    if not nv:
        p("  none")
    # ---- e0
    g = dict(gl, b0=b0, b1=not bad_b1)
    p("\ne0 VERDICT")
    p("   z0-z2, a0-a2, b0, b1 hold => the map is a statement of the validated pieces run together on this binary: %s "
      "(%s)" % ("YES" if all(g.values()) else "NO", ", ".join("%s %s" % (k, "PASS" if v else "FAIL") for k, v in
                                                               g.items())))
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
