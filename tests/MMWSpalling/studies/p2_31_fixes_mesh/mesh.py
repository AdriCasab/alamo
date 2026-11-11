#!/usr/bin/env python3
"""P2-31 stage 4: THE MESH TEST on stage 3's binary (no code change). Legs under S/legs/m/<name>/ (argv.sh, BINARY.txt,
run.log, rc, run/, run_removal_events.csv); a leg already done on this binary with this argv is skipped (resume).
  mesh.py run [names...]   run the legs (all, or the named), serial, -np 4; prune each after it runs (below)
  mesh.py reduce           the tables (one per series), the change on each halving, the adequacy table -> stdout
  mesh.py keys             every series' key differences (only amr.n_cell, timestep, amr.max_grid_size, amr.plot_int
                           and (W) weibull.V0 may differ inside a series; blocking_factor=1 is on all three H legs)
DISK (16 GiB free at z0): legs whose reduction reads only thermo.dat and the event log run with amr.plot_int = -1; the
fixture's plotfile legs keep the plotfiles their windows read (t = 0, 5, 12, 20 s and the last), node plotfiles are not
kept, every file's md5 before the prune is in MANIFEST.md5, the list in PRUNED.txt (P2-25's pattern); an event log over
50 MB is gzipped after the run (its md5 in MANIFEST.md5).
Units: s, K, m (mm where stated), m/s, m3."""
import csv
import glob
import gzip
import hashlib
import importlib.machinery
import importlib.util
import io
import math
import os
import shlex
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
U = os.path.join(T, "unit")
V = os.path.join(T, "validation")
W = os.path.join(T, "studies", "tree_sweep_0922", "witness")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
LEGS = os.path.join(HERE, "legs", "m")
MPI = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4"]
P25D = os.path.join(T, "studies", "p2_25_regime_map")
HOS_R = os.path.join(T, "studies", "p2_16_span_rule", "legs", "HOS_R", "argv.sh")
P29_F0 = {"2": 6.681e-4, "1": 6.132e-4, "05": 6.027e-4}      # P2-9's flat rate, span off (A3_A6_CHECK.txt:89)
B_STOP = 43.2 + 20.0                                          # P2-24 M0's first bc_z_melt (1.5625 mm) + 20 s


def load(path, name):
    ld = importlib.machinery.SourceFileLoader(name, path)
    sp = importlib.util.spec_from_loader(name, ld)
    m = importlib.util.module_from_spec(sp)
    ld.exec_module(m)
    return m


FR = load(os.path.join(U, "face_removal", "test"), "fr_mesh")
FR.OUT = LEGS
FG = load(os.path.join(W, "face_gate.py"), "fg_mesh")
os.environ["P2_25_LEGS"] = LEGS
P25 = load(os.path.join(P25D, "reduce.py"), "p25_mesh")
P25.LEGS = LEGS

P = {"surface_patch.face_list": "1", "spall.buckle_span": "0.012"}
DX = {"2": 0.002, "1": 0.001, "05": 0.0005}
V0C = {"2": "8.0e-9", "1": "1.0e-9", "05": "1.25e-10"}         # the cell's volume (L12)
NOPLOT = {"amr.plot_int": "-1"}


def kv(d):
    return ["%s=%s" % (k, v) for k, v in d.items()]


def fx(dxk, extra, seeds=None):
    """the fixture (unit/face_removal/input) at mesh dxk with `extra` keys and optional seeded void boxes"""
    d = dict(FR.mesh_keys(dxk), **extra)
    if seeds:
        d["spall.seed_void"] = " ".join(str(v) for b in seeds for v in b)
    return [FR.INPUT] + kv(d)


def curved(c, dxk, extra):
    dx = DX[dxk]
    d = dict(FR.MI.curved_keys(c, dx), **extra)
    d["spall.seed_void"] = " ".join(str(v) for b in FR.MI.curved_seeds(c, dx) for v in b)
    return [FR.INPUT] + kv(d)


def p25_tokens(nm):
    toks = shlex.split(open(os.path.join(P25D, "legs", nm, "argv.sh")).read())
    i = next(i for i, t in enumerate(toks) if t.endswith("mmwspalling-3d-g++"))
    return toks[i + 1], [t for t in toks[i + 2:] if not t.startswith(("plot_file=", "stop_time=", "beam.schedule.t_end="))]


def regime(nm, extra):
    deck, ov = p25_tokens(nm)
    d = dict(t.split("=", 1) for t in ov)
    d.update({"stop_time": "20.0", "beam.schedule.t_end": "20.0"}, **NOPLOT)
    d.update(extra)
    return [deck] + kv(d)


def hoeser(n, dt, mgs):
    o = {"amr.n_cell": n, "timestep": dt, "amr.max_grid_size": mgs, "stop_time": "10.0", "amr.plot_int": "-1",
         "amr.thermo.plot_int": "1", "amr.blocking_factor": "1"}
    toks = FG.from_record(HOS_R, os.path.join(LEGS, "_"), o)
    i = next(i for i, t in enumerate(toks) if t.endswith("mmwspalling-3d-g++"))
    return [t for t in toks[i + 1:] if not t.startswith("plot_file=")]


def legs():
    L = {}
    for k in ("2", "1", "05"):
        L["O_" + k] = fx(k, dict(P, **{"spall.enabled": "0", "stop_time": "6.0"}, **NOPLOT))
    for k in ("2", "1", "05"):
        L["F_" + k] = fx(k, P)
        L["F0_" + k] = fx(k, {"surface_patch.face_list": "1"})
    L["F_1_dt"] = fx("1", dict(P, **{"timestep": "0.004", "amr.plot_int": "250"}))
    for k in ("2", "1", "05"):
        L["Wc_" + k] = fx(k, dict(P, **{"weibull.enabled": "1", "weibull.V0": V0C[k], "stop_time": "12.0"}))
    for k in ("2", "05"):
        L["Wf_" + k] = fx(k, dict(P, **{"weibull.enabled": "1", "weibull.V0": "1.0e-9", "stop_time": "12.0"}))
    L["Wc_1_s2"] = fx("1", dict(P, **{"weibull.enabled": "1", "weibull.V0": V0C["1"], "stop_time": "12.0",
                                       "weibull.seed": "54321"}))
    for nm in ("TL35", "TL45"):
        th, ph = FR.TILTS[nm]
        for k in ("2", "1", "05"):
            L["%s_%s" % (nm, k)] = fx(k, dict(P, **{"stop_time": "12.0"}), seeds=FR.tilt_seeds(th, ph, k))
    for k in ("2", "1", "05"):
        L["BORE_" + k] = curved("BORE", k, dict(P, **{"stop_time": "8.0"}, **NOPLOT))
    for k in ("2", "1"):
        L["SPH_" + k] = curved("SPH", k, dict(P, **{"stop_time": "8.0"}, **NOPLOT))
    og = os.path.join(V, "oglesby", "input")
    for nm, n, dt in (("B_3", "32 32 16", "0.2"), ("B_15", "64 64 32", "0.2"), ("B_08", "128 128 64", "0.05")):
        L[nm] = [og] + kv(dict({"amr.n_cell": n, "timestep": dt, "stop_time": repr(B_STOP)}, **NOPLOT))
    for tag, nm in (("R1", "D1_Q2"), ("R2", "D2_Q2"), ("R24", "D2_Q4")):
        L[tag + "_15"] = regime(nm, {})
        L[tag + "_3"] = regime(nm, {"amr.n_cell": "32 32 16", "amr.max_grid_size": "16", "weibull.V0": "3.0517578125e-8"})
    L["R24_08f"] = regime("F_D2_Q4", {"weibull.V0": "3.814697265625e-9"})
    for tag, nm in (("R1", "F_D1_Q2"), ("R2", "F_D2_Q2"), ("R24", "F_D2_Q4")):
        L[tag + "_08"] = regime(nm, {})                      # CONDITIONAL: run only on a failed reproduction
    L["H_2"] = hoeser("25 25 75", "0.008", "25")
    L["H_1"] = hoeser("50 50 150", "0.008", "25")
    L["H_05"] = hoeser("100 100 300", "0.004", "50")
    return L


CONDITIONAL = ("R1_08", "R2_08", "R24_08")
SERIES = {"O": ["O_2", "O_1", "O_05"], "F": ["F_2", "F_1", "F_05"], "F0": ["F0_2", "F0_1", "F0_05"],
          "F_dt": ["F_1", "F_1_dt"], "Wc": ["Wc_2", "Wc_1", "Wc_05"], "Wf": ["Wf_2", "Wc_1", "Wf_05"],
          "W_seed": ["Wc_1", "Wc_1_s2"], "TL35": ["TL35_2", "TL35_1", "TL35_05"], "TL45": ["TL45_2", "TL45_1", "TL45_05"],
          "BORE": ["BORE_2", "BORE_1", "BORE_05"], "SPH": ["SPH_2", "SPH_1"], "B": ["B_3", "B_15", "B_08"],
          "R1": ["R1_3", "R1_15"], "R2": ["R2_3", "R2_15"], "R24": ["R24_3", "R24_15", "R24_08f"],
          "H": ["H_2", "H_1", "H_05"]}
ALLOWED = {"amr.n_cell", "timestep", "amr.max_grid_size", "amr.plot_int", "weibull.V0", "spall.seed_void",
           "geometry.prob_hi", "surface_patch.x0", "surface_patch.y0", "surface_patch.nozzle_z0"}


def sha(p=BIN):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------------- running -----
def prune(d, keep_times):
    """keep the cell plotfiles at the given whole seconds, t = 0 and the last; drop every node plotfile; MANIFEST first"""
    run = os.path.join(d, "run")
    man = []
    for dp, _, fs in os.walk(run):
        for fn in sorted(fs):
            p = os.path.join(dp, fn)
            man.append("%s  %s" % (md5(p), os.path.relpath(p, d)))
    cells = sorted(glob.glob(os.path.join(run, "[0-9]" * 5 + "cell")))
    keep = set()
    if cells:
        keep |= {cells[0], cells[-1]}
        for pf in cells:
            if round(FR.pf_time(pf)) in keep_times and abs(FR.pf_time(pf) - round(FR.pf_time(pf))) < 1e-6:
                keep.add(pf)
    gone = [pf for pf in cells if pf not in keep] + sorted(glob.glob(os.path.join(run, "[0-9]" * 5 + "node")))
    for g in gone:
        shutil.rmtree(g)
    ev = os.path.join(d, "run_removal_events.csv")
    if os.path.exists(ev):
        man.append("%s  %s" % (md5(ev), os.path.relpath(ev, d)))
        if os.path.getsize(ev) > 50 * 2 ** 20:
            with open(ev, "rb") as fi, gzip.open(ev + ".gz", "wb", compresslevel=6) as fo:
                shutil.copyfileobj(fi, fo)
            os.remove(ev)
            gone.append(ev + " (gzipped to .gz)")
    open(os.path.join(d, "MANIFEST.md5"), "w").write("\n".join(man) + "\n")
    open(os.path.join(d, "PRUNED.txt"), "w").write("kept cell plotfiles %s; removed: %s\n" % (
        sorted(os.path.basename(k) for k in keep), [os.path.relpath(g, d) for g in gone]))


def run_leg(name, tail):
    d = os.path.join(LEGS, name)
    argv = MPI + [BIN] + list(tail) + ["plot_file=" + os.path.join(d, "run")]
    a_txt = " ".join(shlex.quote(a) for a in argv) + "\n"
    b = sha()
    if os.path.exists(os.path.join(d, ".done")) and open(os.path.join(d, "argv.sh")).read() == a_txt \
            and open(os.path.join(d, "BINARY.txt")).read().split()[0] == b:
        print("  done  %-9s (kept)" % name, flush=True)
        return 0
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    open(os.path.join(d, "argv.sh"), "w").write(a_txt)
    open(os.path.join(d, "BINARY.txt"), "w").write("%s  bin/mmwspalling-3d-g++\n" % b)
    t0 = time.time()
    with open(os.path.join(d, "run.log"), "w") as lg:
        rc = subprocess.call(argv, cwd=ROOT, stdout=lg, stderr=subprocess.STDOUT)
    open(os.path.join(d, "rc"), "w").write("%d\n" % rc)
    wall = time.time() - t0
    keep = {0, 5, 12, 20} if not name.startswith(("O_", "BORE", "SPH", "B_", "R", "H_")) else set()
    prune(d, keep)
    if rc == 0:
        open(os.path.join(d, ".done"), "w").write("rc 0 wall %.1f\n" % wall)
    for bt in glob.glob(os.path.join(ROOT, "Backtrace.*")):
        if os.path.getmtime(bt) >= t0 - 1:
            shutil.move(bt, os.path.join(d, os.path.basename(bt)))
    print("  ran   %-9s rc=%d wall %6.1f s  (%s)" % (name, rc, wall, time.strftime("%H:%M:%S")), flush=True)
    return rc


def main_run(names):
    L = legs()
    todo = names or [n for n in L if n not in CONDITIONAL]
    print("P2-31 mesh legs on %s: %d to run" % (sha()[:16], len(todo)), flush=True)
    bad = 0
    for n in todo:
        bad += 1 if run_leg(n, L[n]) else 0
    return 1 if bad else 0


def main_keys():
    L = legs()
    bad = 0
    for s, names in SERIES.items():
        D = [dict(t.split("=", 1) for t in L[n][1:]) for n in names]
        keys = set().union(*D)
        diff = sorted(k for k in keys if len({d.get(k) for d in D}) > 1)
        odd = [k for k in diff if k not in ALLOWED and not (s == "W_seed" and k == "weibull.seed")]   # the seed twin
        print("  %-6s %s: differ in %s%s" % (s, " ".join(names), diff, "  ** NOT ALLOWED: %s" % odd if odd else ""))
        bad += 1 if odd else 0
    return 1 if bad else 0


# ------------------------------------------------------------------------------------------------------- reading -----
def events(name):
    d = os.path.join(LEGS, name)
    p = os.path.join(d, "run_removal_events.csv")
    if os.path.exists(p):
        return list(csv.DictReader(open(p)))
    if os.path.exists(p + ".gz"):
        return list(csv.DictReader(io.TextIOWrapper(gzip.open(p + ".gz"), newline="")))
    return []


def thermo(name):
    p = os.path.join(LEGS, name, "run", "thermo.dat")
    rows = [ln.rstrip("\n").split("\t") for ln in open(p) if ln.strip()]
    h = [x.strip() for x in rows[0]]
    return {c: np.array([float(r[i]) if i < len(r) else np.nan for r in rows[1:]]) for i, c in enumerate(h)}


def rc(name):
    p = os.path.join(LEGS, name, "rc")
    return int(open(p).read()) if os.path.exists(p) else None


def rate_win(name, dxk, n, t0, t1):
    """fr.rate() with its window as arguments: the mean normal recession over interior lines between the plotfiles at
    t0 and t1 (lines whose surface is strictly inside at t0 and >= 4 cells off the solid end at t1)"""
    nx, nz = FR.MESH[dxk][0], FR.MESH[dxk][1]
    dx = DX[dxk]
    a = FR.dominant(n)
    mu = abs(n[a])
    sides = [nx, nx, nz]
    Lz = sides[a] * dx
    other = [q for q in range(3) if q != a]
    pos = dict((round(t), p) for t, p in FR.positions(name, dxk, n))
    if t0 not in pos or t1 not in pos:
        return float("nan"), 0
    A, B = pos[t0], pos[t1]
    keep = [k for k in A if all(4 <= k[q] <= sides[other[q]] - 5 for q in range(2))
            and k in B and 0.0 < A[k] / mu < Lz - 1e-12 and B[k] / mu >= 4 * dx]
    if not keep:
        return float("nan"), 0
    return float(np.mean([(A[k] - B[k]) / float(t1 - t0) for k in keep])), len(keep)


def interior(e, nx):
    return 4 <= int(e["col_i"]) <= nx - 5 and 4 <= int(e["col_j"]) <= nx - 5


def med(xs):
    xs = [x for x in xs if x == x]
    return float(np.median(xs)) if xs else float("nan")


def onset(name, dxk):
    ev = events(name)
    if not ev:
        return None
    nx = FR.MESH[dxk][0]
    first = min(ev, key=lambda e: float(e["time"]))
    col = {}
    for e in ev:
        k = (e["col_i"], e["col_j"])
        if k not in col or float(e["time"]) < float(col[k]["time"]):
            col[k] = e
    ci = [e for e in col.values() if interior(e, nx)]
    n_int = (nx - 8) ** 2
    return {"t1": float(first["time"]), "Ts1": float(first["T_surf"]), "l1": float(first["ell_eff"]),
            "t_med": med([float(e["time"]) for e in ci]), "Ts_med": med([float(e["T_surf"]) for e in ci]),
            "l_med": med([float(e["ell_eff"]) for e in ci]), "cols": len(ci), "n_int": n_int, "rows": len(ev)}


def floor_stats(name, dxk, t0=5.0, t1=20.0):
    ev = [e for e in events(name) if t0 <= float(e["time"]) < t1 and interior(e, FR.MESH[dxk][0]) and e["regime"] == "1"]
    by = {}
    for e in ev:
        by.setdefault((e["col_i"], e["col_j"]), []).append(float(e["time"]))
    iv = [b - a for ts in by.values() for a, b in zip(sorted(ts), sorted(ts)[1:])]
    return {"T_s": med([float(e["T_surf"]) for e in ev]), "ell": med([float(e["ell_eff"]) for e in ev]),
            "h_f": med([float(e["h_f"]) for e in ev]), "dt_fire": med(iv), "n": len(ev)}


def first_rows(name, dxk):
    ev = [e for e in events(name) if e["regime"] == "1"]
    if not ev:
        return float("nan"), float("nan")
    nx = FR.MESH[dxk][0]
    col = {}
    for e in ev:
        k = (e["col_i"], e["col_j"])
        col[k] = min(col.get(k, 1e30), float(e["time"]))
    return min(float(e["time"]) for e in ev), med([t for (i, j), t in col.items() if 4 <= int(i) <= nx - 5
                                                     and 4 <= int(j) <= nx - 5])


def ev_volume(name, t, dx):
    return sum(float(e["h_applied"]) for e in events(name) if e["regime"] == "1" and float(e["time"]) < t - 1e-9) * dx * dx


def ledger(name):
    th = thermo(name)
    les = float(np.nanmax(np.abs(th["ledger_err_solid"]))) if "ledger_err_solid" in th else float("nan")
    idr = float(np.nanmax(th["fr_id_remap"])) if "fr_id_remap" in th else float("nan")
    ovl = th.get("fr_overlap")
    return les, idr, (float(np.nansum(ovl)) if ovl is not None else float("nan")), \
        (float(np.nanmax(ovl)) if ovl is not None and len(ovl) else float("nan"))


def hoeser_stats(name, dx, nz):
    ev = [e for e in events(name) if e["regime"] == "1"]
    v5 = sum(float(e["h_applied"]) for e in ev if float(e["time"]) < 5.0 - 1e-9) * dx * dx
    v10 = sum(float(e["h_applied"]) for e in ev) * dx * dx
    vox = set()
    for e in ev:
        nv = int(float(e["n_voided"]))
        if nv <= 0:
            continue
        i, j, k, d = int(e["col_i"]), int(e["col_j"]), int(e["k"]), int(e["dir"])
        a, sg = d // 2, (-1 if d % 2 == 1 else 1)
        for s in range(nv):
            p = [i, j, k]
            p[a] += sg * s
            vox.add(tuple(p))
    if not vox:
        return v5, v10, 0.0, 0.0
    kmin = min(p[2] for p in vox)
    depth = (nz - kmin) * dx
    k_half = nz - 1 - int(round(depth / 2 / dx))
    r_half = max((math.hypot(p[0] + 0.5, p[1] + 0.5) * dx for p in vox if p[2] == k_half), default=0.0)
    return v5, v10, depth, r_half


# ------------------------------------------------------------------------------------------------------- reducing ----
OUT = []


def p(s=""):
    OUT.append(s)


def halving(label, vals, unit, fmt="%.6g"):
    """the three values, the change on each halving, monotone?, the ratio of successive changes, the order"""
    v = list(vals)
    s = "  %-34s %s" % (label, "  ".join((fmt % x if x == x else "nan") for x in v))
    if len(v) >= 2 and all(x == x for x in v):
        ch = ["%+.2f %%" % (100.0 * (b - a) / a if a else float("nan")) for a, b in zip(v, v[1:])]
        s += "   [%s]  changes %s" % (unit, ", ".join(ch))
    p(s)
    return v


ADEQ = []


def adequacy(label, v, cells=("2", "1", "0.5")):
    """e0: from all three cells"""
    if len(v) < 3 or not all(x == x and math.isfinite(x) for x in v):
        ADEQ.append("  %-40s %s" % (label, "not established (fewer than three finite values: %s)" % v))
        return
    d1, d2 = v[1] - v[0], v[2] - v[1]
    last = abs(d2 / v[2]) if v[2] else float("inf")
    mono = (d1 == 0 and d2 == 0) or (d1 * d2 > 0)
    shrink = abs(d2) < abs(d1) or (d1 == 0 and d2 == 0)
    ratio = (d1 / d2) if d2 else float("inf")
    order = math.log2(abs(ratio)) if mono and d2 and d1 else float("nan")
    if last <= 0.03 and mono and shrink:
        verdict = "ADEQUATE at %s mm (last change %.2f %%; error bar +-%.2f %%)" % (cells[1], 100 * last, 100 * last)
    elif last <= 0.03:
        verdict = "within 3 %% on one difference, not confirmed (%s%s)" % ("not monotone" if not mono else "",
                                                                          "" if shrink or not mono else "change not shrinking")
    else:
        verdict = "not established (last change %.2f %% > 3 %%; the finest cell is the only candidate)" % (100 * last)
    ADEQ.append("  %-40s %s | changes %+.3g, %+.3g; ratio %.3g%s" % (
        label, verdict, d1, d2, ratio, "; order %.2f" % order if order == order else ""))


def main_reduce():
    L = legs()
    p("P2-31 stage 4 -- THE MESH TEST (binary %s; legs S/legs/m)" % sha()[:16])
    p("=" * 116)
    # m0
    ran = [n for n in L if os.path.exists(os.path.join(LEGS, n, "rc"))]
    bad = [n for n in ran if rc(n) != 0]
    p("m0 every leg rc 0: %d legs ran, rc != 0 on %s; not run: %s" % (len(ran), bad or "none",
                                                                        [n for n in L if n not in ran]))
    for n in ran:
        if n.startswith(("O_", "F", "W", "TL", "BORE", "SPH")) and rc(n) == 0:
            les, idr, ovs, ovm = ledger(n)
            p("    %-9s max |ledger_err_solid| %.2e, max fr_id_remap %.2e (bar 1e-12), fr_overlap sum %g max %g" % (
                n, les, idr, ovs, ovm))
    p("  m0 F0 (span off, P2-9's form) vs P2-9's 6.681 / 6.132 / 6.027e-4 m/s:")
    z = np.array([0.0, 0.0, 1.0])
    for k in ("2", "1", "05"):
        r, nl = rate_win("F0_" + k, k, z, 5, 20) if rc("F0_" + k) == 0 else (float("nan"), 0)
        p("    F0_%-3s %.4e m/s over %d lines vs P2-9 %.4e (%+.2f %%)" % (k, r, nl, P29_F0[k], 100 * (r / P29_F0[k] - 1)))
    # m1
    p("\nm1 ONSET (evaluate-only, span 0.012): the first event row and the interior columns' first rows (medians)")
    O = [onset("O_" + k, k) if rc("O_" + k) == 0 else None for k in ("2", "1", "05")]
    for key, lab, unit in (("t1", "first row time", "s"), ("Ts1", "first row T_surf", "K"), ("l1", "first row ell", "m"),
                           ("t_med", "column-median first time", "s"), ("Ts_med", "column-median first T_surf", "K"),
                           ("l_med", "column-median first ell", "m")):
        v = halving(lab + " (2 / 1 / 0.5 mm)", [o[key] if o else float("nan") for o in O], unit)
        adequacy("onset " + lab, v)
    p("    interior columns that fired by 6 s: %s" % ", ".join("%s of %s" % (o["cols"], o["n_int"]) if o else "-" for o in O))
    # m2
    p("\nm2 FLAT FLOOR: rate over 5-20 s (rate(), plotfile fr_f) and the event log's interior medians over 5-20 s")
    for s in ("F", "F0"):
        rs = [rate_win("%s_%s" % (s, k), k, z, 5, 20)[0] if rc("%s_%s" % (s, k)) == 0 else float("nan")
              for k in ("2", "1", "05")]
        v = halving("%s rate (2 / 1 / 0.5 mm)" % s, rs, "m/s", "%.5e")
        adequacy("%s flat rate" % s, v)
        st = [floor_stats("%s_%s" % (s, k), k) if rc("%s_%s" % (s, k)) == 0 else None for k in ("2", "1", "05")]
        for key, lab, unit in (("T_s", "median T_surf", "K"), ("ell", "median ell_eff", "m"), ("h_f", "median h_f", "m"),
                               ("dt_fire", "median firing interval", "s")):
            v = halving("%s %s" % (s, lab), [x[key] if x else float("nan") for x in st], unit)
            if s == "F":
                adequacy("F " + lab, v)
    # m3
    p("\nm3 TIME STEP: F_1 (dt 0.008) vs F_1_dt (dt 0.004)")
    r1, r2 = (rate_win(n, "1", z, 5, 20)[0] if rc(n) == 0 else float("nan") for n in ("F_1", "F_1_dt"))
    p("  rate %.5e vs %.5e m/s: %+.3f %% (prediction < 1 %%)" % (r1, r2, 100 * (r2 / r1 - 1)))
    # m4
    p("\nm4 RANDOM FLAWS (weibull.enabled = 1, 12 s): first row overall, column-median first row, rate 5-12 s")
    seed = {}
    for n in ("Wc_1", "Wc_1_s2"):
        if rc(n) == 0:
            f, cm = first_rows(n, "1")
            seed[n] = (f, cm, rate_win(n, "1", z, 5, 12)[0])
    sp = [abs(seed["Wc_1_s2"][i] - seed["Wc_1"][i]) if len(seed) == 2 else float("nan") for i in range(3)]
    p("  the seed twin (Wc_1 seed 12345 vs 54321): |d first| %.4g s, |d column-median| %.4g s, |d rate| %.4g m/s" % tuple(sp))
    for s, names in (("Wc (V0 = the cell)", ("Wc_2", "Wc_1", "Wc_05")), ("Wf (V0 = 1e-9 fixed)", ("Wf_2", "Wc_1", "Wf_05"))):
        vals = []
        for n, k in zip(names, ("2", "1", "05")):
            if rc(n) == 0:
                f, cm = first_rows(n, k)
                vals.append((f, cm, rate_win(n, k, z, 5, 12)[0]))
            else:
                vals.append((float("nan"),) * 3)
        for i, lab, unit in ((0, "first row", "s"), (1, "column-median first row", "s"), (2, "rate 5-12 s", "m/s")):
            v = halving("%s %s" % (s, lab), [x[i] for x in vals], unit)
            res = "resolved" if all(abs(b - a) > sp[i] for a, b in zip(v, v[1:]) if a == a and b == b) else \
                "NOT resolved on >= 1 halving (a change smaller than the seed twin's %.3g)" % sp[i]
            p("      vs the seed spread: %s" % res)
            if s.startswith("Wc"):
                adequacy("Wc " + lab, v)
    if all(rc(n) == 0 for n in ("Wc_2", "Wf_2", "Wc_05", "Wf_05")):
        for k in ("2", "05"):
            a, b = rate_win("Wc_" + k, k, z, 5, 12)[0], rate_win("Wf_" + k, k, z, 5, 12)[0]
            fa, fb = first_rows("Wc_" + k, k)[1], first_rows("Wf_" + k, k)[1]
            p("  at %s mm: Wf / Wc rate %.4f, column-median first row %.3f vs %.3f s (the flaw scale differs by %.4f)" % (
                k, b / a, fb, fa, 8 ** (1 / 13.0) if k == "2" else 8 ** (-1 / 13.0)))
    # m5
    p("\nm5 TILTS AND CURVED against flat AT THE SAME MESH")
    for s in ("TL35", "TL45"):
        ratios = []
        for k in ("2", "1", "05"):
            n = FR.normal(*FR.TILTS[s])
            if rc("%s_%s" % (s, k)) == 0 and rc("F_" + k) == 0:
                rt = rate_win("%s_%s" % (s, k), k, n, 5, 12)[0]
                rf = rate_win("F_" + k, k, z, 5, 12)[0]
                ratios.append(rt / rf)
                les, idr, ovs, ovm = ledger("%s_%s" % (s, k))
                p("    %s_%-3s normal rate %.5e vs flat %.5e m/s, ratio %.4f; fr_overlap sum %g (max %g per step)" % (
                    s, k, rt, rf, rt / rf, ovs, ovm))
            else:
                ratios.append(float("nan"))
        v = halving("%s / flat rate ratio" % s, ratios, "-", "%.4f")
        adequacy("%s / flat ratio" % s, v)
    for c, ks in (("BORE", ("2", "1", "05")), ("SPH", ("2", "1"))):
        rr = []
        for k in ks:
            n = "%s_%s" % (c, k)
            if rc(n) != 0:
                rr.append(float("nan"))
                continue
            dx = DX[k]
            cc = FR.MI.CURVED[c]
            A = 4 * math.pi * cc["R"] ** 2 if c == "SPH" else 2 * math.pi * cc["R"] * cc["depth"] + math.pi * cc["R"] ** 2
            vs = [ev_volume(n, t, dx) for t in (2.0, 4.0, 6.0, 8.0 + 1.0)]
            rv = (vs[3] - vs[1]) / 4.0 / A
            fl = ev_volume("F_" + k, 8.0, dx) - ev_volume("F_" + k, 4.0, dx) if rc("F_" + k) == 0 else float("nan")
            nx = FR.MESH[k][0]
            rf = fl / 4.0 / ((nx * dx) ** 2)
            rr.append(rv / rf)
            les, idr, ovs, ovm = ledger(n)
            p("    %-8s removed volume at 2 / 4 / 6 / 8 s: %s m3; rate per area 4-8 s %.4e m/s vs flat F_%s %.4e (event "
              "log over the top) ratio %.4f; fr_overlap sum %g" % (n, " / ".join("%.4e" % x for x in vs), rv, k, rf,
                                                                    rv / rf, ovs))
        v = halving("%s / flat rate ratio" % c, rr, "-", "%.4f")
        if len(rr) == 3:
            adequacy("%s / flat ratio" % c, v)
    # m6
    p("\nm6 BEAM HEAT-UP (validation/oglesby, #19): melt onset (first finite bc_z_melt, thermo stamp), bc_z_melt and "
      "bc_T_peak at the stop (%.1f s), bc_T_peak at 30 s" % B_STOP)
    Bv = []
    for n in ("B_3", "B_15", "B_08"):
        if rc(n) != 0:
            Bv.append((float("nan"),) * 4)
            continue
        th = thermo(n)
        fin = np.isfinite(th["bc_z_melt"])
        t_on = float(th["time"][fin][0]) if fin.any() else float("nan")
        i30 = int(np.argmin(np.abs(th["time"] - 30.0)))
        Bv.append((t_on, float(th["bc_z_melt"][-1]), float(th["bc_T_peak"][-1]), float(th["bc_T_peak"][i30])))
    for i, lab, unit in ((0, "melt onset", "s"), (1, "bc_z_melt at stop", "m"), (2, "bc_T_peak at stop", "K"),
                         (3, "bc_T_peak at 30 s", "K")):
        v = halving("B %s (3.125 / 1.5625 / 0.78 mm)" % lab, [x[i] for x in Bv], unit)
        adequacy("beam " + lab, v, cells=("3.125", "1.5625", "0.78"))
    # m8
    p("\nm8 THE REGIME BOUNDARY (P2-25's points, 20 s): class, first removal (thermo), melt, N, removed volume")
    for tag, base, fine in (("R1", "D1_Q2", "F_D1_Q2"), ("R2", "D2_Q2", "F_D2_Q2"), ("R24", "D2_Q4", "F_D2_Q4")):
        old = P25.LEGS
        P25.LEGS = os.path.join(P25D, "legs")
        B0, F0 = P25.Leg(base), P25.Leg(fine)
        P25.LEGS = LEGS
        rows = []
        for n in (tag + "_3", tag + "_15", tag + "_08"):
            if rc(n) == 0:
                Lg = P25.Leg(n)
                rows.append((n, Lg.outcome_at(20.0 + 1e-9), Lg.t_rem, Lg.t_melt, Lg.N_at(20.0 + 1e-9), Lg.vol_at(20.0 + 1e-9)))
        rows.append(("P2-25 " + base, B0.outcome_at(20.0), B0.t_rem if B0.t_rem and B0.t_rem <= 20 else None,
                     B0.t_melt if B0.t_melt and B0.t_melt <= 20 else None, B0.N_at(20.0), B0.vol_at(20.0)))
        rows.append(("P2-25 " + fine, F0.outcome_at(20.0 + 1e-9), F0.t_rem, F0.t_melt, F0.N_at(20.0 + 1e-9), F0.vol_at(20.0 + 1e-9)))
        P25.LEGS = old
        for r in rows:
            p("    %-16s %-2s rem %-6s melt %-6s N %5d vol %.4e m3" % (r[0], r[1], P25.fmt(r[2]), P25.fmt(r[3]), r[4], r[5]))
        mine = {r[0]: r for r in rows}
        if tag + "_15" in mine:
            a, b = mine[tag + "_15"], mine["P2-25 " + base]
            same = a[1] == b[1] and (abs(a[4] - b[4]) <= 0.01 * max(b[4], 1))
            p("    %s_15 reproduces P2-25's %s over 20 s: %s (class %s / %s, N %d / %d)%s" % (
                tag, base, "YES" if same else "NO", a[1], b[1], a[4], b[4],
                "" if same else " => run %s_08 (mesh.py run %s_08)" % (tag, tag)))
    if os.path.exists(os.path.join(LEGS, "R24_08f", "rc")) and rc("R24_08f") == 0:
        Lg = P25.Leg("R24_08f")
        p("    R24_08f (F_D2_Q4 with V0 = the 1.5625 mm cell): %s rem %s melt %s N %d vol %.4e" % (
            Lg.outcome_at(20.0 + 1e-9), P25.fmt(Lg.t_rem), P25.fmt(Lg.t_melt), Lg.N_at(20.0 + 1e-9), Lg.vol_at(20.0 + 1e-9)))
    # m7
    p("\nm7 HOESER, SHORT (HOS_R's record, 10 s; quarter domain: volumes x4 for the whole hole)")
    Hs = []
    for n, k, nz in (("H_2", 0.002, 75), ("H_1", 0.001, 150), ("H_05", 0.0005, 300)):
        if rc(n) != 0:
            Hs.append((float("nan"),) * 4)
            continue
        v5, v10, dep, rh = hoeser_stats(n, k, nz)
        les, idr, ovs, ovm = ledger(n)
        Hs.append((v5, v10, dep, rh))
        p("    %-5s V(5 s) %.4e, V(10 s) %.4e m3 (x4: %.4e); depth %.2f mm; radius at half depth %.2f mm; fr_overlap sum %g"
          % (n, v5, v10, 4 * v10, 1e3 * dep, 1e3 * rh, ovs))
    for i, lab, unit in ((0, "V(5 s) quarter", "m3"), (1, "V(10 s) quarter", "m3"), (2, "depth", "m"), (3, "radius at half depth", "m")):
        v = halving("H %s (2 / 1 / 0.5 mm)" % lab, [x[i] for x in Hs], unit)
        adequacy("Hoeser " + lab, v)
    p("\ne0 THE ADEQUACY TABLE (a cell is ADEQUATE when the last halving changes the quantity by <= 3 % AND the series is "
      "monotone with a shrinking change; the last change is its error bar)")
    OUT.extend(ADEQ)
    print("\n".join(OUT))
    return 0


if __name__ == "__main__":
    if sys.argv[1] == "run":
        sys.exit(main_run(sys.argv[2:]))
    if sys.argv[1] == "reduce":
        sys.exit(main_reduce())
    if sys.argv[1] == "keys":
        sys.exit(main_keys())
