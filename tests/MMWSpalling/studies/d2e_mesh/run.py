#!/usr/bin/env python3
"""D2e harness: controlled 2 mm / 1 mm pairs for the mesh-seed diagnosis.

  run.py --cases B0_2mm Pa_2mm ... [--jobs N] [--force]
  run.py --mgs-check          B0 at 2 mm, 16 ms, 30 s: max_grid_size 20 vs 60
  run.py --identity           B0_2mm vs studies/d2c_steady/output/R7b_2mm_dt16
                              (every row with t <= 300 s byte-identical)
  run.py --gate               the 1 mm step gate: B0_1mm vs R7_1mm, 0-150 s
  run.py --kill NAME          stop one run (matches the binary path AND its
                              plot_file, never the bare run name)
  run.py --list               the matrix with dt, K/step and stops

Imports studies/d2c_steady/run.py (D2d's hashed record; not edited) and reuses
its keys(), run_one() and read_thermo(). The D2e policy is applied by
rebinding, as dt_cap.py does:
  * timestep: a fixed step per mesh, DT_MESH = {2 mm: 16 ms, 1 mm: 8 ms}
    (the dt ladder, RESULTS "Timestep ladder"), halved only while the harness
    formula q_max dt/(rho cp dz) exceeds 6 K;
  * no steady-stop (STEADY_LEN = inf); the stall (60 s) and bottom watchdogs
    stay on;
  * stops: P-a 200 s, every other case 300 s.
The geometry is the D2d mesh pair's: 0.40 m full 0.12 x 0.12 quarter domain,
keys = MATRIX["R7b_2mm"] / MATRIX["R7_1mm"] plus one switch per pair (the
switch replaces the scored key in the list, so no key is given twice).
Output: jet profile every 10 s, removal_events.csv on; plotfiles every 100 s
at 2 mm; at 1 mm only at t = 0 and at the stop (a 1 mm plotfile is 1.55 GB;
see RESULTS §2 for why the packet's 150 s interval is not used).
"""
import argparse
import concurrent.futures as cf
import hashlib
import importlib.util
import os
import signal
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D2C = os.path.join(os.path.dirname(HERE), "d2c_steady")
sys.path.insert(0, D2C)
# D2c's run.py under its own module name (this file is also called run.py)
_spec = importlib.util.spec_from_file_location("d2c_run", os.path.join(D2C, "run.py"))
R = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(R)

OUT = os.path.join(HERE, "output")
R.OUT = OUT                        # run_one's makedirs only
R.STEADY_LEN = float("inf")        # no steady-stop (D2e Goal 2)
DT_MESH = {2e-3: 16e-3, 1e-3: 8e-3}
K_MAX = 6.0                        # K/step halving threshold
MGS_FILE = os.path.join(HERE, "mgs.txt")   # written by --mgs-check when 20 == 60
REF_ID = os.path.join(D2C, "output", "R7b_2mm_dt16")
REF_R7 = os.path.join(D2C, "output", "R7_1mm")


def dt_policy(q, dz, cap=None):
    dt = DT_MESH[round(dz, 6)]
    while q * dt / (R.wj.RHOCP * dz) > K_MAX:
        dt *= 0.5
    return dt


R.dt_rule = dt_policy

SWITCH = {
    "B0": ({}, "none (the scored configuration)"),
    "Pa": ({"surface_patch.foot_body_clearance=1": ["surface_patch.foot_body_clearance=0"]},
           "foot_body_clearance = 0"),
    "Pb": ({"surface_patch.jet_T_ent_mode=exhaust": ["surface_patch.jet_T_ent_mode=rec_fixed",
                                                     "surface_patch.jet_T_rec_fixed=1550.0"]},
           "jet_T_ent_mode = rec_fixed, jet_T_rec_fixed = 1550 K"),
    "Pc": ({"surface_patch.foot_rule=pads": ["surface_patch.foot_rule=mean"],
            "surface_patch.foot_body_clearance=1": ["surface_patch.foot_body_clearance=0"]},
           "foot_rule = mean, foot_body_clearance = 0"),
}
STOP = {"B0": 300.0, "Pa": 200.0, "Pb": 300.0, "Pc": 300.0}
MESH = {"2mm": "R7b_2mm", "1mm": "R7_1mm"}
CASES = [f"{p}_{m}" for p in SWITCH for m in MESH]


def mgs():
    try:
        return int(open(MGS_FILE).read().split()[0])
    except (OSError, ValueError, IndexError):
        return 20


def binary_sha():
    h = hashlib.sha256()
    with open(R.BIN, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def job(name, stop=None, grid=None, plot_s=None):
    pair, mesh = name.split("_")
    kw = dict(R.MATRIX[MESH[mesh]])
    kw["stop"] = STOP[pair] if stop is None else stop
    fine = mesh == "1mm"
    kw["plot_s"] = plot_s if plot_s is not None else (kw["stop"] if fine else 100.0)
    k, meta = R.keys(**kw)
    sw, label = SWITCH[pair]
    out = []
    for key in k:
        out += sw.get(key, [key])
    assert all(any(s == key for key in k) for s in sw), f"{name}: a switch key is not in the scored list"
    if grid is None:
        # B0_2mm is the identity check against R7b_2mm_dt16, which ran at 20
        grid = 20 if name == "B0_2mm" else mgs()
    if grid != 20:
        out.append(f"amr.max_grid_size={grid}")
    pf = os.path.join(OUT, name)
    cmd = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", R.BIN,
           os.path.join(R.MEIER, "input_drilling"), f"plot_file={pf}"] + out
    meta.update(name=name, pair=pair, mesh=mesh, switch=label, max_grid_size=grid,
                K_per_step=meta["overshoot"], dt_policy="DT_MESH 16/8 ms, halve above 6 K",
                parent=MESH[mesh], binary_sha256=binary_sha())
    return pf, cmd, meta


def run_all(names, jobs, force):
    js = [job(n) for n in names]
    print(f"{len(js)} runs, {min(jobs, len(js))} in parallel -> {OUT}", flush=True)
    for pf, _, m in js:
        print(f"  {m['name']}: {m['switch']}; dz {m['dz'] * 1e3:g} mm, {m['nxy']}²×{m['nz']}, "
              f"dt {m['dt'] * 1e3:g} ms ({m['K_per_step']:.2f} K/step), stop {m['stop']:g} s, "
              f"mgs {m['max_grid_size']}", flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(jobs, len(js))) as ex:
        for fu in cf.as_completed([ex.submit(R.run_one, pf, c, m, force) for pf, c, m in js]):
            m, status, wall = fu.result()
            print(f"  {m['name']:8s} {status} at t = {m.get('t_end', 0):g} s ({wall / 60:.1f} min)", flush=True)
            ok = ok and not status.startswith("FAILED")
    return ok


def rows_upto(path, tmax, tcol=0, header=1, sep=None):
    """Text rows of a thermo.dat / CSV with time <= tmax (header rows kept)."""
    out = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i < header:
                out.append(line)
                continue
            parts = line.split(sep) if sep else line.split()
            try:
                if float(parts[tcol]) <= tmax + 1e-9:
                    out.append(line)
            except (ValueError, IndexError):
                pass
    return out


def thermo_cols(pa, pb, tmax):
    """{column: max |a - b|} for the thermo columns that differ (t <= tmax)."""
    ra = [l.split() for l in rows_upto(os.path.join(pa, "thermo.dat"), tmax)]
    rb = [l.split() for l in rows_upto(os.path.join(pb, "thermo.dat"), tmax)]
    if len(ra) != len(rb) or ra[0] != rb[0]:
        return None
    a, b = np.array(ra[1:], dtype=float), np.array(rb[1:], dtype=float)
    return {n: float(np.max(np.abs(a[:, i] - b[:, i]))) for i, n in enumerate(ra[0])
            if not np.array_equal(a[:, i], b[:, i])}


def compare(pa, pb, tmax):
    """Byte comparison of thermo / removal events / jet profile rows, t <= tmax.
    Returns [(file, identical, n_rows, max relative difference)]."""
    res = []
    for suf, sep in (("thermo.dat", None), ("_removal_events.csv", ","), ("_jet_profile.csv", ",")):
        fa = os.path.join(pa, suf) if suf == "thermo.dat" else pa + suf
        fb = os.path.join(pb, suf) if suf == "thermo.dat" else pb + suf
        ra, rb = rows_upto(fa, tmax, sep=sep), rows_upto(fb, tmax, sep=sep)
        same = ra == rb
        dmax = 0.0
        if not same:
            if len(ra) != len(rb):
                dmax = float("inf")
            else:
                for la, lb in zip(ra[1:], rb[1:]):
                    xa = np.array((la.split(sep) if sep else la.split()), dtype=float)
                    xb = np.array((lb.split(sep) if sep else lb.split()), dtype=float)
                    if xa.shape != xb.shape:
                        dmax = float("inf")
                        break
                    d = np.abs(xa - xb) / np.maximum(np.maximum(np.abs(xa), np.abs(xb)), 1e-300)
                    dmax = max(dmax, float(np.max(d)))
        res.append((suf.lstrip("_"), same, len(ra) - 1, dmax))
    return res


def mgs_check():
    names = []
    for grid in (20, 60):
        pf, cmd, meta = job("B0_2mm", stop=30.0, grid=grid, plot_s=30.0)
        pf = pf.replace("B0_2mm", f"mgs{grid}_2mm")
        cmd = [c.replace(os.path.join(OUT, "B0_2mm"), pf) for c in cmd]
        meta["name"] = f"mgs{grid}_2mm"
        m, status, wall = R.run_one(pf, cmd, meta, True)
        # run_one polls every 10 s, so time the stepping from the files: the
        # t = 0 plotfile header (written after init) to the last thermo row
        w = (os.path.getmtime(os.path.join(pf, "thermo.dat"))
             - os.path.getmtime(os.path.join(pf, "00000cell", "Header")))
        print(f"  mgs {grid}: {status} at t = {m.get('t_end', 0):g} s, stepping {w:.1f} s "
              f"(run_one {wall:.0f} s)", flush=True)
        names.append((pf, w))
    res = compare(names[0][0], names[1][0], 30.0)
    cols = thermo_cols(names[0][0], names[1][0], 30.0)
    # ledger_err is itself a relative residual (round-off, ~1e-15), so its
    # absolute difference is the relative measure; every other column and
    # both CSVs must be byte-identical or within 1e-12 relative.
    other = {n: d for n, d in (cols or {}).items() if n != "ledger_err"}
    led = (cols or {}).get("ledger_err", 0.0)
    agree = (cols is not None and led <= 1e-12
             and all(same or d <= 1e-12 for f, same, _, d in res if f != "thermo.dat")
             and (not other or res[0][3] <= 1e-12))
    L = ["| file | rows | 20 vs 60 |", "|---|---|---|"]
    for f, same, n, d in res:
        txt = "byte-identical" if same else f"max rel diff {d:.2e}"
        if f == "thermo.dat" and cols:
            txt = "differs only in " + ", ".join(f"{c} (max abs {v:.1e})" for c, v in cols.items())
        L.append(f"| {f} | {n} | {txt} |")
    ratio = names[1][1] / names[0][1]
    L.append(f"\nstepping wall 20: {names[0][1]:.0f} s, 60: {names[1][1]:.0f} s (1 s file-time "
             f"resolution), ratio 60/20 = {ratio:.3f}. Verdict: {'AGREE, use 60' if agree else 'DIFFER, keep 20'}.")
    with open(MGS_FILE, "w") as f:
        f.write(f"{60 if agree else 20}\n")
    print("\n".join(L))
    return agree


def identity():
    pf = os.path.join(OUT, "B0_2mm")
    # compare to B0's last time: stop_time 300 s ends at the last step before
    # it (299.904 s at 16 ms), and the reference ran on past 300 s
    t_end = float(R.read_thermo(pf, ["time"])[-1, 0])
    res = compare(pf, REF_ID, t_end)
    L = [f"| file | rows (t <= {t_end:g} s, B0's end) | B0_2mm vs R7b_2mm_dt16 |", "|---|---|---|"]
    for f, same, n, d in res:
        L.append(f"| {f} | {n} | {'byte-identical' if same else f'DIFFER (max rel {d:.2e})'} |")
    ok = all(s for _, s, _, _ in res)
    L.append(f"\nIdentity: {'PASS' if ok else 'FAIL'}")
    print("\n".join(L))
    return ok


def blocks(pf, edges):
    th = R.read_thermo(pf, ["time", "nozzle_z", "jet_s_c", "jet_T_rec"])
    th = th[th[:, 0] > 0]
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        s = (th[:, 0] >= a) & (th[:, 0] <= b)
        out.append(-np.polyfit(th[s, 0], th[s, 1], 1)[0] * 3600.0 if s.sum() > 2 else np.nan)
    return out, th


def gate():
    pf = os.path.join(OUT, "B0_1mm")
    e = [0.0, 50.0, 100.0, 150.0]
    b, th = blocks(pf, e)
    r, thr = blocks(REF_R7, e)
    L = ["| block | B0_1mm (8 ms) | R7_1mm (4 ms) | diff |", "|---|---|---|---|"]
    ok = True
    for i in range(3):
        d = b[i] / r[i] - 1.0
        ok = ok and abs(d) <= 0.02
        L.append(f"| {e[i]:.0f}–{e[i + 1]:.0f} s | {b[i]:.3f} m/h | {r[i]:.3f} m/h | {d * 100:+.2f} % |")
    for j, name in ((2, "s_c"), (3, "T_rec")):
        vb = float(np.interp(150.0, th[:, 0], th[:, j]))
        vr = float(np.interp(150.0, thr[:, 0], thr[:, j]))
        d = vb / vr - 1.0
        ok = ok and abs(d) <= 0.01
        L.append(f"| {name} at 150 s | {vb:.5g} | {vr:.5g} | {d * 100:+.2f} % |")
    L.append(f"\nStep gate (blocks within 2 %, s_c and T_rec within 1 %): {'PASS' if ok else 'FAIL'}")
    print("\n".join(L))
    return ok


def kill(name):
    """Stop one run: match processes whose argv holds the binary path AND
    plot_file=<its prefix> (a bare-name match also kills shell waiters)."""
    pf = os.path.join(OUT, name)
    tok = f"plot_file={pf}"
    ps = subprocess.run(["ps", "-axo", "pid=,command="], capture_output=True, text=True).stdout
    pids = []
    for line in ps.splitlines():
        pid, _, cmd = line.strip().partition(" ")
        argv = cmd.split()
        if R.BIN in argv and tok in argv and int(pid) != os.getpid():
            pids.append(int(pid))
    for p in pids:
        try:
            os.kill(p, signal.SIGTERM)
        except ProcessLookupError:
            pass
    print(f"SIGTERM to {len(pids)} processes of {name}: {pids}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--mgs-check", action="store_true")
    ap.add_argument("--identity", action="store_true")
    ap.add_argument("--gate", action="store_true")
    ap.add_argument("--kill")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.kill:
        return kill(a.kill)
    if a.list:
        for n in CASES:
            pf, cmd, m = job(n)
            print(f"{n:8s} {m['switch']:52s} dt {m['dt'] * 1e3:g} ms {m['K_per_step']:.2f} K/step "
                  f"stop {m['stop']:g} s plot_int {[c for c in cmd if c.startswith('amr.plot_int')][0]} "
                  f"mgs {m['max_grid_size']}")
        return
    if a.mgs_check:
        sys.exit(0 if mgs_check() else 1)
    if a.identity:
        sys.exit(0 if identity() else 1)
    if a.gate:
        sys.exit(0 if gate() else 1)
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
