#!/usr/bin/env python3
"""D2c deliverable 11: the D2d campaign harness (defined here; D2c runs only
the cost probes).

  run.py --probe 2mm              scored configuration, 2 mm, 30 s (cost probe)
  run.py --probe 1mm              scored configuration, 1 mm, full domain, 10 s
  run.py --project                wall-clock projection from the probes (+ D2d schedule)
  run.py --stage D2d [--cases ..] the D2d matrix (+ the optional runs)
  run.py --write-input            validation/meier/sp_meier_pilot/input_feet_d2c

Common: validation/meier/sp_meier_pilot/input_drilling + CLI overrides (the
D2b keys of d2b_feet_rop/run.py plus the D2c keys); quarter domain 0.12 x
0.12 m, axis at the xlo/ylo corner; beam and bit off; V0 = V_cell; pinned
closure (idle 2); wall treatment A2 (no coherence cap); jet_closure =
enthalpy, decay, entrained mass, exhaust recirculation, T_ent 293.15,
mdot/4, cp 1250, D 7.5 mm.
Domain height (D2d decision 1, before any run): LZ = 0.70 m for the matrix
runs, LZ_MESH = 0.40 m for the mesh pair R7_1mm / R7b_2mm. The hand model
needs 0.65 m for a steady window and runs 13 % fast, so 0.70 m leaves
margin; the mesh pair stays at 0.40 m because it is a matched transient
comparison (decision 2), not a steady one.
Scored configuration (ACTIVE_STEP §Scored configuration):
  Martin h_ref at T_nozzle (walljet.h_anchor), far law power n = 1,
  jet_decay_diameter = momentum with jet_De_ref = nozzle.de_ref(T_nozzle)
  (full jet; not scaled with mdot/4), jet_core_length = 8,
  jet_free_surface = 1 (aspect 1, exponent 1), jet_bin_update = exponential,
  jet_negative_flux = count, feet pads [0.028, 0.040] m, 50 mm, 3 pads,
  quantile 0.9, foot_body_clearance = 1.
Matrix (D2d):
  R1_scored     as pre-registered                       0.70 m, cap 1300 s
  R2_clamp      far = clamp                             0.70 m, cap 1300 s
  R3_n05        n = 0.5: NOT RUN (needs >= 1.55 m; hand model only)
  R4_1600       T_nozzle 1600 K                         0.70 m, cap 1300 s
  R5_1750       T_nozzle 1750 K                         0.70 m, cap 1300 s
  R6_d2bjet_n1  D2b jet closures + n = 1: nozzle D, core 5, free surface off;
                exponential, count and clearance stay on (decision 5)
  R7_1mm        1 mm, full domain                       0.40 m, cap 700 s
  R7b_2mm       mesh partner: same keys, 2 mm           0.40 m, cap 700 s
  Oa_core5 / Ob_ricou / Oc_1450   optional (core 5, core 3.125, T 1450 K)
Watchdogs (poll thermo.dat every 10 s of wall time):
  stalled       nozzle_z has not descended for 60 s of simulated time;
  bottom        centre face (nozzle_z - jet_s_c) <= 40 mm above the domain bottom;
  steady-stop   a steady window (>= 150 s, 50 s descent sub-fits within ±10 %,
                |d jet_s_c/dt| < 0.02 mm/s; score.py's rule on thermo.dat) has
                been found: the run ends 50 s later.
  cap           stop_time per run (CAP_2MM = 1300 s, CAP_1MM = 700 s; the mesh
                pair 700 s, decision 3).
Plotfiles are sparse (PLOT_S_2MM / PLOT_S_1MM): about 0.36 GB per 2 mm
plotfile at 0.70 m and 1.55 GB per 1 mm plotfile, so R7 writes about 4.7 GB.
Time step: the overshoot rule (q dt/(rho Cp dz) <= 5 K, dt = 4 ms halved) at
q = max over s in [50, 90] mm of h(0, s) (T_nozzle - 821) - loss(821), i.e.
the hand transient's maximum T_stag (= T_nozzle, PREDICTIONS.md).
No parameter is adjusted toward Meier's ROP, hole diameter or volume rate.
"""
import argparse
import concurrent.futures as cf
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import handmodel as hm  # noqa: E402
import nozzle as nz  # noqa: E402

fm, wj = hm.fm, hm.wj
ROOT = HERE.split(os.sep + os.path.join("tests", "MMWSpalling") + os.sep)[0]
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
MEIER = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier", "sp_meier_pilot")
OUT = os.path.join(HERE, "output")
D = wj.D
DT_CAP = 4.0e-3
LXY = 0.12
LZ = 0.70                     # D2d decision 1 (matrix runs)
LZ_MESH = 0.40                # D2d decision 2 (the matched-window mesh pair)
CAP_2MM, CAP_1MM = 1300.0, 700.0
PLOT_S_2MM, PLOT_S_1MM = 200.0, 350.0
STALL_S, BOTTOM_M = 60.0, 0.040
STEADY_LEN, SUB, SUB_TOL, DSC_MAX, STEADY_STOP = 150.0, 50.0, 0.10, 2.0e-5, 50.0
WALL_A2 = ["spall.flake_coherence_length=0.0", "spall.surface_normal=0"]
PAR_FACTOR = 0.72             # D2b Stage B: slowest parallel rate / lone probe rate (1.37e7 / 1.91e7)


def dt_rule(q, dz, cap=DT_CAP):
    dt = cap
    while q * dt / (wj.RHOCP * dz) > 5.0:
        dt *= 0.5
    return dt


def q_max(h_ref, T_nozzle, far, n):
    hs = max(float(wj.h_py(h_ref, 0.0, s, far, n)) for s in (0.050, 0.06, 0.07, 0.08, 0.09))
    return hs * (T_nozzle - wj.T_FIRE) - float(wj.loss(wj.T_FIRE))


def keys(T=1900.0, far="power", n=1.0, momentum=True, core=8.0, fs=True, aspect=1.0, exponent=1.0,
         dz=2e-3, lz=LZ, stop=None, plot_s=None, profile_s=10.0):
    fine = dz < 1.5e-3
    if stop is None:
        stop = CAP_1MM if fine else CAP_2MM
    if plot_s is None:
        plot_s = PLOT_S_1MM if fine else PLOT_S_2MM
    h = hm.h_ref_for(T)
    q = q_max(h, T, far, n)
    dt = dt_rule(q, dz)
    nxy, nzc = int(round(LXY / dz)), int(round(lz / dz))
    k = [f"stop_time={stop!r}", f"timestep={dt!r}",
         f"amr.plot_int={int(round(plot_s / dt))}",
         f"amr.thermo.plot_int={max(1, int(round(0.1 / dt)))}", "amr.thermo.int=1",
         f"amr.n_cell={nxy} {nxy} {nzc}", f"geometry.prob_hi={LXY!r} {LXY!r} {lz!r}",
         "bit.enabled=0", "beam.P0=0.0", f"weibull.V0={dz ** 3!r}",
         "surface_patch.enabled=1", "surface_patch.mode=convective_flame",
         "surface_patch.robin_form=pinned", "surface_patch.pinned_idle_cycles=2.0",
         "surface_patch.face=zhi", "surface_patch.x0=0.0", "surface_patch.y0=0.0",
         "surface_patch.radius=0.2",
         f'surface_patch.h_expr="{wj.h_expr(h, far, n)}"',
         "surface_patch.nozzle_descent=feet", "surface_patch.foot_r_inner=0.028",
         "surface_patch.foot_r_outer=0.040", "surface_patch.foot_standoff=0.050",
         "surface_patch.foot_rule=pads", "surface_patch.foot_npads=3",
         "surface_patch.foot_pad_quantile=0.9", "surface_patch.foot_body_clearance=1",
         "surface_patch.jet_closure=enthalpy", "surface_patch.jet_stagnation=decay",
         "surface_patch.jet_entrained_mass=1", "surface_patch.jet_T_ent_mode=exhaust",
         "surface_patch.jet_T_ent=293.15",
         f"surface_patch.jet_mdot={wj.MDOT / 4.0!r}", f"surface_patch.jet_cp={wj.CP!r}",
         f"surface_patch.jet_D={D!r}", f"surface_patch.jet_T_nozzle={T!r}",
         f"surface_patch.jet_core_length={core!r}",
         "surface_patch.jet_bin_update=exponential", "surface_patch.jet_negative_flux=count",
         f"surface_patch.jet_profile_interval={profile_s!r}",
         "surface.follow_mask=1", "energy_ledger.enabled=1", "spall.removal_events_csv=1",
         ] + WALL_A2
    if momentum:
        k += ["surface_patch.jet_decay_diameter=momentum", f"surface_patch.jet_De_ref={nz.de_ref(T)!r}"]
    if fs:
        k += ["surface_patch.jet_free_surface=1", f"surface_patch.jet_fs_aspect={aspect!r}",
              f"surface_patch.jet_fs_exponent={exponent!r}"]
    meta = dict(anchor="M", h_ref=h, T_nozzle=T, mode="exhaust", entrained=1, far=far, n=n,
                momentum=momentum, core=core, fs=fs, aspect=aspect, exponent=exponent, dz=dz,
                nxy=nxy, nz=nzc, lxy=LXY, lz=lz, patch_r=0.2, stop=stop, wall_treatment="A2", dt=dt,
                q_max=q, overshoot=q * dt / (wj.RHOCP * dz), sector=0.25,
                De_ref=nz.de_ref(T) if momentum else None)
    return k, meta


MATRIX = {
    "R1_scored": dict(),
    "R2_clamp": dict(far="clamp"),
    "R4_1600": dict(T=1600.0),
    "R5_1750": dict(T=1750.0),
    # decision 5: only the momentum diameter, the core length and the free
    # surface revert; exponential, count and clearance stay on
    "R6_d2bjet_n1": dict(momentum=False, core=5.0, fs=False),
    "R7_1mm": dict(dz=1e-3, lz=LZ_MESH),
    "R7b_2mm": dict(lz=LZ_MESH, stop=CAP_1MM, plot_s=PLOT_S_1MM),
    "Oa_core5": dict(core=5.0),
    "Ob_ricou": dict(core=3.125),
    "Oc_1450": dict(T=1450.0),
}
# R3_n05 (n = 0.5) is NOT run: PREDICTIONS.md P1 needs a >= 1.55 m domain.
D2D_DEFAULT = ["R1_scored", "R2_clamp", "R4_1600", "R5_1750", "R6_d2bjet_n1",
               "R7b_2mm", "R7_1mm"]


def job(name, **kw):
    k, meta = keys(**kw)
    pf = os.path.join(OUT, name)
    cmd = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN,
           os.path.join(MEIER, "input_drilling"), f"plot_file={pf}"] + k
    meta["name"] = name
    return pf, cmd, meta


def read_thermo(pf, cols):
    p = os.path.join(pf, "thermo.dat")
    try:
        with open(p) as f:
            head = f.readline().split()
            idx = [head.index(c) for c in cols]
            rows = []
            for line in f:
                parts = line.split()
                if len(parts) == len(head):
                    rows.append([float(parts[i]) for i in idx])
        return np.array(rows) if rows else None
    except (OSError, ValueError):
        return None


def steady_found(a):
    """a = [time, nozzle_z, jet_s_c] rows (t > 0). The score's window rule on
    [t0, t_last]: returns t0 of the earliest window of >= 150 s, else None."""
    t, z, sc = a[:, 0], a[:, 1], a[:, 2]
    t_last = t[-1]
    for t0 in np.arange(0.0, t_last - STEADY_LEN + 1e-9, 10.0):
        nsub = int((t_last - t0) // SUB)
        if nsub < 2:
            continue
        subs = []
        for k in range(nsub):
            s = (t >= t_last - SUB * (k + 1)) & (t <= t_last - SUB * k)
            subs.append(-np.polyfit(t[s], z[s], 1)[0])
        mu = np.mean(subs)
        if mu <= 0.0 or np.max(np.abs(np.array(subs) / mu - 1.0)) > SUB_TOL:
            continue
        s = t >= t0
        if abs(np.polyfit(t[s], sc[s], 1)[0]) >= DSC_MAX:
            continue
        return float(t0)
    return None


def run_one(pf, cmd, meta, force):
    done = pf + ".done"
    if os.path.exists(done) and not force:
        return meta, "skip (done)", 0.0
    for p in (pf, pf + "_h_col_events.csv", pf + "_removal_events.csv",
              pf + "_jet_profile.csv", done):
        if os.path.isdir(p):
            shutil.rmtree(p)
        elif os.path.exists(p):
            os.remove(p)
    os.makedirs(OUT, exist_ok=True)
    t0w = time.time()
    status = "ok"
    last_check = 0.0
    steady_at = None
    with open(pf + ".log", "w") as log:
        log.write(" ".join(cmd) + "\n")
        log.flush()
        proc = subprocess.Popen(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
        while proc.poll() is None:
            time.sleep(10.0)
            a = read_thermo(pf, ["time", "nozzle_z", "jet_s_c"])
            if a is None:
                continue
            a = a[a[:, 0] > 0.0]
            if a.shape[0] < 2:
                continue
            t, z, sc = a[-1]
            best = np.minimum.accumulate(a[:, 1])
            last_desc = a[np.nonzero(np.diff(best) < 0)[0][-1] + 1, 0] if np.any(np.diff(best) < 0) else a[0, 0]
            reason = None
            if t - last_desc >= STALL_S:
                reason = "stalled"
            elif sc > 0.0 and z - sc <= BOTTOM_M:
                reason = "bottom"
            elif steady_at is None and t >= STEADY_LEN and t - last_check >= 10.0:
                last_check = t
                w0 = steady_found(a)
                if w0 is not None:
                    steady_at = t
                    log.write(f"\nWATCHDOG: steady window from {w0} s found at t = {t} s; "
                              f"stopping at {t + STEADY_STOP} s\n")
                    log.flush()
            if reason is None and steady_at is not None and t >= steady_at + STEADY_STOP:
                reason = "steady-stop"
            if reason:
                status = reason
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=60)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                log.write(f"\nWATCHDOG: {reason} at t = {t} s\n")
                break
        rc = proc.returncode
    wall_s = time.time() - t0w
    a = read_thermo(pf, ["time"])
    meta = dict(meta, wall_s=wall_s, status=status, t_end=float(a[-1, 0]) if a is not None else 0.0)
    if status != "ok" or rc == 0:
        with open(done, "w") as f:
            f.write(repr(meta) + "\n")
        return meta, status, wall_s
    return meta, f"FAILED rc={rc}", wall_s


def rate_of(meta):
    return meta["nxy"] ** 2 * meta["nz"] * (meta["t_end"] / meta["dt"]) / meta["wall_s"]


def project():
    import ast
    import predict as pr
    probes = {}
    for n in ("probe_2mm", "probe_1mm"):
        p = os.path.join(OUT, n + ".done")
        if os.path.exists(p):
            probes[n] = ast.literal_eval(open(p).read())
    if len(probes) < 2:
        raise SystemExit("run both probes first")
    r2, r1 = rate_of(probes["probe_2mm"]), rate_of(probes["probe_1mm"])
    L = ["| probe | cells | dt | simulated | wall | cell-steps/s | wall s per simulated s |",
         "|---|---|---|---|---|---|---|"]
    for n, m in probes.items():
        L.append(f"| {n} | {m['nxy']}² × {m['nz']} | {m['dt'] * 1e3:g} ms | {m['t_end']:.1f} s | "
                 f"{m['wall_s']:.0f} s | {rate_of(m):.3g} | {m['wall_s'] / m['t_end']:.1f} |")
    L += ["", "Projection: to each run's expected end, the earliest of the hand model's bottom "
          "time at that run's height, its steady-stop and its cap. Alone = the probe rate scaled "
          f"by the cell count; three at a time = × {PAR_FACTOR} (D2b Stage B slowest / lone).", "",
          "| run | dz, height | dt | expected end (hand) | alone | three at a time | to the cap alone |",
          "|---|---|---|---|---|---|---|"]
    rows = {}
    for name, kw in MATRIX.items():
        k, m = keys(**kw)
        c = dict(T=kw.get("T", 1900.0),
                 cfg=hm.cfg_for(kw.get("T", 1900.0), kw.get("far", "power"), kw.get("n", 1.0),
                                kw.get("core", 8.0), kw.get("momentum", True)),
                 fs=hm.FS(kw.get("fs", True), 1.0, 1.0), clear=True)
        res = hm.transient(c["T"], c["cfg"], c["fs"], t_end=m["stop"], clear=True)
        tb = pr.t_reach(res, m["lz"] - BOTTOM_M)          # per-run height (decision 1)
        t_sw = None
        t0, _ = hm.steady_window(res)
        if t0 is not None:
            t_sw = t0 + hm.STEADY_LEN + STEADY_STOP
        ends = [x for x in (tb, t_sw) if x is not None] + [m["stop"]]
        t_exp = min(ends)
        why = "bottom" if t_exp == tb else ("steady-stop" if t_exp == t_sw else "cap")
        rate = r1 if m["dz"] < 1.5e-3 else r2
        cs = m["nxy"] ** 2 * m["nz"] / m["dt"]
        alone = cs * t_exp / rate / 3600.0
        cap = cs * m["stop"] / rate / 3600.0
        rows[name] = (alone, t_exp)
        L.append(f"| {name} | {m['dz'] * 1e3:g} mm, {m['lz']:g} m | {m['dt'] * 1e3:g} ms | "
                 f"{t_exp:.0f} s ({why}) | {alone:.2f} h | {alone / PAR_FACTOR:.2f} h | "
                 f"{cap:.2f} h |")
    k7, m7 = keys(**MATRIX["R7_1mm"])
    cs7 = m7["nxy"] ** 2 * m7["nz"] / m7["dt"]
    t14 = 14.0 * 3600.0 * r1 / cs7
    L += ["", f"Run 7 (1 mm, {m7['lz']:g} m): {t14:.0f} s simulated in 14 h alone; its cap is "
          f"{m7['stop']:g} s (decision 3)."]
    return L, rows, t14


def write_input():
    """validation/meier/sp_meier_pilot/input_feet_d2c = input_drilling with the
    scored D2c keys written in (input_drilling and input_feet untouched)."""
    src = os.path.join(MEIER, "input_drilling")
    dst = os.path.join(MEIER, "input_feet_d2c")
    # input_feet_d2c is the frozen D2c artefact: the scored keys at the D2c
    # height and cap (0.40 m, 1200 s, plot 100 s). The D2d campaign heights and
    # caps live in MATRIX, not here, so regenerating leaves the file unchanged.
    ks, meta = keys(lz=0.40, stop=1200.0, plot_s=100.0)
    ks = ["plot_file=tests/MMWSpalling/validation/meier/sp_meier_pilot/output/feet_d2c"] + ks
    lines = open(src).read().split("\n")
    body = lines[next(i for i, l in enumerate(lines) if l.startswith("alamo.program")):]
    appended = []
    for kv in ks:
        k, v = kv.split("=", 1)
        pat = re.compile(r"^" + re.escape(k) + r"\s*=")
        hit = [i for i, l in enumerate(body) if pat.match(l)]
        new = f"{k} = {v}"
        if hit:
            body[hit[0]] = new
        else:
            appended.append(new)
    hdr = f"""#@
#@ [3d]
#@ dim=3
#@ check=false
#@ exe=mmwspalling
#@

# ============================================================================
# Meier 2017 Grimsel pilot -- D2c scored configuration (for D2d run 1).
# ============================================================================
# Generated by tests/MMWSpalling/studies/d2c_steady/run.py --write-input from
# input_drilling (unchanged; input_feet, the D2b configuration, is untouched).
# Quarter domain {meta['nxy']} x {meta['nxy']} x {meta['nz']} at 2 mm ({LXY} x {LXY} x {meta['lz']} m, axis at
# the xlo/ylo corner), beam and bit off, V0 = V_cell, pinned closure, wall A2.
# Burner on feet (annulus [28, 40] mm, 50 mm, 3 pads at 0.9) with body
# clearance; wall-jet enthalpy closure: Martin h_ref {meta['h_ref']:g} at
# T_nozzle {meta['T_nozzle']:g} K with the far law (12 D/s)^1 beyond 12 D,
# momentum decay diameter jet_De_ref {meta['De_ref']:.6g} m (full jet; nozzle.py)
# with core 8, entrained mass, exhaust recirculation, free-surface dilution
# (aspect 1, exponent 1), exponential bin update, rock-to-gas heat counted.
# stop {meta['stop']:g} s (cap), dt {meta['dt'] * 1e3:g} ms. The D2d harness adds the watchdogs.
# No parameter was adjusted toward Meier's ROP, hole diameter or volume rate.
"""
    with open(dst, "w") as f:
        f.write(hdr + "\n" + "\n".join(body).rstrip("\n") + "\n\n# --- D2b/D2c keys not in input_drilling ---\n"
                + "\n".join(appended) + "\n")
    print(f"wrote {dst}: {len(ks) - len(appended)} keys replaced, {len(appended)} appended")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--probe", choices=("2mm", "1mm"))
    ap.add_argument("--project", action="store_true")
    ap.add_argument("--stage", choices=("D2d",))
    ap.add_argument("--cases", nargs="*", default=None)
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--write-input", action="store_true")
    args = ap.parse_args()
    if args.write_input:
        write_input()
        return
    if args.project:
        L, _, _ = project()
        print("\n".join(L))
        return
    if not os.path.exists(BIN):
        sys.exit(f"missing binary {BIN}; build first")
    if args.probe == "2mm":
        jobs = [job("probe_2mm", stop=30.0)]
    elif args.probe == "1mm":
        jobs = [job("probe_1mm", dz=1e-3, stop=10.0)]
    elif args.stage == "D2d":
        names = args.cases or D2D_DEFAULT
        jobs = [job(n, **MATRIX[n]) for n in names]
    else:
        sys.exit("give --probe, --project, --stage D2d or --write-input")
    print(f"{len(jobs)} runs, {min(args.jobs, len(jobs))} in parallel -> {OUT}", flush=True)
    for pf, _, m in jobs:
        print(f"  {m['name']}: T {m['T_nozzle']:g} K, h_ref {m['h_ref']:g}, far {m['far']} n {m['n']:g}, "
              f"{'momentum D_e ' + format(m['De_ref'], '.4g') if m['momentum'] else 'nozzle D'}, core {m['core']:g}, "
              f"fs {int(m['fs'])}, dz {m['dz'] * 1e3:g} mm, {m['nxy']}²×{m['nz']}, dt {m['dt'] * 1e3:g} ms, "
              f"stop {m['stop']:g} s, q_max {m['q_max'] / 1e6:.2f} MW/m², overshoot {m['overshoot']:.2f} K/step",
              flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(args.jobs, len(jobs))) as ex:
        for fu in cf.as_completed([ex.submit(run_one, pf, c, m, args.force) for pf, c, m in jobs]):
            m, status, wall_s = fu.result()
            extra = f", {rate_of(m):.3g} cell-steps/s" if m.get("wall_s") and m.get("t_end") else ""
            print(f"  {m['name']:14s} {status} at t = {m.get('t_end', 0):g} s ({wall_s / 60:.1f} min{extra})",
                  flush=True)
            ok = ok and not status.startswith("FAILED")
    print("RUNS OK" if ok else "RUNS HAD FAILURES")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
