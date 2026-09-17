#!/usr/bin/env python3
"""D2b study: burner on feet (surface_patch.nozzle_descent = feet) with the
D2a2 wall-jet closure and exhaust recirculation. 3-D quarter domain only.

  run.py --predict            hand model (feetmodel.py) + Stage A rule -> RESULTS.md section 0
  run.py --stage A            wall-treatment decision runs (J-M, 150 s)
  run.py --probe              cost probe (10 s of J-M, 2 mm)
  run.py --stage B --wall A2  scored runs J-M, J-5, J-10 (600 s / watchdog / bottom)
  run.py --probe-1mm --wall A2
  run.py --stage D --wall A2 [--stop-d 300]   mesh check (1 mm + 2 mm, trimmed r <= 60 mm)
  run.py --stage C --wall A2  sensitivities (J-M: fixed 293 K; fixed + nozzle mass; 1436 K)
  run.py --write-input-feet --wall A2   validation/meier/sp_meier_pilot/input_feet

Common: input_drilling + CLI overrides; quarter domain, axis at the xlo/ylo
corner; 2 mm, 0.12 x 0.12 x 0.30 m; beam off, bit off; V0 = V_cell; pinned
closure (idle 2); follow_mask; ledger; removal_events.csv; jet_closure =
enthalpy, decay, entrained mass, jet_T_ent_mode = exhaust, T_ent 293.15,
mdot/4, cp 1250, D 7.5 mm, T_nozzle 1900; feet on [0.028, 0.040] m,
50 mm stand-off, 3 pads, 0.9 quantile.

Watchdog (polls thermo.dat): stop cleanly and mark
  "stalled" if nozzle_z has not descended for 60 s of simulated time;
  "bottom"  if the centre face (nozzle_z - jet_s_c) is <= 40 mm above the
            domain bottom.
Time step: q dt/(rho Cp dz) <= 5 K with q = h(0, s >= SOD) (T_stag,max - 821),
T_stag,max = T_rec + 0.75 (T_nozzle - T_rec) with the hand model's T_rec
(fixed mode: T_rec = 293.15); dt = 4 ms halved until it holds.
No parameter is adjusted toward Meier's ROP, hole diameter or volume rate.
"""
import argparse
import concurrent.futures as cf
import math
import os
import shutil
import signal
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import feetmodel as fm  # noqa: E402

wj = fm.wj
ROOT = HERE.split(os.sep + os.path.join("tests", "MMWSpalling") + os.sep)[0]
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
MEIER = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier", "sp_meier_pilot")
OUT = os.path.join(HERE, "output")
D = wj.D
DT_CAP = 4.0e-3
T_REC_HAND = {("M", 1900.0): 1679.0, ("5", 1900.0): 1397.0, ("10", 1900.0): 1212.0,
              ("M", 1436.0): 1314.0}
WALLS = {
    "A1": ["spall.flake_coherence_length=0.004", "spall.surface_normal=0"],
    "A2": ["spall.flake_coherence_length=0.0", "spall.surface_normal=0"],
    "A3": ["spall.flake_coherence_length=0.004", "spall.surface_normal=1"],
}
STALL_S, BOTTOM_M = 60.0, 0.040


def dt_rule(q, dz, cap=DT_CAP):
    dt = cap
    while q * dt / (wj.RHOCP * dz) > 5.0:
        dt *= 0.5
    return dt


def q_max(hk, Tn, mode):
    T_rec = T_REC_HAND.get((hk, Tn), 293.15) if mode == "exhaust" else 293.15
    T_st = T_rec + 0.75 * (Tn - T_rec)
    h = float(wj.h_anchor(hk, Tn))
    hs = max(float(wj.h_py(h, 0.0, s)) for s in (0.050, 0.06, 0.07, 0.08, 0.09))
    return hs * (T_st - wj.T_FIRE) - float(wj.loss(wj.T_FIRE))


def keys(hk, Tn, mode="exhaust", entrained=1, dz=2e-3, nxy=60, nz=150, lxy=0.12, lz=0.30,
         patch_r=0.2, stop=600.0, wall="A2", plot_s=100.0):
    h = float(wj.h_anchor(hk, Tn))
    q = q_max(hk, Tn, mode)
    dt = dt_rule(q, dz)
    k = [f"stop_time={stop!r}", f"timestep={dt!r}",
         f"amr.plot_int={int(round(plot_s / dt))}",
         f"amr.thermo.plot_int={max(1, int(round(0.1 / dt)))}", "amr.thermo.int=1",
         f"amr.n_cell={nxy} {nxy} {nz}", f"geometry.prob_hi={lxy!r} {lxy!r} {lz!r}",
         "bit.enabled=0", "beam.P0=0.0", f"weibull.V0={dz ** 3!r}",
         "surface_patch.enabled=1", "surface_patch.mode=convective_flame",
         "surface_patch.robin_form=pinned", "surface_patch.pinned_idle_cycles=2.0",
         "surface_patch.face=zhi", "surface_patch.x0=0.0", "surface_patch.y0=0.0",
         f"surface_patch.radius={patch_r!r}",
         f'surface_patch.h_expr="{wj.h_expr(h)}"',
         "surface_patch.nozzle_descent=feet", "surface_patch.foot_r_inner=0.028",
         "surface_patch.foot_r_outer=0.040", "surface_patch.foot_standoff=0.050",
         "surface_patch.foot_rule=pads", "surface_patch.foot_npads=3",
         "surface_patch.foot_pad_quantile=0.9",
         "surface_patch.jet_closure=enthalpy", "surface_patch.jet_stagnation=decay",
         f"surface_patch.jet_entrained_mass={entrained}",
         f"surface_patch.jet_T_ent_mode={mode}", "surface_patch.jet_T_ent=293.15",
         f"surface_patch.jet_mdot={wj.MDOT / 4.0!r}", f"surface_patch.jet_cp={wj.CP!r}",
         f"surface_patch.jet_D={D!r}", f"surface_patch.jet_T_nozzle={Tn!r}",
         "surface_patch.jet_profile_interval=10.0",
         "surface.follow_mask=1", "energy_ledger.enabled=1", "spall.removal_events_csv=1",
         ] + WALLS[wall]
    meta = dict(anchor=hk, h_ref=h, T_nozzle=Tn, mode=mode, entrained=entrained, dz=dz,
                nxy=nxy, nz=nz, lxy=lxy, lz=lz, patch_r=patch_r, stop=stop, wall=wall, dt=dt,
                q_max=q, overshoot=q * dt / (wj.RHOCP * dz), sector=0.25)
    return k, meta


def job(name, **kw):
    k, meta = keys(**kw)
    pf = os.path.join(OUT, name)
    cmd = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN,
           os.path.join(MEIER, "input_drilling"), f"plot_file={pf}"] + k
    meta["name"] = name
    return pf, cmd, meta


def stage_jobs(stage, wall, stop_d=None):
    if stage == "A":
        return [job(f"A_{w}", hk="M", Tn=1900.0, stop=150.0, wall=w) for w in ("A1", "A2", "A3")]
    if stage == "B":
        return [job(f"B_J{hk}_{wall}", hk=hk, Tn=1900.0, stop=600.0, wall=wall) for hk in ("M", "5", "10")]
    if stage == "C":
        return [job(f"C_fixed_{wall}", hk="M", Tn=1900.0, mode="fixed", stop=600.0, wall=wall),
                job(f"C_fixed_nozmass_{wall}", hk="M", Tn=1900.0, mode="fixed", entrained=0,
                    stop=600.0, wall=wall),
                job(f"C_1436_{wall}", hk="M", Tn=1436.0, stop=600.0, wall=wall)]
    if stage == "D":
        s = stop_d or 300.0
        lz = 0.30
        return [job(f"D_1mm_{wall}", hk="M", Tn=1900.0, dz=1e-3, nxy=60, nz=int(round(lz / 1e-3)),
                    lxy=0.06, lz=lz, patch_r=0.06, stop=s, wall=wall, plot_s=100.0),
                job(f"D_2mm_{wall}", hk="M", Tn=1900.0, dz=2e-3, nxy=30, nz=int(round(lz / 2e-3)),
                    lxy=0.06, lz=lz, patch_r=0.06, stop=s, wall=wall, plot_s=100.0)]
    raise SystemExit(f"unknown stage {stage}")


def last_row(pf):
    try:
        with open(os.path.join(pf, "thermo.dat")) as f:
            head = f.readline().split()
            last = None
            for line in f:
                if line.strip():
                    last = line
        return head, (last.split() if last else None)
    except Exception:
        return None, None


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
    t0 = time.time()
    status = "ok"
    with open(pf + ".log", "w") as log:
        log.write(" ".join(cmd) + "\n")
        log.flush()
        proc = subprocess.Popen(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
        best_z, best_t = math.inf, 0.0
        while proc.poll() is None:
            time.sleep(10.0)
            head, row = last_row(pf)
            if not row or head is None or "nozzle_z" not in head:
                continue
            try:
                t = float(row[0])
                z = float(row[head.index("nozzle_z")])
                sc = float(row[head.index("jet_s_c")])
            except (ValueError, IndexError):
                continue
            if t <= 0.0:
                continue
            if z < best_z - 1e-9:
                best_z, best_t = z, t
            reason = None
            if t - best_t >= STALL_S:
                reason = "stalled"
            elif sc > 0.0 and z - sc <= BOTTOM_M:
                reason = "bottom"
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
    wall = time.time() - t0
    head, row = last_row(pf)
    meta = dict(meta, wall=wall, status=status, t_end=float(row[0]) if row else 0.0)
    if status != "ok" or rc == 0:
        with open(done, "w") as f:
            f.write(repr(meta) + "\n")
        return meta, status, wall
    return meta, f"FAILED rc={rc}", wall


STAGE_A_RULE = """### Stage A decision rule (written before Stage A runs)

Runs: J-M, 150 s, with A1 `flake_coherence_length = 0.004`, A2
`flake_coherence_length = 0`, A3 `spall.surface_normal = 1`.
Scoring window: 50–150 s.

1. **Energy ratio** R = ρCp·ΔT_fire·(removed volume rate) / (absorbed face power).
   - The removed volume rate is Σ h_applied·dx·dy over the window.
   - ΔT_fire is the window-mean regime-1 T_top − 293.15 K.
   - The absorbed face power is the window mean of the jet profile's P_bin
     sum (= jet_P_face when there is no floor supply).
   - Whole hole: R ≤ 1.05. R per 1 D ring is reported; a rim ring may
     exceed 1.
2. **Overheating:** regime-1 T_top p99 ≤ 873 K, i.e. the pinned T_fire band
   (C1 / S1b: 821–823 K at V0 = V_cell) + 50 K.
3. **Artefacts, none allowed.**
   - Needles: a column whose end depth exceeds every 4-neighbour's by more
     than 2 dz (4 mm).
   - Never-fired columns inside r_wall(2 mm).
4. **Choice among treatments passing 1–3:** A2 > A3 > A1.
5. **Tie-break, if none passes 1–3.** Accept the treatments whose violations
   do not touch the foot annulus:
   - ring-bin R ≤ 1.05;
   - no ring column with T_top > 873 K;
   - no never-fired ring column;
   - no needle in the ring.

   Choose among them in the order A2 > A3 > A1, and carry the violation as a
   caveat.
6. **If every treatment violates inside the ring, stop and report.**
7. **ROP proximity to Meier is not looked at.**
"""


def predict():
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        fm.main()
    text = buf.getvalue()
    path = os.path.join(HERE, "RESULTS.md")
    body = f"""# D2b results: burner on feet, scored Meier ROP / hole diameter

## 0. Pre-registration (written 2026-09-16, before any D2b run)

### Hand prediction (deliverable 5; `feetmodel.py`, output verbatim)

The model is described in the `feetmodel.py` docstring. In short:
- the ring is the 0.9 area-quantile radius of the annulus (38.97 mm) at
  s = 50 mm;
- the centre stand-off s_c is set by centre ROP = v;
- s is linear from s_c (r ≤ r_core) to 50 mm at the ring;
- the face is at T_fire = 821 K;
- the march ends at the ring in steady state (columns beyond have frozen),
  and T_rec = T_exhaust is iterated to its fixed point.

```
{text.rstrip()}
```

**Reading, before running:**
- **Exhaust recirculation has no centre equilibrium.** The entrained gas is
  the jet's own exhaust, so T_stag → T_rec as s_c grows, and m = mdot/φ grows
  without bound. The centre is faster than the ring at every stand-off, and
  the pit is bounded only by the wall treatment (Stage A) or the domain
  bottom.
- **In the deep-pit limit the gas over the face is uniform at T_rec,** set by
  the global balance P_face = mdot·cp·(T_nozzle − T_rec). The ring ROP then
  follows from h at the ring:

  | case | ring ROP | T_rec | P_face |
  |---|---|---|---|
  | **J-M, 1900 K (scored)** | **1.82 m/h** | 1679 K | 4.1 kW |
  | J-5, 1900 K | 4.3 m/h | 1397 K | 9.5 kW |
  | J-10, 1900 K | 5.9 m/h | 1212 K | 13.0 kW |
  | J-M, 1436 K | 0.97 m/h | 1314 K | 2.3 kW |

- **Fixed 293 K entrainment (either mass treatment): stall predicted.** The
  centre drills ahead until its T_stag reaches T_fire (s_c = 114 mm), and the
  ring ROP falls to zero on the way (0.07 m/h at s_c = 100 mm). Only a wall
  treatment that holds the centre near the ring (s_c ≲ 60 mm) would keep it
  drilling.
- **Hole.**
  - The wall is at the column's own freeze depth
    D(r) = 50 mm · v(r)/(v − v(r)). The packet's d(r) = 50 mm · v/(v − v(r))
    is the feet's depth at that moment, exactly 50 mm deeper.
  - Depth-mean Ø over 0–0.5 m: 110–113 mm for all four exhaust cases.
  - Minimum Ø at 0.5 m: 85 mm (still narrowing toward 2·r_q = 78 mm).
  - The shallow wall reaches the domain edge: v(r) > 0 out to 200 mm in the
    hand model.
- **D2a2 cross-check.** With D2a2's measured s_c the model gives ring ROPs of
  0.39–0.95 m/h against D2a2's measured 0.54–1.27 m/h: J-M +20 %, J-5 −16 %
  (nozzle mass) and +21 % (entrained), J-10 −45 % and −25 %. The model uses a
  firing face everywhere, while D2a2's face was mostly unfired, so this
  agreement is order-of-magnitude, not closure-level.

{STAGE_A_RULE}
<!-- RESULTS -->

<!-- DISCUSSION -->
"""
    if os.path.exists(path):
        raise SystemExit(f"{path} exists; refusing to overwrite the pre-registration")
    with open(path, "w") as f:
        f.write(body)
    print(body)


def write_input_feet(wall, stop):
    """validation/meier/sp_meier_pilot/input_feet = input_drilling with the
    Stage B J-M keys written in (existing keys replaced in place, new keys
    appended). input_drilling itself is not touched."""
    import re
    src = os.path.join(MEIER, "input_drilling")
    dst = os.path.join(MEIER, "input_feet")
    ks, meta = keys(hk="M", Tn=1900.0, stop=stop, wall=wall)
    ks = [f"plot_file=tests/MMWSpalling/validation/meier/sp_meier_pilot/output/feet"] + ks
    lines = open(src).read().split("\n")
    body_start = next(i for i, l in enumerate(lines) if l.startswith("alamo.program"))
    body = lines[body_start:]
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
# Meier 2017 Grimsel pilot -- D2b burner-on-feet scored configuration (J-M).
# ============================================================================
# Generated by tests/MMWSpalling/studies/d2b_feet_rop/run.py --write-input-feet
# from input_drilling (unchanged) with the D2b Stage B J-M keys written in:
# quarter domain 60 x 60 x 150 at 2 mm (axis at the xlo/ylo corner), beam and
# bit off, V0 = V_cell, pinned closure, wall-jet enthalpy closure with the
# stagnation decay, entrained mass and exhaust recirculation (Martin h at the
# measured mdot, T_nozzle = 1900 K), burner on feet (annulus [28, 40] mm,
# 50 mm stand-off, 3 pads at the 0.9 quantile), wall treatment {wall}
# (the D2b Stage A choice), stop {stop:g} s (dt {meta['dt'] * 1e3:g} ms).
# Scored by test_feet (check-only over the study output; see test_feet).
# No parameter was adjusted toward Meier's ROP, hole diameter or volume rate.
"""
    with open(dst, "w") as f:
        f.write(hdr + "\n" + "\n".join(body).rstrip("\n") + "\n\n# --- D2b keys not in input_drilling ---\n"
                + "\n".join(appended) + "\n")
    print(f"wrote {dst}: {len(ks) - len(appended)} keys replaced, {len(appended)} appended")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predict", action="store_true")
    ap.add_argument("--stage", choices=("A", "B", "C", "D"))
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--probe-1mm", action="store_true")
    ap.add_argument("--wall", choices=tuple(WALLS), default=None)
    ap.add_argument("--stop-d", type=float, default=None)
    ap.add_argument("--cases", nargs="*", default=None)
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--write-input-feet", action="store_true")
    ap.add_argument("--stop-feet", type=float, default=600.0)
    args = ap.parse_args()
    if args.write_input_feet:
        if not args.wall:
            sys.exit("--write-input-feet needs --wall")
        write_input_feet(args.wall, args.stop_feet)
        return
    if args.predict:
        predict()
        return
    if not os.path.exists(BIN):
        sys.exit(f"missing binary {BIN}; build first")
    if args.probe:
        jobs = [job("probe", hk="M", Tn=1900.0, stop=10.0, wall=args.wall or "A1")]
    elif args.probe_1mm:
        if not args.wall:
            sys.exit("--probe-1mm needs --wall")
        pf, cmd, m = stage_jobs("D", args.wall)[0]
        jobs = [job("probe_1mm", hk="M", Tn=1900.0, dz=1e-3, nxy=60, nz=300, lxy=0.06, lz=0.30,
                    patch_r=0.06, stop=10.0, wall=args.wall)]
    elif args.stage:
        if args.stage != "A" and not args.wall:
            sys.exit("stages B, C, D need --wall (the Stage A choice)")
        jobs = stage_jobs(args.stage, args.wall, args.stop_d)
        if args.cases:
            jobs = [j for j in jobs if j[2]["name"] in args.cases]
    else:
        sys.exit("give --predict, --stage, --probe or --probe-1mm")
    print(f"{len(jobs)} runs, {min(args.jobs, len(jobs))} in parallel -> {OUT}", flush=True)
    for pf, _, m in jobs:
        print(f"  {m['name']}: h_ref {m['h_ref']:g}, T_nozzle {m['T_nozzle']:g}, {m['mode']}"
              f"{'' if m['entrained'] else ' nozzle-mass'}, wall {m['wall']}, dz {m['dz'] * 1e3:g} mm, "
              f"dt {m['dt'] * 1e3:g} ms, stop {m['stop']} s, q_max {m['q_max'] / 1e6:.2f} MW/m^2, "
              f"overshoot {m['overshoot']:.2f} K/step", flush=True)
    ok = True
    with cf.ThreadPoolExecutor(max_workers=min(args.jobs, len(jobs))) as ex:
        for fu in cf.as_completed([ex.submit(run_one, pf, c, m, args.force) for pf, c, m in jobs]):
            m, status, wall = fu.result()
            print(f"  {m['name']:22s} {status} at t = {m.get('t_end', 0):g} s  ({wall / 60:.1f} min)",
                  flush=True)
            if (args.probe or args.probe_1mm) and not status.startswith("FAILED"):
                cells = m["nxy"] * m["nxy"] * m["nz"]
                steps = m["t_end"] / m["dt"]
                rate = cells * steps / wall
                print(f"    {rate:.3g} cell-steps/s", flush=True)
                for st in ("B", "C", "D") if args.probe else ("D",):
                    for _, _, mm in stage_jobs(st, args.wall or "A1"):
                        est = mm["nxy"] ** 2 * mm["nz"] * mm["stop"] / mm["dt"] / rate
                        print(f"    estimate {mm['name']}: {est / 3600:.2f} h alone (to its full stop)",
                              flush=True)
            ok = ok and not status.startswith("FAILED")
    print("RUNS OK" if ok else "RUNS HAD FAILURES")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
