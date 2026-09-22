#!/usr/bin/env python3
"""D2f harness: the contact support rule (pads, no clip) at q = 1.0 / 0.9 / 0.8.

  run.py --cases Q10_2mm Q09_2mm Q08_2mm [--jobs 3] [--force]
  run.py --identity        Q09_2mm vs studies/d2e_mesh/output/Pa_2mm, rows
                           t <= 200 s (Pa's end) byte-identical
  run.py --kill NAME       binary path + plot_file match (never pkill -f)
  run.py --list

Loads studies/d2e_mesh/run.py by path (its SWITCH / STOP / job / dt policy /
mgs 60 / kill, unchanged) and registers three switches, each replacing the
scored keys foot_pad_quantile and foot_body_clearance in d2c's keys() list:
  Q10: foot_pad_quantile = 1.0 (the strict per-pad maximum), clearance 0
  Q09: 0.9, clearance 0 (= D2e P-a)
  Q08: 0.8, clearance 0
Stop 300 s for all; stall (60 s) and bottom watchdogs on; no steady-stop;
dt 16 ms at 2 mm and 8 ms at 1 mm; 1 mm plotfiles only at 0 and the stop.
Output: studies/d2f_support/output/.
"""
import argparse
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D2E = os.path.join(os.path.dirname(HERE), "d2e_mesh")
_spec = importlib.util.spec_from_file_location("d2e_run", os.path.join(D2E, "run.py"))
E = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E)
R = E.R

OUT = os.path.join(HERE, "output")
E.OUT = OUT                      # job() and kill() read E.OUT
R.OUT = OUT                      # run_one's makedirs
QS = {"Q10": 1.0, "Q09": 0.9, "Q08": 0.8}
for name, q in QS.items():
    E.SWITCH[name] = ({"surface_patch.foot_pad_quantile=0.9": [f"surface_patch.foot_pad_quantile={q!r}"],
                       "surface_patch.foot_body_clearance=1": ["surface_patch.foot_body_clearance=0"]},
                      f"foot_pad_quantile = {q}, foot_body_clearance = 0")
    E.STOP[name] = 300.0
CASES = [f"{p}_{m}" for p in QS for m in ("2mm", "1mm")]
REF_PA = os.path.join(D2E, "output", "Pa_2mm")


def identity():
    pf = os.path.join(OUT, "Q09_2mm")
    t_end = float(R.read_thermo(REF_PA, ["time"])[-1, 0])
    res = E.compare(pf, REF_PA, t_end)
    L = [f"| file | rows (t <= {t_end:g} s, Pa_2mm's end) | Q09_2mm vs d2e Pa_2mm |", "|---|---|---|"]
    for f, same, n, d in res:
        L.append(f"| {f} | {n} | {'byte-identical' if same else f'DIFFER (max rel {d:.2e})'} |")
    ok = all(s for _, s, _, _ in res)
    L.append(f"\nIdentity: {'PASS' if ok else 'FAIL'}")
    print("\n".join(L))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--identity", action="store_true")
    ap.add_argument("--kill")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.kill:
        return E.kill(a.kill)
    if a.list:
        for n in CASES:
            pf, cmd, m = E.job(n)
            print(f"{n:8s} {m['switch']:48s} dt {m['dt'] * 1e3:g} ms {m['K_per_step']:.2f} K/step "
                  f"stop {m['stop']:g} s {[c for c in cmd if c.startswith('amr.plot_int')][0]} "
                  f"mgs {m['max_grid_size']}")
        return
    if a.identity:
        sys.exit(0 if identity() else 1)
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if E.run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
