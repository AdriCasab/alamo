#!/usr/bin/env python3
"""D2g Stage 1 harness: the D2f support rule (Q09: pads, q 0.9, no clip) run
long, to the bottom watchdog of the 0.40 m domain.

  run.py --cases L09_2mm L09_1mm [--jobs 2] [--force]
  run.py --identity        each L09 thermo.dat, rows t <= 299.9 s, vs the
                           matching d2f Q09 run (byte-identical)
  run.py --kill NAME       binary path + plot_file match (never pkill -f)
  run.py --list

Loads studies/d2f_support/run.py by path (which loads d2e's; neither is
edited) and adds switch L09 = Q09's key replacements with a 700 s stop cap.
The bottom watchdog (centre <= 40 mm above the domain bottom) ends the runs
(~550 s); stall (60 s) is on; no steady-stop. dt 16 ms at 2 mm, 8 ms at
1 mm; mgs 60. Plotfiles: 2 mm every 100 s; 1 mm only at t = 0 (plot_int =
stop/dt, and a watchdog stop writes no final plotfile). Output:
studies/d2g_mesh/output/.
"""
import argparse
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D2F = os.path.join(os.path.dirname(HERE), "d2f_support")
_spec = importlib.util.spec_from_file_location("d2f_run", os.path.join(D2F, "run.py"))
F = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(F)
E, R = F.E, F.R

OUT = os.path.join(HERE, "output")
E.OUT = OUT
R.OUT = OUT
E.SWITCH["L09"] = (E.SWITCH["Q09"][0], E.SWITCH["Q09"][1] + " (long: stop cap 700 s, bottom watchdog)")
E.STOP["L09"] = 700.0
CASES = ["L09_2mm", "L09_1mm"]
REF = {"L09_2mm": os.path.join(D2F, "output", "Q09_2mm"), "L09_1mm": os.path.join(D2F, "output", "Q09_1mm")}
T_ID = 299.9


def identity():
    ok = True
    L = ["| run | rows (t <= 299.9 s) | vs d2f Q09 |", "|---|---|---|"]
    for n, ref in REF.items():
        a = E.rows_upto(os.path.join(OUT, n, "thermo.dat"), T_ID)
        b = E.rows_upto(os.path.join(ref, "thermo.dat"), T_ID)
        same = a == b
        ok = ok and same
        L.append(f"| {n} | {len(a) - 1} | {'byte-identical' if same else 'DIFFER'} |")
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
            print(f"{n:8s} {m['switch']}; dt {m['dt'] * 1e3:g} ms {m['K_per_step']:.2f} K/step "
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
