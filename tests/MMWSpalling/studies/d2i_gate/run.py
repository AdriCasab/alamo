#!/usr/bin/env python3
"""D2i harness: the D2f support rule with the idle clock off (D2h's I0), run
long, to the bottom watchdog of the 0.40 m domain — the decisive mesh gate.

  run.py --cases LI0_2mm LI0_1mm I20_1mm [--jobs 3] [--force]
  run.py --identity        each LI0 thermo.dat, rows t <= 349.9 s, vs the
                           matching d2h I0 run (byte-identical)
  run.py --kill NAME       binary path + plot_file match (never pkill -f)
  run.py --list

Loads studies/d2h_openloop/run.py by path (which loads d2f -> d2e -> d2c;
none of them is edited) and registers two switches on top of D2f's Q09 keys
(pads, q 0.9, no clip):

  LI0  = d2h I0 (surface_patch.pinned_idle_cycles 2.0 -> 1.0e9, idle clock
         off), stop cap 700 s; the bottom watchdog (centre <= 40 mm above the
         domain floor) ends the runs (~550-620 s). Stall (60 s) on, no
         steady-stop. dt 16 ms at 2 mm, 8 ms at 1 mm; mgs 60.
  I20  = the same keys with pinned_idle_cycles = 20, at 1 mm, stop 350 s. The
         third point of the idle sensitivity trio (2 -> d2g L09_1mm, 1e9 ->
         d2h I0_1mm, both existing data). At 1 mm because that is where the
         silent front appeared; at 2 mm the clock is already immaterial.
         Its key differs from I0's, so it has no identity partner.

`pinned_idle_cycles` is set per input only. The source default (2) is
untouched: changing it would rebase every pinned-closure result since C1.

Plotfiles every 50 s at both meshes with the full field list, as D2h wrote
them (a 1 mm plotfile is 1.58 GB: ~20 GB for the long 1 mm run, ~12 GB for
I20_1mm, against ~475 GB free). Output: studies/d2i_gate/output/.
"""
import argparse
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D2H = os.path.join(os.path.dirname(HERE), "d2h_openloop")
_spec = importlib.util.spec_from_file_location("d2h_run", os.path.join(D2H, "run.py"))
H = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(H)
E, R = H.E, H.R

OUT = os.path.join(HERE, "output")
E.OUT = OUT
R.OUT = OUT

IDLE_KEY = "surface_patch.pinned_idle_cycles=2.0"
I0 = dict(E.SWITCH["I0"][0])                       # Q09 keys + idle clock off
E.SWITCH["LI0"] = (I0, "Q09 + pinned_idle_cycles = 1e9 (idle off), long")
E.STOP["LI0"] = 700.0

I20 = dict(E.SWITCH["Q09"][0])
I20[IDLE_KEY] = ["surface_patch.pinned_idle_cycles=20.0"]
E.SWITCH["I20"] = (I20, "Q09 + pinned_idle_cycles = 20")
E.STOP["I20"] = 350.0

CASES = ["LI0_2mm", "LI0_1mm", "I20_1mm"]
REF = {"LI0_2mm": os.path.join(D2H, "output", "I0_2mm"),
       "LI0_1mm": os.path.join(D2H, "output", "I0_1mm")}
T_ID = 349.9

job = H.job                                        # d2h's wrapper: plotfiles every 50 s


def identity():
    ok = True
    L = [f"| run | rows (t <= {T_ID:g} s) | vs d2h I0 |", "|---|---|---|"]
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
            pf, cmd, m = job(n)
            print(f"{n:8s} {m['switch']}; dt {m['dt'] * 1e3:g} ms {m['K_per_step']:.2f} K/step "
                  f"stop {m['stop']:g} s {[c for c in cmd if c.startswith('amr.plot_int')][0]} "
                  f"mgs {m['max_grid_size']}")
            print("   idle:", " ".join(c for c in cmd if "idle" in c))
        return
    if a.identity:
        sys.exit(0 if identity() else 1)
    bad = [n for n in (a.cases or []) if n not in CASES]
    if not a.cases or bad:
        sys.exit(f"give --cases from {CASES}" + (f" (unknown: {bad})" if bad else ""))
    sys.exit(0 if E.run_all(a.cases, a.jobs, a.force) else 1)


if __name__ == "__main__":
    main()
