#!/usr/bin/env python3
"""P2-31 c1 diagnosis: where does the overlap rule make an OBSERVABLE difference? `Upd::h_line` was read at two places
only -- a +z line's own cell (the column's phi shift) and a voided cell (the plotted thickness / rate) -- so an overlap
in an interior cell crossed by two lines changes nothing. The event log is written in the lines' processing order, so
the OLD rule (last writer) can be replayed offline against the new one (the +z line's own h; the largest h).
  ovl_search.py replay <events.csv> <nx> <ny> <nz> <dx>    the steps where old != new
  ovl_search.py run                                         the fixture's tilts, SPH and BORE at 1 mm to 20 s (event
                                                            logs only) on bin/, into S/legs/ovl_search/, then replay"""
import collections
import csv
import importlib.machinery
import importlib.util
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
H_SPALL_MAX = 0.015


def replay(path, nx, ny, nz, dx):
    span = max(nx, ny, nz)
    steps = collections.OrderedDict()
    for e in csv.DictReader(open(path)):
        if e["regime"] == "1":
            steps.setdefault(e["time"], []).append(e)
    out = []
    for t, rows in steps.items():
        last, hmax, nl, own, vox = {}, {}, collections.Counter(), {}, set()
        for e in rows:
            i, j, k, d = int(e["col_i"]), int(e["col_j"]), int(e["k"]), int(e["dir"])
            run_, nv, mu, h = int(e["run"]), int(e["n_voided"]), float(e["mu"]), float(e["h_line"])
            ax, sg = d // 2, (-1 if d % 2 else 1)
            Dmax = max(1, min(span, int(math.floor(H_SPALL_MAX / (mu * dx)))))
            Wr = min(min(Dmax, run_ - 1), run_ - 1 - nv)
            for jd in range(0, min(run_ - 1, Wr + nv) + 1):
                p = [i, j, k]
                p[ax] += sg * jd
                key = tuple(p)
                last[key] = h
                hmax[key] = max(hmax.get(key, 0.0), h)
                nl[key] += 1
                if jd < nv:
                    vox.add(key)
            if d == 5:
                last[(i, j, k)] = h
                own[(i, j, k)] = h
        a = [k for k, h in own.items() if last[k] != h]
        b = [k for k in vox if last[k] != hmax[k]]
        if a or b:
            out.append((t, len(rows), sum(1 for c in nl.values() if c > 1), len(a), len(b)))
    return out


def load(path, name):
    ld = importlib.machinery.SourceFileLoader(name, path)
    sp = importlib.util.spec_from_loader(name, ld)
    m = importlib.util.module_from_spec(sp)
    ld.exec_module(m)
    return m


def main_run():
    fr = load(os.path.join(ROOT, "tests", "MMWSpalling", "unit", "face_removal", "test"), "fr_ovl_search")
    fr.OUT = os.path.join(HERE, "legs", "ovl_search")
    os.makedirs(fr.OUT, exist_ok=True)
    k = {"surface_patch.face_list": "1", "stop_time": "20.0", "amr.plot_int": "-1"}
    legs = []
    for nm, (th, ph) in fr.TILTS.items():
        lg = nm + "_1"
        fr.run_leg(lg, dict(fr.mesh_keys("1"), **k), seeds=fr.tilt_seeds(th, ph, "1"))
        legs.append((lg, (24, 24, 64), 0.024 / 24))
    for c in ("SPH", "BORE"):
        lg = c + "_1"
        fr.run_leg(lg, dict(fr.MI.curved_keys(c, 0.001), **k), seeds=fr.MI.curved_seeds(c, 0.001))
        L = fr.MI.CURVED[c]["L"]
        n = int(round(L / 0.001))
        legs.append((lg, (n, n, n), L / n))
    for lg, (nx, ny, nz), dx in legs:
        p = os.path.join(fr.OUT, lg, "run_removal_events.csv")
        f = replay(p, nx, ny, nz, dx) if os.path.exists(p) else None
        print("%-8s %s" % (lg, "no events" if f is None else
                           "%d step(s) with old != new; first %s" % (len(f), f[:3])), flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "replay":
        for f in replay(sys.argv[2], *[int(x) for x in sys.argv[3:6]], float(sys.argv[6])):
            print(f)
    else:
        main_run()
