#!/usr/bin/env python3
"""P2-29 stage 1 -- THE HU DECKS on the production kernel (bin_1, S/legs/hu_new/) against the OLD path's 09-25 outputs.

The user's trim (2026-11-10): no old side is re-run where an output of the unchanged code is on disk -- the baselines are
the 09-25 outputs of hu_conduction/output/ and hu_thermoelastic/output/, read from their archived copies
(S/archive/<test>_output, SR 9); hu_thermoelastic_amr is RETIRED (not ported, not run); hu_breakage_index has no old
output and is scored against Hu's data only. The three old legs the chain completed before the trim (S/legs/hu_old/:
hu_conduction x 2 on bin_pre) are reported as a cross-check that the 09-25 record IS the unchanged code's output.

  a8  [gate] the breakage decks (hu_breakage_index x 3), the NEW runs: at every cell plotfile where
      level 1 exists, every level-0 top-row cell inside the patch (r <= surface_patch.radius about (x0, y0)) lies
      inside a level-1 grid -- SurfaceOwner(lev) is `lev == finest_level`, so an uncovered patch cell would run
      UNHELD. While level 0 is the only level it owns the patch (held), and that is reported, not failed.
  m4  hu_conduction, m5 hu_thermoelastic: each test on the 09-25 output and on the new run, side by side;
  m6  hu_breakage_index: the test on the new runs (Hu's f_b bars). Logs: S/K1_<test>_<old|new|binpre>.log.
Usage: hu_compare.py a8 | report   (a8 exits 1 on a failed gate)
"""
import glob
import os
import re
import subprocess
import sys

import yt

yt.set_log_level(50)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
PY = "/Users/tzetze20/Desktop/code/.venv/bin/python"
HU = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "hu")
LEGS = {"new": os.path.join(HERE, "legs", "hu_new"), "binpre": os.path.join(HERE, "legs", "hu_old")}
ARCH = {t: os.path.join(HERE, "archive", t + "_output") for t in ("hu_conduction", "hu_thermoelastic")}   # 09-25
AMR = [("hu_breakage_index", "granite2_37"), ("hu_breakage_index", "sandstone2_37"),
       ("hu_breakage_index", "sandstone2_89")]


def meta(run):
    d = {}
    for ln in open(os.path.join(run, "metadata")):
        if " = " in ln:
            k, v = ln.rstrip("\n").split(" = ", 1)
            d[k.strip()] = v.strip()
    return d


def coverage(run):
    m = meta(run)
    x0, y0, r = (float(m["surface_patch." + k]) for k in ("x0", "y0", "radius"))
    rows = []
    for pf in sorted(glob.glob(os.path.join(run, "*cell"))):
        ds = yt.load(pf)
        g1 = [(g.LeftEdge.to_value(), g.RightEdge.to_value()) for g in ds.index.grids if g.Level == 1]
        lo0 = ds.domain_left_edge.to_value()
        dx0 = (ds.domain_right_edge - ds.domain_left_edge).to_value() / ds.domain_dimensions
        ztop = float(ds.domain_right_edge[2]) - 0.5 * dx0[2]
        n = cov = 0
        for i in range(int(ds.domain_dimensions[0])):      # the level-0 TOP ROW, cell centres inside the patch
            for j in range(int(ds.domain_dimensions[1])):
                c = (lo0[0] + (i + 0.5) * dx0[0], lo0[1] + (j + 0.5) * dx0[1], ztop)
                if (c[0] - x0) ** 2 + (c[1] - y0) ** 2 > r * r:
                    continue
                n += 1
                cov += any(all(lo[d] < c[d] < hi[d] for d in range(3)) for lo, hi in g1)
        rows.append((float(ds.current_time), int(ds.index.max_level), n, cov))
    return rows


def a8():
    bad = []
    for side in ("new",):
        print("a8 -- level-1 coverage of the patch's top cells, %s runs" % side)
        for d, rock in AMR:
            run = os.path.join(LEGS[side], d, rock)
            if not os.path.isdir(run):
                print("    %-22s %-14s (no run)" % (d, rock))
                if side == "new":
                    bad.append("%s/%s: no run" % (d, rock))
                continue
            rows = coverage(run)
            with1 = [x for x in rows if x[1] >= 1]
            unc = [x for x in with1 if x[3] < x[2]]
            only0 = [x[0] for x in rows if x[1] == 0]
            print("    %-22s %-14s %d plotfiles: level 1 present on %d, every patch top cell covered on %d%s; level 0 alone "
                  "(owns the patch) at t = %s; %d patch top cells per level-0 plotfile" % (
                      d, rock, len(rows), len(with1), len(with1) - len(unc),
                      ("; UNCOVERED at t = %s" % [(x[0], x[2] - x[3]) for x in unc]) if unc else "",
                      only0, rows[0][2] if rows else 0))
            if side == "new" and (unc or not with1):
                bad.append("%s/%s: %s" % (d, rock, "uncovered" if unc else "level 1 never present"))
    print("RESULT: %s" % ("FAIL: " + "; ".join(bad) if bad else "PASS"))
    return 1 if bad else 0


# (test, side) -> the test's run-directory arguments
RUNS = {
    ("hu_conduction", "old"): [os.path.join(ARCH["hu_conduction"], r) for r in ("granite2", "sandstone2")],
    ("hu_conduction", "binpre"): [os.path.join(LEGS["binpre"], "hu_conduction", r) for r in ("granite2", "sandstone2")],
    ("hu_conduction", "new"): [os.path.join(LEGS["new"], "hu_conduction", r) for r in ("granite2", "sandstone2")],
    ("hu_thermoelastic", "old"): [os.path.join(ARCH["hu_thermoelastic"], r) for r in ("granite2", "sandstone2")],
    ("hu_thermoelastic", "new"): [os.path.join(LEGS["new"], "hu_thermoelastic", r) for r in ("granite2", "sandstone2")],
    ("hu_breakage_index", "new"): [os.path.join(LEGS["new"], "hu_breakage_index")],
}
ORDER = [("hu_conduction", "m4", ("old", "binpre", "new")), ("hu_thermoelastic", "m5", ("old", "new")),
         ("hu_breakage_index", "m6", ("new",))]
KEEP = re.compile(r"L[1-6]|peak|f_b|fb|onset|PASS|FAIL|baseline|pct|%", re.I)


def report():
    rc = 0
    names = {"old": "the 09-25 output (archived)", "binpre": "bin_pre re-run before the trim (a cross-check)",
             "new": "the production kernel, bin_1"}
    for t, row, sides in ORDER:
        print("=" * 110)
        print("%s  (%s)" % (t, row))
        for side in sides:
            r = subprocess.run([PY, os.path.join(HU, t, "test")] + RUNS[(t, side)], cwd=ROOT, capture_output=True,
                               text=True)
            out = (r.stdout or "") + (r.stderr or "")
            open(os.path.join(HERE, "K1_%s_%s.log" % (t, side)), "w").write(out)
            print("  --- %s: %s (rc %d) ---" % (side, names[side], r.returncode))
            for ln in out.splitlines():
                if KEEP.search(ln) and "Loading" not in ln and "Wrote plot" not in ln:
                    print("    " + ln[:180])
            if side == "new" and r.returncode != 0:
                rc = 1
    print("=" * 110)
    print("RESULT: the new runs' tests %s" % ("all PASS" if rc == 0 else "NOT all PASS (see above)"))
    return rc


if __name__ == "__main__":
    sys.exit({"a8": a8, "report": report}[sys.argv[1]]())
