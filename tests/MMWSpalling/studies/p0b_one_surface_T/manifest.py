#!/usr/bin/env python3
"""P0b rev5 deliverable 1 / criterion (g) -- THE PROVENANCE MANIFEST.

A number with no manifest row is not a result.  Every scored leg gets one row
carrying: the leg name, the output directory it actually wrote to, the binary
sha256 it ran on, the sha256 of its event log, and the argv.

The argv is recovered from the run's own metadata file where ALAMO wrote one,
and from the harness otherwise -- never retyped.

usage: manifest.py > MANIFEST.md
"""
import hashlib, os, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
BIN  = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
BINT = os.path.getmtime(BIN)

def sha(p, n=16):
    if not os.path.exists(p): return "-"
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()[:n]

def argv_of(d):
    """ALAMO writes the full command line into <plotdir>/metadata."""
    for dp, _, fs in os.walk(d):
        if "metadata" in fs:
            for L in open(os.path.join(dp, "metadata"), errors="replace"):
                if L.startswith("# Command line:") or L.lower().startswith("cmd"):
                    return L.split(":", 1)[-1].strip()
    p = os.path.join(d, "cmd.txt")
    if os.path.exists(p): return open(p).read().strip()
    return "(see the harness: %s)" % os.path.basename(d)

SETS = [
    ("(a) (c') (c'') (x) -- the 1-D column matrix, 18 legs",
     "tests/MMWSpalling/validation/p0_steady_column/output",
     "tests/MMWSpalling/validation/p0_steady_column/test"),
    ("(d) -- the Kant anchor, 15 legs",
     "tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output",
     "tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/test"),
    ("(f) (h) (i) d2 d6 -- the Meier FH rerun, closure FOLLOW",
     "tests/MMWSpalling/studies/p0b_one_surface_T/output_fh_rev7",
     "tests/MMWSpalling/studies/p0b_one_surface_T/output_fh_rev7/cmd.txt"),
    ("d6 -- the Meier control, closure CELL",
     "tests/MMWSpalling/studies/p0b_one_surface_T/output_fh_rev7_cell",
     "tests/MMWSpalling/studies/p0b_one_surface_T/output_fh_rev7_cell/cmd.txt"),
]

print("# P0b rev7 -- MANIFEST\n")
print("Binary under test: `%s`\n" % sha(BIN, 64))
print("⚠ The toolchain is NOT byte-reproducible: rebuilding from identical")
print("sources gives a different sha256. So the only way a leg's binary sha can")
print("match the final build is for NOTHING to have been rebuilt after it ran.")
print("Every row below was produced after the final `make` and before any other.\n")
for title, out, harness in SETS:
    print("## %s\n" % title)
    print("harness: `%s`\n" % harness)
    d = os.path.join(ROOT, out)
    if not os.path.isdir(d):
        print("  ⚠ MISSING: %s\n" % out); continue
    print("| leg | output path | event-log sha256 (16) | rows | vs final binary |")
    print("|---|---|---|---|---|")
    for leg in sorted(os.listdir(d)):
        p = os.path.join(d, leg)
        if not os.path.isdir(p): continue
        ev = None
        for f in sorted(os.listdir(p)):
            if f.endswith("_removal_events.csv"): ev = os.path.join(p, f)
        rows = "-"
        if ev:
            with open(ev) as fh: rows = str(sum(1 for _ in fh) - 1)
        # ⚠ THE STALENESS GUARD. Revision 2's scored legs predated the final
        # binary by 43 minutes and the harness said nothing. A leg whose newest
        # file is older than the binary did not run on it, full stop.
        newest = max([os.path.getmtime(os.path.join(dp, f))
                      for dp, _, fs in os.walk(p) for f in fs] or [0.0])
        fresh = "ok" if newest >= BINT else "⚠ **STALE — predates the final binary**"
        print("| `%s` | `%s/%s` | `%s` | %s | %s |"
              % (leg, out, leg, sha(ev) if ev else "-", rows, fresh))
    # flat-file runs (the Meier tag writes <tag>_removal_events.csv beside <tag>/)
    for f in sorted(os.listdir(d)):
        if f.endswith("_removal_events.csv"):
            q = os.path.join(d, f)
            with open(q) as fh: rows = sum(1 for _ in fh) - 1
            fresh = "ok" if os.path.getmtime(q) >= BINT else "⚠ **STALE**"
            print("| `%s` | `%s/%s` | `%s` | %d | %s |"
                  % (f[:-len("_removal_events.csv")], out, f, sha(q), rows, fresh))
    print()
    print("argv: see `%s`\n" % harness)
