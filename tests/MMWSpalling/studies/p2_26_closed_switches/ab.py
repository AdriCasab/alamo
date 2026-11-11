#!/usr/bin/env python3
"""P2-26 m2 / m3 -- every deck stage 2 edited OUTSIDE the gate, run BEFORE and AFTER, same argv, same horizon.

  before: the pre-image deck (src_pre_2/tests/...) on studies/p2_24_melt_vapour/bin_p2_24 (98a261d3; P2-25 changed no
          src/, so this is P2-25's code; stage 1 is proven 0/2639 + face IDENTICAL against it and no edited deck sets a
          B3 switch, so before -> after isolates C2 + C3)
  after : the edited deck on bin/ (BINARY_2.txt)
Each pair runs -np 4 to the horizon in HORIZON (the deck's own stop_time where that is cheap, else a stated prefix),
plot_file in a scratch directory, amr.plot_int large (the first and last plotfile only). Compared: every
*_events.csv byte for byte (row counts, first differing line), thermo.dat's shared columns row by row, the last
plotfile's Level_0 fabs by md5. A deck whose pair is IDENTICAL over its full stop_time keeps its test verdict by
identity (m2: "passes unchanged"); a pair that differs is LIVE and its record changes BY DESIGN (m2: the old and new
value of what moved is the line printed here). m3 = the AFTER leg exits rc 0.
Usage: ab.py [DECK ...]  > AB.txt   (default: all; exit 1 if an AFTER leg fails to run)
"""
import glob
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
OLD_BIN = os.path.join(T, "studies", "p2_24_melt_vapour", "bin_p2_24", "mmwspalling-3d-g++")
NEW_BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
PRE = os.path.join(HERE, "src_pre_2", "tests", "MMWSpalling")
SCR = os.path.join(HERE, "ab_legs")
# deck -> horizon [s] (None = the deck's own stop_time); a prefix where the full run costs hours
HORIZON = {
    "unit/gradient_blind/input": None,
    "unit/surface_T_recon/input": None,
    "unit/seed_void/input": None,
    "validation/p0_steady_column/input": None,
    "validation/kant/sp_kant_onset_meiercfg/input_model": None,
    "validation/kant/sp_kant_onset_meiercfg/input_weibull": None,
    "validation/meier/sp_meier_pilot/input_energy": 4.0,
    "validation/meier/sp_meier_pilot/input_drilling": 1.0,
    "validation/meier/sp_meier_pilot/input_drilling_25mm": 2.0,
    "validation/meier/sp_meier_pilot/input_2d_dev_melt": 15.0,
    "extensions/hu_spall_onset_mmwbeam_lefm/input_pulsed_brittle_400w": None,
    "extensions/hu_spall_onset_mmwbeam_lefm/input_pulsed_brittle_400w_hires": 8.0,
    "extensions/hu_spall_onset_mmwbeam_lefm/input_wellbore_3d": 8.0,
    "extensions/hu_spall_onset_mmwbeam_lefm/input_wellbore_3d_fine": 0.5,
    "extensions/hu_spall_onset_mmwbeam_lefm/input_wellbore_3d_smallgrain": 0.5,
}
TIMEOUT = 1800   # s per leg; a leg past it is reported, never compared


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def stop_time(deck):
    m = re.search(r"^stop_time\s*=\s*([^\s#]+)", open(deck).read(), re.M)
    return m.group(1) if m else "?"


def leg(binary, deck, tag, horizon):
    d = os.path.join(SCR, tag)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    argv = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", binary, deck,
            "plot_file=%s" % os.path.join(d, "run"), "amr.plot_int=1000000000"]
    if horizon is not None:
        argv.append("stop_time=%r" % horizon)
    open(os.path.join(d, "argv.sh"), "w").write(" ".join(argv) + "\n")
    t0 = time.time()
    try:
        r = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True, timeout=TIMEOUT)
        rc = r.returncode
        open(os.path.join(d, "run.log"), "w").write((r.stdout or "") + (r.stderr or ""))
    except subprocess.TimeoutExpired:
        rc = "timeout"
    return rc, time.time() - t0, d


def compare(da, db):
    out, same = [], True
    ca = sorted(os.path.basename(p) for p in glob.glob(os.path.join(da, "run*_events.csv")))
    cb = sorted(os.path.basename(p) for p in glob.glob(os.path.join(db, "run*_events.csv")))
    for n in sorted(set(ca) | set(cb)):
        if n not in ca or n not in cb:
            out.append("%s only %s" % (n, "before" if n in ca else "after"))
            same = False
            continue
        A = open(os.path.join(da, n)).read().split("\n")
        B = open(os.path.join(db, n)).read().split("\n")
        if A == B:
            out.append("%s identical (%d rows)" % (n.replace("run_", ""), len(A) - 2))
        else:
            k = next((i for i in range(min(len(A), len(B))) if A[i] != B[i]), min(len(A), len(B)))
            out.append("%s rows %d -> %d, first differing line %d" % (n.replace("run_", ""), len(A) - 2, len(B) - 2, k))
            same = False
    ta, tb = os.path.join(da, "run", "thermo.dat"), os.path.join(db, "run", "thermo.dat")
    if os.path.exists(ta) and os.path.exists(tb):
        A, B = open(ta).read().split("\n"), open(tb).read().split("\n")
        ha, hb = A[0].split(), B[0].split()
        timers = {"fl_t_path", "fr_t_removal", "fr_t_step"}
        idx = [(ha.index(c), hb.index(c)) for c in hb if c in ha and c not in timers]
        gone = [c for c in ha if c not in hb]
        first = None
        for k in range(1, min(len(A), len(B))):
            a, b = A[k].split(), B[k].split()
            if a and b and any(a[i] != b[j] for i, j in idx):
                first = k
                break
        out.append("thermo: %d -> %d rows%s, shared columns %s" % (
            len(A) - 2, len(B) - 2, (", columns gone %s" % gone) if gone else "",
            "identical" if first is None else "first differ at row %d (%s)" % (
                first, ", ".join(hb[j] for i, j in idx if A[first].split()[i] != B[first].split()[j])[:120])))
        same = same and first is None and len(A) == len(B)
    pa = sorted(p for p in glob.glob(os.path.join(da, "run", "*cell")))
    pb = sorted(p for p in glob.glob(os.path.join(db, "run", "*cell")))
    if pa and pb:
        fa = sorted(glob.glob(os.path.join(pa[-1], "Level_0", "Cell_D_*")))
        fb = sorted(glob.glob(os.path.join(pb[-1], "Level_0", "Cell_D_*")))
        eq = [md5(x) == md5(y) for x, y in zip(fa, fb)] if len(fa) == len(fb) else [False]
        out.append("last plotfile %s vs %s: %s" % (os.path.basename(pa[-1]), os.path.basename(pb[-1]),
                                                   "fabs identical" if all(eq) else "fabs differ"))
        same = same and all(eq) and os.path.basename(pa[-1]) == os.path.basename(pb[-1])
    return same, out


def main():
    decks = sys.argv[1:] or list(HORIZON)
    bad = 0
    print("P2-26 m2/m3 A/B: before = pre-image deck on %s, after = edited deck on bin/ (%s)" % (
        os.path.relpath(OLD_BIN, ROOT), hashlib.sha256(open(NEW_BIN, "rb").read()).hexdigest()[:16]))
    for rel in decks:
        hz = HORIZON[rel]
        old_deck, new_deck = os.path.join(PRE, rel), os.path.join(T, rel)
        tag = rel.replace("/", "__")
        ra, sa, da = leg(OLD_BIN, old_deck, tag + "__before", hz)
        rb, sb, db = leg(NEW_BIN, new_deck, tag + "__after", hz)
        full = hz is None
        print("\n%s  (horizon %s s%s; before rc %s %.0f s, after rc %s %.0f s)" % (
            rel, stop_time(new_deck) if full else hz, " = the deck's stop_time" if full else
            " of %s, a PREFIX" % stop_time(new_deck), ra, sa, rb, sb))
        if rb != 0:
            bad += 1
            print("  ** m3: the edited deck does not run on the new binary")
        if ra == 0 and rb == 0:
            same, lines = compare(da, db)
            for ln in lines:
                print("  " + ln)
            print("  => %s" % ("IDENTICAL%s" % (" over the full run: the test verdict is unchanged by identity"
                                               if full else " over the prefix") if same else
                               "LIVE: the deleted line moved this deck (its record changes BY DESIGN)"))
        else:
            print("  => not compared (a leg did not complete)")
    print("\nm3 %s (every edited deck runs on the new binary)" % ("PASS" if not bad else "FAIL"))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
