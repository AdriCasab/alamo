#!/usr/bin/env python3
"""P2-26 a4 -- the two B1 parse fixes, probed on BOTH binaries (one step each).

  L1 the flake-cap check reads the DECK's spall.h_spall_max (it read the member initialiser 0.010, because the parse
     ran after the check). ⚠ The packet's probe (flake_h_max 0.02, h_spall_max 0.015) aborts on BOTH binaries
     (0.02 >= 0.010 and >= 0.015), so it cannot show the order moved; it is kept as a record. The two discriminating
     pairs (z0 section 5) are the bar:
       P1a  0.012 / 0.015 : before ABORTS (0.012 >= 0.010), after RUNS   (0.012 < 0.015)
       P1b  0.008 / 0.005 : before RUNS   (0.008 < 0.010),  after ABORTS (0.008 >= 0.005), naming both keys
     on unit/gradient_blind/input, -np 1, stop_time 0.008 (one step).
  L2 surface_patch.side_power_interval is parsed under every jet_closure: the CAL deck (jet_closure = none) with
     side_power_interval = 1.0 aborts before (the strict scope: "specified but not used") and runs after;
     stop_time = one step (timestep 0.008).
"before" = studies/p2_24_melt_vapour/bin_p2_24 (98a261d3): P2-25 changed no src/, so its parse order is P2-25's
(c052c348 itself was not archived); "after" = bin/ (BINARY_1.txt). Writes PROBES.txt.
Usage: probes.py   (exit 1 on a failed bar)
"""
import hashlib
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = os.path.join(ROOT, "tests", "MMWSpalling")
BEFORE = os.path.join(T, "studies", "p2_24_melt_vapour", "bin_p2_24", "mmwspalling-3d-g++")
AFTER = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
GB = os.path.join(T, "unit", "gradient_blind", "input")
CAL = os.path.join(T, "validation", "kant", "grimsel_calibration", "input")
OUT = []


def say(s):
    print(s)
    OUT.append(s)


def leg(binary, deck, extra):
    d = tempfile.mkdtemp(prefix="p226probe_")
    argv = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "1", binary, deck] + extra + \
           ["stop_time=0.008", "plot_file=%s" % os.path.join(d, "run")]
    r = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True)
    txt = re.sub(r"\x1b\[[0-9;]*m", "", (r.stdout or "") + (r.stderr or ""))
    return r.returncode, txt


def main():
    for b in (BEFORE, AFTER):
        say("binary %-80s sha256 %s" % (os.path.relpath(b, ROOT), hashlib.sha256(open(b, "rb").read()).hexdigest()[:16]))
    bad = 0
    rows = [
        ("P1a", GB, ["spall.flake_h_max=0.012", "spall.h_spall_max=0.015"], True, False),
        ("P1b", GB, ["spall.flake_h_max=0.008", "spall.h_spall_max=0.005"], False, True),
        ("P1c (packet's, record)", GB, ["spall.flake_h_max=0.02", "spall.h_spall_max=0.015"], None, None),
        ("P2", CAL, ["surface_patch.side_power_interval=1.0"], True, False),
    ]
    for tag, deck, extra, want_before, want_after in rows:
        res = []
        for b, want in ((BEFORE, want_before), (AFTER, want_after)):
            rc, txt = leg(b, deck, extra)
            ab = rc != 0
            names = ("spall.flake_h_max" in txt and "spall.h_spall_max" in txt) if tag.startswith("P1") \
                else ("side_power_interval" in txt and "specified but not used" in txt)
            ok = want is None or (ab == want and (not ab or names))
            bad += 0 if ok else 1
            res.append("%s rc %d%s%s" % ("ABORTS" if ab else "runs  ", rc,
                                         " (names %s)" % ("both keys" if tag.startswith("P1") else
                                                          "the key, unused") if ab and names else "",
                                         "" if ok else " ** WRONG"))
        say("  %-24s %-62s before: %-40s after: %s" % (tag, " ".join(extra) + "  on " + os.path.relpath(deck, T),
                                                       res[0], res[1]))
    say("a4 %s (L1: P1a old aborts/new runs, P1b old runs/new aborts naming both keys; L2: old refuses the key, "
        "new runs)" % ("PASS" if not bad else "FAIL"))
    open(os.path.join(HERE, "PROBES.txt"), "w").write("\n".join(OUT) + "\n")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
