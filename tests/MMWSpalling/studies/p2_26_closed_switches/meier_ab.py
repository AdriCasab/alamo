#!/usr/bin/env python3
"""P2-26 m2 -- the Meier pilot test (validation/meier/sp_meier_pilot/test) scored BEFORE and AFTER, so its record's
change is stated with the old and the new value of each scored quantity.
  before: the pre-image decks (src_pre_2/.../input_energy, input_drilling; spall.flake_coherence_length = 0.004) on
          studies/p2_24_melt_vapour/bin_p2_24 (P2-25's code)
  after : the edited decks (the coherence line deleted, C2) on bin/ (BINARY_2.txt)
The test module is loaded as a library and its own scorers run unchanged (series, score_energy, score_drilling);
only the binary, the deck and the output directory differ (plot_file on the command line, a scratch tree here --
the test's output/ is not touched). Writes MEIER_AB.txt.
"""
import contextlib
import importlib.machinery
import importlib.util
import io
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
TD = os.path.join(ROOT, "tests", "MMWSpalling", "validation", "meier", "sp_meier_pilot")
PRE = os.path.join(HERE, "src_pre_2", "tests", "MMWSpalling", "validation", "meier", "sp_meier_pilot")
OLD_BIN = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_24_melt_vapour", "bin_p2_24", "mmwspalling-3d-g++")
NEW_BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
SCR = os.path.join(HERE, "meier_ab")


def load():
    l = importlib.machinery.SourceFileLoader("meier_test", os.path.join(TD, "test"))
    s = importlib.util.spec_from_loader("meier_test", l)
    m = importlib.util.module_from_spec(s)
    l.exec_module(m)
    return m


def main():
    m = load()
    out = []
    for side, binary, deckdir in (("before", OLD_BIN, PRE), ("after", NEW_BIN, TD)):
        res = {}
        for name in ("energy", "drilling"):
            case = dict(m.CASES[name])
            case["input"] = os.path.join(deckdir, os.path.basename(case["input"]))
            case["outdir"] = os.path.join(SCR, side, name)
            os.makedirs(os.path.dirname(case["outdir"]), exist_ok=True)
            cmd = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", binary, case["input"],
                   "plot_file=%s" % case["outdir"]]
            with open(case["outdir"] + ".log", "w") as lg:
                lg.write(" ".join(cmd) + "\n")
                lg.flush()
                rc = subprocess.call(cmd, cwd=ROOT, stdout=lg, stderr=subprocess.STDOUT)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                ok = (m.score_energy if name == "energy" else m.score_drilling)(m.series(case)) if rc == 0 else False
            res[name] = (rc, ok, buf.getvalue())
        out.append("=== %s (%s, decks %s)" % (side, os.path.relpath(binary, ROOT), os.path.relpath(deckdir, ROOT)))
        for name, (rc, ok, txt) in res.items():
            out.append("  case %s: rc %d, %s" % (name, rc, "PASS" if ok else "FAIL"))
            out += ["    " + ln for ln in txt.strip().split("\n") if ln.strip()]
    open(os.path.join(HERE, "MEIER_AB.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
