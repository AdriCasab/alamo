#!/usr/bin/env python3
"""P2-29 smoke: every new deck runs a few steps on bin/ (-np 4) into S/smoke/<name>; prints rc, each level's dt message
(dt, dt_lim, the ratio) and the Temp range of the last plotfile. Usage: smoke.py [name ...] (default: every deck)."""
import glob
import os
import re
import subprocess
import sys

import yt

yt.set_log_level(50)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
T = "tests/MMWSpalling/"
DECKS = {
    "ZQ": (T + "validation/zhang/zhang_oglesby/input_phase", ["stop_time=0.5", "amr.plot_int=10"]),
    "beam_tau": (T + "unit/beam_tau/input_phase", []),
    "beam": (T + "unit/beam/input_phase", ["stop_time=0.002", "amr.plot_int=10"]),
    "equilibrium": (T + "unit/equilibrium/input_phase", ["stop_time=20", "amr.plot_int=10"]),
    "stefan": (T + "unit/stefan/input_phase", ["stop_time=0.002", "amr.plot_int=10"]),
    "convective_patch": (T + "unit/convective_patch/input_phase", ["stop_time=0.1", "amr.plot_int=5"]),
    "hu_cond_g": (T + "validation/hu/hu_conduction/input_granite2_phase", ["stop_time=1.0", "amr.plot_int=5"]),
    "hu_cond_s": (T + "validation/hu/hu_conduction/input_sandstone2_phase", ["stop_time=1.0", "amr.plot_int=5"]),
    "hu_te_g": (T + "validation/hu/hu_thermoelastic/input_granite2_phase", ["stop_time=0.4", "amr.plot_int=1"]),
    "hu_te_s": (T + "validation/hu/hu_thermoelastic/input_sandstone2_phase", ["stop_time=0.2", "amr.plot_int=1"]),
    "hu_bi_g37": (T + "validation/hu/hu_breakage_index/input_granite2_37_phase", ["stop_time=0.4", "amr.plot_int=1"]),
    "hu_bi_s37": (T + "validation/hu/hu_breakage_index/input_sandstone2_37_phase", ["stop_time=0.2", "amr.plot_int=1"]),
    "hu_bi_s89": (T + "validation/hu/hu_breakage_index/input_sandstone2_89_phase", ["stop_time=0.2", "amr.plot_int=1"]),
}
MSG = re.compile(r"explicit time-step limit, level (\d+): dt = (\S+) s, dt_lim = dx\^2/\((\d+) alpha\) = (\S+) s")


def main(argv):
    names = argv or list(DECKS)
    bad = 0
    for n in names:
        deck, extra = DECKS[n]
        out = os.path.join(HERE, "smoke", n)
        argv_ = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", "bin/mmwspalling-3d-g++", deck] + \
            extra + ["plot_file=" + out]
        r = subprocess.run(argv_, cwd=ROOT, capture_output=True, text=True)
        txt = (r.stdout or "") + (r.stderr or "")
        os.makedirs(os.path.join(HERE, "smoke"), exist_ok=True)
        open(out + ".log", "w").write(" ".join(argv_) + "\n" + txt)
        msgs = ["L%s dt %s lim(%s) %s ratio %.3f" % (m[0], m[1], m[2], m[3], float(m[1]) / float(m[3])
                                                    if float(m[3]) not in (0.0,) else 0.0)
                for m in MSG.findall(txt)]
        tr = "-"
        pfs = sorted(glob.glob(os.path.join(out, "*cell")))
        if r.returncode == 0 and pfs:
            ad = yt.load(pfs[-1]).all_data()
            tr = "T %.3f .. %.3f K" % (float(ad["Temp"].min()), float(ad["Temp"].max()))
        warn = "WARN" if "exceeds the explicit limit" in txt else ""
        print("%-16s rc %d  %s  %s %s" % (n, r.returncode, "; ".join(msgs) or "(no dt message)", tr, warn))
        if r.returncode != 0:
            bad += 1
            print("    " + "\n    ".join([ln for ln in txt.splitlines() if "bort" in ln or "rror" in ln][:4]))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
