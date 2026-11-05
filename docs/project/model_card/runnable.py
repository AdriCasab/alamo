#!/usr/bin/env python3
"""Model card, pass 2: does each deck still RUN on today's binary? Each deck's own recorded argv (argv.sh of
its leg) or its input file is started with stop_time = 1e-9 and plot_file pointed at the scratchpad; the rc
and the first abort / unused-input message are recorded. Nothing in the repository is written except
RUNNABLE.txt next to this script. Run from the repo root."""
import os
import re
import shlex
import subprocess
import hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
SCR = "/private/tmp/claude-501/-Users-tzetze20-amr-tools-alamo/8ac6af7a-7af4-46dd-a6e4-ee2b54287af0/scratchpad/runnable"
S = "tests/MMWSpalling/studies"
V = "tests/MMWSpalling/validation"
# label, kind, path: kind argv = a leg's argv.sh; kind input = an input file run with no overrides
DECKS = [
    ("KANT_onset input_model", "input", f"{V}/kant/sp_kant_onset_meiercfg/input_model"),
    ("KANT_onset input_weibull", "input", f"{V}/kant/sp_kant_onset_meiercfg/input_weibull"),
    ("KANT_dref G2.0_130 (P2-2e)", "argv", f"{S}/p2_2e_cref/legs/G2.0_130/argv.sh"),
    ("KANT_onset V80 (P2-2e)", "argv", f"{S}/p2_2e_cref/legs/V80/argv.sh"),
    ("KANT_cal grimsel (P2-13)", "input", f"{V}/kant/grimsel_calibration/input"),
    ("HOSER_47 G1 (P2-11)", "argv", f"{S}/p2_11_heat/legs/G1/argv.sh"),
    ("HOSER input_lab", "input", f"{V}/hoeser/input_lab"),
    ("HOSER_0 HG (P2-7)", "argv", f"{S}/p2_7_hoeser_heat/legs/HG/argv.sh"),
    ("HOSER_0_S s100 (whatif)", "argv", f"{S}/whatif_uptake/legs/s100/argv.sh"),
    ("MEIER ctrl8/m13 (P2-2f)", "argv", f"{S}/p2_2c_meier/ctrl8/m13/argv.sh"),
    ("MEIER_dref SD (P2-2e)", "argv", f"{S}/p2_2e_cref/legs/SD/argv.sh"),
    ("MEIER input_feet", "input", f"{V}/meier/sp_meier_pilot/input_feet"),
    ("MEIER input_feet_d2c", "input", f"{V}/meier/sp_meier_pilot/input_feet_d2c"),
    ("MEIER input_drilling", "input", f"{V}/meier/sp_meier_pilot/input_drilling"),
    ("MMW_ZHANG input", "input", f"{V}/zhang/zhang_oglesby/input"),
    ("HU onset granite2", "input", f"{V}/hu/hu_spall_onset/input_granite2"),
    ("HU end_to_end granite2", "input", f"{V}/hu/hu_end_to_end/input_granite2"),
]


def strip_ansi(s):
    return re.sub(r"\x1b\[[0-9;]*m", "", s)


def argv_for(kind, path, out):
    if kind == "input":
        return ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN, os.path.join(ROOT, path),
                "stop_time=1e-9", "plot_file=" + out]
    toks = shlex.split(open(os.path.join(ROOT, path)).read())
    if toks and toks[0] == "caffeinate":
        toks = toks[2:] if toks[1].startswith("-") else toks[1:]
    res = []
    for t in toks:
        k = t.split("=", 1)[0]
        if k == "plot_file":
            res.append("plot_file=" + out)
        elif k == "stop_time":
            res.append("stop_time=1e-9")
        elif t.endswith("mmwspalling-3d-g++"):
            res.append(BIN)
        else:
            res.append(t)
    if not any(t.startswith("stop_time=") for t in res):
        res.append("stop_time=1e-9")
    return res


def main():
    os.makedirs(SCR, exist_ok=True)
    sha = hashlib.sha256(open(BIN, "rb").read()).hexdigest()[:16]
    L = ["Model card -- does each deck run on today's binary (sha %s)? stop_time 1e-9, plot_file -> scratchpad" % sha,
         "=" * 110]
    for lab, kind, path in DECKS:
        out = os.path.join(SCR, re.sub(r"[^A-Za-z0-9_.-]", "_", lab))
        os.makedirs(out, exist_ok=True)
        if not os.path.exists(os.path.join(ROOT, path)):
            L.append("  %-30s MISSING %s" % (lab, path))
            continue
        cmd = argv_for(kind, path, os.path.join(out, "run"))
        p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, errors="replace")
        log = strip_ansi(p.stdout + p.stderr)
        open(os.path.join(out, "run.log"), "w").write(log)
        msgs = []
        lines = log.splitlines()
        for i, ln in enumerate(lines):
            if ("ABORT" in ln or "EXCEPTION" in ln or "not used" in ln or "retired" in ln
                    or "not found in database" in ln) and "Backtrace" not in ln and "MPI_ABORT" not in ln \
                    and "invoking MPI_ABORT" not in ln:
                msgs.append(re.sub(r"\s+", " ", ln.strip())[:230])
        status = "RUNS" if p.returncode == 0 else "ABORTS (rc %d)" % p.returncode
        if p.returncode != 0 and "not found in database" in log:
            status = "TEMPLATE"        # a base input whose harness supplies required keys (amr.n_cell ...)
        L.append("  %-30s %-15s %s" % (lab, status, path))
        seen = set()
        for m in msgs:
            if m not in seen:
                L.append("        " + m)
                seen.add(m)
            if len(seen) >= 8:
                break
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(HERE, "RUNNABLE.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
