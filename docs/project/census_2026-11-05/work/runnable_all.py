#!/usr/bin/env python3
"""Census part 4 helper (read-only w.r.t. src/ and tests/): does each validation / extension / unit input deck
PARSE and take its first step on today's binary? Same method as docs/project/model_card/runnable.py: the deck
with stop_time = 1e-9 and plot_file pointed at the scratchpad, -np 4; rc and first ABORT line recorded.
Writes RUNNABLE_ALL.txt next to this script. Run from the repo root."""
import glob
import hashlib
import os
import re
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
SCR = os.environ.get("SCR", "/tmp/census_runnable")
T = "tests/MMWSpalling"
decks = []
for pat in ("validation/*/*/input*", "validation/*/input*", "validation/hu/hu_spall_onset/diag/input*",
            "validation/kant/grimsel_calibration/pre/input*", "extensions/*/input*", "unit/*/input*"):
    decks += sorted(glob.glob(os.path.join(ROOT, T, pat)))
decks = [d for d in decks if os.path.isfile(d) and not d.endswith(".py")]


def strip(s):
    return re.sub(r"\x1b\[[0-9;]*m", "", s)


sha = hashlib.sha256(open(BIN, "rb").read()).hexdigest()[:16]
L = ["census part 4 -- each deck on today's binary (sha %s), stop_time 1e-9, -np 4, plot_file -> scratch" % sha,
     "=" * 110]
for d in decks:
    rel = os.path.relpath(d, ROOT)
    out = os.path.join(SCR, re.sub(r"[^A-Za-z0-9_.-]", "_", rel))
    os.makedirs(out, exist_ok=True)
    cmd = ["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN, d, "stop_time=1e-9",
           "plot_file=" + os.path.join(out, "run")]
    try:
        p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, errors="replace", timeout=300)
        log, rc = strip(p.stdout + p.stderr), p.returncode
    except subprocess.TimeoutExpired:
        log, rc = "", "TIMEOUT"
    msg = ""
    for ln in log.splitlines():
        if "ABORT" in ln or "not found in database" in ln or "Error" in ln:
            msg = ln.strip()[:220]
            break
    verdict = "RUNS" if rc == 0 else ("TIMEOUT" if rc == "TIMEOUT" else "FAILS (rc %s)" % rc)
    L.append("  %-14s %s" % (verdict, rel))
    if msg and rc != 0:
        L.append("        " + msg)
    print(L[-1] if not msg else L[-2] + "\n" + L[-1], flush=True)
open(os.path.join(HERE, "RUNNABLE_ALL.txt"), "w").write("\n".join(L) + "\n")
