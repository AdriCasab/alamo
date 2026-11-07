#!/usr/bin/env python3
"""P2-18: does every kept deck still START on the final binary after the A2 cleanup? Each deck
(the same set a2_apply.py edited: unit/, validation/, extensions/ input* and *.inputs-including
decks, + the s1 gate deck) runs -np 4 with stop_time = 1e-9; rc and the abort line are recorded
and compared with the census's pre-change list (work/RUNNABLE_ALL.txt: RUNS / ABORT per deck)."""
import os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import a2_apply as A
ROOT = A.ROOT; BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")
SCR = "/private/tmp/claude-501/-Users-tzetze20-amr-tools-alamo/79a0a434-469e-4979-a234-ccc76ca8c68d/scratchpad/p2_18/runnable"
os.makedirs(SCR, exist_ok=True)
census = {}
for ln in open(os.path.join(ROOT, "docs/project/census_2026-11-05/work/RUNNABLE_ALL.txt")):
    m = re.match(r"\s+(RUNS|ABORT\S*|\S+)\s+(tests/\S+)", ln)
    if m: census[m.group(2)] = m.group(1)
rows, bad = [], 0
for p in A.deck_files():
    if p.endswith(".inputs"): continue
    rel = os.path.relpath(p, ROOT)
    out = os.path.join(SCR, rel.replace("/", "_"))
    r = subprocess.run(["mpirun", "--oversubscribe", "--bind-to", "none", "-np", "4", BIN, p, "stop_time=1e-9", "plot_file=" + out],
                       cwd=ROOT, capture_output=True, text=True)
    log = r.stdout + r.stderr
    why = ""
    if r.returncode != 0:
        m = re.search(r"(specified but not used[\s\S]{0,120}|required value for \S+ missing|retired[^\n]{0,100}|Abort[^\n]{0,100})", log)
        why = re.sub(r"\s+", " ", m.group(1))[:140] if m else log.strip().split("\n")[-1][:140]
    was = census.get(rel, "?")
    verdict = "RUNS" if r.returncode == 0 else "ABORT"
    flag = "" if (was == "?" or was.startswith(verdict)) else "  <-- CHANGED vs census"
    if flag: bad += 1
    rows.append("  %-6s (census %-6s) %s%s %s" % (verdict, was, rel, flag, why))
    print(rows[-1], flush=True)
print("%d decks; %d changed verdict vs the census" % (len(rows), bad))
