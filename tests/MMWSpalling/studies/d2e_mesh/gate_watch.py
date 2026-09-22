#!/usr/bin/env python3
"""Wait until B0_1mm passes 150 s, read the 1 mm step gate, and on a FAIL stop
the 8 ms 1 mm runs at once (D2e §4 decision rule). Writes gate_watch.log."""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
TH = os.path.join(HERE, "output", "B0_1mm", "thermo.dat")
log = open(os.path.join(HERE, "gate_watch.log"), "a")
while True:
    try:
        last = open(TH).read().rstrip().splitlines()[-1].split()[0]
        if float(last) >= 150.5:
            break
    except (OSError, IndexError, ValueError):
        pass
    time.sleep(30)
r = subprocess.run([PY, os.path.join(HERE, "run.py"), "--gate"], capture_output=True, text=True)
log.write(time.strftime("%F %T") + "\n" + r.stdout + r.stderr)
if r.returncode != 0:
    for n in ("Pa_1mm", "B0_1mm"):
        k = subprocess.run([PY, os.path.join(HERE, "run.py"), "--kill", n], capture_output=True, text=True)
        log.write(k.stdout + k.stderr)
    log.write("GATE FAILED: 8 ms 1 mm runs stopped\n")
else:
    log.write("GATE PASSED\n")
log.close()
