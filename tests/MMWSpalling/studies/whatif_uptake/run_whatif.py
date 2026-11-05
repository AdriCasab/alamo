#!/usr/bin/env python3
"""What-if legs: P2-7 HG with the h_expr leading constant x S (see README.md). Serial, -np 4, .done sentinel."""
import os, shlex, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BASE = os.path.join(HERE, "..", "p2_7_hoeser_heat", "legs", "HG", "argv.sh")
H0G = "7272.3958"
SCALES = {"s0581": 0.581, "s150": 0.15, "s087": 0.15 * 0.581, "s100": 0.10,
          "s058": 0.0581, "s034": 0.0581 * 0.581}

def argv(name):
    toks = shlex.split(open(BASE).read())
    out, hit = [], 0
    for t in toks:
        if t.startswith("surface_patch.h_expr="):
            assert t.count(H0G + "*") == 1
            t = t.replace(H0G + "*", f"{float(H0G) * SCALES[name]:.6f}*"); hit += 1
        elif t.startswith("amr.plot_int="):
            t = "amr.plot_int=3750"; hit += 1
        elif t.startswith("plot_file="):
            t = "plot_file=" + os.path.join(HERE, "legs", name, "run"); hit += 1
        out.append(t)
    assert hit == 3, hit
    return out

def main():
    subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])
    for name in sys.argv[1:] or list(SCALES):
        d = os.path.join(HERE, "legs", name); os.makedirs(d, exist_ok=True)
        cmd = argv(name)
        open(os.path.join(d, "argv.sh"), "w").write(" ".join(shlex.quote(c) for c in cmd) + "\n")
        open(os.path.join(d, "DROPPED.txt"), "w").write(
            f"what-if leg {name}: base = p2_7_hoeser_heat/legs/HG/argv.sh token for token;\n"
            f"changed: h_expr leading constant {H0G} x {SCALES[name]:.6f}; amr.plot_int 625 -> 3750; plot_file.\n")
        p = os.path.join(d, ".done")
        if os.path.exists(p): os.remove(p)
        t0 = time.time()
        with open(os.path.join(d, "run.log"), "w") as f:
            rc = subprocess.call(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=ROOT)
        w = time.time() - t0
        print(f"  {name:6s} S {SCALES[name]:.5f} rc={rc} wall {w:6.1f} s", flush=True)
        if rc == 0: open(p, "w").write(f"rc 0 wall {w:.1f}\n")
if __name__ == "__main__":
    main()
