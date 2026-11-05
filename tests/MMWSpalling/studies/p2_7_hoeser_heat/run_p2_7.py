#!/usr/bin/env python3
"""P2-7 deliverable 2 -- the Hoeser legs at her MEASURED heat.

Every leg is P2-5 leg A's argv (`p2_5_hoeser/legs/A/argv.sh`, read token by
token, never rebuilt) with ONLY the keys in the table below replaced, plus
`plot_file`. The argv and the list of replaced/added keys are written to each
leg directory (standing rule 8: a leg's argv is part of its provenance; A's own
drops from P2-3 W1 are copied verbatim into DROPPED.txt).

LEGS (ACTIVE_STEP.md goal 2)
  A0   none -- A byte-for-byte on the 100 mm block (the b1 identity gate).
  HA   her Ø100 x 150 mm length: prob_hi z 0.15, n_cell z 150, nozzle_z0 0.2132
       (rock top 0.15 + her 63.2 mm, so s = 0.0632 and s_lid 0.0637 unchanged);
       h_expr unchanged (uniform plateau, her 2018 model's assumption).
  HG   HA's block + the Goldstein 1986 radial shape at her Ø9 disc mean.
  HGd  HG + spallation.sp.constraint_depth = D_REF_M from p2_2e_cref/D_REF.txt.

MODES
  run_p2_7.py A0 HA HG HGd        the legs (serial, -np 4, sentinel .done)
  run_p2_7.py --probe <legs>      ONE step (stop_time 1e-9) into probe/<leg>:
                                  rc + the unused-ParmParse block  (a3)
  run_p2_7.py --smoke <legs>      stop_time 2 into smoke/<leg>: witnesses
  run_p2_7.py --parseback         HG's shipped h_expr vs an independently
                                  built equivalent, 0.05 s each into parse/  (a3)
  run_p2_7.py --argv <leg>        print the argv and exit
"""
import argparse
import hashlib
import os
import shlex
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(STUD, "..", "..", ".."))
BASE = os.path.join(STUD, "p2_5_hoeser", "legs", "A", "argv.sh")
BASE_DROPPED = os.path.join(STUD, "p2_5_hoeser", "legs", "A", "DROPPED.txt")
D_REF = os.path.join(STUD, "p2_2e_cref", "D_REF.txt")
BIN = os.path.join(ROOT, "bin", "mmwspalling-3d-g++")

# ---- the packet's strings, LITERAL (goal 2) -------------------------------
H_A = ("if((r>0.02)and(s<0.0637),0.001,"
       "7100*pow(0.0632/max(0.0632,s),1.3185305455))")
H_G = ("if((r>0.02)and(s<0.0637),0.001,"
       "7272.3958*533/(533+44*pow(r/0.0079,1.285))"
       "*pow(0.0632/max(0.0632,s),1.3185305455))")
# The parse-back's INDEPENDENT construction of H_G: no `/` anywhere, the
# Goldstein factor as pow(1 + (44/533) (r/D)^1.285, -1), r/D as r*(1/D), and
# the standoff factor as s0^n * max(s0,s)^-n (the PARSER_DEFECT.txt rule).
H_G_EQV = ("if((r>0.02)and(s<0.0637),0.001,"
           "7272.3958*pow(1+{c}*pow(r*{iD},1.285),-1)"
           "*{s0n}*pow(max(0.0632,s),-1.3185305455))").format(
               c=repr(44.0 / 533.0), iD=repr(1.0 / 0.0079),
               s0n=repr(0.0632 ** 1.3185305455))
BLOCK150 = {"geometry.prob_hi": "0.05 0.05 0.15",
            "amr.n_cell": "50 50 150",
            "surface_patch.nozzle_z0": "0.2132"}


def d_ref():
    for ln in open(D_REF):
        if ln.startswith("D_REF_M="):
            return ln.split("=", 1)[1].strip()
    sys.exit(f"{D_REF}: no D_REF_M= line")


def legs():
    dr = d_ref()
    # ⚠ the packet writes the key as 0.002; D_REF.txt carries 0.0020. Same
    # number -- asserted, and the packet's literal is what ships.
    assert abs(float(dr) - 0.002) < 1e-15, f"D_REF_M {dr} != the packet's 0.002"
    return {
        "A0":  {},
        "HA":  dict(BLOCK150),
        "HG":  dict(BLOCK150, **{"surface_patch.h_expr": H_G}),
        "HGd": dict(BLOCK150, **{"surface_patch.h_expr": H_G,
                                 "spallation.sp.constraint_depth": "0.002"}),
    }


def base_tokens():
    """A's argv as launched: (prefix tokens up to the deck path, [key=val])."""
    toks = shlex.split(open(BASE).read())
    for n, t in enumerate(toks):
        if t.endswith("input_lab"):
            return toks[:n + 1], toks[n + 1:]
    sys.exit("A's argv: no input_lab deck path")


def argv_of(name, outdir, extra=None):
    """A's argv with the leg's keys replaced IN PLACE (order kept), new keys
    appended, plot_file pointed at outdir. Returns (cmd, replaced, added)."""
    pre, kv = base_tokens()
    want = dict(legs()[name])
    if extra:
        want.update(extra)
    keys = [t.split("=", 1)[0] for t in kv]
    assert all("=" in t for t in kv), "A's argv has a bare token"
    assert len(keys) == len(set(keys)), "A's argv has a duplicate key"
    out, replaced = [], []
    for t in kv:
        k = t.split("=", 1)[0]
        if k == "plot_file":
            out.append("plot_file=" + os.path.join(outdir, "run"))
        elif k in want:
            out.append(f"{k}={want[k]}")
            if t != out[-1]:
                replaced.append(k)
        else:
            out.append(t)
    added = [k for k in want if k not in keys]
    out += [f"{k}={want[k]}" for k in added]
    return pre + out, replaced, added


def sha(path, n=24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:n]


def write_prov(d, name, cmd, replaced, added):
    with open(os.path.join(d, "argv.sh"), "w") as f:
        f.write(" ".join(shlex.quote(c) for c in cmd) + "\n")
    with open(os.path.join(d, "DROPPED.txt"), "w") as f:
        f.write(f"P2-7 leg {name}: base = p2_5_hoeser/legs/A/argv.sh "
                f"(sha256 {sha(BASE)}), token for token\n")
        f.write(f"binary bin/mmwspalling-3d-g++ sha256 {sha(BIN)}\n\n")
        f.write("-- keys REPLACED in place (P2-7 goal 2)\n")
        for k in replaced or ["(none)"]:
            f.write(f"   {k}\n")
        f.write("\n-- keys ADDED\n")
        for k in added or ["(none)"]:
            f.write(f"   {k}\n")
        f.write("\n-- keys DROPPED by P2-7: none. plot_file re-pointed.\n")
        f.write("\n-- A's own drops from P2-3 W1, verbatim:\n\n")
        f.write(open(BASE_DROPPED).read())


def unused_block(log):
    L = open(log).read().splitlines()
    out, on = [], False
    for l in L:
        if l.startswith("Unused ParmParse Variables"):
            on = True
            continue
        if on:
            if not l.strip():
                break
            out.append(l.strip())
    return out


def run_one(name, d, cmd):
    os.makedirs(d, exist_ok=True)
    for f in (".done",):
        p = os.path.join(d, f)
        if os.path.exists(p):
            os.remove(p)
    t0 = time.time()
    log = os.path.join(d, "run.log")
    with open(log, "w") as f:
        rc = subprocess.call(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=ROOT)
    wall = time.time() - t0
    if rc == 0:
        open(os.path.join(d, ".done"), "w").write(f"{{'rc': 0, 'wall': {wall:.1f}}}\n")
    return rc, wall, log


def first_token(d, col):
    L = open(os.path.join(d, "run", "thermo.dat")).read().split("\n")
    h = L[0].split()
    rows = [l.split() for l in L[1:] if l.strip()]
    body = [r for r in rows if float(r[h.index("time")]) > 0]
    # ⚠ the PRINTED token, not a recomputed double (SR 7 iv: thermo is %g).
    return body[0][h.index(col)], body[0][h.index("time")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("legs", nargs="*")
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--parseback", action="store_true")
    ap.add_argument("--argv", action="store_true")
    a = ap.parse_args()
    # standing rule 11: the harness holds the host awake itself.
    if not a.argv:
        subprocess.Popen(["caffeinate", "-dis", "-w", str(os.getpid())])
    if a.parseback:
        print("P2-7 (a3) -- HG's shipped h_expr vs an independent construction,"
              " through the BINARY")
        print(f"  shipped    {H_G}")
        print(f"  equivalent {H_G_EQV}")
        res = {}
        for tag, expr in (("ship", H_G), ("eqv", H_G_EQV), ("HA", H_A)):
            d = os.path.join(HERE, "parse", tag)
            nm = "HA" if tag == "HA" else "HG"
            cmd, rep, add = argv_of(nm, d, {"surface_patch.h_expr": expr,
                                            "stop_time": "0.05",
                                            "amr.thermo.plot_int": "1"})
            os.makedirs(d, exist_ok=True)
            write_prov(d, f"parse/{tag}", cmd, rep, add)
            rc, wall, log = run_one(tag, d, cmd)
            tok = first_token(d, "patch_P_robin") if rc == 0 else (None, None)
            res[tag] = (rc, tok)
            print(f"  {tag:5s} rc={rc}  patch_P_robin at t={tok[1]}: {tok[0]}"
                  f"   unused: {unused_block(log)}")
        ok = (res["ship"][0] == 0 and res["eqv"][0] == 0
              and res["ship"][1][0] == res["eqv"][1][0])
        print(f"  shipped == equivalent (printed token): {'IDENTICAL' if ok else 'MISMATCH'}")
        if res["HA"][0] == 0 and res["ship"][0] == 0:
            print(f"  HG/HA step-1 patch_P_robin = "
                  f"{float(res['ship'][1][0])/float(res['HA'][1][0]):.5f}")
        sys.exit(0 if ok else 1)
    L = legs()
    for n in a.legs:
        if n not in L:
            sys.exit(f"unknown leg {n}; have {sorted(L)}")
    for n in a.legs:
        if a.probe:
            d = os.path.join(HERE, "probe", n)
            cmd, rep, add = argv_of(n, d, {"stop_time": "1e-9"})
        elif a.smoke:
            d = os.path.join(HERE, "smoke", n)
            cmd, rep, add = argv_of(n, d, {"stop_time": "2.0"})
        else:
            d = os.path.join(HERE, "legs", n)
            cmd, rep, add = argv_of(n, d)
        if a.argv:
            print(" \\\n  ".join(shlex.quote(c) for c in cmd))
            print(f"# replaced {rep}  added {add}")
            continue
        os.makedirs(d, exist_ok=True)
        write_prov(d, n, cmd, rep, add)
        rc, wall, log = run_one(n, d, cmd)
        print(f"  {n:4s} rc={rc} np={cmd[cmd.index('-np')+1]} wall {wall:7.1f} s"
              f"  replaced {rep} added {add}")
        if a.probe or a.smoke:
            print(f"       unused ParmParse: {unused_block(log)}")


if __name__ == "__main__":
    main()
