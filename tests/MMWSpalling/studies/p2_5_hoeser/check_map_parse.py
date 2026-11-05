#!/usr/bin/env python3
"""P2-5 rev 2 (a1 vii) -- END-TO-END PARSER CHECK of every leg's h_expr.

⚠⚠ WHY THIS EXISTS. The packet's a1 "branch witness" evaluates the derivation's
own mirror of each string -- my arithmetic on both sides of the comparison -- so
it cannot detect the PARSER disagreeing with the formula. It did not: rev 2's
first draft shipped `max(s,0.0632)` in a divisor, which this binary's expression
parser mis-evaluates (the max call becomes a multiplier), and the map delivered
0.26 W where the design asked ~240 W. Only the smoke's POWER reading caught it.

This closes that gap without trusting either side's arithmetic: for every leg it
runs the SHIPPED string and an INDEPENDENTLY CONSTRUCTED mathematically
identical string through the binary for one step, and requires the two
`patch_P_robin` values to agree to the last printed digit. A parser pathology in
either construction shows up as a disagreement.

  shipped:    h0*exp(-(r*r)/(2 sig^2))      * pow(s0/max(s0,s), n)
  equivalent: h0*exp(-(r*r)*(1/(2 sig^2)))  * pow(s0,n)*pow(max(s0,s),-n)

Both forms were spot-verified against a constant control before being used here.
"""
import os
import shlex
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from p2_5_legs import LEGS, ORDER, SHAPES, S0, R_OPEN, LID_H, s_lid
import derive_lab_map as D

TMP = os.environ.get("TMPDIR", "/tmp") + "/p2_5_parsecheck"
BASE = os.path.join(HERE, "legs", "SMOKE", "argv.sh")


def equiv(leg, n_of_K):
    """The same h, built a different way: no `/` by a max call at all."""
    dz, _, K, z0, shape = LEGS[leg]
    sig = D.sigma_m(shape)
    h0s = D.h0_of(shape)
    n = n_of_K[K]
    dec = f"pow({D.g(S0)},{D.g(n)})*pow(max({D.g(S0)},s),{D.g(-n)})"
    core = (f"{D.g(h0s)}*{dec}" if sig is None
            else f"{D.g(h0s)}*exp(-(r*r)*{D.g(1.0/(2.0*sig*sig))})*{dec}")
    return f"if((r>{D.g(R_OPEN)})and(s<{D.g(s_lid(leg))}),{D.g(LID_H)},{core})"


def run(expr, tag):
    base = shlex.split(open(BASE).read())
    out = os.path.join(TMP, tag)
    os.makedirs(out, exist_ok=True)
    cmd = []
    for t in base:
        if t.startswith("surface_patch.h_expr="):
            cmd.append("surface_patch.h_expr=" + expr)
        elif t.startswith("stop_time="):           cmd.append("stop_time=0.05")
        elif t.startswith("amr.thermo.plot_int="): cmd.append("amr.thermo.plot_int=1")
        elif t.startswith("plot_file="):           cmd.append("plot_file=" + out + "/run")
        else:                                      cmd.append(t)
    rc = subprocess.call(cmd, stdout=open(out + ".log", "w"),
                         stderr=subprocess.STDOUT)
    if rc != 0:
        return rc, None
    L = open(os.path.join(out, "run", "thermo.dat")).read().split("\n")
    h = L[0].split()
    rows = [l.split() for l in L[1:] if l.strip()]
    body = [r for r in rows if float(r[h.index("time")]) > 0]
    # ⚠ the PRINTED token, not a recomputed double (SR 7 iv: thermo is %g).
    return rc, body[0][h.index("patch_P_robin")]


def main():
    De = D.read_De()
    her = D.her_table()
    n_of_K = {K: D.fit_n(K, De, her)[0] for K in (30, 20)}
    print("P2-5 rev 2 (a1 vii) -- the shipped string vs an independent "
          "construction, through the BINARY")
    print("=" * 92)
    print(f"  {'leg':>9s} {'rc':>3s} {'P_robin shipped':>18s}"
          f" {'P_robin equivalent':>20s}  verdict")
    bad = 0
    for leg in ORDER:
        rc1, p1 = run(D.h_expr(leg, n_of_K), leg + "_ship")
        rc2, p2 = run(equiv(leg, n_of_K), leg + "_eqv")
        ok = (rc1 == 0 and rc2 == 0 and p1 == p2)
        bad += 0 if ok else 1
        print(f"  {leg:>9s} {max(rc1, rc2):>3d} {str(p1):>18s} {str(p2):>20s}"
              f"  {'ok' if ok else 'MISMATCH'}")
    print("=" * 92)
    if bad:
        print(f"FAIL (a1 vii): {bad} leg(s) disagree -- the parser does not "
              f"evaluate the shipped formula.")
        sys.exit(1)
    print("PASS (a1 vii): every leg's shipped h_expr and an independently")
    print("  constructed equivalent produce the SAME printed patch_P_robin, so")
    print("  the parser evaluates the formula the derivation intends.")


if __name__ == "__main__":
    main()
