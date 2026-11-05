#!/usr/bin/env python3
"""P2-2b (b2): the dense tables against an INDEPENDENT twin.

WHAT IS COMPARED, AND WHY IT IS NOT CIRCULAR. `b6_dump` prints
`strain_of` / `E_of` / `sigma_thermal` from the REAL `src/Numeric/Material/
Grimsel.H`. This script recomputes all three from `grimsel.py`, which reads
Table 3 out of the reference file and interpolates it by a direct scan -- no
dense table anywhere. So the two implementations share no code, and agreement
is evidence rather than tautology.

⚠⚠ WHY THE `OFF` SWEEP IS THE ONE THAT MATTERS. P2-2b replaced `E_of`'s 10-row
scan with a dense 1 K lookup. The claim is that linear interpolation between
nodes is EXACT because E/E_0 is piecewise linear and every Table 3 kink lands
on an integer degree. At integer theta -- the packet's 708-point sweep -- every
interpolation scheme agrees by construction, so that sweep has no power over
this change. The `OFF` sweep steps by 980/707 K, landing between nodes almost
everywhere, and is what could actually fail.

Usage: python3 score_b2.py > B2_TABLES.txt   (exit 1 if any bar is missed)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import grimsel as G  # noqa: E402

BAR = 1e-12


def err(got, want):
    """Relative where the reference is non-zero, absolute at an exact zero."""
    return abs(got - want) / abs(want) if want != 0.0 else abs(got - want)


def main():
    path = os.path.join(HERE, "B2_CPP.txt")
    if not os.path.exists(path):
        print("FAIL: B2_CPP.txt missing -- run `sh standalone/build.sh && "
              "./b6_dump > B2_CPP.txt` first")
        return 1

    rows = {"NODE": [], "OFF": []}
    with open(path) as fh:
        next(fh)
        for line in fh:
            p = line.split()
            if len(p) == 5 and p[0] in rows:
                rows[p[0]].append(tuple(float(x) for x in p[1:]))

    print("P2-2b (b2) THE DENSE TABLES vs THE INDEPENDENT TWIN")
    print("bar: <= %.0e relative, on both sweeps\n" % BAR)
    fails = []
    for sweep in ("NODE", "OFF"):
        pts = rows[sweep]
        if not pts:
            fails.append("%s sweep is empty" % sweep)
            continue
        worst = {"strain": (0.0, 0.0), "E": (0.0, 0.0), "sigma": (0.0, 0.0)}
        for th, strain, E, sigma in pts:
            for key, got, want in (
                    ("strain", strain, G.thermal_strain(th)),
                    ("E",      E,      G.E_of(th)),
                    ("sigma",  sigma,  G.sigma_formula(th + 273.15)[0])):
                e = err(got, want)
                if e > worst[key][0]:
                    worst[key] = (e, th)
        note = ("integer theta -- grid NODES, no power over the dense lookup"
                if sweep == "NODE" else
                "theta = 20 + i*980/707 -- BETWEEN nodes, this is the real test")
        print("  --- %s sweep, %d points (%s) ---" % (sweep, len(pts), note))
        for key in ("strain", "E", "sigma"):
            e, th = worst[key]
            ok = e <= BAR
            print("      %-7s worst %.3e at theta %8.3f   %s"
                  % (key, e, th, "ok" if ok else "FAIL"))
            if not ok:
                fails.append("%s/%s worst %.3e at theta %.3f (bar %.0e)"
                             % (sweep, key, e, th, BAR))
        print()

    # Non-vacuity: the OFF sweep must actually land off-node, or it silently
    # degenerates into the NODE sweep and the gate tests nothing.
    off_th = [th for th, _, _, _ in rows["OFF"]]
    n_on_node = sum(1 for th in off_th if abs(th - round(th)) < 1e-9)
    print("  non-vacuity: %d of %d OFF points sit on a 1 K node"
          % (n_on_node, len(off_th)))
    # The bar is "this sweep probes BETWEEN nodes", so it is a fraction, not a
    # count: a step of 980/707 K lands on an integer occasionally and that is
    # harmless. Requiring >= 90 % off-node states the intent; a count of 2 was
    # an arbitrary number that would have failed a sweep doing its job.
    if len(off_th) and n_on_node > 0.10 * len(off_th):
        fails.append("the OFF sweep lands on %d of %d nodes (> 10 %%) -- it is "
                     "not probing between them" % (n_on_node, len(off_th)))

    if fails:
        print("\nFAILURES:")
        for f in fails:
            print("  " + f)
        print("\nFAIL")
        return 1
    print("\nstrain, stiffness and the one stress formula all reproduce an")
    print("independent twin to %.0e BETWEEN grid nodes as well as on them." % BAR)
    print("\nPASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
