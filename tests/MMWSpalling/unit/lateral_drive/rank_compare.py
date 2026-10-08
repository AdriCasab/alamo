#!/usr/bin/env python3
"""P1g-2 deliverable 14 / (b7): decomposition independence, done correctly.

usage:  rank_compare.py <legA> <legB>        # names under output/
        rank_compare.py --path <dirA> <dirB>

⚠⚠ WHY `diff -r` CANNOT BE USED, AND WHY THIS FILE EXISTS. Two packets in a row
specified the rank check as `diff -r` across rank counts. That can never pass,
for two structural reasons that have nothing to do with the physics:

  1. `Checkpoint` and `Header` record the BOX LAYOUT. At np = 2 the domain is cut
     into 2 boxes, at np = 4 into 4. The files therefore differ by construction.
  2. Plotfile cell data is written one `Cell_D_NNNNN` file PER BOX. np = 4 writes
     more of them than np = 2, so the file NAMED `Cell_D_00000` holds a DIFFERENT
     REGION of the domain in each run. Comparing them by name is meaningless --
     it is not a weaker test, it is a test of nothing.

So the comparison is made on **decomposition-invariant output** only:

  * `thermo.dat`   -- globally reduced scalars, one row per step
  * the event CSVs -- per-column records, globally reduced before they are written

⚠⚠ TWO COLUMNS ARE EXCLUDED BY NAME, NOT BY TOLERANCE. `ledger_err` and
`ledger_err_solid` are *differences of large rank-reduced sums*, so their
floating-point value depends on summation order at the last bit. Measured across
np = 2 vs 4 they differ by 2.22e-16 -- one double epsilon -- while all other
columns agree bit for bit. Excluding two named columns and printing their
residual is honest; applying a blanket tolerance to all 60+ columns would hide a
real regression somewhere else. ⚠ The exclusion is REPORTED on every run, so a
reader always sees what was not compared.
"""
import os, sys

RESIDUAL_COLS = ("ledger_err", "ledger_err_solid")
HERE = os.path.dirname(os.path.abspath(__file__))


def thermo_path(d):
    for c in (os.path.join(d, "run", "thermo.dat"), os.path.join(d, "thermo.dat")):
        if os.path.exists(c):
            return c
    return None


def compare_thermo(A, B):
    pa, pb = thermo_path(A), thermo_path(B)
    if not pa or not pb:
        print("    thermo.dat MISSING (%s / %s)" % (bool(pa), bool(pb)))
        return False, 0, 0.0
    la = open(pa).read().split("\n")
    lb = open(pb).read().split("\n")
    ha, hb = la[0].split(), lb[0].split()
    if ha != hb:
        onlya = [c for c in ha if c not in hb]
        onlyb = [c for c in hb if c not in ha]
        print("    COLUMN SETS DIFFER -- only in A: %s | only in B: %s"
              % (onlya, onlyb))
        return False, 0, 0.0
    ra = [r.split() for r in la[1:] if r.strip()]
    rb = [r.split() for r in lb[1:] if r.strip()]
    if len(ra) != len(rb):
        print("    ROW COUNTS DIFFER: %d vs %d" % (len(ra), len(rb)))
        return False, 0, 0.0
    bad, resid = {}, 0.0
    for x, y in zip(ra, rb):
        for i, (u, v) in enumerate(zip(x, y)):
            if u == v:
                continue
            try:
                d = abs(float(u) - float(v))
            except ValueError:
                d = float("inf")
            if ha[i] in RESIDUAL_COLS:
                resid = max(resid, d)
            else:
                bad[ha[i]] = max(bad.get(ha[i], 0.0), d)
    npc = len(ha) - sum(1 for c in ha if c in RESIDUAL_COLS)
    print("    thermo.dat: %d physics columns compared, %d differ"
          % (npc, len(bad)))
    for k, v in sorted(bad.items(), key=lambda kv: -kv[1]):
        print("        %-28s max abs %.3e" % (k, v))
    print("    EXCLUDED BY NAME (differences of rank-reduced sums): %s"
          % ", ".join(c for c in RESIDUAL_COLS if c in ha))
    print("        their max discrepancy here: %.3e" % resid)
    return not bad, npc, resid


def compare_events(A, B):
    """the per-column event logs, byte for byte -- these carry the lateral state."""
    ok = True
    found = 0
    # ⚠ The event CSVs are named from the plot_file PREFIX, so depending on how
    # the leg was launched they sit either INSIDE the plot directory (prefix
    # ".../LEG/run") or BESIDE it as a sibling (prefix ".../LEG"). Check both, or
    # the comparison silently finds nothing and the row passes on thermo.dat
    # alone -- the per-column lateral state is exactly what this row is for.
    cands = []
    for suffix in ("removal_events.csv", "h_col_events.csv"):
        cands.append((os.path.join(A, "run_" + suffix),
                      os.path.join(B, "run_" + suffix), "run_" + suffix))
        cands.append((A + "_" + suffix, B + "_" + suffix,
                      os.path.basename(A) + "_" + suffix))
    for pa, pb, suffix in cands:
        if not (os.path.exists(pa) and os.path.exists(pb)):
            continue
        found += 1
        import filecmp
        same = filecmp.cmp(pa, pb, shallow=False)
        print("    %-26s byte-identical: %s" % (suffix, same))
        ok = ok and same
    if found == 0:
        print("    (no event CSVs found beside either leg)")
    return ok, found


def main():
    av = sys.argv[1:]
    if len(av) >= 3 and av[0] == "--path":
        A, B = av[1], av[2]
    elif len(av) >= 2:
        A = os.path.join(HERE, "output", av[0])
        B = os.path.join(HERE, "output", av[1])
    else:
        print(__doc__)
        return 2
    print("=" * 74)
    print("(b7) DECOMPOSITION INDEPENDENCE  --  %s  vs  %s"
          % (os.path.basename(A), os.path.basename(B)))
    print("=" * 74)
    for d in (A, B):
        if not os.path.isdir(d):
            print("  FAIL: missing %s" % d)
            return 1
    t_ok, npc, resid = compare_thermo(A, B)
    e_ok, nfound = compare_events(A, B)
    print("-" * 74)
    # ⚠ non-vacuity: if neither comparand existed the run must not read as a pass
    if npc == 0 and nfound == 0:
        print("  UNSCOREABLE: nothing comparable was found. NOT a pass.")
        return 1
    ok = t_ok and e_ok
    print("  -> %s" % ("RANK-INDEPENDENT on every compared quantity"
                       if ok else "DIFFERS -- see above"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
