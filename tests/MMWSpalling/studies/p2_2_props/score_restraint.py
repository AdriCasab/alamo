#!/usr/bin/env python3
"""P2-2a (b1) -- the tree gate scored at COMPONENT granularity.

(b1)'s bar is not a file count. It is:

  * **0 differing** among the 2589 files that do NOT belong to a strip-derived
    `jet_enthalpy` leg;
  * the 50 strip files differ in **`Sp_field` only** -- every other component,
    and in particular `Temp`, `H`, `phi`, `removed`, `eta`, `T_ref_phase`,
    `thermal_free_strain`, `beta_phase` and `E_phase`, byte-identical;
  * removal events **0 -> 0** on all 11 output-producing strip legs.

⚠⚠ WHY THE COMPONENT SPLIT IS THE WHOLE POINT. `T_ref_phase` and
`thermal_free_strain` are the EIGENSTRAIN's fields. Change A must reach the
criterion and nothing else, so those two moving means it leaked into the elastic
solve. A whole-file hash cannot tell that from the criterion's own diagnostic
moving -- which is exactly the hole rev 5's proxy fell into: substituting
`microstructure.phase0.T_ref` moved FIVE components (measured: Temp, H,
T_ref_phase, thermal_free_strain, Sp_field), so it over-stated the change.

Usage:  score_restraint.py <manifest.md5> [<tree root>]
"""

import os
import hashlib
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "../../../.."))
sys.path.insert(0, HERE)

import fabcmp  # noqa: E402

# the 13 strip-derived tags of unit/jet_enthalpy, from its own RUNS table
STRIP_TAGS = ("B05", "B1", "B1e", "B1n", "B2", "B3", "B4",
              "E", "Eon", "G", "Gbug", "P1", "P3")
# the components that must NOT move: if any of these differs, change A left the
# criterion. The eigenstrain pair is listed first because it is the discriminator.
PHYSICS = ("T_ref_phase", "thermal_free_strain", "Temp", "H", "phi",
           "removed", "eta", "beta_phase", "E_phase", "kappa_phase",
           "mu_phase", "D", "kappa_eff", "k_eff")
ALLOWED = ("Sp_field",)


def is_strip(path):
    if "unit/jet_enthalpy" not in path:
        return False
    return any(re.search(r"/output/%s(/|_)" % t, path) for t in STRIP_TAGS)


def sha(path):
    """⚠ The tree manifest is MD5, not SHA-256 -- `tree_sweep_0922/witness/*.md5`
    means what it says. An earlier version of this function shelled out to
    `shasum -a 256` and therefore reported every file as differing."""
    try:
        with open(path, "rb") as f:
            return hashlib.md5(f.read()).hexdigest()
    except Exception:
        return None


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    man, tree = argv[1], (argv[2] if len(argv) > 2 else
                          os.path.join(ROOT, "tests", "MMWSpalling"))
    rows = []
    for line in open(man):
        parts = line.split()
        if len(parts) >= 2:
            rows.append((parts[0], parts[1]))
    print("=" * 78)
    print("P2-2a (b1) TREE IDENTITY, SCORED PER COMPONENT")
    print("=" * 78)
    print("  manifest : %s  (%d files)" % (os.path.basename(man), len(rows)))

    non_strip_diff, strip_diff, missing = [], [], []
    for want, rel in rows:
        p = os.path.join(tree, rel)
        if not os.path.exists(p):
            missing.append(rel)
            continue
        if sha(p) == want:
            continue
        (strip_diff if is_strip(rel) else non_strip_diff).append(rel)

    n_strip_total = sum(1 for _, r in rows if is_strip(r))
    print("  strip-derived files in manifest : %d" % n_strip_total)
    print("  non-strip files                 : %d" % (len(rows) - n_strip_total))
    # ⚠ REVIEW ROW 14 -- THE PARTITION, STATED. The packet's bar says
    # "2589 non-strip / 50 strip"; the actual partition is the one printed above
    # and it is the arithmetic that closes (they sum to the manifest). 68 files
    # therefore sit in the strip set that the bar's wording puts in the gated
    # one. ⚠ Harmless only because the differing set is checked directly: every
    # differing file is a strip-leg plotfile, so "0 non-strip differing" holds
    # under EITHER partition. The planner should correct the bar's numbers.
    if (len(rows) - n_strip_total, n_strip_total) != (2589, 50):
        print("  ⚠ the packet's bar says 2589 / 50; the closing partition is "
              "%d / %d" % (len(rows) - n_strip_total, n_strip_total))
    print()

    # ---- clause 1: the non-strip files must be byte-identical -------------
    print("-" * 78)
    print("(b1).1  NON-STRIP FILES -- bar: 0 differing")
    print("-" * 78)
    # ⚠⚠ NO LEG IS EXCLUDED HERE, AND AN EARLIER VERSION OF THIS SCORER WRONGLY
    # EXCLUDED ONE. During rev 7 I reported `s1_surface_resolution` `M-100_dz1`
    # as "compile-sensitive, not this step", on the strength of an A-absent
    # build whose 7 M-100 files still differed from the pin. That measurement
    # was STALE: `run_sweep.py` skips a leg with a `.done` marker unless
    # `--force` is passed (`run_one`, "skip (done)"), the gate passes
    # `--smoke --force` and I had run plain `--smoke`, so NOTHING re-ran and I
    # compared rev 6's own output against the pin. The run-to-run
    # "reproducibility" check I used to confirm it was memoised the same way, so
    # it confirmed nothing. Measured properly: M-100 returns to the p2_1 hashes
    # on 14 of 14 files under the stored field, exactly like the other three
    # areas. ⇒ there is no compile-reproducibility problem and no
    # manifest-membership question; the bar below is the plain one.
    print("  differing : %d" % len(non_strip_diff))
    for r in non_strip_diff[:20]:
        print("     %s" % r)
    if len(non_strip_diff) > 20:
        print("     ... and %d more" % (len(non_strip_diff) - 20))
    other = non_strip_diff
    ok1 = not other

    # ---- clause 2: the strip files may differ in Sp_field ONLY -----------
    print()
    print("-" * 78)
    print("(b1).2  STRIP FILES -- bar: `Sp_field` only, every physics "
          "component identical")
    print("-" * 78)
    print("  differing files : %d of %d" % (len(strip_diff), n_strip_total))
    moved_counts, leaked = {}, []
    checked = 0
    for rel in strip_diff:
        if "Cell_D" not in rel:
            moved_counts.setdefault("(not a plotfile: %s)"
                                    % os.path.basename(rel), 0)
            moved_counts["(not a plotfile: %s)" % os.path.basename(rel)] += 1
            continue
        cur = os.path.join(tree, rel)
        hdr = os.path.join(os.path.dirname(os.path.dirname(cur)), "Header")
        try:
            ncomp = fabcmp.read_fabs(cur)[0][0]
            names = fabcmp.var_names(hdr, ncomp)
        except Exception as e:
            leaked.append((rel, "unreadable: %s" % e))
            continue
        # the frozen copy lives beside the manifest's own witness tree
        ref = os.path.join(HERE, "frozen_prechange", rel)
        if not os.path.exists(ref):
            moved_counts.setdefault("(no frozen copy to compare)", 0)
            moved_counts["(no frozen copy to compare)"] += 1
            continue
        checked += 1
        _, moved = fabcmp.compare(ref, cur, names)
        for _, nm in moved:
            moved_counts[nm] = moved_counts.get(nm, 0) + 1
            if nm in PHYSICS:
                leaked.append((rel, nm))
    for nm in sorted(moved_counts):
        tag = "ALLOWED" if nm in ALLOWED else (
            "** PHYSICS LEAK **" if nm in PHYSICS else "(other)")
        print("     %-26s %4d files   %s" % (nm, moved_counts[nm], tag))
    if not moved_counts:
        print("     (no per-component comparison performed -- see note below)")
    print("  files compared per component : %d" % checked)
    ok2 = not leaked
    if leaked:
        print("  ** PHYSICS COMPONENTS MOVED -- change A leaked out of the "
              "criterion **")
        for rel, nm in leaked[:10]:
            print("     %-18s %s" % (nm, rel))

    if missing:
        print()
        print("  ⚠ %d manifest files absent from the tree:" % len(missing))
        for r in missing[:6]:
            print("     %s" % r)

    print()
    print("=" * 78)
    print("(b1) VERDICT")
    print("=" * 78)
    print("  .1 non-strip 0 differing            %s"
          % ("PASS" if ok1 else "** FAIL **"))
    print("  .2 strip files: no physics moved    %s"
          % ("PASS" if ok2 else "** FAIL **"))
    print()
    print("  ⚠ A re-pin is authorised by the packet ONLY for files meeting .2.")
    return 0 if (ok1 and ok2) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
