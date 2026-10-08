#!/usr/bin/env python3
"""P1g-0 (j): the Kant and p0_steady_column anchors, BYTE for BYTE.

⚠⚠ `--check` CANNOT validate an anchor after a source change -- it re-reduces
output written by the OLD binary (standing rule 9). So the legs are re-run on
the shipping binary and the DATA FILES are compared to the pre-change archive.

⚠ Only the file classes the identity gate hashes are compared: plotfile cell
data, the event CSVs and thermo.dat. Logs are excluded -- they carry wall-clock
timings and absolute paths and were never byte-stable.

⚠ Files present in the ARCHIVE but not the re-run are reported separately and
are NOT counted as matches. P1f's count was inflated exactly this way: a
retired leg (`D150_2.0`, no longer in any `test`) was compared with itself,
turning 1,512 comparable files into a claimed 840.
"""
import filecmp, os, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "..", ".."))
PAIRS = [
    ("kant/sp_kant_onset",
     "tests/MMWSpalling/validation/kant/sp_kant_onset/output"),
    ("kant/sp_kant_onset_meiercfg",
     "tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/output"),
    ("kant/sp_kant_pressure_sweep",
     "tests/MMWSpalling/validation/kant/sp_kant_pressure_sweep/output"),
    ("p0_steady_column",
     "tests/MMWSpalling/validation/p0_steady_column/output"),
]
SUFFIX = ("thermo.dat",)
CONTAINS = ("Cell_D_", "_removal_events.csv", "_h_col_events.csv",
            "_clusters.csv", "_jet_profile.csv")


def hashable(name):
    return name.endswith(SUFFIX) or any(c in name for c in CONTAINS)


def walk(base):
    got = {}
    for d, _sub, files in os.walk(base):
        for f in files:
            if not hashable(f):
                continue
            p = os.path.join(d, f)
            got[os.path.relpath(p, base)] = p
    return got


def main():
    tot = same = diff = only_a = only_b = 0
    bad = []
    print("=" * 78)
    print("(j) ANCHOR BYTE-IDENTITY -- re-run on the shipping binary vs archive")
    print("=" * 78)
    for tag, rel in PAIRS:
        live = os.path.join(ROOT, rel)
        # ⚠ accept either step's archive, preferring the most recent: P1g-0's
        # anchors were re-run and proven byte-identical on that binary, so the
        # older archive is transitively a valid "before" for this step too.
        arch = next((live + sfx for sfx in (".p1g0b_archive", ".p1g0_archive")
                     if os.path.isdir(live + sfx)), live + ".p1g0_archive")
        if not (os.path.isdir(live) and os.path.isdir(arch)):
            print("  %-32s SKIP (live %s, archive %s)"
                  % (tag, os.path.isdir(live), os.path.isdir(arch)))
            continue
        A, B = walk(arch), walk(live)
        ka, kb = set(A), set(B)
        common = sorted(ka & kb)
        s = d = 0
        for k in common:
            if filecmp.cmp(A[k], B[k], shallow=False):
                s += 1
            else:
                d += 1
                if len(bad) < 15:
                    bad.append("%s :: %s" % (tag, k))
        oa, ob = sorted(ka - kb), sorted(kb - ka)
        tot += len(common); same += s; diff += d
        only_a += len(oa); only_b += len(ob)
        print("  %-32s %5d compared  %5d identical  %5d DIFFER   "
              "archive-only %3d  rerun-only %3d"
              % (tag, len(common), s, d, len(oa), len(ob)))
        for k in oa[:4]:
            print("        archive-only: %s" % k)
        for k in ob[:4]:
            print("        rerun-only  : %s" % k)
    print("-" * 78)
    print("  TOTAL %d comparable files: %d byte-identical, %d differ"
          % (tot, same, diff))
    print("  archive-only %d   rerun-only %d  (neither counted as a match)"
          % (only_a, only_b))
    for b in bad:
        print("    DIFFERS: %s" % b)
    ok = tot > 0 and diff == 0
    print("  (j) -> %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
