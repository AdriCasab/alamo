#!/usr/bin/env python3
"""P1g-2c (c1)/(c1b): is the fix BIT-EXACT where one row touches a cell?

  usage: replay_anchors.py --tag diag
         replay_anchors.py --tag fixed --against diag

Hashes the deterministic outputs of every leg the committed fixture runs, on the
build currently in `bin/`, and stores them under `anchors_<tag>.json`. With
`--against` it diffs two such manifests file by file.

⚠⚠ WHY THE BASELINE IS THE ORD-2 BINARY AND NOT P1g-2's ARCHIVED OUTPUTS. The
toolchain is **not byte-reproducible** (`studies/p1e_staircase/MANIFEST.md:108`),
so a cross-build byte comparison can fail for reasons that have nothing to do
with this fix -- and a build-induced mismatch is indistinguishable from a real
break. The only sound comparison is diagnostic-build vs fixed-build, where the
ONLY difference is the write-back.

⚠⚠ AND WHY THE `evidence_p1g2` LEGS ARE NOT REPLAYED HERE. The packet asks for
them to be "replayed from their `metadata`". They cannot be: `metadata` is a
RESOLVED parameter dump (every key including defaults), not an argv, and the
runner that launched them was a session scratch script that no longer exists --
`evidence_p1g2/PROVENANCE.txt` records binaries and compilation stamps but no
command line. Rebuilding an argv from a resolved dump is RECONSTRUCTION, not
replay, and a wrong reconstruction would produce exactly the spurious mismatch
this script exists to avoid. So (c1b) is reported as unscoreable, with the cause
named, and (c1) is scored on a STRICTLY LARGER replayable anchor set: every leg
the committed test launches, which includes the new single-touch shift leg.

⚠ Only deterministic files are hashed -- the per-event CSVs and `thermo.dat`.
`metadata` carries a compilation timestamp and would differ on every build for a
reporting reason; plotfiles are not committed and no scored row reads them
(`evidence_p1g2/PROVENANCE.txt`).
"""
import argparse
import hashlib
import json
import os
import sys

ROOT = "/Users/tzetze20/amr_tools/alamo"
OUT = os.path.join(ROOT, "tests/MMWSpalling/unit/lateral_drive/output")
HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(ROOT, "bin/mmwspalling-3d-g++")

# the legs whose single-row write-back must not move. G1_FLAT_SHIFT is the one
# that matters most: it is the only committed leg that RUNS the shift on a wall
# where every open face is +x, i.e. exactly one row per cell.
SINGLE_TOUCH = ("G1_FLAT_SHIFT", "G1_OFF", "G1_BANK", "G1_NEG", "G1_FLAKEE",
                "STEP1D", "STEP1R", "RISER", "DIAG")
# ⚠ G1_DIAG45 is EXPECTED to change: it is the two-axis leg and the fix is
# precisely a change there. Listed so it is never silently counted as an anchor.
EXPECTED_TO_CHANGE = ("G1_DIAG45",)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def collect():
    """{relpath: sha256} over deterministic files under the fixture output."""
    out = {}
    if not os.path.isdir(OUT):
        return out
    for dp, _, fns in os.walk(OUT):
        for fn in fns:
            if not (fn.endswith("_removal_events.csv") or fn == "thermo.dat"
                    or fn.endswith("_h_col_events.csv")):
                continue
            p = os.path.join(dp, fn)
            out[os.path.relpath(p, OUT)] = sha(p)
    return out


def leg_of(rel):
    return rel.split(os.sep)[0].replace("_removal_events.csv", "") \
             .replace("_h_col_events.csv", "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--against")
    a = ap.parse_args()

    man = collect()
    bsha = sha(BIN)
    path = os.path.join(HERE, "anchors_%s.json" % a.tag)
    json.dump({"binary": bsha, "files": man}, open(path, "w"), indent=1,
              sort_keys=True)

    P = print
    P("=" * 78)
    P("P1g-2c ANCHOR MANIFEST  tag=%s" % a.tag)
    P("=" * 78)
    P("  binary %s…" % bsha[:16])
    P("  %d deterministic files hashed under" % len(man))
    P("    %s" % OUT)
    legs = sorted({leg_of(r) for r in man})
    P("  legs present: %s" % ", ".join(legs))
    P("")

    if not a.against:
        P("  baseline written to %s" % os.path.basename(path))
        P("")
        P("  ⚠ (c1b) the `evidence_p1g2` legs are NOT in this manifest and")
        P("    cannot be: no argv was recorded for them (PROVENANCE.txt carries")
        P("    binaries and compilation stamps only) and their launcher was a")
        P("    session scratch script. (c1b) is UNSCOREABLE -- the cause is a")
        P("    provenance gap in P1g-2, not this fix. Same class as Leg B's")
        P("    missing binary in P1g-2b.")
        return 0

    old = json.load(open(os.path.join(HERE, "anchors_%s.json" % a.against)))
    P("  comparing against tag=%s (binary %s…)"
      % (a.against, old["binary"][:16]))
    if old["binary"] == bsha:
        P("  ⚠⚠ THE TWO BINARIES ARE THE SAME SHA -- the comparison is vacuous.")
        P("     Did the fix build actually happen?")
    P("")
    of = old["files"]
    allk = sorted(set(man) | set(of))
    changed, missing, same = [], [], []
    for k in allk:
        if k not in man or k not in of:
            missing.append(k)
        elif man[k] != of[k]:
            changed.append(k)
        else:
            same.append(k)

    P("  %-44s %s" % ("file", "verdict"))
    for k in allk:
        v = ("ONLY IN ONE MANIFEST" if k in missing
             else ("**DIFFERS**" if k in changed else "identical"))
        P("  %-44s %s" % (k[:44], v))
    P("")
    P("  identical %d   differing %d   one-sided %d"
      % (len(same), len(changed), len(missing)))
    P("")

    # ---- the verdict, by leg class -------------------------------------
    bad = sorted({leg_of(k) for k in changed} & set(SINGLE_TOUCH))
    exp = sorted({leg_of(k) for k in changed} & set(EXPECTED_TO_CHANGE))
    P("  SINGLE-TOUCH legs that changed : %s" % (", ".join(bad) or "NONE"))
    P("  legs EXPECTED to change        : %s" % (", ".join(exp) or "none"))
    P("")
    ok = not bad
    if ok:
        P("  (c1) **PASS** -- the fix is BIT-EXACT on every single-touch leg,")
        P("       including G1_FLAT_SHIFT, which runs the shift on a one-row")
        P("       wall. So a cell touched once is still ASSIGNED, not")
        P("       incremented, and no frozen flat-wall anchor drifts.")
    else:
        P("  (c1) **FAIL** -- a single-touch leg moved. That is the signature of")
        P("       a blanket `Hw += (rNew - rH)`: algebraically identity-")
        P("       preserving, but `old + (new - old) != new` in floating point.")
    if not exp:
        P("")
        P("  ⚠⚠ AND G1_DIAG45 DID NOT CHANGE, which would mean the fix did not")
        P("     reach the two-axis path at all -- check it is not behind a key.")
        ok = False
    P("")
    P("  ⚠ (c1b) unscoreable: the `evidence_p1g2` legs have no recorded argv.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
