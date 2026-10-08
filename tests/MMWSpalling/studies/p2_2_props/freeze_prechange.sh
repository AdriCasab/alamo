#!/bin/bash
# P2-2a deliverable 1 / (a0) -- FREEZE THE PRE-CHANGE PLOTFILES.
#
# ⚠⚠ WHAT "PRE-CHANGE" MEANS HERE, AND WHY IT IS NOT THE p2_1 TREE.
# rev 6 ran the full gate with change A live, and `jet_enthalpy` rmtree's and
# re-runs its own legs, so the genuine pre-P2-2a plotfiles were overwritten on
# 2026-10-07 at 16:36. There is no archive. "Freeze before editing source"
# arrived one revision too late.
#
# What replaces it is BETTER than a snapshot, and measured rather than assumed:
# an **A-ABSENT** build -- changes B, C and D landed, change A reverted -- whose
# output was verified against the frozen `p2_1.md5` pin:
#
#     unit/jet_enthalpy      288 of 288 files MATCH   (so QF returns to zero)
#     unit/beam_void_closure 169 of 169 files MATCH
#     validation/.../dev2d    64 of  64 files MATCH
#     studies/s1_...         135 of 142 MATCH, 7 DIFFER  <- see below
#
# ⇒ changes B, C and D are **byte-inert**, which rev 6 inferred but never
# measured, and the A-absent tree reproduces the pin everywhere the step can
# reach. So these copies are a valid pre-change reference for attributing
# change A at component granularity, and unlike a snapshot they are
# REGENERABLE.
#
# ⚠⚠ THE M-100 ROW ABOVE WAS A STALE READ, AND THE CORRECTION MATTERS.
# I first reported `s1_surface_resolution` `M-100_dz1` as differing on the
# A-absent build and concluded it was compile-sensitive rather than this step's
# doing. It was not. `run_sweep.py` SKIPS any leg carrying a `.done` marker
# unless `--force` is passed (`run_one` returns "skip (done)"); the gate invokes
# it as `--smoke --force` and I had run plain `--smoke`. So no leg re-ran and I
# compared rev 6's own leftover output against the pin.
# ⚠ The run-to-run "reproducibility" check I ran to validate that conclusion was
# memoised in exactly the same way, so it confirmed nothing -- the non-vacuity
# trap, applied to my own verification rather than to the mechanism.
# ⇒ MEASURED PROPERLY (gate, --force): under the stored field M-100 returns to
# the p2_1 hashes on 14 of 14 files, like the other three areas. All four legs
# that moved under rev 6 were change A, and the toolchain's byte-reproducibility
# is not in question here.
#
# Usage: bash freeze_prechange.sh            (copies; no compute)

set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
T="$(cd "$HERE/../.." && pwd)"          # tests/MMWSpalling
DEST="$HERE/frozen_prechange"

mkdir -p "$DEST"
copy () {   # copy <relative path under tests/MMWSpalling>
    local rel="$1"
    if [ ! -e "$T/$rel" ]; then echo "  MISSING  $rel"; return; fi
    mkdir -p "$DEST/$(dirname "$rel")"
    cp -R "$T/$rel" "$DEST/$(dirname "$rel")/" 2>/dev/null \
        && echo "  frozen   $rel" || echo "  FAILED   $rel"
}

echo "=== the 13 strip-derived jet_enthalpy tags (the legs (b1) authorises) ==="
for t in B05 B1 B1e B1n B2 B3 B4 E Eon G Gbug P1 P3; do
    copy "unit/jet_enthalpy/output/$t"
done

echo "=== QF -- a QUARTER leg that moved under rev 6 and must now read zero ==="
copy "unit/jet_enthalpy/output/QF"

echo "=== the three other areas rev 6 moved, plus the M-2000 control ==="
copy "studies/s1_surface_resolution/output_smoke/M-100_dz1"
copy "studies/s1_surface_resolution/output_smoke/M-2000_dz1"
copy "unit/beam_void_closure/output"
copy "validation/meier/sp_meier_pilot/output/dev2d"

# the binary that produced them, so the provenance is not an argument
if [ -f "$T/../../bin/mmwspalling-3d-g++" ]; then
    shasum -a 256 "$T/../../bin/mmwspalling-3d-g++" > "$DEST/BINARY_A_absent.sha256"
    echo "=== A-absent binary $(cut -c1-24 "$DEST/BINARY_A_absent.sha256") ==="
fi

echo
echo "frozen under: $DEST"
du -sh "$DEST" 2>/dev/null
