#!/bin/sh
# P2-2b (b5): bound the cost of the two hot paths -- the crack integral and
# the H->T inversion.
#
# ⚠⚠ WHAT THIS MEASURES AND WHAT IT CANNOT. There is NO per-kernel timer in
# either source file (zero BL_PROFILE scopes), so a per-path "share of
# runtime" is not obtainable without instrumenting the binary -- which would
# change the binary that ships and force a second tree gate for the sake of a
# report row. What is measured instead is the TOTAL wall time of a leg whose
# cost is almost entirely those two paths, with the law OFF and ON.
#
# ⚠ OFF vs ON ON THE SAME BINARY, not before-binary vs after-binary. The
# toolchain is not byte-reproducible, so a two-binary comparison would
# confound the law's cost with compiler variation. Under `law = constant` the
# arithmetic is the pre-P2-2b arithmetic bit for bit ((b1)/(b7)), so the OFF
# leg IS the "before" measurement, taken with the shipping compiler.
#
# THE DECK: unit/scalar_flaw -- 1 x 1 x 125 cells, 20,000 steps, spall on.
# A single column, so there is almost no cost outside the depth scan and the
# enthalpy inversion; a 3-D deck would dilute both behind the conduction
# sweep and measure something else.
set -e
ROOT=$(cd "$(dirname "$0")/../../../.." && pwd)
BIN=$ROOT/bin/mmwspalling-3d-g++
IN=$ROOT/tests/MMWSpalling/unit/scalar_flaw/input
TMP=${TMPDIR:-/tmp}/b5_$$
mkdir -p "$TMP"
echo "P2-2b (b5) THE COST OF THE TWO HOT PATHS"
echo "binary $(shasum -a 256 "$BIN" | cut -c1-24)"
echo "deck   unit/scalar_flaw  1x1x125, 20000 steps, spall.enabled = 1"
echo
printf "%-10s %8s %8s %8s   %s\n" leg rep1 rep2 rep3 mean
for leg in off on; do
  if [ "$leg" = on ]; then EXTRA=microstructure.phase0.law=grimsel; else EXTRA=; fi
  tot=0; out=""
  for r in 1 2 3; do
    s=$(date +%s.%N)
    ( cd "$ROOT" && mpirun --oversubscribe --bind-to none -np 1 "$BIN" "$IN" \
        $EXTRA plot_file="$TMP/${leg}_$r/run" >"$TMP/${leg}_$r.log" 2>&1 )
    e=$(date +%s.%N)
    d=$(echo "$e - $s" | bc)
    out="$out $(printf '%8.2f' "$d")"
    tot=$(echo "$tot + $d" | bc)
  done
  printf "%-10s%s   %8.2f s\n" "$leg" "$out" "$(echo "$tot / 3" | bc -l)"
done
rm -rf "$TMP"
