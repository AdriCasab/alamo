#!/bin/bash
# P1g-2c deliverable 7 / (d1): re-run P1g-2b's CTRL4 on the FIXED binary.
#
# ⚠⚠ THE ARGV IS P1g-2b's CTRL4 VERBATIM. It is read out of that leg's own
# recorded `CTRL4.argv.sh` -- which itself carried Leg B's argv verbatim -- and
# only the `plot_file` is re-pointed. Retyping it would be a silent chance to
# differ from the leg whose 8.94e-3 ledger failure this step exists to close, and
# (d2)'s bore comparison is against that same leg's sibling CTRL2.
#
# ⚠ The recorded file also carries a trailing `sha  path` line; it is stripped.
#
# ⚠ caffeinate: not checkpointed; a host sleep loses the whole 36-minute leg.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo
S=$ROOT/tests/MMWSpalling/studies/p1g2c_shiftfix
SRC=$ROOT/tests/MMWSpalling/studies/p1g2b_subcell/CTRL4.argv.sh
BIN=$ROOT/bin/mmwspalling-3d-g++
NAME=CTRL4_FIXED

[ -f "$SRC" ] || { echo "FATAL: no archived argv at $SRC"; exit 1; }
[ -x "$BIN" ] || { echo "FATAL: no binary at $BIN"; exit 1; }
mkdir -p "$S/$NAME"

CMD=$(python3 - "$SRC" "$S/$NAME" <<'PY'
import re, sys
src, out = sys.argv[1], sys.argv[2]
lines = [ln for ln in open(src).read().split("\n")
         if ln.strip() and not re.match(r"^[0-9a-f]{64}\s+\S+$", ln.strip())]
cmd = " ".join(lines).strip()
# re-point ONLY the plot_file; every other token is carried through untouched
cmd = re.sub(r"plot_file=\S+", "plot_file=%s/FH" % out, cmd)
sys.stdout.write(cmd)
PY
)
[ -n "$CMD" ] || { echo "FATAL: could not rebuild the command"; exit 1; }
case "$CMD" in
  *spall.lateral_thermal=1*) ;;
  *) echo "FATAL: the carried argv lost spall.lateral_thermal=1"; exit 1 ;;
esac

shasum -a 256 "$BIN" > "$S/$NAME.sha_before"
printf '%s\n' "$CMD" > "$S/$NAME.argv.sh"
cat "$S/$NAME.sha_before" >> "$S/$NAME.argv.sh"

echo "=== $NAME start $(date +%T)"
eval "$CMD" > "$S/$NAME.log" 2>&1
RC=$?
shasum -a 256 "$BIN" > "$S/$NAME.sha_after"
echo "=== $NAME rc=$RC  $(date +%T)"
if ! cmp -s "$S/$NAME.sha_before" "$S/$NAME.sha_after"; then
  echo "GUARD FAILURE: the binary CHANGED during $NAME -- the leg is VOID."
  exit 2
fi
echo "=== $NAME binary unchanged across the leg"
exit $RC
