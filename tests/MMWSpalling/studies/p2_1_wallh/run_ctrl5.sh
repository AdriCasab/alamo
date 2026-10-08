#!/bin/bash
# P2-1 deliverables 9 and 10: CTRL5 = CTRL4_FIXED's argv + ONE key.
#
# ⚠⚠ THE ARGV IS CTRL4_FIXED's VERBATIM. It is read out of that leg's own
# recorded `CTRL4_FIXED.argv.sh` -- which itself carried P1g-2b's CTRL4, which
# carried Leg B -- and only the `plot_file` is re-pointed before exactly ONE key
# is appended. Retyping it would be a silent chance to differ from the leg (c2)
# and (d2)/(c3) compare against, which is the whole value of the comparison.
#
# ⚠ DELIVERABLE 9, the provenance lesson P1g-2c paid for: the argv is written
# NEXT TO the binary hash, before the run, and the hash is re-checked after.
# P1g-2c's (c1b) was unscoreable because the `evidence_p1g2` legs recorded their
# binary and their resolved parameters but never their COMMAND, so they could be
# reconstructed but not replayed. A resolved-parameter dump is not an argv.
#
# ⚠ The recorded file also carries a trailing `sha  path` line; it is stripped.
# ⚠ caffeinate: runs are not checkpointed; a host sleep loses the whole leg.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo
S=$ROOT/tests/MMWSpalling/studies/p2_1_wallh
SRC=$ROOT/tests/MMWSpalling/studies/p1g2c_shiftfix/CTRL4_FIXED.argv.sh
BIN=$ROOT/bin/mmwspalling-3d-g++
NAME=CTRL5
THE_ONE_KEY="surface_patch.side_face_h_law=derived"

[ -f "$SRC" ] || { echo "FATAL: no archived argv at $SRC"; exit 1; }
[ -x "$BIN" ] || { echo "FATAL: no binary at $BIN"; exit 1; }
mkdir -p "$S/$NAME"

CMD=$(python3 - "$SRC" "$S/$NAME" "$THE_ONE_KEY" <<'PY'
import re, sys
src, out, key = sys.argv[1], sys.argv[2], sys.argv[3]
lines = [ln for ln in open(src).read().split("\n")
         if ln.strip() and not re.match(r"^[0-9a-f]{64}\s+\S+$", ln.strip())]
cmd = " ".join(lines).strip()
# re-point ONLY the plot_file; every other token is carried through untouched
cmd = re.sub(r"plot_file=\S+", "plot_file=%s/FH" % out, cmd)
sys.stdout.write(cmd + " " + key)
PY
)
[ -n "$CMD" ] || { echo "FATAL: could not rebuild the command"; exit 1; }

# ⚠ ONE key, asserted rather than trusted: the packet's scope is one key and a
# second would make (c2) unattributable, which is the mistake D2v made and D2w
# had to repair. Count the tokens that differ from the carried argv.
NDIFF=$(python3 - "$SRC" "$CMD" <<'PY'
import re, sys
src = [ln for ln in open(sys.argv[1]).read().split("\n")
       if ln.strip() and not re.match(r"^[0-9a-f]{64}\s+\S+$", ln.strip())]
a = set(" ".join(src).split()); b = set(sys.argv[2].split())
extra = [t for t in (b - a) if not t.startswith("plot_file=")]
print(len(extra)); sys.stderr.write("added: %s\n" % " ".join(sorted(extra)))
PY
)
echo "=== $NAME adds $NDIFF key(s) to the carried argv"
[ "$NDIFF" = "1" ] || { echo "FATAL: CTRL5 must add exactly ONE key, not $NDIFF"; exit 1; }
case "$CMD" in
  *"$THE_ONE_KEY"*) ;;
  *) echo "FATAL: the one key is missing"; exit 1 ;;
esac
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
