#!/bin/bash
# D2q-2c pre-flight 3: run post_all.sh with the binary pinned.
#
# WHY THIS EXISTS. At the D2q-2b checkpoint the implementer rebuilt while an
# identity sweep was running. The sweep's later cases then used a DIFFERENT
# binary from its earlier ones, so the 2639-file comparison was meaningless --
# and nothing in the tooling noticed. The result was caught only because the
# implementer happened to check `pgrep` afterwards.
#
# Byte-identity is the only oracle this line has. An oracle that can be
# silently corrupted by a background `make` is not an oracle. So: record the
# binary's fingerprint before the sweep, verify it after, and FAIL LOUDLY if it
# moved. Guardrails enforced by memory have already failed once.
#
#   usage: guarded_post_all.sh <label>
#
# Exit codes: 0 sweep ran and the binary never moved
#             1 usage / missing binary
#             2 THE BINARY CHANGED MID-SWEEP -- result is invalid, discard it
set -u
ROOT=/Users/tzetze20/amr_tools/alamo
W="$ROOT/tests/MMWSpalling/studies/tree_sweep_0922/witness"
BIN="$ROOT/bin/mmwspalling-3d-g++"
LABEL="${1:?usage: guarded_post_all.sh <label>}"

cd "$ROOT" || exit 1
[ -x "$BIN" ] || { echo "GUARD: no binary at $BIN"; exit 1; }

# mtime+size+hash: mtime alone can be preserved by a toolchain, a hash alone
# costs a full read but is definitive. Take both.
fingerprint() { stat -f "%m %z" "$BIN"; shasum -a 256 "$BIN" | awk '{print $1}'; }

BEFORE=$(fingerprint)
echo "GUARD: binary pinned at $(echo "$BEFORE" | tail -1 | cut -c1-16)…  label=$LABEL"
echo "GUARD: do NOT run 'make' until this sweep finishes."

bash "$W/post_all.sh" "$LABEL"
RC=$?

AFTER=$(fingerprint)
if [ "$BEFORE" != "$AFTER" ]; then
    echo "================================================================"
    echo "GUARD FAILURE: the binary changed WHILE the sweep was running."
    echo "  before: $(echo "$BEFORE" | tr '\n' ' ')"
    echo "  after : $(echo "$AFTER"  | tr '\n' ' ')"
    echo "The $LABEL sweep mixed two binaries and its identity result is"
    echo "MEANINGLESS. Discard $W/$LABEL.md5 and $W/${LABEL}_wit, then re-run"
    echo "on a settled binary. Do not report this sweep."
    echo "================================================================"
    exit 2
fi
echo "GUARD: binary unchanged across the sweep; $LABEL identity result is valid."
exit $RC
