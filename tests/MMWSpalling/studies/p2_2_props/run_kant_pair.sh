#!/bin/bash
# P2-2a (c1) -- THE SHORT KANT PAIR, ON THE SHIPPING BINARY.
#
# ⚠⚠ WHY THIS SCRIPT EXISTS. The rev-7 Commands block invoked it and it did not
# exist (review row 10), and the (c1) evidence in the repository was from a
# 10/07 16:12:50 build while the shipping binary was 10/08 10:53:57 (review
# must-fix 1). The physics argument -- the Kant IC equals the typed `T_ref`, so
# the stored field IS the old constant there -- is sound but it is an argument,
# and (c1) is a gate.
#
# ⚠⚠ SO THIS SCRIPT ASSERTS ITS OWN PROVENANCE. Every leg's
# `metadata:Compilation_Time` is compared against the binary actually invoked.
# That is the check whose absence let stale output be reported as a result for a
# whole revision: a directory full of plotfiles looks like evidence whatever
# produced it.
#
# ⚠ `P1_KEYS` ARE NOT PASSED, and that is the point. Until P2-2a the harness
# passed `spall.driving_stress=differential` plus a flaw refit x10 to 80 um,
# sized so the differential drive's 0.330 and K_I ~ sqrt(a) cancel. `differential`
# is retired, so the refit lost its partner; these legs run the SHIPPED flaw
# under the one surviving form.
#
# THE BAR, and which numbers are NOT it:
#   463.3 / 465.8 K onset SUPERHEAT on R150 / R500 at 2 mm -- the deck's own
#   recorded ABSOLUTE onset with the shipped flaw.
#   ⚠ NOT 754.57 K: that is the WITHDRAWN `differential` + 80 um figure, which
#     rev 7's bar text still names (review must-fix 2).
#   ⚠ NOT 520.68 K: that is `absolute` + the 80 um refit, i.e. the broken
#     intermediate state with the keys unconsumed. The review proposes it as the
#     corrected bar; it is not the bar for the configuration that now ships.
#   ⚠ PROVENANCE CAVEAT on the pair: `463.3` is independently recorded (an
#     S150-type leg reproduced it to -1.5 K in the Oct-2 P1f archive), but the
#     `465.8` in this deck's comment traces to an earlier packet's `R500`, not
#     to the same leg family. The pairing is therefore less clean than it looks,
#     even though the shipping binary reproduces BOTH to 0.01 K.
#
# Usage: bash run_kant_pair.sh            (4 ranks, ~2 min, no Meier leg)

set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../../.." && pwd)"
DECK="$ROOT/tests/MMWSpalling/validation/kant/sp_kant_onset_meiercfg/input_model"
BIN="$ROOT/bin/mmwspalling-3d-g++"
OUT="$HERE/C1_KANT"
LOG="$HERE/C1_KANT.txt"

[ -x "$BIN" ] || { echo "no binary at $BIN"; exit 1; }
BSHA="$(shasum -a 256 "$BIN" | cut -d' ' -f1)"

# runs are not checkpointed; hold the host awake even for a short pair
CAF=""
command -v caffeinate >/dev/null && CAF="caffeinate -dis"

rm -rf "$OUT"; mkdir -p "$OUT"
{
  echo "=============================================================================="
  echo "P2-2a (c1) -- KANT ONSET ON THE SHIPPING BINARY"
  echo "=============================================================================="
  echo "  binary        : ${BSHA:0:24}"
  echo "  deck          : validation/kant/sp_kant_onset_meiercfg/input_model"
  echo "  P1_KEYS       : NOT PASSED (the differential-matched refit is withdrawn)"
  echo "  bar           : 463.3 / 465.8 K superheat on R150 / R500 at 2 mm"
  echo
  printf "  %-12s %-10s %-12s %-12s %-22s\n" leg mesh "dT [K]" "T_surf [K]" "built"
  echo "  ----------------------------------------------------------------------------"
} > "$LOG"

fail=0
for mesh in 12 24; do
  [ "$mesh" = 12 ] && mm=2mm || mm=1mm
  for h in 150 500; do
    tag="R${h}_${mm}"
    $CAF mpirun --oversubscribe --bind-to none -np 4 "$BIN" "$DECK" \
         "amr.n_cell=$mesh $mesh $mesh" "amr.max_grid_size=$mesh" \
         "surface_patch.h_conv=$h" "spall.surface_T=reconstructed" \
         "stop_time=45.0" "plot_file=$OUT/$tag/run" > "$OUT/$tag.log" 2>&1
    rc=$?
    # ⚠ THE PROVENANCE ASSERTION: the leg must have been produced by the binary
    # whose sha is recorded above. This is the check whose absence let a 10/07
    # build be reported against a 10/08 binary for a whole revision.
    built="$(grep -h 'Compilation_Time' "$OUT/$tag"/run*/metadata 2>/dev/null \
             | head -1 | sed 's/.*= *//' | tr -d '\r')"
    bdate="$(grep -h 'Compilation_Date' "$OUT/$tag"/run*/metadata 2>/dev/null \
             | head -1 | sed 's/.*= *//' | tr -d '\r')"
    read -r dT Ts <<EOF
$(awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i}
           NR>1 && $h["regime"]==1 {printf "%.2f %.3f", $h["T_surf"]-293.15, $h["T_surf"]; exit}' \
   "$OUT/$tag/run_removal_events.csv" 2>/dev/null \
   || echo "NA NA")
EOF
    printf "  %-12s %-10s %-12s %-12s %-22s\n" \
           "R$h" "$mm" "${dT:-NA}" "${Ts:-NA}" "${bdate:-?} ${built:-?}" >> "$LOG"
    [ "$rc" = 0 ] || fail=1
  done
done

{
  echo
  echo "  ⚠ the 'built' column is each leg's OWN metadata, so a stale directory" 
  echo "    cannot be reported as a result. Compare it against the sha at the top."
  echo
  echo "  scored at 2 mm, against the deck's recorded absolute onset 463.3 / 465.8:"
} >> "$LOG"

score () {   # score <leg tag> <bar>
    local tag="$1" bar="$2"
    local dT
    dT="$(awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i}
                   NR>1 && $h["regime"]==1 {printf "%.2f", $h["T_surf"]-293.15; exit}' \
          "$OUT/$tag/run_removal_events.csv" 2>/dev/null)"
    local blt
    blt="$(grep -h 'Compilation_Time' "$OUT/$tag"/run*/metadata 2>/dev/null \
           | head -1 | sed 's/.*= *//' | tr -d '\r')"
    awk -v d="$dT" -v b="$bar" -v t="$tag" -v blt="$blt" 'BEGIN{
        diff = d - b; if (diff < 0) diff = -diff;
        printf "    %-12s dT %-9s bar %-9s |diff| %-8.2f %s   built %s\n",
               t, d, b, diff, (diff <= 0.01 ? "PASS" : "** FAIL **"), blt;
        exit (diff <= 0.01 ? 0 : 1) }' >> "$LOG"
}
score R150_2mm 463.3 || fail=1
score R500_2mm 465.8 || fail=1
{
  echo
  echo "  ⚠ the 1 mm rows are the MESH, not a bar: the recorded pair was taken at"
  echo "    2 mm, so scoring them against it would measure the mesh and call it a"
  echo "    regression."
  if [ "$fail" = 0 ]; then
    echo "  (c1) PASS -- onset unchanged to 0.01 K on the shipping binary."
  else
    echo "  (c1) ** FAIL **"
  fi
} >> "$LOG"

cat "$LOG"
exit $fail
