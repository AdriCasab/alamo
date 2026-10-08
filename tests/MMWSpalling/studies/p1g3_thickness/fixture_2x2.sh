#!/bin/bash
# P1g-3 closure, ord 0 (ii)+(p2): the four 2x2 ledger legs, IN-REPO, on the
# SHIPPING binary.
#
#   usage: fixture_2x2.sh            (runs all four, then scores)
#          fixture_2x2.sh --score    (scores what is already on disk)
#
# ⚠⚠ WHY THIS SCRIPT EXISTS AT ALL. P1g-3's (b2) 2x2 was run in the session
# scratchpad, which is session-scoped and is not the repository. Rule 8 of this
# programme: a number is backed only if BOTH the file it was read from AND the
# code that reduced it are in the repository. Four legs in /tmp back nothing.
#
# ⚠⚠ WHY THE FIXTURE IS p0_steady_column AND NOT unit/lateral_drive. P1g-3's
# first (b2) read 8.9348e-16 on ALL FOUR legs and looked like a clean
# conservation pass. It was a dead comparand: `unit/lateral_drive` has ZERO
# regime-1 (vertical) events, so the four legs were four identical runs being
# compared with each other. The flake-thickness keys only act on a vertical
# firing, so the fixture has to HAVE vertical firings. p0_steady_column does.
# => (p2)'s bar has two halves and the second one is the real one: the four
# values must be <= 5e-3 AND they must DIFFER.
#
# ⚠ p0_steady_column exits rc=1 BY DESIGN -- it is a standing red carried in
# ACTIVE_STEP.md's guardrails. rc=1 is therefore NOT a failure here and the
# script records it rather than aborting. What would be a failure is a missing
# or short thermo.dat, which is checked explicitly.
#
# ⚠⚠ THREE keys are added, not two, and the third is the one that bites.
# `amr.thermo.int=1` sets how often the ledger is COMPUTED; `energy_ledger.enabled=1`
# turns it on; but `amr.thermo.plot_int` sets how often it is WRITTEN and it
# DEFAULTS TO -1 = never (src/Integrator/Integrator.cpp:126, :1149). Without it the
# leg runs all 7500 steps, exits rc=0, writes every plotfile -- and produces NO
# thermo.dat at all. That failure looks exactly like a clean run. All three are
# identical on all four legs, so none of them can carry the signal.
# ⚠ plot_int=1 (not 10) because the bar is a MAX over steps: sampling would
# understate the worst residual.
set -u
ROOT=/Users/tzetze20/amr_tools/alamo
S=$ROOT/tests/MMWSpalling/studies/p1g3_thickness
INP=$ROOT/tests/MMWSpalling/validation/p0_steady_column/input
BIN=$ROOT/bin/mmwspalling-3d-g++
OUT=$S/ledger_2x2

LEGS="V_OFF V_LEN V_FRAC V_BOTH"
keys_for() {
  case "$1" in
    V_OFF)  echo "" ;;
    V_LEN)  echo "spall.flake_length=flux" ;;
    V_FRAC) echo "spall.flake_fraction=weibull" ;;
    V_BOTH) echo "spall.flake_length=flux spall.flake_fraction=weibull" ;;
  esac
}

[ -x "$BIN" ] || { echo "FATAL: no binary at $BIN"; exit 1; }
SHA=$(shasum -a 256 "$BIN" | awk '{print $1}')

if [ "${1:-}" != "--score" ]; then
  mkdir -p "$OUT"
  cat > "$OUT/PROVENANCE.txt" <<EOF
P1g-3 (b2)/(p2) -- the 2x2 energy-ledger legs for the flake-thickness keys.

binary   : $SHA
          (the SHIPPING binary; P1g-3's own (b2) ran on a180a191..., which was
           superseded 16:31 that day when the LATERAL thickness site was wired.
           ⚠ The toolchain is NOT byte-reproducible, so the earlier figures
           cannot be carried across the rebuild -- hence this re-run.)
date     : $(date '+%F %T')
fixture  : $INP
added    : amr.thermo.int=1 amr.thermo.plot_int=1 energy_ledger.enabled=1
           (all three on all four legs; plot_int is what WRITES thermo.dat
            and defaults to never -- without it a leg looks clean and
            produces no ledger)
2x2      : spall.flake_length = twocell|flux  x  spall.flake_fraction = cf|weibull
legs     : V_OFF (neither)  V_LEN (length)  V_FRAC (fraction)  V_BOTH (both)
exit code: p0_steady_column returns rc=1 by design -- a STANDING RED carried in
           ACTIVE_STEP.md. Recorded per leg below, not treated as failure.

⚠⚠ This fixture replaces unit/lateral_drive, on which the same 2x2 was VACUOUS:
   that input has 0 regime-1 events, so all four legs returned an identical
   8.9348e-16 and the row proved nothing. p0_steady_column fires vertically.
EOF
  for L in $LEGS; do
    K=$(keys_for "$L")
    rm -rf "$OUT/$L"
    CMD="caffeinate -dis mpirun --oversubscribe --bind-to none -np 4 $BIN $INP \
amr.thermo.int=1 amr.thermo.plot_int=1 energy_ledger.enabled=1 $K plot_file=$OUT/$L/FH"
    printf '%s\n' "$CMD" > "$OUT/$L.argv.sh"
    printf '%s  %s\n' "$SHA" "$BIN" >> "$OUT/$L.argv.sh"
    echo "=== $L start $(date +%T)   keys: ${K:-<none>}"
    eval "$CMD" > "$OUT/$L.log" 2>&1
    RC=$?
    echo "=== $L rc=$RC  $(date +%T)   (rc=1 is the standing red)"
    echo "$L rc=$RC" >> "$OUT/PROVENANCE.txt"
    T="$OUT/$L/FH/thermo.dat"
    if [ ! -s "$T" ]; then
      echo "FATAL: $L produced no thermo.dat -- this IS a failure, unlike rc=1"
      exit 1
    fi
  done
fi

python3 - "$OUT" <<'PY'
import os, statistics as st, sys
OUT = sys.argv[1]
LEGS = ("V_OFF", "V_LEN", "V_FRAC", "V_BOTH")
BAR = 5.0e-3
print()
print("=" * 78)
print("P1g-3 (p2)  the 2x2 energy ledger, in-repo, on the shipping binary")
print("=" * 78)
res = {}
for L in LEGS:
    t = os.path.join(OUT, L, "FH/thermo.dat")
    if not os.path.exists(t):
        print("  %-8s MISSING" % L); continue
    lines = open(t).read().split("\n")
    hdr = lines[0].split()
    ix = {n: i for i, n in enumerate(hdr)}
    need = ("ledger_err_solid", "ledger_err", "ledger_E_flake")
    miss = [n for n in need if n not in ix]
    if miss:
        print("  %-8s no column(s): %s" % (L, miss)); continue
    es, er, ef, n = [], [], [], 0
    for ln in lines[1:]:
        f = ln.split()
        if len(f) < len(hdr):
            continue
        n += 1
        es.append(abs(float(f[ix["ledger_err_solid"]])))
        er.append(abs(float(f[ix["ledger_err"]])))
        ef.append(float(f[ix["ledger_E_flake"]]))
    res[L] = (max(es), max(er), ef[-1] if ef else float("nan"), n)
print("  %-8s %16s %16s %16s %7s" % ("leg", "max|err_solid|", "max|err|",
                                     "E_flake end [J]", "rows"))
for L in LEGS:
    if L in res:
        a, b, c, n = res[L]
        print("  %-8s %16.6e %16.6e %16.6e %7d" % (L, a, b, c, n))
print()
if len(res) < 4:
    print("  (p2) FAIL -- fewer than four legs scored")
    sys.exit(1)
vals = [res[L][0] for L in LEGS]
worst = max(vals)
print("  half 1 -- CONSERVATION: max over all four = %.6e   [bar <= %.1e]"
      % (worst, BAR))
ok1 = worst <= BAR
print("           %s" % ("PASS" if ok1 else "FAIL"))
uniq = len({"%.12e" % v for v in vals})
print("  half 2 -- NON-VACUITY: %d distinct values of ledger_err_solid among 4"
      % uniq)
# the flake energy is the quantity the keys actually move; show it too
fl = [res[L][2] for L in LEGS]
uniqf = len({"%.12e" % v for v in fl})
print("           %d distinct values of ledger_E_flake among 4" % uniqf)
ok2 = uniq == 4
print("           %s" % ("PASS" if ok2 else
                         "FAIL -- identical values mean a DEAD COMPARAND, "
                         "not conservation"))
if ok2:
    print("           spread %.3e to %.3e (%.1fx)"
          % (min(vals), max(vals), max(vals) / min(vals) if min(vals) else 0))
print()
print("  (p2) %s" % ("PASS" if (ok1 and ok2) else "FAIL"))
sys.exit(0 if (ok1 and ok2) else 1)
PY
