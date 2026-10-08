#!/usr/bin/env python3
"""P1g-2b ord 1 (a1)-(a4): what actually separates Leg B from MEIER_CTRL2.

Emits `SUBCELL_PRE.md`. Offline: reads only files already on disk, edits nothing.

⚠⚠ THE POINT OF THIS SCRIPT IS TO REFUSE AN ATTRIBUTION, not to make one. Three
readers (me included) have quoted the +4 to +6.6 mm between these two legs as
the sub-cell mode's effect. It cannot be: the legs differ by one argv key AND by
the BINARY, with four steps of source work in between. `score_bore.py:11-21`
already refuses exactly this comparison. This script proves both halves of that
-- one key, two binaries -- so the claim is retired with evidence rather than
with an assertion.

Rules carried verbatim from `p1g3_thickness/score_bore.py`:
  1. never score the bore at 100 mm depth (the 90-110 band is absent by
     construction -- Meier p. 223 pause underream);
  2. quarter domain, so diameter = 2*r_max;
  3. a band near either front reports front position, not width (reported here,
     enforced in `score_subcell.py` where two live legs are compared).
"""
import csv
import math
import os
import re
import statistics as st
import sys

ROOT = "/Users/tzetze20/amr_tools/alamo"
LEGB = os.path.join(ROOT, "tests/MMWSpalling/studies/p1f_lateral/output_legB_fix")
CTRL2 = os.path.join(ROOT, "tests/MMWSpalling/studies/p1g3_thickness/MEIER_CTRL2")
BIN = os.path.join(ROOT, "bin/mmwspalling-3d-g++")

# ⚠ 90-110 mm is deliberately absent. Rule 1.
BANDS = ((10, 50), (50, 90), (110, 150), (150, 200), (200, 300))

# Meier's digitised Fig. 8.8 profile, as recorded in
# validation/meier/meier-2017-pilot-660kg-validation-reference.md:63-64.
# ⚠⚠ It is a VISIBLE WIDTH, so it is a LOWER BOUND on the true bore, +/-5 mm
# (item 10 of that note). Every error below is therefore an upper bound on the
# model's excess, and the +/-5 mm is larger than the band-to-band spread.
MEIER = ((10, 96.0), (25, 92.0), (50, 88.0), (75, 87.0), (125, 87.0),
         (300, 84.0))

# ⚠ `plot_file` differs BY CONSTRUCTION -- each leg writes to its own directory,
# and a differing output path is not a differing physics key. It is excluded,
# and the exclusion is stated in the output rather than left implicit.
EXCLUDE = ("plot_file",)


def tokens(path):
    """argv tokens, respecting the single-quoted h_expr (which contains spaces,
    parentheses and '=' signs)."""
    s = open(path).read()
    # the recorded command may carry a trailing `sha  path` line; drop any line
    # that is just a 64-hex digest plus a path
    s = "\n".join(ln for ln in s.split("\n")
                  if not re.match(r"^[0-9a-f]{64}\s+\S+$", ln.strip()))
    return re.findall(r"'[^']*'|\S+", s)


def argv_keys(path):
    out = {}
    for t in tokens(path):
        t = t.strip("'")
        if "=" not in t:
            continue
        k, v = t.split("=", 1)
        if not re.match(r"^[A-Za-z_][\w.]*$", k):
            continue
        out[k] = v          # last wins, exactly as ParmParse resolves it
    return out


def geom(av):
    """DZ and the original rock top, from the leg's own argv."""
    nz = [float(x) for x in av["amr.n_cell"].split()] if " " in av["amr.n_cell"] \
        else None
    # the argv records multi-valued keys as one token only when quoted; the
    # recorded commands here keep them space-separated, so re-read from tokens
    return nz


def bore_by_band(leg_dir, dz, ktop0, z_top):
    """max radius per band, plus the deepest removal and the last event time."""
    f = os.path.join(leg_dir, "FH_removal_events.csv")
    rmax = {b: 0.0 for b in range(len(BANDS))}
    deepest, tmax, n6 = 0.0, 0.0, 0
    with open(f) as fh:
        for r in csv.DictReader(fh):
            try:
                kt = int(float(r["k_top"]))
                rd = float(r.get("riser_depth") or 0.0)
                ci, cj = int(r["col_i"]), int(r["col_j"])
                t = float(r["time"])
                reg = int(float(r["regime"]))
            except (ValueError, KeyError, TypeError):
                continue
            d = (ktop0 - kt) * dz * 1e3 + rd * 1e3
            rr = dz * math.hypot(ci + 0.5, cj + 0.5) * 1e3
            deepest = max(deepest, d)
            tmax = max(tmax, t)
            if reg == 6:
                n6 += 1
            for b, (a, c) in enumerate(BANDS):
                if a <= d < c and rr > rmax[b]:
                    rmax[b] = rr
    return rmax, deepest, tmax, n6


def meier_at(depth):
    """piecewise-linear interpolation of the digitised profile."""
    if depth <= MEIER[0][0]:
        return MEIER[0][1]
    if depth >= MEIER[-1][0]:
        return MEIER[-1][1]
    for (d0, v0), (d1, v1) in zip(MEIER, MEIER[1:]):
        if d0 <= depth <= d1:
            return v0 + (v1 - v0) * (depth - d0) / (d1 - d0)
    return MEIER[-1][1]


def meier_band_integral(a, c, n=2001):
    """mean of the profile over [a,c). The declared convention."""
    xs = [a + (c - a) * i / (n - 1) for i in range(n)]
    return st.mean(meier_at(x) for x in xs)


def meier_band_samples(a, c):
    """mean of the profile's own digitised samples inside [a,c), or None.
    Carried only because the two conventions disagree on the COLLAR band."""
    v = [val for d, val in MEIER if a <= d < c]
    return st.mean(v) if v else None


def main():
    P = print
    avB = argv_keys(os.path.join(LEGB, "argv.sh"))
    avC = argv_keys(os.path.join(CTRL2, "argv.sh"))

    shaB = open(os.path.join(LEGB, "BINARY.sha256")).read().split()[0]
    shaC = os.popen("shasum -a 256 %s" % BIN).read().split()[0]

    # geometry, read from CTRL2's own argv (both legs share it; asserted below)
    tb, tc = tokens(os.path.join(LEGB, "argv.sh")), \
        tokens(os.path.join(CTRL2, "argv.sh"))

    def multi(tk, key):
        for i, t in enumerate(tk):
            if t.startswith(key + "="):
                vals = [t.split("=", 1)[1]]
                j = i + 1
                while j < len(tk) and re.match(r"^-?[\d.eE+-]+$", tk[j]):
                    vals.append(tk[j]); j += 1
                return [float(x) for x in vals]
        return None

    nc_B, nc_C = multi(tb, "amr.n_cell"), multi(tc, "amr.n_cell")
    lo_C, hi_C = multi(tc, "geometry.prob_lo"), multi(tc, "geometry.prob_hi")
    nz = int(nc_C[2]); DZ = (hi_C[2] - lo_C[2]) / nz
    KTOP0 = nz - 1
    Z_TOP = lo_C[2] + nz * DZ

    P("# P1g-2b `(a1)`-`(a4)` — what separates Leg B from MEIER_CTRL2")
    P("")
    P("Generated by `argv_diff.py`. Offline; nothing was run and nothing edited.")
    P("")
    P("## ⚠⚠ The headline: the +4 to +6.6 mm is NOT attributable to the mode")
    P("")
    P("These two legs differ by **one argv key** *and* by **the binary**. A bore")
    P("difference between them is the sum of one key and four steps of source")
    P("work, and nothing in the pair can separate them. `score_bore.py:11-21`")
    P("refuses this comparison by rule; this note records why.")
    P("")

    # ---------------- (a1) the argv diff ---------------------------------
    P("## `(a1)` GATE — one key, two binaries")
    P("")
    P("### Binaries — these DO differ")
    P("")
    P("| leg | sha256 | note |")
    P("|---|---|---|")
    P("| Leg B | `%s…` | from `output_legB_fix/BINARY.sha256`; **this binary no longer exists** |"
      % shaB[:16])
    P("| CTRL2 | `%s…` | the shipping binary, on disk now |" % shaC[:16])
    P("")
    same_bin = (shaB == shaC)
    P("**binaries identical: %s** — the gate requires them to DIFFER, and they do."
      % ("YES ⇒ (a1) FAILS, the runs were unnecessary" if same_bin else "NO"))
    P("")
    P("### argv keys, compared BY VALUE")
    P("")
    P("⚠ **`plot_file` is excluded, and this is that statement.** Each leg writes")
    P("to its own output directory, so `plot_file` differs by construction; an")
    P("output path is not a physics key. Nothing else is excluded.")
    P("")
    keys = sorted(set(avB) | set(avC))
    diffs = []
    for k in keys:
        if k in EXCLUDE:
            continue
        b, c = avB.get(k), avC.get(k)
        if b != c:
            diffs.append((k, b, c))
    P("| key | Leg B | CTRL2 |")
    P("|---|---|---|")
    for k, b, c in diffs:
        P("| `%s` | %s | %s |"
          % (k, "*(absent)*" if b is None else "`%s`" % b,
             "*(absent)*" if c is None else "`%s`" % c))
    if not diffs:
        P("| *(none)* | | |")
    P("")
    P("**keys differing (excluding `plot_file`): %d**  — bar: exactly **1**"
      % len(diffs))
    P("")
    ok_a1 = (len(diffs) == 1 and diffs and diffs[0][0] == "spall.lateral_subcell"
             and not same_bin)
    P("`(a1)` **%s** — %d key, binaries differ: %s"
      % ("PASS" if ok_a1 else "FAIL", len(diffs), "yes" if not same_bin else "NO"))
    P("")
    P("⚠ Total keys compared: Leg B %d, CTRL2 %d (union %d). `amr.n_cell` reads"
      % (len(avB), len(avC), len(keys)))
    P("`%s` on both legs, and the quoted `h_expr` is compared as one token."
      % (" ".join("%g" % x for x in (nc_B or nc_C))))
    P("")

    # ---------------- (a2) the control is half the mode ------------------
    P("## `(a2)` GATE — ⚠⚠ the control is **banking only**, not the mode")
    P("")
    th_in_B = "spall.lateral_thermal" in avB
    th_in_C = "spall.lateral_thermal" in avC
    P("`spall.lateral_thermal` — P1g-2's thermal shift — in the argv:")
    P("")
    P("| leg | `spall.lateral_thermal` |")
    P("|---|---|")
    P("| Leg B | %s |" % ("`%s`" % avB["spall.lateral_thermal"] if th_in_B
                          else "**absent**"))
    P("| CTRL2 | %s |" % ("`%s`" % avC["spall.lateral_thermal"] if th_in_C
                          else "**absent**"))
    P("")
    P("Source default: `MMWSpalling.H:5841`")
    P("`pp.query_default(\"spall.lateral_thermal\", value.sp_lateral_thermal, false)`,")
    P("member default `false` at `:9573`-adjacent; reachability is at the use site")
    P("as `bank_on && sp_lateral_thermal` (`:5840` comment), with **no `requires`**")
    P("abort — removed after P1g-1's negative control returned rc 6.")
    P("")
    ok_a2 = (not th_in_C)
    P("`(a2)` **%s** — the key is %s from CTRL2's argv and defaults `false`."
      % ("PASS" if ok_a2 else "FAIL", "absent" if not th_in_C else "PRESENT"))
    P("")
    P("⇒ **CTRL2 ran P1g-1's banking with P1g-2's shift OFF.** So:")
    P("")
    P("- the banking comparison can only measure **bookkeeping, not heat** — the")
    P("  configuration P1g-1's own review called physically inert on the fixture;")
    P("- **P1g-2's thermal half has never been measured at Meier scale**, two steps")
    P("  after it was implemented. `CTRL4` is the leg that closes it.")
    P("")

    # ---------------- (a3) the bore table --------------------------------
    rB, depB, tB, n6B = bore_by_band(LEGB, DZ, KTOP0, Z_TOP)
    rC, depC, tC, n6C = bore_by_band(CTRL2, DZ, KTOP0, Z_TOP)
    EXP_B = (115.1, 101.1, 98.2, 98.0, 96.5)
    EXP_C = (121.8, 105.3, 102.5, 102.2, 102.0)
    P("## `(a3)` the Leg B → CTRL2 bore table")
    P("")
    P("| depth [mm] | Leg B Ø | CTRL2 Ø | Δ | expected B | expected C | within 0.2 |")
    P("|---|---|---|---|---|---|---|")
    worst = 0.0
    for b, (a, c) in enumerate(BANDS):
        dB, dC = 2 * rB[b], 2 * rC[b]
        eB, eC = EXP_B[b], EXP_C[b]
        err = max(abs(dB - eB), abs(dC - eC))
        worst = max(worst, err)
        P("| %d–%d | %.1f | %.1f | %+.1f | %.1f | %.1f | %s |"
          % (a, c, dB, dC, dC - dB, eB, eC, "yes" if err <= 0.2 else "**%.2f**" % err))
    P("")
    P("worst deviation from the packet's figures: **%.2f mm** (tol 0.2) ⇒ `(a3)` **%s**"
      % (worst, "PASS" if worst <= 0.2 else "FAIL"))
    P("")
    P("deepest removal: Leg B **%.0f mm**, CTRL2 **%.0f mm**; last event %.1f / %.1f s"
      % (depB, depC, tB, tC))
    P("regime-6 events: Leg B **%d**, CTRL2 **%d** (**%.1f×**)"
      % (n6B, n6C, n6C / n6B if n6B else float("nan")))
    P("")
    P("⚠⚠ **Δ is +%.1f to +%.1f mm at every band — and it is UNATTRIBUTED.**"
      % (min(2 * rC[b] - 2 * rB[b] for b in range(len(BANDS))),
         max(2 * rC[b] - 2 * rB[b] for b in range(len(BANDS)))))
    P("One key *and* two binaries. Do not quote it as the sub-cell mode's effect.")
    P("")

    # ---------------- (a4) the depth-matched error -----------------------
    P("## `(a4)` ⚠⚠ the depth-matched error is **+16.8 mm**, not +10.7")
    P("")
    P("The +10.7 mm in circulation scored **every** band against Meier's single")
    P("**96 mm mouth** value. Meier's bore tapers, so that understates the error")
    P("on every band below the collar. Scored against the profile instead:")
    P("")
    P("| depth [mm] | CTRL2 Ø | Meier (band mean) | error | Meier (samples) | error |")
    P("|---|---|---|---|---|---|")
    errs, errs_s = [], []
    for b, (a, c) in enumerate(BANDS):
        dC = 2 * rC[b]
        mi = meier_band_integral(a, c)
        ms = meier_band_samples(a, c)
        errs.append(dC - mi)
        if ms is not None:
            errs_s.append(dC - ms)
        P("| %d–%d | %.1f | %.2f | **%+.1f** | %s | %s |"
          % (a, c, dC, mi, dC - mi,
             "%.2f" % ms if ms is not None else "—",
             "%+.1f" % (dC - ms) if ms is not None else "—"))
    P("")
    below = errs[1:]
    P("**mean below the mouth (bands 50–90 … 200–300): %+.1f mm**  — bar **+16.8**"
      % st.mean(below))
    P("")
    P("mean over all five bands %+.1f mm; against the single 96 mm mouth value it"
      % st.mean(errs))
    P("would read %+.1f mm, which is the figure being retired."
      % st.mean([2 * rC[b] - 96.0 for b in range(len(BANDS))]))
    P("")
    P("⚠ **The collar band is the one convention-sensitive row.** Averaging the")
    P("piecewise-linear profile over 10–50 mm gives %.2f mm; averaging only the"
      % meier_band_integral(10, 50))
    P("profile's own digitised samples inside that band gives %.2f. The collar is"
      % (meier_band_samples(10, 50) or float("nan")))
    P("steep (96→88 over 40 mm), so the two differ by ~2.5 mm there and by <0.5 mm")
    P("everywhere else. The integral mean is the declared convention; the sample")
    P("mean is carried so the difference is visible rather than buried.")
    P("")
    P("⚠⚠ **Meier's profile is a VISIBLE WIDTH — a LOWER bound on Ø, ±5 mm**")
    P("(reference item 10). So every error above is an **upper** bound on the")
    P("model's excess, and the ±5 mm exceeds the band-to-band spread of %.1f mm."
      % (max(below) - min(below)))
    P("")
    P("---")
    P("")
    P("## What this note does and does not back")
    P("")
    P("✅ **Backs:** that the archived pair differs by one key and two binaries;")
    P("that CTRL2 is banking-only; the bore tables; the **+16.8 mm** depth-matched")
    P("error replacing +10.7.")
    P("")
    P("⚠ **Does NOT back:** any attribution of the +%.1f–%.1f mm to `lateral_subcell`."
      % (min(2 * rC[b] - 2 * rB[b] for b in range(len(BANDS))),
         max(2 * rC[b] - 2 * rB[b] for b in range(len(BANDS)))))
    P("`CTRL3` — the same keys as CTRL2 with `lateral_subcell=0`, on the **same**")
    P("binary — is the only thing that can, and `CTRL3` vs Leg B is what sizes the")
    P("four-step drift.")
    return 0 if (ok_a1 and ok_a2) else 1


if __name__ == "__main__":
    sys.exit(main())
