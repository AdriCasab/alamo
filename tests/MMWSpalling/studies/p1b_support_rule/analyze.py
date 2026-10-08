#!/usr/bin/env python3
"""P1b criteria (c), (e), (f) and deliverable 5 -- a reduction, starts no run.

  (c) the support-firing count per 25 s window must DECAY, monotonically after
      an initial transient, to < 20 % of its peak by 400-500 s. ⚠ If it does
      not decay the rule has become a removal law and the packet FAILS; that is
      the result, not a tuning target.
  (e) frozen columns at P1's cut (t <= 177.6 s, r < 40 mm, last fired < 120 s):
      152 -> < 30.
  (f) adjacent-column depth step inside r < 40 mm: median 21 mm -> < 6 mm.
  d5  the exposed-side histogram over time, and the support-firing T_surf
      median (⚠ witness 8: it must sit at the criterion's ~768 K, NOT ~240 K
      below it, which is where the absolute form at the same flaw would put it).

Usage:  analyze.py [<event-log>]     (default: output_fh/FH_removal_events.csv)
"""
import csv, math, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DX = DZ = 0.002
KTOP0 = 299
R_BURNER = 0.040
P1_T_END, P1_FREEZE = 177.616, 120.0
RING = (0.044, 0.050)   # the wall band the isotropy witness reads
WIN = 25.0


def load(path):
    rows = []
    for d in csv.DictReader(open(path)):
        rows.append({
            "t": float(d["time"]),
            "i": int(d["col_i"]), "j": int(d["col_j"]),
            "reg": int(float(d["regime"])),
            "nv": int(float(d["n_voided"])),
            "ne": int(float(d.get("n_exposed", 0) or 0)),
            "pil": int(float(d.get("pillar", 0) or 0)),
            "Ts": float(d.get("T_surf_sup", 0) or 0.0),
        })
    rows.sort(key=lambda r: r["t"])
    return rows


def r_of(i, j):
    return math.hypot((i + 0.5) * DX, (j + 0.5) * DX)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        HERE, "output_fh", "FH_removal_events.csv")
    if not os.path.exists(path):
        print("missing: %s" % path); return 1
    ev = load(path)
    t_end = ev[-1]["t"] if ev else 0.0
    print("%s\n  %d events, t_max %.3f s\n" % (path, len(ev), t_end))

    # ---- (c) the decay ---------------------------------------------------
    print("(c) support firings per %.0f s window -- MUST DECAY" % WIN)
    nwin = int(math.ceil(t_end / WIN))
    hist = [0] * max(1, nwin)
    crit = [0] * max(1, nwin)
    for e in ev:
        w = min(len(hist) - 1, int(e["t"] / WIN))
        if e["reg"] == 5: hist[w] += 1
        elif e["reg"] == 1: crit[w] += 1
    peak = max(hist) if hist else 0
    for w, (n, c) in enumerate(zip(hist, crit)):
        bar = "#" * int(40.0 * n / peak) if peak else ""
        print("    %6.0f-%-6.0f  support %5d  criterion %6d  %s"
              % (w * WIN, (w + 1) * WIN, n, c, bar))
    tail = [n for w, n in enumerate(hist) if w * WIN >= 400.0]
    print("    peak %d; tail (400-500 s) %s" % (peak, tail))
    if peak and tail:
        print("    tail/peak = %.3f   (bar: < 0.20)  ==> %s"
              % (max(tail) / peak,
                 "PASS" if max(tail) / peak < 0.20 else "FAIL"))
    else:
        print("    ⚠ not scorable: the run did not reach 400 s")

    # ---- d5 the exposed-side histogram and the firing temperature --------
    print("\nd5  exposed sides at support firings, and T_surf")
    sup = [e for e in ev if e["reg"] == 5]
    if sup:
        h = [0] * 5
        for e in sup: h[min(4, e["ne"])] += 1
        tot = len(sup)
        print("    n_exposed  " + "  ".join(
            "%d:%4.0f %%" % (k, 100.0 * h[k] / tot) for k in range(5)))
        print("    pillars    %d of %d (%.1f %%)"
              % (sum(e["pil"] for e in sup), tot,
                 100.0 * sum(e["pil"] for e in sup) / tot))
        Ts = sorted(e["Ts"] for e in sup if e["Ts"] > 0.0)
        if Ts:
            print("    T_surf at a support firing: median %.1f K  "
                  "(p10 %.1f, p90 %.1f, n = %d)"
                  % (statistics.median(Ts), Ts[len(Ts) // 10],
                     Ts[9 * len(Ts) // 10], len(Ts)))
            print("    ⚠ witness 8: this must sit at the criterion's ~768 K, "
                  "NOT at the absolute\n       form's ~525 K.")
        # time-resolved
        print("    n_exposed >= 3 share by window:")
        for w in range(len(hist)):
            s = [e for e in sup if int(e["t"] / WIN) == w]
            if not s: continue
            print("      %6.0f-%-6.0f  %4.0f %%  (n = %d)"
                  % (w * WIN, (w + 1) * WIN,
                     100.0 * sum(1 for e in s if e["ne"] >= 3) / len(s),
                     len(s)))
    else:
        print("    ⚠⚠ NO SUPPORT FIRINGS AT ALL -- the rule is inert "
              "(packet witness 5).")

    # ---- ⚠⚠ THE AZIMUTHAL-ISOTROPY WITNESS -------------------------------
    # The check this step needed and did not have. A face-neighbour count is NOT
    # grid-isotropic: a ROUND inward-descending step has TWO inward-lower face
    # neighbours along its DIAGONALS (both i-1 and j-1 sit at smaller radius,
    # hence deeper) and only ONE along its AXES (there the tangential neighbours
    # are at the same radius and the same depth). So support_min_sides = 2 fires
    # on the diagonals of a SMOOTH round wall and spares its axes -- which is a
    # removal law with a square-grid signature, not a repair.
    # ⚠ A MEDIAN DEPTH STEP CANNOT SEE THIS. Criterion (f) passed at 2.0 mm while
    # the 44-50 mm ring carried an 18x azimuthal spread just outside its radius.
    print("\n⚠⚠ azimuthal isotropy -- final depth of the %.0f-%.0f mm ring by "
          "15 deg sector\n    (a round hole is flat here; rev 7 reads 1.1x. "
          "bar: spread <= 1.03x)" % (1000 * RING[0], 1000 * RING[1]))
    dep_all = {}
    for e in ev:
        k = (e["i"], e["j"])
        dep_all[k] = dep_all.get(k, KTOP0) - e["nv"]
    bins = {}
    for (i, j), kt in dep_all.items():
        x, y = (i + 0.5) * DX, (j + 0.5) * DX
        r = math.hypot(x, y)
        if not (RING[0] <= r < RING[1]):
            continue
        b = int(math.degrees(math.atan2(y, x)) // 15) * 15
        bins.setdefault(b, []).append((KTOP0 - kt) * DZ)
    meds = []
    for b in sorted(bins):
        m = statistics.median(bins[b])
        meds.append(m)
        print("    %3d-%3d deg  n=%2d  median %6.1f mm  max %6.1f mm"
              % (b, b + 15, len(bins[b]), 1000 * m, 1000 * max(bins[b])))
    if meds and min(meds) > 0:
        sp = max(meds) / min(meds)
        print("    spread %.1fx  ==> %s" % (sp, "PASS" if sp <= 1.03 else
              "FAIL -- the rule is drilling the grid diagonals"))
    # support firings by azimuth, over the band they concentrate in
    cnt = {}
    for e in ev:
        if e["reg"] != 5:
            continue
        x, y = (e["i"] + 0.5) * DX, (e["j"] + 0.5) * DX
        r = math.hypot(x, y)
        if not (0.025 <= r < 0.050):
            continue
        b = int(math.degrees(math.atan2(y, x)) // 15) * 15
        cnt[b] = cnt.get(b, 0) + 1
    if cnt:
        ax = cnt.get(0, 0) + cnt.get(75, 0)
        dg = cnt.get(30, 0) + cnt.get(45, 0)
        print("    support firings 25-50 mm: near-axis %d, near-diagonal %d  "
              "=> %.1f:1" % (ax, dg, dg / max(1, ax)))

    # ---- (e) frozen columns, P1's cut ------------------------------------
    print("\n(e) frozen columns at P1's cut (t <= %.3f, r < %.0f mm, "
          "last fired < %.0f s)" % (P1_T_END, 1000 * R_BURNER, P1_FREEZE))
    last, final = {}, {}
    for e in ev:
        if e["t"] > P1_T_END: break
        k = (e["i"], e["j"])
        last[k] = max(last.get(k, -1.0), e["t"])
        final[k] = final.get(k, KTOP0) - e["nv"]
    under = [k for k in final if r_of(*k) < R_BURNER]
    frozen = [k for k in under if last[k] < P1_FREEZE]
    # ⚠ NOT SCORABLE before the run passes P1's cut. On a shorter log every
    # column is trivially "frozen" -- its last firing is before t = 120 s
    # because the run never reached 120 s -- and the criterion reads 316/316,
    # which looks like a catastrophic failure and is an artefact of the cut.
    if t_end < P1_T_END:
        print("    under the burner %d, frozen %d  -- ⚠ NOT SCORABLE: the log "
              "ends at %.1f s,\n        before the %.1f s cut, so every column "
              "is trivially frozen."
              % (len(under), len(frozen), t_end, P1_T_END))
    else:
        print("    under the burner %d, frozen %d  (P1: 316 / 152; bar < 30)  "
              "==> %s" % (len(under), len(frozen),
                          "PASS" if len(frozen) < 30 else "FAIL"))

    # ---- (f) roughness ---------------------------------------------------
    def depth(k): return (KTOP0 - final[k]) * DZ
    steps = []
    for (i, j) in under:
        for di, dj in ((1, 0), (0, 1)):
            nb = (i + di, j + dj)
            if nb in final and r_of(*nb) < R_BURNER:
                steps.append(abs(depth((i, j)) - depth(nb)))
    if steps:
        med = 1000 * statistics.median(steps)
        print("\n(f) adjacent-column depth step, r < %.0f mm: median %.1f mm, "
              "p90 %.1f mm\n    (P1: median 21 mm; bar < 6 mm)  ==> %s"
              % (1000 * R_BURNER, med,
                 1000 * sorted(steps)[9 * len(steps) // 10],
                 "PASS" if med < 6.0 else "FAIL"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
