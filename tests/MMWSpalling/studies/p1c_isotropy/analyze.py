#!/usr/bin/env python3
"""P1c criteria (d) (e) (g) (h) (i) (m) -- a reduction, starts no run.

⚠ Firings are counted from the EVENT LOG, never by summing a thermo.dat counter:
`amr.thermo.plot_int` is not 1, and summing a per-step counter over the sampled
rows undercounts by that factor (P1b takeaway 11).

Roundness ((f), (j), (k), (l)) lives in roundness.py.
"""
import csv, math, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)
DZ = DX = 0.002
KTOP0 = 299
R_BURNER = 0.040
P1_T_END, P1_FREEZE = 177.616, 120.0
WIN = 25.0
P1B = os.path.join(STUD, "p1b_support_rule", "output_fh")


def load(path):
    rows = []
    for d in csv.DictReader(open(path)):
        rows.append({"t": float(d["time"]),
                     "i": int(d["col_i"]), "j": int(d["col_j"]),
                     "reg": int(float(d["regime"])),
                     "nv": int(float(d["n_voided"])),
                     "ne": int(float(d.get("n_exposed", 0) or 0)),
                     "pil": int(float(d.get("pillar", 0) or 0)),
                     "Ts": float(d.get("T_surf_sup", 0) or 0.0)})
    rows.sort(key=lambda r: r["t"])
    return rows


def r_of(i, j):
    return math.hypot((i + 0.5) * DX, (j + 0.5) * DX)


def sector_of(i, j):
    return int(math.degrees(math.atan2(j + 0.5, i + 0.5)) // 15) * 15


def thermo_max(path, key):
    if not os.path.exists(path):
        return None
    L = open(path).read().strip().split("\n")
    h = L[0].split()
    if key not in h:
        return None
    i = h.index(key)
    v = [abs(float(r.split()[i])) for r in L[1:] if len(r.split()) > i]
    return max(v) if v else None


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "output_fh")
    ev_path = os.path.join(out, "FH_removal_events.csv")
    if not os.path.exists(ev_path):
        print("missing: %s" % ev_path); return 1
    ev = load(ev_path)
    t_end = ev[-1]["t"] if ev else 0.0
    sup = [e for e in ev if e["reg"] == 5]
    fails, reports = [], []

    def gate(ok, msg):
        print(("  %s  " % ("ok  " if ok else "FAIL")) + msg)
        if not ok: fails.append(msg)

    print("P1c leg: %s\n  %d events, t_max %.3f s\n" % (ev_path, len(ev), t_end))

    # ---- (d) the run completes -- PRIMARY ---------------------------------
    print("(d) the run reaches stop_time -- PRIMARY GATE")
    log = os.path.join(out, "FH.log")
    txt = open(log).read() if os.path.exists(log) else ""
    coll = ("collision" in txt.lower()) or ("RAN OUT OF ROCK" in txt)
    gate(t_end >= 499.9 and not coll,
         "t_max %.3f s (bar 500.0), collision/abort in log: %s"
         % (t_end, "YES" if coll else "no"))

    # ---- (e) the stall stays cured -- PRIMARY -----------------------------
    print("(e) the stall stays cured -- PRIMARY GATE")
    last, final = {}, {}
    for e in ev:
        if e["t"] > P1_T_END: break
        k = (e["i"], e["j"])
        last[k] = max(last.get(k, -1.0), e["t"])
        final[k] = final.get(k, KTOP0) - e["nv"]
    under = [k for k in final if r_of(*k) < R_BURNER]
    frozen = [k for k in under if last[k] < P1_FREEZE]
    gate(len(frozen) == 0, "frozen columns %d of %d (bar 0; P1 had 152)"
         % (len(frozen), len(under)))
    dep = {}
    for e in ev:
        dep[(e["i"], e["j"])] = dep.get((e["i"], e["j"]), KTOP0) - e["nv"]
    def d_of(c): return (KTOP0 - dep[c]) * DZ
    # ⚠ THE STEP IS SCORED AT THE SAME CUT AS THE FROZEN COUNT (t <= 177.616 s),
    # because that is the half of (e) whose cut the packet states, and because it
    # is what makes P1b read the 2.0 / 8.0 mm the packet quotes. Scored at the
    # FINAL state instead, P1b reads 4.0 / 16.0 and would fail its own p90 bar --
    # so the cut is load-bearing and is reported both ways below.
    def step_stats(state):
        st = []
        for (i, j) in [c for c in state if r_of(*c) < R_BURNER]:
            for di, dj in ((1, 0), (0, 1)):
                nb = (i + di, j + dj)
                if nb in state and r_of(*nb) < R_BURNER:
                    st.append(abs((state[(i, j)] - state[nb]) * DZ))
        if not st:
            return None, None
        return (1000 * statistics.median(st),
                1000 * sorted(st)[9 * len(st) // 10])
    med, p90 = step_stats(final)                 # t <= 177.616 s -- GATED
    fmed, fp90 = step_stats(dep)                 # t <= 500 s -- reported
    if med is not None:
        gate(med <= 6.0 and p90 <= 8.0,
             "step at P1's cut: median %.1f mm (bar <= 6), p90 %.1f mm "
             "(bar <= 8)   [P1b 2.0 / 8.0; P1 21 / 70]" % (med, p90))
        print("        reported, FINAL state at 500 s: median %.1f mm, p90 "
              "%.1f mm   [P1b 4.0 / 16.0]" % (fmed, fp90))

    # ---- (g) the ledger ---------------------------------------------------
    print("(g) ledger_err_solid")
    m = thermo_max(os.path.join(out, "FH", "thermo.dat"), "ledger_err_solid")
    gate(m is not None and m <= 5e-3,
         "worst |ledger_err_solid| %s (bar 5e-3; P1b 7.60e-4)"
         % ("n/a" if m is None else "%.4e" % m))

    # ---- (h) the firing temperature --------------------------------------
    print("(h) support-firing temperature")
    T = [e["Ts"] for e in sup if e["Ts"] > 0.0]
    if T:
        gate(min(T) >= 700.0,
             "min %.1f K, median %.1f K (bar: none below 700; P1b min 738.7). "
             "The absolute form at this flaw would fire at ~525 K."
             % (min(T), statistics.median(T)))
    else:
        print("  ⚠ no support firing carried a temperature -- see (m)")

    # ---- (i) the decay -- REPORT, NOT A GATE ------------------------------
    print("\n(i) the decay -- ⚠⚠ REPORTED, NOT A GATE. At P1b's 0.49 %% static\n"
          "    reach the rule is very nearly switched off, so this passes\n"
          "    because there is almost nothing to decay. NOT evidence of repair.")
    nwin = int(math.ceil(500.0 / WIN))
    hist = [0] * nwin
    crit = [0] * nwin
    for e in ev:
        w = min(nwin - 1, int(e["t"] / WIN))
        if e["reg"] == 5: hist[w] += 1
        elif e["reg"] == 1: crit[w] += 1
    peak = max(hist) if hist else 0
    tail = sum(hist[-4:]) / 4.0
    for w in range(nwin):
        if hist[w] or w % 4 == 0:
            print("    %6.0f-%-6.0f  support %5d  criterion %6d"
                  % (w * WIN, (w + 1) * WIN, hist[w], crit[w]))
    print("    peak %d, tail (mean of last four bins) %.1f  ->  tail/peak %s"
          % (peak, tail, ("%.3f" % (tail / peak)) if peak else "n/a (peak 0)"))
    print("    P1b for comparison: peak 1716, tail 1062.5, tail/peak 0.619")

    # ---- (m) the firings themselves --------------------------------------
    print("\n(m) support firings -- REPORTED, and witness 7")
    print("    total %d   (P1b 22,800; static replay floor for min_sides=3 "
          "was 111)" % len(sup))
    if sup:
        h = {}
        for e in sup: h[e["ne"]] = h.get(e["ne"], 0) + 1
        print("    n_exposed histogram: %s" % dict(sorted(h.items())))
        po = [e for e in sup if e["ne"] < 3 and e["pil"]]
        print("    pillar-only firings: %d of %d (%.1f %%)"
              % (len(po), len(sup), 100.0 * len(po) / len(sup)))
        if po:
            cnt = {}
            for e in po: cnt[sector_of(e["i"], e["j"])] = \
                cnt.get(sector_of(e["i"], e["j"]), 0) + 1
            print("    pillar-only by sector: %s" % dict(sorted(cnt.items())))
            print("    pillar-only max radius: %.1f mm"
                  % (1000 * max(r_of(e["i"], e["j"]) for e in po)))
        # ⚠⚠ WITNESS 7 -- was spall.support_min_sides actually READ?
        # ⚠ The test is NOT "no firing at n_exposed = 2". The `|| pillar` branch
        # fires an opposing-void column at ANY min_sides (Removal.H:1511), so
        # n_exposed = 2 firings are LEGITIMATE when every one of them is a
        # pillar -- and the packet predicts that is precisely what will fire
        # here. The discriminating statement is: every sub-threshold firing must
        # be a pillar, AND the total must differ from P1b's 22,800.
        sub = [e for e in sup if e["ne"] < 3]
        sub_nonpillar = [e for e in sub if not e["pil"]]
        print("    ⚠⚠ witness 7 (was min_sides=3 honoured?)")
        print("        firings below the 3-side threshold: %d, of which "
              "NON-pillar: %d (must be 0)" % (len(sub), len(sub_nonpillar)))
        print("        total %d vs P1b's 22800 (must differ)" % len(sup))
        if sub_nonpillar:
            fails.append("(witness 7) %d firings with n_exposed < 3 and NOT a "
                         "pillar -- spall.support_min_sides=3 was not honoured"
                         % len(sub_nonpillar))
        if len(sup) == 22800:
            fails.append("(witness 7) the support-firing count equals P1b's "
                         "22800 exactly -- the key was almost certainly ignored")
    else:
        print("    ⚠⚠ ZERO support firings. (m) cannot separate 'the rule did "
              "nothing'\n        from 'the key was never read' on its own -- "
              "check sup_cols in thermo.dat.")

    # P6: floor depth vs P1b
    print("\nP6 floor rate and volume vs P1b (reported)")
    axis = [c for c in dep if r_of(*c) < 0.010]
    if axis:
        print("    floor depth at r < 10 mm: median %.1f mm, max %.1f mm"
              % (1000 * statistics.median([d_of(c) for c in axis]),
                 1000 * max(d_of(c) for c in axis)))
    nv_here = sum(e["nv"] for e in ev)
    p1b_ev = os.path.join(P1B, "FH_removal_events.csv")
    nv_p1b = "n/a"
    if os.path.exists(p1b_ev):
        nv_p1b = str(sum(int(float(d["n_voided"]))
                         for d in csv.DictReader(open(p1b_ev))))
    print("    cells removed (ONE QUADRANT; x4 for any Meier comparison): %d"
          "   -- P1b: %s" % (nv_here, nv_p1b))
    print("    ⚠ P5 predicts this is LOWER than P1b's. A volume drop is the "
          "signature of\n      removing an artefact, not a regression.")
    print()
    if fails:
        for f in fails: print("FAIL: %s" % f)
        return 1
    print("all gated rows in this script PASS (roundness (f) is in roundness.py)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
