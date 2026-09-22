#!/usr/bin/env python3
"""D2k arbiter analysis: the three CRITERION.md conditions, the band table in
the D2j-0b format, the side-power fraction and the f sensitivity. Writes
analysis.md and d2k_bands.png.

Reuses studies/d2j0b_plane/analyze.py (loaded by path; frozen, not edited) for
col_rate / band_timing / r_front / metrics, so every number is computed the
same way as the C0/C1/C2 reference it is compared against.
"""
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUD = os.path.dirname(HERE)


def load(name, path):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


A = load("d2j0b_analyze", os.path.join(STUD, "d2j0b_plane", "analyze.py"))
K = load("d2k_run", os.path.join(HERE, "run.py"))
OUT = K.OUT
D2J0B_OUT = os.path.join(STUD, "d2j0b_plane", "output")
TOL = 0.05
BLOCKS = [(50.0, 100.0), (100.0, 150.0), (150.0, 200.0)]   # the deciding blocks
# D2j-0b, key off (CRITERION.md "numbers to beat")
REF_GAP = {"C0": (-10.4, -15.3, -20.1), "C1": (-13.5, -28.7, -45.5)}


def metrics_in(out, name):
    """A.metrics against a chosen output directory."""
    old = A.OUT
    A.OUT = out
    A.R.OUT = out
    A.R.J.OUT = out
    try:
        return A.metrics(name)
    finally:
        A.OUT = old
        A.R.OUT = old
        A.R.J.OUT = old


def thermo(out, name):
    p = os.path.join(out, name, "thermo.dat")
    head = open(p).readline().split()
    a = np.loadtxt(p, skiprows=1, ndmin=2)
    return {c: a[:, i] for i, c in enumerate(head)}


def gaps(d2, d1):
    out = []
    for b in BLOCKS:
        if b in d2["blocks"] and b in d1["blocks"]:
            r2, r1 = d2["blocks"][b], d1["blocks"][b]
            out.append((b, r2, r1, r1 / r2 - 1.0))
    return out


def block_mean(th, a, b, col):
    s = (th["time"] >= a) & (th["time"] <= b)
    return float(np.mean(th[col][s])) if s.any() else float("nan")


def main():
    L = ["# D2k arbiter analysis (side-face flux on the D2j-0b C1 cascade)", ""]
    d2 = metrics_in(OUT, "C1_2mm")
    d1 = metrics_in(OUT, "C1_1mm")
    d07 = metrics_in(OUT, "C1f07_2mm") if os.path.exists(os.path.join(OUT, "C1f07_2mm.done")) else None

    L += ["## 1. Disk-mean rate (r < 30 mm) by block, 2 mm / 1 mm", "",
          "| block | 2 mm | 1 mm | gap | D2j-0b C1 (key off) | D2j-0b C0 |",
          "|---|---|---|---|---|---|"]
    g = gaps(d2, d1)
    for i, (b, r2, r1, gp) in enumerate(g):
        L.append(f"| {b[0]:.0f}–{b[1]:.0f} s | {r2:.4f} | {r1:.4f} | **{gp*100:+.1f} %** | "
                 f"{REF_GAP['C1'][i]:+.1f} % | {REF_GAP['C0'][i]:+.1f} % |")
    ok_a = all(abs(gp) <= TOL for _, _, _, gp in g) and len(g) == len(BLOCKS)
    L += ["", f"**(a) {'PASS' if ok_a else 'FAIL'}** — every deciding block within 5 %."]

    L += ["", "## 2. Band timing (D2j-0b format): t_stop / t_cross / Δ_own / Δ_outer", "",
          "A band **stalls** only if it STOPS FIRING before its own face crosses the plane.",
          "D2j-0b marks a band still firing at the end with t_stop = t_end (`250 (firing)`),",
          "so t_stop >= t_end - 5 s is NOT a stall however large t_cross - t_stop looks.",
          "",
          "| run | band [mm] | n | t_stop | t_cross | Δ_own | below plane at stop [mm] | Δ_outer |",
          "|---|---|---|---|---|---|---|---|"]
    d_outer = {}
    stalled = {}
    for d in (d2, d1):
        bands = sorted(d["bands"])
        for idx, key in enumerate(bands):
            ts, tc, db, n = d["bands"][key]
            d_own = tc - ts
            nxt = bands[idx + 1] if idx + 1 < len(bands) else None
            do = (ts - d["bands"][nxt][1]) if nxt else float("nan")
            d_outer.setdefault(d["name"], {})[key] = do
            firing = ts >= d["t_end"] - 5.0
            if (not firing) and d_own > 10.0 and key[1] <= 28:
                stalled.setdefault(d["name"], []).append(key)
            L.append(f"| {d['name']} | {key[0]}–{key[1]} | {n} | "
                     f"{A.fmt_t(ts)}{' (firing)' if firing else ''} | {A.fmt_t(tc)} | "
                     f"{'—' if firing else A.fmt_t(d_own)} | {db*1e3:.1f} | {A.fmt_t(do)} |")
    no_stall = not stalled
    common = [k for k in d_outer[d2["name"]] if k in d_outer[d1["name"]]]
    rel = []
    for k in common:
        a2, a1 = d_outer[d2["name"]][k], d_outer[d1["name"]][k]
        if np.isfinite(a2) and np.isfinite(a1) and abs(a2) > 1e-9:
            rel.append(abs(a1 / a2 - 1.0))
    ok_b = no_stall or (len(rel) > 0 and max(rel) <= 0.20)
    L += ["", f"**(b) {'PASS' if ok_b else 'FAIL'}** — "
          + ("no band inside 28 mm stalls with Δ_own > 10 s at either mesh (the better outcome)."
             if no_stall else
             f"bands stalling: {stalled}; Δ_outer mesh spread max "
             f"{max(rel)*100:.0f} %" if rel else "Δ_outer not comparable (no finite pairs).")]

    L += ["", "## 3. r_front(t) at 1 mm [mm] (inner edge of the stalled run reaching 30 mm)", "",
          "| run | " + " | ".join(f"{t:.0f} s" for t in A.FRONT_T) + " |",
          "|---|" + "---|" * len(A.FRONT_T)]
    for d in (d2, d1):
        L.append(f"| {d['name']} | " + " | ".join(
            "—" if t not in d["front"] or not np.isfinite(d["front"][t])
            else f"{d['front'][t]*1e3:.0f}" for t in A.FRONT_T) + " |")
    fr = [d1["front"][t] for t in A.FRONT_T if t in d1["front"] and np.isfinite(d1["front"][t])]
    ok_c = len(fr) <= 1 or all(fr[i + 1] >= fr[i] - 1e-9 for i in range(len(fr) - 1))
    L += ["", f"**(c) {'PASS' if ok_c else 'FAIL'}** — "
          + ("no inward-walking front at 1 mm." if ok_c
             else f"r_front walks inward: {[f'{x*1e3:.0f}' for x in fr]} mm. "
                  "D2j-0b C1 walked 26 -> 22 -> 18 mm.")]

    L += ["", "## 4. Side power as a fraction of top-face power (reported, not deciding)", "",
          "| run | block | side_P_face [W] | patch_P_robin [W] | side / top | side_faces | max abs ledger_err |",
          "|---|---|---|---|---|---|---|"]
    for n in ("C1_2mm", "C1_1mm") + (("C1f07_2mm",) if d07 else ()):
        th = thermo(OUT, n)
        for b in BLOCKS:
            sp = block_mean(th, b[0], b[1], "side_P_face")
            tp = block_mean(th, b[0], b[1], "patch_P_robin")
            sf = block_mean(th, b[0], b[1], "side_faces")
            s = (th["time"] >= b[0]) & (th["time"] <= b[1])
            le = float(np.max(np.abs(th["ledger_err"][s]))) if s.any() else float("nan")
            L.append(f"| {n} | {b[0]:.0f}–{b[1]:.0f} s | {sp:.1f} | {tp:.1f} | "
                     f"{sp/tp*100 if tp else float('nan'):.1f} % | {sf:.0f} | {le:.1e} |")

    if d07:
        L += ["", "## 5. f sensitivity at 2 mm (pre-registered, never fitted)", "",
              "| f | run | " + " | ".join(f"{b[0]:.0f}–{b[1]:.0f} s" for b in BLOCKS) + " |",
              "|---|---|" + "---|" * len(BLOCKS)]
        for lab, d in (("1.0", d2), ("0.7", d07)):
            L.append(f"| {lab} | {d['name']} | " + " | ".join(
                f"{d['blocks'][b]:.4f}" for b in BLOCKS if b in d["blocks"]) + " |")
        rows = [d2["blocks"][b] for b in BLOCKS if b in d2["blocks"]]
        r07 = [d07["blocks"][b] for b in BLOCKS if b in d07["blocks"]]
        if len(rows) == len(r07):
            L.append("")
            L.append("f = 0.7 vs f = 1.0: " + ", ".join(
                f"{(y/x-1)*100:+.1f} %" for x, y in zip(rows, r07)))

    verdict = "PASS" if (ok_a and ok_b and ok_c) else "FAIL"
    L += ["", f"## 6. Verdict: {verdict}", "",
          f"(a) blocks " + ", ".join(f"{gp*100:+.1f}" for _, _, _, gp in g)
          + f" % -> {'ok' if ok_a else 'FAILS'}; (b) {'ok' if ok_b else 'FAILS'}; "
          f"(c) {'ok' if ok_c else 'FAILS'}.",
          "",
          "Per CRITERION.md: " + ("run the conditional Meier pair under a separately hashed "
                                  "D2i criterion." if verdict == "PASS" else
                                  "**stop** — do not run the Meier pair, do not tune f, and "
                                  "name what the residual looks like.")]
    txt = "\n".join(L)
    open(os.path.join(HERE, "analysis.md"), "w").write(txt + "\n")
    print(txt)
    fig(d2, d1, d07)


def fig(d2, d1, d07):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    mids = [0.5 * (b[0] + b[1]) for b in BLOCKS]
    g = gaps(d2, d1)
    ax[0].plot(mids[:len(g)], [x[3] * 100 for x in g], "s-", lw=2, ms=8,
               label="D2k: side-face flux on")
    for k, st in (("C1", "x--"), ("C0", "^:")):
        ax[0].plot(mids, REF_GAP[k], st, color="gray" if k == "C1" else "C2",
                   label=f"D2j-0b {k} (key off)")
    ax[0].axhspan(-5, 5, color="g", alpha=0.1, label="±5 %")
    ax[0].axhline(0, color="k", lw=0.6)
    ax[0].set(xlabel="block centre [s]", ylabel="1 mm / 2 mm disk-mean rate − 1 [%]",
              title="Mesh gap, disk mean r < 30 mm")
    ax[0].grid(alpha=0.3)
    ax[0].legend(fontsize=8)
    for d, st in ((d2, "o-"), (d1, "s-")):
        ks = sorted(d["bands"])
        ax[1].plot([0.5 * (k[0] + k[1]) for k in ks],
                   [d["bands"][k][0] for k in ks], st, label=f"{d['name']} t_stop")
        ax[1].plot([0.5 * (k[0] + k[1]) for k in ks],
                   [d["bands"][k][1] for k in ks], st, alpha=0.4,
                   label=f"{d['name']} t_cross")
    ax[1].set(xlabel="band centre [mm]", ylabel="time [s]",
              title="Band stop vs plane crossing")
    ax[1].grid(alpha=0.3)
    ax[1].legend(fontsize=7)
    f.suptitle("D2k: flux on exposed vertical faces, against the D2j-0b C1 cascade", y=0.99)
    f.tight_layout()
    f.savefig(os.path.join(HERE, "d2k_bands.png"), dpi=110)
    plt.close(f)
    print("wrote d2k_bands.png")


if __name__ == "__main__":
    main()
