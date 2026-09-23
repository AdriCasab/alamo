#!/usr/bin/env python3
"""D2t scoring. Criteria are frozen in CRITERION.md (sha256 in CRITERION.sha256),
hashed before the first campaign leg launched.

  analyze.py            score everything that has finished
  analyze.py --quiet    tables only

THE DELIVERABLE IS A SHAPE METRIC, NOT A RATE. Every metric in this project so
far has been a disk-mean rate, which is exactly how a mask-edge artefact
survived four packets: an averaged number cannot say WHERE two meshes disagree.
Here:

  z(r)  mean face height per FIXED 2 mm physical radial bin (asserted: the bin
        width never moves with the mesh, or the metric would be self-referential)
  D(z)  diameter versus depth on a FIXED 10-100 mm ladder

A disk-mean rate appears once, in (d), purely as a regression check against a
reference recomputed on the same mask. It is never the result.
"""
import argparse
import ast
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_sp = importlib.util.spec_from_file_location("d2t_run", os.path.join(HERE, "run.py"))
R = importlib.util.module_from_spec(_sp)
_sp.loader.exec_module(R)

OUT = R.OUT
FROZ = os.path.join(os.path.dirname(HERE), "d2r1_aboveplane_sides", "output")
LZ, Z0, FEED = 0.15, 0.170, 4.3056e-4
BIN = 2.0e-3                                  # FIXED PHYSICAL radial bin
LADDER = np.arange(0.010, 0.1001, 0.010)      # fixed depth ladder, 10-100 mm
LADDER_EXTRA = 0.110                          # reported, ungated
BLOCKS = [(50.0, 100.0), (100.0, 150.0), (150.0, 200.0), (200.0, 250.0)]
BANDS = [(0.0, R.R_FLAT, "core r<18"), (R.R_FLAT, R.R_TAPER, "band 18-34"),
         (R.R_TAPER, R.R_PATCH, "floor 34-50")]
REF_D = {"2mm": 1.6667, "1mm": 1.6349}        # T_2mm / T_1mm, r<18 mm, 150-200 s
R_HALF = 0.026                                # h >= 0.5 x peak
FLOOR_CLEAR = 0.010                           # domain watch


def have(n, d=OUT):
    return os.path.exists(os.path.join(d, n + ".done"))


def load(n, d=OUT):
    m = ast.literal_eval(open(os.path.join(d, n + ".done")).read())
    dz, nxy = m["dz"], m["nxy"]
    c = (np.arange(nxy) + 0.5) * dz
    r = np.hypot(*np.meshgrid(c, c, indexing="ij"))
    a = np.loadtxt(os.path.join(d, n + "_removal_events.csv"), delimiter=",",
                   skiprows=1, ndmin=2)
    return m, r, a


def z_face(m, r, a):
    """Final face height per column: LZ - (cells voided)*dz. Columns with no
    event keep the initial top, which is what a zero sum gives."""
    nv = np.zeros_like(r)
    np.add.at(nv, (a[:, 1].astype(int), a[:, 2].astype(int)), a[:, 10])
    return LZ - nv * m["dz"], nv


def profile(m, r, zf, rmax=None):
    """z(r) on FIXED 2 mm physical bins. The bin width is a constant of this
    module, never derived from dz -- N4."""
    assert BIN == 2.0e-3, "the radial bin must be a fixed physical width"
    rmax = rmax if rmax is not None else R.R_PATCH
    edges = np.arange(0.0, rmax + 0.5 * BIN, BIN)
    rc, zc, nc = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (r >= lo) & (r < hi)
        if not s.any():
            continue
        rc.append(0.5 * (lo + hi))
        zc.append(float(zf[s].mean()))
        nc.append(int(s.sum()))
    return np.array(rc), np.array(zc), np.array(nc)


def diameter(rc, zc, depth):
    """D(z) = 2*max{ r : z(r) <= z }, z = LZ - depth. nan if nothing reached it."""
    z = LZ - depth
    s = np.nonzero(zc <= z + 1e-12)[0]
    return 2.0 * rc[s].max() if len(s) else np.nan


def rate(r, a, win, mask):
    s = (a[:, 0] >= win[0]) & (a[:, 0] <= win[1])
    v = np.zeros_like(r)
    np.add.at(v, (a[s, 1].astype(int), a[s, 2].astype(int)), a[s, 9])
    return float(v[mask].mean()) / (win[1] - win[0]) * 3600.0


def latch(m, r, a, dt=0.5):
    """Above-plane count at each block end, and plane crossings per block.
    s(t) = z_n(t) - z_face(t); z_n falls continuously, z_face only at events."""
    nxy = m["nxy"]
    t = np.arange(0.0, m["stop"] + dt, dt)
    ncol = nxy * nxy
    nv = np.zeros((ncol, len(t)))
    ci = a[:, 1].astype(int) * nxy + a[:, 2].astype(int)
    ti = np.clip(np.searchsorted(t, a[:, 0]), 0, len(t) - 1)
    np.add.at(nv, (ci, ti), a[:, 10])
    zf = LZ - np.cumsum(nv, axis=1) * m["dz"]
    s = (Z0 - FEED * t)[None, :] - zf
    inp = (r <= R.R_PATCH).ravel()
    above = s <= 0.0
    out = []
    for lo, hi in BLOCKS:
        k = np.searchsorted(t, hi) - 1
        w = (t >= lo) & (t <= hi)
        cross = int((np.diff(above[inp][:, w].astype(int), axis=1) != 0).sum())
        out.append((lo, hi, int(above[inp, k].sum()), cross))
    return out


def score(tag, legs, rows, quiet):
    """(a)-style shape convergence between a 2 mm and a 1 mm leg."""
    (m2, r2, a2), (m1, r1, a1) = legs
    zf2, _ = z_face(m2, r2, a2)
    zf1, _ = z_face(m1, r1, a1)
    rc2, zc2, nc2 = profile(m2, r2, zf2)
    rc1, zc1, nc1 = profile(m1, r1, zf1)
    assert np.allclose(rc2, rc1), "bins differ between meshes"
    if not quiet:
        print(f"\n--- {tag}: z(r) on fixed 2 mm bins, recession = LZ - z [mm] ---")
        print(" r[mm]   n(2mm) n(1mm)   rec_2mm  rec_1mm     diff   |  step2  step1")
        d2 = np.r_[0.0, np.diff(zc2)]
        d1 = np.r_[0.0, np.diff(zc1)]
        for i in range(len(rc2)):
            print(f"{rc2[i]*1e3:6.1f} {nc2[i]:7d} {nc1[i]:6d}   "
                  f"{(LZ-zc2[i])*1e3:8.2f} {(LZ-zc1[i])*1e3:8.2f} "
                  f"{(zc1[i]-zc2[i])*1e3:+8.2f}   | {d2[i]*1e3:6.1f} {d1[i]*1e3:6.1f}")
    core = rc2 < R.R_FLAT
    a_core = float(np.abs(zc1[core] - zc2[core]).max()) * 1e3
    lad = []
    for dep in list(LADDER) + [LADDER_EXTRA]:
        D2, D1 = diameter(rc2, zc2, dep), diameter(rc1, zc1, dep)
        lad.append((dep, D2, D1))
    gated = [x for x in lad if x[0] <= LADDER[-1] + 1e-12]
    # In BINS, not mm. D lives on an exact 2*r_centre lattice whose spacing is
    # 2*BIN = 4 mm, so every nonzero difference is a multiple of one bin and a
    # comparison in mm fails on 3.6e-15 of float representation. The bar
    # ("<= 1 bin") is unchanged; only the arithmetic is exact.
    dd = [int(round(abs(D1 - D2) / (2.0 * BIN))) for _, D2, D1 in gated
          if np.isfinite(D1) and np.isfinite(D2)]
    a_band = max(dd) if dd else float("nan")
    miss = [d * 1e3 for d, D2, D1 in gated if not (np.isfinite(D1) and np.isfinite(D2))]
    if not quiet:
        print(f"\n--- {tag}: D(z), fixed 10-100 mm ladder (110 mm ungated) ---")
        print(" depth[mm]   D_2mm[mm]  D_1mm[mm]   |diff|[mm]")
        for dep, D2, D1 in lad:
            f = "" if dep <= LADDER[-1] + 1e-12 else "   (ungated)"
            print(f"{dep*1e3:9.0f}   {D2*1e3:9.1f}  {D1*1e3:9.1f}   "
                  f"{abs(D1-D2)*1e3:9.1f}{f}")
    rows.append((tag, a_core, a_band, miss, rc2, zc2, zc1, nc2, nc1))
    return a_core, a_band, miss


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    q = ap.parse_args().quiet

    print("D2t — retire the mask edge with a smooth radial taper, and score SHAPE")
    print(f"criteria  CRITERION.md   sha256 {open(os.path.join(HERE,'CRITERION.sha256')).read().split()[0][:16]}...")
    print(f"item 2    TAPER_TRACE.md sha256 {open(os.path.join(HERE,'TAPER_TRACE.sha256')).read().split()[0][:16]}... (PASSED T1-T8)")
    print(f"taper flat to {R.R_FLAT*1e3:g} mm, floor by {R.R_TAPER*1e3:g} mm, "
          f"patch radius {R.R_PATCH*1e3:g} mm; bins {BIN*1e3:g} mm FIXED PHYSICAL\n")

    need = ["P_2mm", "P_1mm", "PA_2mm", "PA_1mm"]
    missing = [n for n in need if not have(n)]
    if missing:
        print(f"not finished: {missing}")
        return 1
    L = {n: load(n) for n in need}

    verdict = []

    # ---------- domain watch ----------
    print("--- domain watch: deepest column vs the domain floor ---")
    ok_dom = True
    for n in need:
        m, r, a = L[n]
        zf, _ = z_face(m, r, a)
        inp = r <= R.R_PATCH
        cl = float(zf[inp].min())
        bad = cl < FLOOR_CLEAR
        ok_dom &= not bad
        print(f"  {n:7s} deepest face {cl*1e3:6.1f} mm, recession "
              f"{(LZ-cl)*1e3:6.1f} mm, clearance {cl*1e3:5.1f} mm "
              f"{'*** BELOW 10 mm — DEPTH RESULT VOID ***' if bad else 'ok'}")
    verdict.append(("domain watch", ok_dom, ""))

    # ---------- (a) ----------
    rows = []
    a_core, a_band, miss = score("(a)  P_1mm vs P_2mm", (L["P_2mm"], L["P_1mm"]), rows, q)
    # ---------- (c) ----------
    c_core, c_band, cmiss = score("(c)  PA_1mm vs PA_2mm", (L["PA_2mm"], L["PA_1mm"]), rows, q)

    # ---------- PA - P on D(z) ----------
    print("\n--- (c) reported: PA - P on D(z), the first measurement of what "
          "annulus heating does to a wall slope ---")
    print(" depth[mm]    2 mm: P    PA     diff  |   1 mm: P    PA     diff   [mm]")
    pa_p = []
    prof = {}
    for n in need:
        m, r, a = L[n]
        zf, _ = z_face(m, r, a)
        prof[n] = profile(m, r, zf)
    for dep in list(LADDER) + [LADDER_EXTRA]:
        row = []
        for mesh in ("2mm", "1mm"):
            Dp = diameter(*prof[f"P_{mesh}"][:2], dep)
            Da = diameter(*prof[f"PA_{mesh}"][:2], dep)
            row += [Dp, Da, Da - Dp]
            if dep <= LADDER[-1] + 1e-12 and np.isfinite(Da - Dp):
                pa_p.append(int(round(abs(Da - Dp) / (2.0 * BIN))))
        print(f"{dep*1e3:9.0f}  {row[0]*1e3:8.1f} {row[1]*1e3:6.1f} {row[2]*1e3:+7.1f}  |"
              f" {row[3]*1e3:9.1f} {row[4]*1e3:6.1f} {row[5]*1e3:+7.1f}"
              f"{'   (ungated)' if dep > LADDER[-1]+1e-12 else ''}")
    pa_p_max = max(pa_p) if pa_p else float("nan")

    # ---------- (b) ----------
    print("\n--- (b) smoothness witness ---")
    ok_b = True
    for n in need:
        m, r, a = L[n]
        zf, nv = z_face(m, r, a)
        rc, zc, nc = prof[n]
        step = float(np.abs(np.diff(zc)).max()) * 1e3
        fired = np.zeros(r.shape, bool)
        fired[a[:, 1].astype(int), a[:, 2].astype(int)] = True
        never = int(((r <= R_HALF) & ~fired).sum())
        ok_b &= never == 0
        print(f"  {n:7s} largest step in z(r) between adjacent 2 mm bins = "
              f"{step:6.2f} mm ({step/(m['dz']*1e3):.1f} cells); "
              f"never-spalled columns with r <= {R_HALF*1e3:g} mm: {never} "
              f"{'ok' if never == 0 else '*** FAIL ***'}")
    verdict.append(("(b) no never-spalled column inside 26 mm", ok_b, ""))

    # ---------- (d) ----------
    print("\n--- (d) the core is not broken: r < 18 mm, block 150-200 s ---")
    ok_d = True
    for n in need:
        m, r, a = L[n]
        v = rate(r, a, (150.0, 200.0), r < R.R_FLAT)
        ref = REF_D[m["mesh"]]
        dev = (v / ref - 1.0) * 100.0
        bad = abs(dev) > 5.0
        ok_d &= not bad
        print(f"  {n:7s} {v:.4f} m/h vs recomputed reference {ref:.4f} "
              f"({dev:+.2f} %) {'*** FAIL ***' if bad else 'ok'}")
    verdict.append(("(d) core within 5 % of the recomputed reference", ok_d, ""))

    # ---------- radial rate decomposition ----------
    print("\n--- radial rate decomposition [m/h], and the mesh gap per band ---")
    for pair, lab in ((("P_2mm", "P_1mm"), "P"), (("PA_2mm", "PA_1mm"), "PA")):
        print(f"  {lab}:")
        for lo, hi, nm in BANDS:
            g = []
            for w in BLOCKS:
                v = {}
                for n in pair:
                    m, r, a = L[n]
                    v[m["mesh"]] = rate(r, a, w, (r >= lo) & (r < hi))
                g.append((v["2mm"], v["1mm"],
                          (v["1mm"] / v["2mm"] - 1.0) * 100.0 if v["2mm"] > 0 else np.nan))
            print(f"    {nm:12s} " + "  ".join(
                f"{w[0]:.0f}-{w[1]:.0f}: {a_:.3f}/{b:.3f} ({c:+5.1f}%)"
                for w, (a_, b, c) in zip(BLOCKS, g)))

    # ---------- (e) ----------
    print("\n--- (e) the latch (REPORTED, NOT GATED): above-plane columns at "
          "each block end, and plane crossings in the block ---")
    for n in need:
        m, r, a = L[n]
        print(f"  {n:7s} " + "  ".join(
            f"{lo:.0f}-{hi:.0f}: {ab} above, {cr} crossings"
            for lo, hi, ab, cr in latch(m, r, a)))

    # ---------- verdict ----------
    print("\n=== VERDICT ===")
    ok_a = (a_core <= 4.0) and (a_band <= 1) and not miss
    ok_c = (c_core <= 4.0) and (c_band <= 1) and not cmiss
    print(f"(a) SCORED  core max |z_1mm - z_2mm| over r < 18 mm = {a_core:.2f} mm "
          f"(bar 4)  {'PASS' if a_core <= 4.0 else 'FAIL'}")
    print(f"            band max |D_1mm - D_2mm| on the 10-100 mm ladder = "
          f"{a_band:d} bin ({a_band * 2.0 * BIN * 1e3:.0f} mm in diameter) "
          f"(bar 1 bin)  {'PASS' if a_band <= 1 else 'FAIL'}"
          + (f"   [{len(miss)} ladder depths unreachable at one mesh: "
             f"{[f'{d:.0f}' for d in miss]}]" if miss else ""))
    print(f"(c) WALL    core {c_core:.2f} mm (bar 4), band {c_band:d} bin "
          f"(bar 1)  {'PASS' if ok_c else 'FAIL'}")
    print(f"    PA - P  max |D_PA - D_P| on the ladder = {pa_p_max:d} bin "
          f"({pa_p_max * 2.0 * BIN * 1e3:.0f} mm) -> P2 "
          f"{'HELD' if pa_p_max >= 1 else 'REFUTED'} (refuted under 1 bin)")
    print("    *** the D(z) lattice spacing IS 1 bin (4 mm in diameter), so this "
          "metric\n        cannot resolve a difference finer than its own bar. "
          "'1 bin' means\n        ADJACENT BIN, not '4.0 mm measured'. Carried as a "
          "caveat, not a result.")
    for nm, ok, _ in verdict:
        print(f"{nm}: {'PASS' if ok else 'FAIL'}")
    print(f"\n(a) overall: {'PASS' if ok_a else 'FAIL'}")
    if not ok_a and ok_d and ok_b:
        print("  -> P6 reading: the core passes and (b) is clean, so the band's "
              "failure says h STEPPING AT THE NOZZLE PLANE CANNOT REPRESENT A WALL.\n"
              "     The successor smooths h in s over several cells. The taper is "
              "NOT re-tuned to make this pass.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
