#!/usr/bin/env python3
"""P2-20 -- the report rows c1-c4 and the verdict e0, from the legs and the gate logs. Usage: reduce.py > RESULTS.txt

The axial column under the beam is read from each cell plotfile with P2-8's reader (studies/p2_8_face_list/
reduce_p2_8.py read_pf over fabcmp), the geometry from each run's RESOLVED metadata (no LZ constant, SR 5); the
column is the cell pair nearest (beam.x0, beam.y0) -- 0.05 m is a cell face at 128 cells, so the cell just past it.

  c1  Z1 and Z1w: raw and brightness (1 - R) T at 25 / 50 / 100 / 150 s vs #28c (the test's own table, its log)
  c2  Z2 (Zhang's recorded inputs: alpha 14, R 0) raw T at 50 / 100 / 150 s vs Zhang's curve 989 / 1400 / 1678 K
  c3  Z3a-c (test #19 geometry, 2 kW): the first time the PEAK (anywhere on the axis) reaches T_m = 1510.5 K,
      linear in t between the 10 s plotfiles; the depth of the peak and the surface (top-cell) T then
  c4  Z1 at 150 s: the peak and its depth; the FACE-referenced e-fold depth of the axial rise (the P2-19 helper's
      definition, re-implemented here for a report: the top cell's rise at z = 0 as the face sample, ln-linear
      crossing down the cell centres) and the PEAK-referenced one, vs 1/alpha
"""
import importlib.machinery
import importlib.util
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
STUD = os.path.join(ROOT, "tests", "MMWSpalling", "studies")
LEGS = os.path.join(HERE, "legs")
T_M, T0 = 1510.5, 297.5
ZHANG = {50: 989.0, 100: 1400.0, 150: 1678.0}


def load(name, path):
    keep = sys.argv
    sys.argv = [path]
    sys.path.insert(0, os.path.dirname(path))
    ld = importlib.machinery.SourceFileLoader(name, path)
    m = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, ld))
    ld.exec_module(m)
    sys.argv = keep
    return m


P8 = load("r8_p220", os.path.join(STUD, "p2_8_face_list", "reduce_p2_8.py"))


def meta(d):
    kv = {}
    for ln in open(os.path.join(d, "run", "metadata")):
        if " = " in ln:
            k, v = ln.rstrip("\n").split(" = ", 1)
            kv[k.strip()] = v.strip()
    return kv


def geom(kv):
    nx, ny, nz = (int(v) for v in kv["amr.n_cell"].split())
    lo = [float(v) for v in kv["geometry.prob_lo"].split()]
    hi = [float(v) for v in kv["geometry.prob_hi"].split()]
    return dict(nx=nx, ny=ny, nz=nz, dx=(hi[0] - lo[0]) / nx, dz=(hi[2] - lo[2]) / nz, lo=lo, hi=hi)


def plot_time(d, pf):
    hdr = open(os.path.join(d, "run", pf, "Header")).read().split("\n")
    nvar = int(hdr[1])
    return float(hdr[2 + nvar + 1])            # Header: version, nvar, names..., dim, time


def columns(leg):
    """[(t, z from the top face [m] per cell centre (top first), T)] for every cell plotfile"""
    d = os.path.join(LEGS, leg)
    kv = meta(d)
    g = geom(kv)
    i = min(g["nx"] - 1, int((float(kv["beam.x0"]) - g["lo"][0]) / g["dx"]))
    j = min(g["ny"] - 1, int((float(kv["beam.y0"]) - g["lo"][1]) / g["dx"]))
    out = []
    for pf in sorted(p for p in os.listdir(os.path.join(d, "run")) if re.fullmatch(r"\d{5}cell", p)):
        Tf = P8.read_pf(d, pf, g, comps=("Temp",))["Temp"]
        col = Tf[i, j, ::-1]                       # top cell first
        z = (np.arange(g["nz"]) + 0.5) * g["dz"]
        out.append((plot_time(d, pf), z, col))
    return out, kv, g


def efold(z, rise, start, r0):
    """first depth below `start` (index of the reference sample at depth z0) where the rise falls to r0/e;
    ln-linear between bracketing samples, linear when the lower one is <= 0 (the P2-19 definition)"""
    tgt = r0 / math.e
    za, ra = (0.0, r0) if start < 0 else (z[start], r0)
    for k in range(start + 1, len(z)):
        zb, rb = z[k], rise[k]
        if rb <= tgt:
            if rb > 0:
                return za + (zb - za) * math.log(ra / tgt) / math.log(ra / rb)
            return za + (zb - za) * (ra - tgt) / (ra - rb)
        za, ra = zb, rb
    return float("nan")


def log_table(path, rows=(25, 50, 100, 150)):
    if not os.path.exists(path):
        return {}
    out = {}
    for ln in open(path).read().splitlines():
        t = ln.split()
        if len(t) >= 3 and t[0].isdigit() and int(t[0]) in rows:
            out[int(t[0])] = (float(t[1]), float(t[2]))
    return out


def C1C2():
    L = ["c1/c2 -- surface-centre T [K] (raw | brightness (1 - R) T) vs #28c uncorrected (x0.7) / as printed (/0.7), "
         "and Z2 vs Zhang's curve", "=" * 100]
    ref = {50: 1410.0, 100: 1630.0, 150: 1716.0}
    for leg in ("Z1", "Z1w", "Z2"):
        tab = log_table(os.path.join(HERE, "TEST_%s.log" % leg))
        if not tab:
            L.append("   %-4s (no test log)" % leg)
            continue
        cells = []
        for t in (25, 50, 100, 150):
            raw, br = tab.get(t, (float("nan"), float("nan")))
            extra = ""
            if t in ref:
                extra = " vs %.0f/%.0f (raw/uncorr %.2f)" % (0.7 * ref[t], ref[t], raw / (0.7 * ref[t]))
            if leg == "Z2" and t in ZHANG:
                extra += " vs Zhang %.0f (raw/Zhang %.2f)" % (ZHANG[t], raw / ZHANG[t])
            cells.append("t %3d: %6.1f | %6.1f%s" % (t, raw, br, extra))
        L.append("   %-4s %s" % (leg, ";  ".join(cells)))
    return "\n".join(L)


def C3():
    L = ["c3 -- test #19 (virgin Barre, 2 kW, 50 mm): the first time the axial PEAK reaches T_m = %.1f K" % T_M,
         "=" * 100]
    res = {}
    for leg in ("Z3a", "Z3b", "Z3c"):
        if not os.path.exists(os.path.join(LEGS, leg, ".done")):
            L.append("   %-4s (not run)" % leg)
            continue
        cols, kv, g = columns(leg)
        prev = None
        on = None
        for t, z, T in cols:
            ip = int(np.argmax(T))
            if T[ip] >= T_M:
                if prev is None:
                    on = (t, z[ip], T[0], T[ip])
                else:
                    tp, Tp = prev
                    f = (T_M - Tp) / (T[ip] - Tp)
                    on = (tp + f * (t - tp), z[ip], T[0], T[ip])
                break
            prev = (t, T[ip])
        last = cols[-1]
        res[leg] = on
        L.append("   %-4s w0 %s m: %s;  at %.0f s peak %.1f K at %.2f mm, top cell %.1f K" % (
            leg, kv["beam.omega0"],
            "onset t %.1f s, the peak %.2f mm below the top face (top cell %.1f K)" % (on[0], 1e3 * on[1], on[2])
            if on else "the peak never reaches T_m", last[0], float(np.max(last[2])),
            1e3 * float(last[1][int(np.argmax(last[2]))]), float(last[2][0])))
    L.append("   #19 (DOE p. 54 Fig. 34): melt onset ~300 +- 50 s after beam-on at ~2 kW; 'the melt initially starts "
             "2-4 mm below the surface' (p. 53). P5: onset in [100, 400] s across the bracket, peak 1-5 mm deep")
    return "\n".join(L), res


def C4():
    cols, kv, g = columns("Z1")
    t, z, T = cols[-1]
    er, td, lam = float(kv["beam.eps_r"]), float(kv["beam.tan_delta"]), float(kv["beam.lambda"])
    alpha = 2.0 * math.pi * math.sqrt(er) * (td / lam)
    rise = T - T0
    ip = int(np.argmax(T))
    e_face = efold(z, rise, -1, rise[0])
    e_pk = efold(z, rise, ip, rise[ip]) if ip > 0 else float("nan")
    return ("c4 -- Z1 at %.0f s: peak %.1f K at %.2f mm (top cell %.1f K); face-referenced e-fold depth %s; "
            "peak-referenced (absolute depth) %s; 1/alpha = %.1f mm, block 25 mm" % (
                t, T[ip], 1e3 * z[ip], T[0],
                "%.2f mm" % (1e3 * e_face) if math.isfinite(e_face) else "UNDEFINED inside the block",
                "%.2f mm" % (1e3 * e_pk) if math.isfinite(e_pk) else "UNDEFINED inside the block", 1e3 / alpha))


def gate(path, pat):
    if not os.path.exists(path):
        return False
    return re.search(pat, open(path).read(), re.M) is not None


def main():
    out = []
    for f in ("Z0.txt", "Z1_retired_keys.log", "Z2_beam.log", "Z3_check_shared.txt", "GATE_tree.txt",
              "GATE_face_cmp.txt", "A2_beam_void_closure.log", "TEST_Z1.log", "TEST_Z1x.log"):
        p = os.path.join(HERE, f)
        out += ["== %s" % f, open(p).read().rstrip() if os.path.exists(p) else "NOT ON DISK", ""]
    out += [C1C2(), ""]
    c3, r3 = C3()
    out += [c3, "", C4(), ""]
    g = {
        "z1": gate(os.path.join(HERE, "Z1_retired_keys.log"), r"^PASS"),
        "z2": gate(os.path.join(HERE, "Z2_beam.log"), r"PASS"),
        "z3": gate(os.path.join(HERE, "Z3_check_shared.txt"), r"^z3 PASS"),
        "a0": gate(os.path.join(HERE, "GATE_tree.txt"), r"vs p2_18: 0 of 2639 changed, 0 only in new, 0 only in ref"),
        "a1": gate(os.path.join(HERE, "GATE_face_cmp.txt"), r"VERDICT: IDENTICAL"),
        "a2": gate(os.path.join(HERE, "A2_beam_void_closure.txt"), r"^a2 PASS"),
        "b1b2": gate(os.path.join(HERE, "TEST_Z1.log"), r"^PASS$"),
        "b3": gate(os.path.join(HERE, "TEST_Z1x.log"), r"PASS  b3"),
    }
    out.append("e0 VERDICT")
    out.append("   the deposition is permittivity-derived, closed and gated on every path: %s (%s)" % (
        "YES" if all(g.values()) else "NO", ", ".join("%s %s" % (k, "PASS" if v else "FAIL") for k, v in g.items())))
    out.append("   c3 onsets: %s" % ", ".join("%s %s" % (k, "%.0f s" % v[0] if v else "none") for k, v in r3.items()))
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
