#!/usr/bin/env python3
"""AMR band gate scoring (CRITERION.md).

  gate.py identity A B      row-by-row _removal_events.csv comparison (G1)
  gate.py band  REF TEST    D2j-0b metrics on both runs, block gaps, cost (G2)
"""
import ast
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
STUD = os.path.dirname(HERE)


def load_csv(n):
    p = os.path.join(OUT, n + "_removal_events.csv")
    a = np.loadtxt(p, delimiter=",", skiprows=1, ndmin=2)
    return a if a.size else np.zeros((0, 11))


def meta(n):
    return ast.literal_eval(open(os.path.join(OUT, n + ".done")).read())


def identity(a_name, b_name):
    A, B = load_csv(a_name), load_csv(b_name)
    ma, mb = meta(a_name), meta(b_name)
    print(f"G1 identity: {a_name} ({A.shape[0]} rows, t_end {ma['t_end']:g}, {ma['wall_s']/60:.1f} min) "
          f"vs {b_name} ({B.shape[0]} rows, t_end {mb['t_end']:g}, {mb['wall_s']/60:.1f} min)")
    same_bytes = open(os.path.join(OUT, a_name + "_removal_events.csv"), "rb").read() == \
        open(os.path.join(OUT, b_name + "_removal_events.csv"), "rb").read()
    print(f"  byte-identical CSV: {same_bytes}")
    if A.shape != B.shape:
        print(f"  FAIL: row counts differ {A.shape[0]} vs {B.shape[0]}")
        # first divergence
        n = min(A.shape[0], B.shape[0])
        d = np.where(np.any(A[:n] != B[:n], axis=1))[0]
        if d.size:
            k = d[0]
            print(f"  first differing row {k}: A={A[k].tolist()}\n                        B={B[k].tolist()}")
        else:
            print(f"  first {n} rows identical; extra rows start at t = "
                  f"{(A if A.shape[0] > n else B)[n, 0]:g} s")
        return False
    exact_cols = [0, 1, 2, 3, 4, 9, 10]      # time col_i col_j regime k_top h_applied n_voided
    rel_cols = [5, 7, 8]                     # T_top Sp_top h_scan
    ok = True
    for c in exact_cols:
        bad = np.where(A[:, c] != B[:, c])[0]
        if bad.size:
            ok = False
            k = bad[0]
            print(f"  col {c}: {bad.size} rows differ, first row {k}: {A[k, c]!r} vs {B[k, c]!r}")
    for c in rel_cols:
        den = np.maximum(np.abs(A[:, c]), 1e-300)
        rel = np.abs(A[:, c] - B[:, c]) / den
        if rel.max() > 1e-10:
            ok = False
            k = int(rel.argmax())
            print(f"  col {c}: max rel diff {rel.max():.3e} at row {k}")
    print(f"  total removed volume: {A[:, 9].sum() if A.size else 0:.6e} vs "
          f"{B[:, 9].sum() if B.size else 0:.6e} (sum of h_applied)")
    print("  G1 " + ("PASS" if ok else "FAIL"))
    return ok


def band(ref, test):
    _sp = importlib.util.spec_from_file_location(
        "d2j0b_an", os.path.join(STUD, "d2j0b_plane", "analyze.py"))
    AN = importlib.util.module_from_spec(_sp)
    _sp.loader.exec_module(AN)
    AN.OUT = OUT
    mr, mt = meta(ref), meta(test)
    dr, dt_ = AN.metrics(ref), AN.metrics(test)
    print(f"G2 band: ref {ref} ({mr['wall_s']/60:.1f} min) vs test {test} ({mt['wall_s']/60:.1f} min), "
          f"wall ratio {mt['wall_s']/mr['wall_s']:.2f}")
    print(f"  v_c (r < 6 mm, 50 s..end): {dr['v_c']:.3f} vs {dt_['v_c']:.3f} m/h")
    ok = True
    print("  block          ref [m/h]   test [m/h]   gap")
    for b in sorted(set(dr["blocks"]) & set(dt_["blocks"])):
        g = dt_["blocks"][b] / dr["blocks"][b] - 1.0 if dr["blocks"][b] else float("nan")
        flag = "" if abs(g) <= 0.05 else "  <-- > 5 %"
        if b[1] <= 200.0 and abs(g) > 0.05:
            ok = False
        print(f"  {b[0]:>4.0f}-{b[1]:<4.0f}    {dr['blocks'][b]:9.3f}    {dt_['blocks'][b]:9.3f}   {g:+6.1%}{flag}")
    print("  band [mm]    t_stop ref/test    t_cross ref/test")
    for k in dr["bands"]:
        br, bt = dr["bands"][k], dt_["bands"][k]
        print(f"  {k[0]:>2}-{k[1]:<2}        {AN.fmt_t(br[0]):>4}/{AN.fmt_t(bt[0]):<4}"
              f"        {AN.fmt_t(br[1]):>4}/{AN.fmt_t(bt[1]):<4}")
    print("  r_front [mm]: " + ", ".join(
        f"t={t:g}: {AN.fmt_t(dr['front'].get(t, np.nan) * 1e3 if np.isfinite(dr['front'].get(t, np.nan)) else np.nan)}"
        f"/{AN.fmt_t(dt_['front'].get(t, np.nan) * 1e3 if np.isfinite(dt_['front'].get(t, np.nan)) else np.nan)}"
        for t in AN.FRONT_T if t in dr["front"] and t in dt_["front"]))
    print("  G2 " + ("PASS" if ok else "FAIL"))
    return ok


if __name__ == "__main__":
    mode = sys.argv[1]
    ok = identity(sys.argv[2], sys.argv[3]) if mode == "identity" else band(sys.argv[2], sys.argv[3])
    sys.exit(0 if ok else 1)
