#!/usr/bin/env python3
"""P2-31: the face gate compared COMPONENT BY COMPONENT (face_<a> vs face_<b>, by the md5 lists).
For every path whose md5 differs: a plotfile file -> the plotfile's components, FAB by FAB (bit-identical or the max
|d|), components only on one side; thermo.dat -> the columns that moved (the wall-clock timers excluded), gone, new;
a CSV -> the rows that differ (count) and the columns where they differ.
  usage: components.py <a> <b> [--legs L1,L2,...]      prints a per-leg table; exit 0
`diff(a, b)` returns {leg: {"plot": {comp: maxabs}, "plot_gone": set, "plot_new": set, "thermo": {col: maxabs},
"thermo_gone", "thermo_new", "csv": {file: (rows_a, rows_b, n_rows_differing, cols)}, "other": [paths]}}."""
import csv
import importlib.machinery
import importlib.util
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
W = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "tree_sweep_0922", "witness")
TIMERS = {"fr_t_removal", "fr_t_step", "fl_t_path"}


def _load(path, name):
    ld = importlib.machinery.SourceFileLoader(name, path)
    sp = importlib.util.spec_from_loader(name, ld)
    m = importlib.util.module_from_spec(sp)
    ld.exec_module(m)
    return m


AC = _load(os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_28_column_path", "audit_components.py"), "ac31")


def md5s(label):
    return AC.load_list(os.path.join(W, "face_%s.md5" % label))


def leg_of(rel):
    parts = rel.split("/")
    return "/".join(parts[:2]) if parts[0] in ("fr", "fl", "fb", "fh") else parts[0]


def plot_components(old_pf, new_pf):
    """{comp: max|d|} over the components present on both sides (0.0 = bit-identical), gone, new"""
    _, no = AC.read_header(old_pf)
    _, nn = AC.read_header(new_pf)
    common = [c for c in no if c in nn]
    out = {c: 0.0 for c in common}
    fo, fn = AC.read_cell_h(old_pf), AC.read_cell_h(new_pf)
    co, cn = {}, {}
    for (a_f, a_o), (b_f, b_o) in zip(fo, fn):
        co.setdefault(a_f, AC.read_fabs(old_pf, a_f))
        cn.setdefault(b_f, AC.read_fabs(new_pf, b_f))
        a, b = co[a_f][a_o], cn[b_f][b_o]
        for c in common:
            x, y = a[2][no.index(c)], b[2][nn.index(c)]
            if not np.array_equal(x.view(np.uint64), y.view(np.uint64)):
                both = np.isfinite(x) & np.isfinite(y)
                d = float(np.max(np.abs(x[both] - y[both]))) if both.any() else float("inf")
                out[c] = max(out[c], d if d > 0 else 1e-300)     # a NaN pattern change counts as a difference
    return out, set(no) - set(nn), set(nn) - set(no)


def thermo_cols(path):
    rows = [ln.rstrip("\n").split("\t") for ln in open(path) if ln.strip()]
    h = [x.strip() for x in rows[0]]
    cols = {c: [] for c in h}
    for r in rows[1:]:
        for c, v in zip(h, r):
            try:
                cols[c].append(float(v))
            except ValueError:
                cols[c].append(float("nan"))
    return cols


def csv_diff(pa, pb):
    A = list(csv.reader(open(pa)))
    B = list(csv.reader(open(pb)))
    cols = set()
    nd = 0
    for ra, rb in zip(A[1:], B[1:]):
        if ra != rb:
            nd += 1
            cols |= {A[0][i] for i, (x, y) in enumerate(zip(ra, rb)) if x != y and i < len(A[0])}
    return len(A) - 1, len(B) - 1, nd, cols


def diff(a, b, legs=None):
    A, B = md5s(a), md5s(b)
    Wa, Wb = os.path.join(W, "face_" + a), os.path.join(W, "face_" + b)
    out = {}

    def L(rel):
        lg = leg_of(rel)
        return out.setdefault(lg, {"plot": {}, "plot_gone": set(), "plot_new": set(), "thermo": {}, "thermo_gone": set(),
                                   "thermo_new": set(), "csv": {}, "other": [], "only_a": [], "only_b": []})
    for rel in sorted(set(A) - set(B)):
        L(rel)["only_a"].append(rel)
    for rel in sorted(set(B) - set(A)):
        L(rel)["only_b"].append(rel)
    done_pf = set()
    for rel in sorted(set(A) & set(B)):
        if A[rel] == B[rel]:
            continue
        lg = leg_of(rel)
        if legs and lg not in legs:
            continue
        e = L(rel)
        pf = AC.pf_root(rel)
        if pf:
            if pf in done_pf:
                continue
            done_pf.add(pf)
            comps, gone, new = plot_components(os.path.join(Wa, pf), os.path.join(Wb, pf))
            for c, d in comps.items():
                if d > 0:
                    e["plot"][c] = max(e["plot"].get(c, 0.0), d)
            e["plot_gone"] |= gone
            e["plot_new"] |= new
        elif rel.endswith("thermo.dat"):
            ca, cb = thermo_cols(os.path.join(Wa, rel)), thermo_cols(os.path.join(Wb, rel))
            for c in ca:
                if c in TIMERS or c not in cb:
                    continue
                xa, xb = ca[c], cb[c]
                if len(xa) != len(xb) or any(not (x == y or (x != x and y != y)) for x, y in zip(xa, xb)):
                    fin = [abs(x - y) for x, y in zip(xa, xb) if math.isfinite(x) and math.isfinite(y)]
                    e["thermo"][c] = max(fin) if fin else float("inf")
            e["thermo_gone"] |= set(ca) - set(cb)
            e["thermo_new"] |= set(cb) - set(ca)
        elif rel.endswith(".csv"):
            e["csv"][rel.split("/")[-1]] = csv_diff(os.path.join(Wa, rel), os.path.join(Wb, rel))
        else:
            e["other"].append(rel)
    return out


def substantive(e):
    return bool(e["plot"] or e["plot_gone"] or e["plot_new"] or e["thermo"] or e["thermo_gone"] or e["thermo_new"]
                or e["csv"] or e["other"] or e["only_a"] or e["only_b"])


def describe(lg, e):
    s = []
    if e["plot"]:
        s.append("plot components moved: " + ", ".join("%s(%.3g)" % kv for kv in sorted(e["plot"].items())))
    if e["plot_gone"] or e["plot_new"]:
        s.append("plot components gone %s, new %s" % (sorted(e["plot_gone"]), sorted(e["plot_new"])))
    if e["thermo"]:
        s.append("thermo columns moved: " + ", ".join("%s(%.3g)" % kv for kv in sorted(e["thermo"].items())))
    if e["thermo_gone"] or e["thermo_new"]:
        s.append("thermo columns gone %s, new %s" % (sorted(e["thermo_gone"]), sorted(e["thermo_new"])))
    for f, (ra, rb, nd, cols) in sorted(e["csv"].items()):
        s.append("%s: rows %d / %d, %d differ in %s" % (f, ra, rb, nd, sorted(cols)[:12]))
    if e["other"]:
        s.append("other %s" % e["other"][:4])
    if e["only_a"] or e["only_b"]:
        s.append("only-in-a %d, only-in-b %d" % (len(e["only_a"]), len(e["only_b"])))
    return "  %-12s %s" % (lg, "; ".join(s) if s else "timers only")


def main():
    a, b = sys.argv[1], sys.argv[2]
    legs = set(sys.argv[sys.argv.index("--legs") + 1].split(",")) if "--legs" in sys.argv else None
    D = diff(a, b, legs)
    print("face_%s vs face_%s, component by component (legs with any md5 difference):" % (a, b))
    for lg in sorted(D):
        print(describe(lg, D[lg]))
    print("substantive on: %s" % (sorted(lg for lg, e in D.items() if substantive(e)) or "none"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
