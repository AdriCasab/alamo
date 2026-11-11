#!/usr/bin/env python3
"""P2-28 b1 -- THE FACE GATE'S COMPONENT + THERMO AUDIT across the column path's deletion (face_<old> -> face_<new>).

  usage: audit_components.py face <old> <new> [--mutate]
PASS iff
  * plotfiles: each differing plotfile lost EXACTLY the declared fields it carried (FIELDS, the ones present in its old
    Header), every KEPT component bit-identical FAB by FAB, the Header identical outside the component list, the Cell_H
    box list identical and every kept component's printed min / max unchanged; a plotfile that carried none of them is
    byte-identical;
  * thermo.dat: the columns that left are EXACTLY the declared ones the old file carried (COLS), `ledger_err` is renamed
    `ledger_err_raw` with every value identical, every other column identical row by row (the wall-clock timers
    excluded), nothing new;
  * every other path (the CSVs) byte-identical; only-in-old = EXACTLY the declared CSV (FEET_S's h_col log); nothing
    only-in-new.
--mutate: the planted mutation -- the declared field set with one name dropped (Sp_field); the audit MUST then FAIL
(an audit that passes a wrong declaration is tolerant by accident).
"""
import hashlib
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
T = os.path.abspath(os.path.join(HERE, "..", ".."))
W = os.path.join(T, "studies", "tree_sweep_0922", "witness")
FIELDS = ["Sp_field", "Sp_cluster_id", "h_spall_field", "crack_event_count", "T_surf_field", "ell_eff_field"]
COLS_RX = re.compile(r"^(recon_(cols|flat|ell_clamped|T_clamped|fallback|hf_capped|hf_bounded|dT_mean)|patch_P_robin"
                     r"|surf_cols_robin|surf_q_mean|pinned_\w+)$")
RENAME = {"ledger_err": "ledger_err_raw"}
# ⚠ DECLARED DEVIATION: the packet declared FEET_S's h_col events CSV as the one CSV that leaves. It never existed on the
# face list -- only the COLUMN path wrote it (the face gate's p2_28a list has no such path) -- so nothing can leave:
# "no path only on either side" is the (stronger) bar.
GONE_PATHS = set()
TIMERS = {"fr_t_removal", "fr_t_step", "fl_t_path"}
FAB = re.compile(rb"FAB \(\(8, \(([\d ]+)\)\),\(8, \(([\d ]+)\)\)\)\(\(([-\d,]+)\) \(([-\d,]+)\) \(([\d,]+)\)\) (\d+)\n")


def load_list(path):
    out = {}
    for ln in open(path):
        if ln.strip():
            h, rel = ln.rstrip("\n").split("  ", 1)
            out[rel] = h
    return out


def pf_root(rel):
    m = re.match(r"(.*/\d+(?:cell|node))/(?:Level_0/Cell_D_\d+|Header|Level_0/Cell_H|Checkpoint)$", rel)
    return m.group(1) if m else None


def read_header(pf):
    L = open(os.path.join(pf, "Header")).read().split("\n")
    n = int(L[1])
    return L, L[2:2 + n]


def read_cell_h(pf):
    txt = open(os.path.join(pf, "Level_0", "Cell_H")).read()
    return [(m.group(1), int(m.group(2))) for m in re.finditer(r"FabOnDisk: (\S+) (\d+)", txt)]


def read_fabs(pf, fname):
    raw = open(os.path.join(pf, "Level_0", fname), "rb").read()
    out, pos = {}, 0
    while pos < len(raw):
        m = FAB.match(raw, pos)
        order = m.group(2).split()
        lo = [int(x) for x in m.group(3).split(b",")]
        hi = [int(x) for x in m.group(4).split(b",")]
        nc = int(m.group(6))
        npts = int(np.prod([h - l + 1 for l, h in zip(lo, hi)]))
        dt = np.dtype("<f8") if order[0] == b"8" else np.dtype(">f8")
        data = np.frombuffer(raw, dtype=dt, count=nc * npts, offset=m.end()).reshape(nc, npts)
        out[pos] = (npts, nc, data, (lo, hi))
        pos = m.end() + nc * npts * 8
    return out


def mm_blocks(pf):
    L = open(os.path.join(pf, "Level_0", "Cell_H")).read().split("\n")
    last = max(i for i, ln in enumerate(L) if ln.startswith("FabOnDisk:"))
    boxes = "\n".join(ln for n, ln in enumerate(L[:last + 1]) if not ln.startswith("FabOnDisk:") and n != 2)
    rest = [ln for ln in L[last + 1:] if ln.strip()]
    nb = int(rest[0].split(",")[0])
    return boxes, [r.rstrip(",").split(",") for r in rest[1:1 + nb]], \
        [r.rstrip(",").split(",") for r in rest[2 + nb:2 + 2 * nb]]


def audit_plotfile(old_pf, new_pf, fields):
    Lo, no = read_header(old_pf)
    Ln, nn = read_header(new_pf)
    want = [c for c in no if c in fields]
    gone = [c for c in no if c not in nn]
    kept = [c for c in no if c in nn]
    if gone != want or kept != nn:
        return "FAIL: lost %s, declared %s (kept order %s)" % (gone, want, kept == nn)
    if not gone:
        return "no-declared-field"
    fo, fn = read_cell_h(old_pf), read_cell_h(new_pf)
    if len(fo) != len(fn):
        return "FAIL: FAB count %d vs %d" % (len(fo), len(fn))
    co, cn = {}, {}
    idx = [no.index(c) for c in kept]
    for (a_f, a_o), (b_f, b_o) in zip(fo, fn):
        co.setdefault(a_f, read_fabs(old_pf, a_f))
        cn.setdefault(b_f, read_fabs(new_pf, b_f))
        a, b = co[a_f][a_o], cn[b_f][b_o]
        if a[3] != b[3] or a[1] != len(no) or b[1] != len(nn):
            return "FAIL: FAB box / ncomp"
        if not np.array_equal(a[2][idx].view(np.uint64), b[2].view(np.uint64)):
            return "FAIL: kept components differ: %s" % [kept[k] for k in range(len(kept))
                                                        if not np.array_equal(a[2][idx[k]].view(np.uint64),
                                                                              b[2][k].view(np.uint64))]
    if Lo[2 + len(no):] != Ln[2 + len(nn):] or Lo[0] != Ln[0]:
        return "FAIL: Header differs outside the component list"
    ho, hn = mm_blocks(old_pf), mm_blocks(new_pf)
    if ho[0] != hn[0]:
        return "FAIL: Cell_H box list"
    for bo, bn in ((ho[1], hn[1]), (ho[2], hn[2])):
        if len(bo) != len(bn) or any([ro[i] for i in idx] != rn for ro, rn in zip(bo, bn)):
            return "FAIL: Cell_H min/max of a kept component"
    return "removed-ok"


def thermo_audit(a, b):
    ra = [ln.rstrip("\n").split("\t") for ln in open(a) if ln.strip()]
    rb = [ln.rstrip("\n").split("\t") for ln in open(b) if ln.strip()]
    ha, hb = [h.strip() for h in ra[0]], [h.strip() for h in rb[0]]
    if len(ra) != len(rb):
        return "FAIL: %d vs %d rows" % (len(ra), len(rb)), []
    gone = [h for h in ha if h not in hb and RENAME.get(h) not in hb]
    new = [h for h in hb if h not in ha and h not in RENAME.values()]
    want = [h for h in ha if COLS_RX.match(h)]
    if gone != want or new:
        return "FAIL: columns gone %s (declared %s), new %s" % (gone, want, new), gone
    for h in ha:
        if h in gone or h in TIMERS:
            continue
        hn = RENAME.get(h, h) if RENAME.get(h, h) in hb else h
        i, j = ha.index(h), hb.index(hn)
        if any(x[i] != y[j] for x, y in zip(ra[1:], rb[1:])):
            return "FAIL: column %s -> %s differs" % (h, hn), gone
    return "ok", gone


def audit(old, new, fields):
    A = load_list(os.path.join(W, "face_%s.md5" % old))
    B = load_list(os.path.join(W, "face_%s.md5" % new))
    ra, rb = os.path.join(W, "face_" + old), os.path.join(W, "face_" + new)
    bad = []
    only_a = sorted(set(A) - set(B))
    only_b = sorted(set(B) - set(A))
    if set(only_a) != GONE_PATHS:
        bad.append("only in %s: %s (declared %s)" % (old, only_a[:6], sorted(GONE_PATHS)))
    if only_b:
        bad.append("only in %s: %s" % (new, only_b[:6]))
    pfs = sorted({pf_root(p) for p in A if pf_root(p)} & {pf_root(p) for p in B if pf_root(p)})
    st = {}
    for r in pfs:
        s = audit_plotfile(os.path.join(ra, r), os.path.join(rb, r), fields)
        if s == "no-declared-field":
            files = [p for p in A if pf_root(p) == r]
            s = "identical" if all(A[p] == B.get(p) for p in files) else "FAIL: no declared field yet differs"
        st[r] = s
        if s.startswith("FAIL"):
            bad.append("%s: %s" % (r, s))
    thermos = sorted(p for p in A if p.endswith("thermo.dat") and p in B)
    tcols = {}
    for p in thermos:
        s, gone = thermo_audit(os.path.join(ra, p), os.path.join(rb, p))
        tcols[p] = gone
        if s != "ok":
            bad.append("%s: %s" % (p, s))
    other = sorted(p for p in set(A) & set(B) if not pf_root(p) and not p.endswith("thermo.dat") and A[p] != B[p])
    for p in other:
        bad.append("non-plotfile path differs: %s" % p)
    legs = sorted({p.split("/run")[0] for p in A})
    print("component + thermo audit face_%s -> face_%s; declared fields %s" % (old, new, fields))
    print("  %d / %d paths; only-in-%s %s; only-in-%s %d" % (len(A), len(B), old, only_a, new, len(only_b)))
    print("  plotfiles %d: removed-ok %d, identical %d, FAIL %d" % (
        len(pfs), sum(1 for v in st.values() if v == "removed-ok"), sum(1 for v in st.values() if v == "identical"),
        sum(1 for v in st.values() if v.startswith("FAIL"))))
    for leg in legs:
        lp = [r for r in pfs if r.startswith(leg + "/") and r.endswith("cell")]
        if not lp:
            continue
        _, no = read_header(os.path.join(ra, lp[-1]))
        th = [p for p in thermos if p.startswith(leg + "/")]
        print("    %-14s fields -%s; thermo -%s" % (leg, [c for c in no if c in FIELDS],
                                                    tcols.get(th[0], []) if th else "-"))
    print("  thermo files %d (ledger_err -> ledger_err_raw where present); other paths %d differing" % (
        len(thermos), len(other)))
    for b in bad[:30]:
        print("  FAIL  %s" % b)
    print("AUDIT: %s" % ("PASS" if not bad else "FAIL (%d)" % len(bad)))
    return 0 if not bad else 1


def main():
    a = sys.argv[1:]
    if len(a) >= 3 and a[0] == "face":
        fields = [f for f in FIELDS if f != "Sp_field"] if "--mutate" in a else FIELDS
        return audit(a[1], a[2], fields)
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main())
