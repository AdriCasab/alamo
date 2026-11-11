#!/usr/bin/env python3
"""P2-27 L7 / c1 / c2 -- THE COMPONENT AUDIT of a deregistration: two gate lists, every differing plotfile read
component by component.

  usage: audit_components.py archive <tree-label>          copy every plotfile the tree list hashes (its Cell_D files
                                                           plus each plotfile's Header and Level_0/Cell_H) into
                                                           tree_archive_<label>/ BEFORE the next sweep overwrites them
         audit_components.py tree <old> <new> [--removed D,kappa_eff,k_eff]
         audit_components.py face <old> <new> [--removed ...]
PASS iff
  * every non-plotfile path (thermo.dat, the CSVs) is byte-identical, and nothing is only on one side;
  * every plotfile Cell_D that differs belongs to a plotfile whose component list lost EXACTLY the removed set (in the
    order kept), and every remaining component's data is bit-identical (each FAB, each value);
  * every plotfile that did NOT lose the set (a node plotfile: D / kappa_eff / k_eff are cell fields) is identical;
  * a moved plotfile's Header differs only in the component count and the removed names; its Cell_H keeps the box
    list, and every KEPT component's per-FAB min and max (as printed) are unchanged (the FabOnDisk offsets move: the
    FABs are shorter).
Old tree plotfiles come from tree_archive_<old>/ (the sweep rewrites the live outputs); new ones are the live outputs.
Face plotfiles are witness/face_<label>/ on both sides.
"""
import hashlib
import os
import re
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
T = os.path.abspath(os.path.join(HERE, "..", ".."))
W = os.path.join(T, "studies", "tree_sweep_0922", "witness")
REMOVED = ["D", "kappa_eff", "k_eff"]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_list(path):
    out = {}
    for ln in open(path):
        h, p = ln.rstrip("\n").split(None, 1)
        out[p.strip()] = h
    return out


def pf_root(rel):
    """'.../00012cell/Level_0/Cell_D_00003' -> '.../00012cell' (None for a non-plotfile path)"""
    m = re.match(r"(.*/\d+(?:cell|node))/Level_0/Cell_D_\d+$", rel)
    return m.group(1) if m else None


def read_header(pf):
    L = open(os.path.join(pf, "Header")).read().split("\n")
    n = int(L[1])
    return L, L[2:2 + n]


def read_cell_h(pf):
    """[(Cell_D file, offset)] in FAB order, from Level_0/Cell_H"""
    txt = open(os.path.join(pf, "Level_0", "Cell_H")).read()
    return [(m.group(1), int(m.group(2))) for m in re.finditer(r"FabOnDisk: (\S+) (\d+)", txt)]


FAB = re.compile(rb"FAB \(\(8, \(([\d ]+)\)\),\(8, \(([\d ]+)\)\)\)\(\(([-\d,]+)\) \(([-\d,]+)\) \(([\d,]+)\)\) (\d+)\n")


def read_fabs(pf, fname):
    """{offset: (npts, ncomp, data[ncomp, npts])} for every FAB in one Cell_D file"""
    raw = open(os.path.join(pf, "Level_0", fname), "rb").read()
    out, pos = {}, 0
    while pos < len(raw):
        m = FAB.match(raw, pos)
        if not m:
            raise ValueError("%s/%s: no FAB header at byte %d" % (pf, fname, pos))
        order = m.group(2).split()
        lo = [int(x) for x in m.group(3).split(b",")]
        hi = [int(x) for x in m.group(4).split(b",")]
        nc = int(m.group(6))
        npts = int(np.prod([h - l + 1 for l, h in zip(lo, hi)]))
        dt = np.dtype("<f8") if order == [b"8", b"7", b"6", b"5", b"4", b"3", b"2", b"1"] else np.dtype(">f8")
        start = m.end()
        data = np.frombuffer(raw, dtype=dt, count=nc * npts, offset=start).reshape(nc, npts)
        out[pos] = (npts, nc, data, (lo, hi))
        pos = start + nc * npts * 8
    return out


def audit_plotfile(old_pf, new_pf, removed):
    """('same' | 'removed-ok' | 'FAIL: why', n_fabs)"""
    Lo, no = read_header(old_pf)
    Ln, nn = read_header(new_pf)
    if no == nn:
        return "unchanged-names", 0
    gone = [c for c in no if c not in nn]
    kept = [c for c in no if c in nn]
    if sorted(gone) != sorted(removed) or kept != nn:
        return "FAIL: components %s -> %s (lost %s, want exactly %s)" % (len(no), len(nn), gone, removed), 0
    fo, fn = read_cell_h(old_pf), read_cell_h(new_pf)
    if len(fo) != len(fn):
        return "FAIL: FAB count %d vs %d" % (len(fo), len(fn)), 0
    cache_o, cache_n = {}, {}
    idx_o = [no.index(c) for c in kept]
    for (fo_f, fo_off), (fn_f, fn_off) in zip(fo, fn):
        if fo_f not in cache_o:
            cache_o[fo_f] = read_fabs(old_pf, fo_f)
        if fn_f not in cache_n:
            cache_n[fn_f] = read_fabs(new_pf, fn_f)
        a, b = cache_o[fo_f][fo_off], cache_n[fn_f][fn_off]
        if a[0] != b[0] or a[3] != b[3]:
            return "FAIL: FAB box differs (%s vs %s)" % (a[3], b[3]), 0
        if b[1] != len(nn) or a[1] != len(no):
            return "FAIL: FAB ncomp %d/%d vs header %d/%d" % (a[1], b[1], len(no), len(nn)), 0
        if not np.array_equal(a[2][idx_o].view(np.uint64), b[2].view(np.uint64)):
            bad = [kept[k] for k in range(len(kept))
                   if not np.array_equal(a[2][idx_o[k]].view(np.uint64), b[2][k].view(np.uint64))]
            return "FAIL: kept components differ bit-wise: %s" % bad, 0
    # the Header beyond the component block must be identical
    if Lo[2 + len(no):] != Ln[2 + len(nn):] or Lo[0] != Ln[0]:
        return "FAIL: Header differs outside the component list", 0
    # Cell_H: the box list identical; the per-FAB min and max blocks equal on every KEPT component (as printed)
    ho, hn = mm_blocks(old_pf), mm_blocks(new_pf)
    if ho[0] != hn[0]:
        return "FAIL: Cell_H box list differs", 0
    for blk_o, blk_n, tag in ((ho[1], hn[1], "min"), (ho[2], hn[2], "max")):
        if len(blk_o) != len(blk_n):
            return "FAIL: Cell_H %s block length" % tag, 0
        for ro, rn in zip(blk_o, blk_n):
            if [ro[i] for i in idx_o] != rn:
                return "FAIL: Cell_H %s of a kept component differs" % tag, 0
    return "removed-ok", len(fo)


def mm_blocks(pf):
    """(box-list text, [min rows], [max rows]) of Level_0/Cell_H; each row is the list of printed values"""
    L = open(os.path.join(pf, "Level_0", "Cell_H")).read().split("\n")
    last = max(i for i, ln in enumerate(L) if ln.startswith("FabOnDisk:"))
    # line 3 is the component count, which a deregistration changes by design: not part of the box list
    boxes = "\n".join(ln for n, ln in enumerate(L[:last + 1]) if not ln.startswith("FabOnDisk:") and n != 2)
    rest = [ln for ln in L[last + 1:] if ln.strip()]
    nb = int(rest[0].split(",")[0])
    mins = [r.rstrip(",").split(",") for r in rest[1:1 + nb]]
    maxs = [r.rstrip(",").split(",") for r in rest[2 + nb:2 + 2 * nb]]
    return boxes, mins, maxs


TIMERS = {"fr_t_removal", "fr_t_step", "fl_t_path"}


def thermo_same(a, b):
    ra = [ln.rstrip("\n").split("\t") for ln in open(a) if ln.strip()]
    rb = [ln.rstrip("\n").split("\t") for ln in open(b) if ln.strip()]
    if len(ra) != len(rb) or [h.strip() for h in ra[0]] != [h.strip() for h in rb[0]]:
        return False
    keep = [i for i, h in enumerate(ra[0]) if h.strip() not in TIMERS]
    return all([x[i] for i in keep] == [y[i] for i in keep] for x, y in zip(ra, rb))


def archive(label):
    L = load_list(os.path.join(W, label + ".md5"))
    dst = os.path.join(HERE, "tree_archive_" + label)
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    pfs, n = set(), 0
    for rel in L:
        src = os.path.join(T, rel)
        os.makedirs(os.path.dirname(os.path.join(dst, rel)), exist_ok=True)
        shutil.copy2(src, os.path.join(dst, rel))
        n += 1
        r = pf_root(rel)
        if r:
            pfs.add(r)
    for r in sorted(pfs):
        for f in ("Header", os.path.join("Level_0", "Cell_H")):
            os.makedirs(os.path.dirname(os.path.join(dst, r, f)), exist_ok=True)
            shutil.copy2(os.path.join(T, r, f), os.path.join(dst, r, f))
    # the archive must reproduce the list it came from
    bad = [rel for rel in L if md5(os.path.join(dst, rel)) != L[rel]]
    print("archived %d listed files + the Header and Cell_H of %d plotfiles into %s; %d md5 mismatch(es)"
          % (n, len(pfs), os.path.relpath(dst, T), len(bad)))
    return 0 if not bad else 1


def audit(kind, old, new, removed):
    A = load_list(os.path.join(W, ("face_" if kind == "face" else "") + old + ".md5"))
    B = load_list(os.path.join(W, ("face_" if kind == "face" else "") + new + ".md5"))
    root_a = os.path.join(W, "face_" + old) if kind == "face" else os.path.join(HERE, "tree_archive_" + old)
    root_b = os.path.join(W, "face_" + new) if kind == "face" else T
    bad = []
    only = sorted(set(A) ^ set(B))
    for p in only[:10]:
        bad.append("only in %s: %s" % (old if p in A else new, p))
    if len(only) > 10:
        bad.append("... %d more only-on-one-side paths" % (len(only) - 10))
    diff = sorted(p for p in set(A) & set(B) if A[p] != B[p])
    pfs, other = {}, []
    for p in diff:
        r = pf_root(p) or (re.match(r"(.*/\d+(?:cell|node))/(?:Header|Level_0/Cell_H)$", p).group(1)
                           if re.match(r".*/\d+(?:cell|node)/(?:Header|Level_0/Cell_H)$", p) else None)
        if r:
            pfs.setdefault(r, []).append(p)
        else:
            other.append(p)
    for p in other:
        # thermo.dat: the face legs log WALL-CLOCK timers (fr_t_removal, fr_t_step, fl_t_path), never equal between
        # runs; every other column must be identical (face_gate.py cmp's rule)
        if p.endswith("thermo.dat") and thermo_same(os.path.join(root_a, p), os.path.join(root_b, p)):
            continue
        bad.append("non-plotfile path differs: %s" % p)
    # every plotfile in the lists: a differing one must be removed-ok; a removed-ok one must differ in every Cell_D
    status = {}
    for r in sorted(pfs):
        st, nf = audit_plotfile(os.path.join(root_a, r), os.path.join(root_b, r), removed)
        status[r] = st
        if st != "removed-ok":
            bad.append("%s: %s" % (r, st))
    n_cd = sum(1 for ps in pfs.values() for p in ps if "Cell_D" in p)
    cell_pf = sorted({pf_root(p) for p in A if pf_root(p) and pf_root(p).endswith("cell")})
    node_pf = sorted({pf_root(p) for p in A if pf_root(p) and pf_root(p).endswith("node")})
    moved_node = [r for r in pfs if r.endswith("node")]
    unmoved_cell = [r for r in cell_pf if r not in pfs]
    print("component audit %s %s -> %s (removed %s)" % (kind, old, new, removed))
    print("  %d / %d paths; %d differ: %d Cell_D of %d plotfiles, %d Header/Cell_H, %d other"
          % (len(A), len(B), len(diff), n_cd, len(pfs), len(diff) - n_cd - len(other), len(other)))
    print("  cell plotfiles: %d listed, %d moved (all removed-ok: %s); node plotfiles: %d listed, %d moved"
          % (len(cell_pf), len([r for r in pfs if r.endswith("cell")]),
             all(status[r] == "removed-ok" for r in pfs if r.endswith("cell")), len(node_pf), len(moved_node)))
    if unmoved_cell:
        print("  cell plotfiles NOT moved (no D/kappa_eff/k_eff before? they are listed below and must be checked):")
        for r in unmoved_cell[:10]:
            _, names = read_header(os.path.join(root_b, r))
            print("    %s (has %s)" % (r, [c for c in removed if c in names]))
            if any(c in names for c in removed):
                bad.append("%s still carries a removed component" % r)
    for b in bad[:30]:
        print("  FAIL  %s" % b)
    print("AUDIT: %s" % ("PASS" if not bad else "FAIL (%d)" % len(bad)))
    return 0 if not bad else 1


def main():
    a = sys.argv[1:]
    removed = REMOVED
    if "--removed" in a:
        removed = a[a.index("--removed") + 1].split(",")
    if a and a[0] == "archive":
        return archive(a[1])
    if a and a[0] in ("tree", "face"):
        return audit(a[0], a[1], a[2], removed)
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main())
