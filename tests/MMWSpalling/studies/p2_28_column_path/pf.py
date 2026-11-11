#!/usr/bin/env python3
"""P2-28: a single-level AMReX cell plotfile read into whole-domain numpy arrays, no yt (the FAB reader is P2-27's
audit_components.py's: little-endian doubles, component-major, Fortran order within a box).

  load(pf, names=None) -> (time, {name: array[nz, ny, nx]}, (lo, hi))
  plotfiles(run_prefix_dir) -> [(step, path)] of the *cell plotfiles under a run directory, step order
  top_cells(arr3, k_top) -> the k_top slab
"""
import glob
import os
import re

import numpy as np

FAB = re.compile(rb"FAB \(\(8, \(([\d ]+)\)\),\(8, \(([\d ]+)\)\)\)\(\(([-\d,]+)\) \(([-\d,]+)\) \(([\d,]+)\)\) (\d+)\n")


def header(pf):
    L = open(os.path.join(pf, "Header")).read().split("\n")
    n = int(L[1])
    return L, L[2:2 + n], float(L[2 + n + 1])


def fabs(pf):
    """[(lo, hi, data[ncomp, npts])] over every Cell_D file of level 0"""
    out = []
    for fn in sorted(glob.glob(os.path.join(pf, "Level_0", "Cell_D_*"))):
        raw = open(fn, "rb").read()
        pos = 0
        while pos < len(raw):
            m = FAB.match(raw, pos)
            if not m:
                raise ValueError("%s: no FAB header at byte %d" % (fn, pos))
            order = m.group(2).split()
            lo = [int(x) for x in m.group(3).split(b",")]
            hi = [int(x) for x in m.group(4).split(b",")]
            nc = int(m.group(6))
            npts = int(np.prod([h - l + 1 for l, h in zip(lo, hi)]))
            dt = np.dtype("<f8") if order == [b"8", b"7", b"6", b"5", b"4", b"3", b"2", b"1"] else np.dtype(">f8")
            start = m.end()
            data = np.frombuffer(raw, dtype=dt, count=nc * npts, offset=start).reshape(nc, npts)
            out.append((lo, hi, data))
            pos = start + nc * npts * 8
    return out


def load(pf, names=None):
    _, comps, t = header(pf)
    fs = fabs(pf)
    lo = [min(f[0][d] for f in fs) for d in range(3)]
    hi = [max(f[1][d] for f in fs) for d in range(3)]
    shape = (hi[2] - lo[2] + 1, hi[1] - lo[1] + 1, hi[0] - lo[0] + 1)
    want = comps if names is None else [n for n in names if n in comps]
    out = {n: np.full(shape, np.nan) for n in want}
    for (blo, bhi, data) in fs:
        bs = (bhi[2] - blo[2] + 1, bhi[1] - blo[1] + 1, bhi[0] - blo[0] + 1)
        sl = (slice(blo[2] - lo[2], bhi[2] - lo[2] + 1), slice(blo[1] - lo[1], bhi[1] - lo[1] + 1),
              slice(blo[0] - lo[0], bhi[0] - lo[0] + 1))
        for n in want:
            c = comps.index(n)
            out[n][sl] = data[c].reshape(bs)   # Fortran order (x fastest) == C order of [z, y, x]
    return t, out, (lo, hi)


def plotfiles(run_dir_prefix):
    """the *cell plotfiles written under plot_file = <run_dir_prefix> (AMReX appends the step: <prefix>00012cell)"""
    pat = re.compile(r"(\d+)cell$")
    out = []
    for p in glob.glob(run_dir_prefix + "*cell"):
        m = pat.search(os.path.basename(p))
        if m and os.path.isdir(p):
            out.append((int(m.group(1)), p))
    return sorted(out)
