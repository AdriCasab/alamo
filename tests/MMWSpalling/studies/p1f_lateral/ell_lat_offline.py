#!/usr/bin/env python3
"""P1f -- CAN THE LATERAL MODE FIRE AT ALL? Offline, zero compute.

The packet's own arithmetic says the lateral mode needs ell_lat <= ~8.0 mm to
reach the vertical mode's measured 77-103 K firing drive, because
drive_frac = 1 - exp(-flake_h_max/ell_lat) collapses for long ell_lat -- and it
puts diffusive penetration on a non-receding wall at 16.6 mm by 100 s. If that is
right, a soaked wall CANNOT fire and only a fresh corner band can, which is the
packet's own inertness case 8: the bore would be a statement about
spall.flake_h_max rather than about flux.

That is decidable NOW, from P1e's shipped plotfiles, before deliverables 1-3 are
written. This script measures the real lateral decay length on P1e's wall.

  ell_lat = d_normal / ln(dT_face / dT_behind),  d_normal = dx*|n_x|

with |n_x| from the riser heights P1e already uses for the staircase factor, the
face cell being the solid cell adjacent to the void and `behind` the next solid
cell inward along the same axis.
"""
import math, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import yt
yt.set_log_level(50)

T_REF = 293.15
DX = DZ = 0.002
FLAKE_CAP = 1.5e-3
C_F = 0.4
SCOPE_R = 0.040          # surface_patch.side_face_h_r_min: the 7a wall branch
# ⚠⚠ THE VIABILITY THRESHOLD IS FLAKE_CAP / C_F, NOT 8 mm.
# h_f = min(c_f*ell, flake_h_max), so drive_frac reaches its full self-similar
# 1 - e^{-c_f} = 0.3297 only while the cap is NOT binding, i.e. ell <= 1.5/0.4 =
# 3.75 mm. (8 mm was the laxer figure: the ell at which the CAPPED drive still
# matches the vertical mode's measured 77 K firing drive.) Both are reported --
# VIABLE is the full-drive threshold the packet gates on; FIRES is the laxer one.
ELL_VIABLE = FLAKE_CAP / C_F     # 3.75 mm at the shipped cap
ELL_FIRES  = 0.008               # the laxer "still reaches 77 K" figure
# P1f deliverable 1: bands BELOW 35 mm, because the scope cut is at 40 mm and a
# viable population living inside it would make Leg A inert BY SCOPE, not by soak.
BANDS = ((0.020, 0.028), (0.028, 0.035), (0.035, 0.040),
         (0.040, 0.044), (0.044, 0.050), (0.050, 0.060))


def main(path):
    ds = yt.load(path)
    cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                          dims=ds.domain_dimensions)
    T = np.asarray(cg["Temp"])
    rem = np.asarray(cg["removed"])
    nx, ny, nz = T.shape
    solid = rem <= 0.5
    print("%s   t = %.1f s   grid %dx%dx%d"
          % (os.path.basename(path), float(ds.current_time), nx, ny, nz))

    # k_top per column: highest solid k
    ktop = np.full((nx, ny), -1, dtype=int)
    for i in range(nx):
        for j in range(ny):
            s = np.nonzero(solid[i, j, :])[0]
            if s.size:
                ktop[i, j] = s.max()

    def r_of(i, j):
        return math.hypot((i + 0.5) * DX, (j + 0.5) * DX)

    rows = []
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            if ktop[i, j] < 0:
                continue
            r = r_of(i, j)
            if not (BANDS[0][0] <= r < BANDS[-1][1]):
                continue
            # every solid cell of this column that has an open lateral face
            for k in range(0, ktop[i, j] + 1):
                if not solid[i, j, k]:
                    continue
                # riser heights per axis, from the staircase inputs
                hx = hy = 0
                openf = {}
                for d, (di, dj) in enumerate(((-1, 0), (1, 0), (0, -1), (0, 1))):
                    ni, nj = i + di, j + dj
                    if ni < 0 or nj < 0 or ni >= nx or nj >= ny:
                        continue
                    n = ktop[i, j] - ktop[ni, nj]
                    if n < 1:
                        continue
                    if d < 2: hx += n
                    else:     hy += n
                    # the face is open at THIS k only if the neighbour is void here
                    if not solid[ni, nj, k]:
                        openf[d] = (ni, nj)
                if not openf:
                    continue
                hyp = math.hypot(hx, hy)
                if hyp <= 0:
                    continue
                n_x, n_y = hx / hyp, hy / hyp
                for d, (ni, nj) in openf.items():
                    # the solid cell BEHIND the face, opposite the void
                    bi, bj = (i - (ni - i), j - (nj - j))
                    if bi < 0 or bj < 0 or bi >= nx or bj >= ny:
                        continue
                    if not solid[bi, bj, k]:
                        continue
                    comp = abs(n_x) if d < 2 else abs(n_y)
                    if comp <= 1e-12:
                        continue
                    dn = DX * comp
                    dTf = T[i, j, k] - T_REF
                    dTb = T[bi, bj, k] - T_REF
                    if dTf <= 0.0 or dTb <= 0.0 or dTf <= dTb:
                        # no inward decay: a soaked or inverted face. Record as
                        # infinite ell -- this is the cap-limited population and
                        # dropping it would bias the answer optimistic.
                        rows.append((r, float("inf"), dTf, 0.0))
                        continue
                    ell = dn / math.log(dTf / dTb)
                    frac = 1.0 - math.exp(-FLAKE_CAP / ell)
                    rows.append((r, ell, dTf, frac))

    if not rows:
        print("no open lateral faces found"); return 1

    print("\n%d (cell, face) pairs with an open lateral face\n" % len(rows))
    hdr = ("%-14s %6s %8s %8s %9s %8s %8s %8s %9s"
           % ("band [mm]", "n", "ell p50", "ell p90", "VIABLE", "fires", "<=2mm",
              "<dz/2", "drive p50"))
    print(hdr)
    print("-" * len(hdr))
    for lo, hi in BANDS:
        sub = [x for x in rows if lo <= x[0] < hi]
        if not sub:
            print("%-14s %6d" % ("%.0f-%.0f" % (1000*lo, 1000*hi), 0)); continue
        ells = sorted(x[1] for x in sub)
        fin = [e for e in ells if math.isfinite(e)]
        n = len(sub)
        def pct(thr):
            return 100.0 * sum(1 for e in ells if e <= thr) / n
        p50 = ells[n // 2]
        p90 = ells[min(n - 1, 9 * n // 10)]
        dr = sorted(x[3] for x in sub)
        print("%-14s %6d %8s %8s %8.1f%% %7.1f%% %7.1f%% %7.1f%% %9.3f"
              % ("%.0f-%.0f" % (1000*lo, 1000*hi), n,
                 ("%.2f" % (1000*p50)) if math.isfinite(p50) else "inf",
                 ("%.2f" % (1000*p90)) if math.isfinite(p90) else "inf",
                 pct(ELL_VIABLE), pct(ELL_FIRES), pct(0.002), pct(0.001),
                 dr[n // 2]))

    # ---- (q): WHERE does the viable population live, relative to the scope cut?
    viable = [x for x in rows if x[1] <= ELL_VIABLE]
    print("\n⚠⚠ (q) THE RUN-ORDER GATE -- is Leg A inert by SOAK or by SCOPE?")
    print("    viability threshold ell_lat <= %.2f mm (= flake_h_max / c_f, the "
          "full-drive\n    condition); scope cut r >= %.0f mm."
          % (1000*ELL_VIABLE, 1000*SCOPE_R))
    if not viable:
        print("    viable faces anywhere: 0  ⇒ the mode cannot fire on this leg "
              "at all.")
        print("    ⇒ INERT BY SOAK. Neither scope nor criterion is the binding "
              "term.")
        return 0
    inside  = [x for x in viable if x[0] <  SCOPE_R]
    outside = [x for x in viable if x[0] >= SCOPE_R]
    nv = len(viable)
    print("    viable faces: %d   inside r < %.0f mm: %d (%.1f %%)   "
          "in scope r >= %.0f mm: %d (%.1f %%)"
          % (nv, 1000*SCOPE_R, len(inside), 100.0*len(inside)/nv,
             1000*SCOPE_R, len(outside), 100.0*len(outside)/nv))
    rs = sorted(x[0] for x in viable)
    print("    their radius: p10 %.1f  p50 %.1f  p90 %.1f mm"
          % (1000*rs[len(rs)//10], 1000*rs[len(rs)//2],
             1000*rs[9*len(rs)//10]))
    frac_out = 100.0*len(outside)/nv
    if frac_out < 20.0:
        print("    ⇒ ⚠⚠ INERT BY SCOPE: the viable population lives INSIDE the "
              "40 mm cut.\n       The scope restriction, not the criterion, is "
              "the finding. Leg A would\n       fire almost nothing for a "
              "reason that has nothing to do with the thesis.")
    else:
        print("    ⇒ the scope cut keeps %.1f %% of the viable population: Leg A "
              "has faces to fire on." % frac_out)
    return 0


if __name__ == "__main__":
    p = sys.argv[1] if len(sys.argv) > 1 else (
        "/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/studies/"
        "p1e_staircase/output_fh/FH/28125cell")
    sys.exit(main(p))
