#!/usr/bin/env python3
"""Reduce the relief diagnostic: Sp_field (K_I/K_Ic at the surface flaw) on the top cell layer of legs L
(local 1-D law with R = 0.531294) and M (the elastic operator, no relief), same temperature history.
  ratio(col) = sigma_op / sigma_1D = R * Sp_M / Sp_L      (same flaw, same K_Ic(T), same T)
Reports: per plotfile step, the ratio at the spot centre and its radial profile; the first step/temperature
at which any column reaches Sp >= 1 in each leg; the most tensile Sp in M (rim tension)."""
import os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_2_props"))
import fabcmp
R = 0.531294
SUF = sys.argv[1] if len(sys.argv) > 1 else ""      # "" = the 2 mm legs (Lc, M), "1" = the 1 mm twins (Lc1, M1)
DX = 0.1 / (96 if SUF == "1" else 48)
X0 = 0.05

def top_layer(plot, name, down=0):
    """{(i,j): value} on the cell layer k = kmax - down, and the step time"""
    hdr = os.path.join(plot, "Header")
    lines = open(hdr).read().split("\n")
    n = int(lines[1]); names = lines[2:2 + n]; t = float(lines[2 + n + 1])
    c = names.index(name)
    out = {}; kmax = -1
    lvl = os.path.join(plot, "Level_0")
    fabs = []
    for fn in sorted(os.listdir(lvl)):
        if fn.startswith("Cell_D_"):
            fabs += fabcmp.read_fabs_boxed(os.path.join(lvl, fn))
    kmax = max(f["hi"][2] for f in fabs)
    for f in fabs:
        lo, hi = f["lo"], f["hi"]
        if hi[2] < kmax - down or lo[2] > kmax - down: continue
        nx = hi[0] - lo[0] + 1; ny = hi[1] - lo[1] + 1
        comp = f["comp"][c]
        kk = kmax - down - lo[2]
        for j in range(ny):
            for i in range(nx):
                out[(lo[0] + i, lo[1] + j)] = comp[kk * nx * ny + j * nx + i]
    return out, t

def steps(leg):
    d = os.path.join(HERE, leg, "run")
    return sorted(p for p in os.listdir(d) if p.endswith("cell"))

def main():
    # the ratio is read Lc (local law on the CELL temperature, cell Robin form) vs M (the operator, same form);
    # L (the production form: reconstructed surface, surface Robin) gives the production onset for reference
    L = {s: top_layer(os.path.join(HERE, "Lc" + SUF, "run", s), "Sp_field") for s in steps("Lc" + SUF)}
    M = {s: top_layer(os.path.join(HERE, "M" + SUF, "run", s), "Sp_field") for s in steps("M" + SUF)}
    TL = {s: top_layer(os.path.join(HERE, "Lc" + SUF, "run", s), "Temp")[0] for s in steps("Lc" + SUF)}   # the cell T (no reconstruction on these legs)
    T1 = {s: top_layer(os.path.join(HERE, "Lc" + SUF, "run", s), "Temp", 1)[0] for s in steps("Lc" + SUF)}  # one cell down
    if os.path.isdir(os.path.join(HERE, "L", "run")):
        for s in steps("L"):
            f, t = top_layer(os.path.join(HERE, "L", "run", s), "Sp_field")
            if max(f.values()) >= 1.0:
                p = max(f, key=lambda q: f[q]); Tp = top_layer(os.path.join(HERE, "L", "run", s), "T_surf_field")[0][p]
                print("production form L: first Sp >= 1 at step %s t %.3f s, Sp %.3f, T_surf %.1f K" % (s, t, f[p], Tp)); break
    common = sorted(set(L) & set(M))
    out = ["relief diagnostic -- Hoeser flame on a full 100 mm block, spot centred, %.2f mm cells, cell Robin form, cell T" % (DX * 1000),
           "ratio = sigma_operator / sigma_1D-law = R * Sp_M / Sp_L, R = %.6f" % R, ""]
    out.append("%-10s %7s | %-9s %-9s %-9s %-7s | %-8s %-8s %-8s %-8s | %-9s %-7s" % (
        "step", "t [s]", "Sp_L ctr", "Sp_M ctr", "T_s ctr", "ell/dx", "ratio r0", "r=10mm", "r=20mm", "r=30mm", "min Sp_M", "at r"))
    first = {"L": None, "M": None}
    for s in common:
        l, t = L[s]; m, _ = M[s]; T = TL[s]
        ij = sorted(l, key=lambda p: (abs((p[0] + .5) * DX - X0), abs((p[1] + .5) * DX - X0)))
        c = ij[0]
        def ratio_at(r):
            vals = []
            for p in l:
                rr = math.hypot((p[0] + .5) * DX - X0, (p[1] + .5) * DX - X0)
                if abs(rr - r) < DX * 0.75 and abs(l[p]) > 1e-12:
                    vals.append(R * m[p] / l[p])
            return sum(vals) / len(vals) if vals else float("nan")
        rat0 = R * m[c] / l[c] if abs(l[c]) > 1e-12 else float("nan")
        pm = min(m, key=lambda p: m[p])
        Tt, Tb = T[c], T1[s][c]
        ellx = 1.0 / math.log((Tt - 293.15) / (Tb - 293.15)) if (Tt > 293.15 and Tb > 293.15 and Tt > Tb) else float("nan")
        out.append("%-10s %7.3f | %9.4f %9.4f %9.2f %7.3f | %8.3f %8.3f %8.3f %8.3f | %9.4f %7.1f" % (
            s, t, l[c], m[c], T[c], ellx, rat0, ratio_at(0.010), ratio_at(0.020), ratio_at(0.030), m[pm],
            1000 * math.hypot((pm[0] + .5) * DX - X0, (pm[1] + .5) * DX - X0)))
        for nm, f in (("L", l), ("M", m)):
            if first[nm] is None and max(f.values()) >= 1.0:
                p = max(f, key=lambda q: f[q]); first[nm] = (s, t, f[p], T[p])
    out.append("")
    for nm in ("L", "M"):
        out.append("first Sp >= 1 in %s: %s" % (nm, ("step %s t %.3f s, Sp %.3f, T_surf %.1f K" % first[nm]) if first[nm] else "never within the run"))
    txt = "\n".join(out); print(txt)
    open(os.path.join(HERE, "RESULTS%s.txt" % SUF), "w").write(txt + "\n")

if __name__ == "__main__":
    main()
