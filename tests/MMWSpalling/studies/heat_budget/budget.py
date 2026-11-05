#!/usr/bin/env python3
"""Heat budget (README.md): rock-side hot surface + literature gas side vs the measured uptake gaps. Offline."""
import math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p2_7_hoeser_heat"))
import energy_bound as EB

T0, T_ON = EB.T0, 850.0                       # model onset (inside her IR 430-600 degC); floor for T_s (generous)
T_MELT = (1173.15, 1323.15)                   # Hoeser ">900 degC"; Wilkinson "~1050 degC"
# Wilkinson & Tester 1993 Table 3 (V_dr m/h, steady spalling-surface T degC; IR), mean per speed
W3 = {"Barre": [(2.2, (475, 520)), (2.5, (520, 520, 508, 520)), (2.8, (710, 700)), (3.1, (879, 850, 930, 890, 785)),
                (3.3, (937, 937, 837)), (4.3, (925,))],
      "Westerly": [(2.5, (475, 520)), (2.9, (600,)), (5.6, (738,))]}
GAS = {"cavity vortex (Terekhov 2013)": (1.0, 2.0), "hot-wall gas properties": (1.0, 1.3),
       "blowing by spalls": (1.0, 1.18), "emissivity/catalysis (Baukal 1997)": (1.0, 1.12),
       "particles, roughness (mostly enhance)": (1.0, 1.0)}
def cool_factor(shape, mcp=13.6e-3 * 1250.0, Ts=850.0, n=4000):
    """Delivered over Ø40 with the wall-jet gas COOLING as it gives up heat (all 13.6 g/s flowing outward,
    dT_g/dr = -q 2 pi r / (m c_p)), over delivered at a uniform T_f (P2-7's assumption). <= 1."""
    g_disc = EB.disc_mean(EB.goldstein, EB.R_DISC)
    h = (lambda r: EB.H_DISC / g_disc * EB.goldstein(r)) if shape == "goldstein" else (lambda r: EB.H_DISC)
    Tg, P, dr = EB.T_F, 0.0, EB.R_OPEN / n
    for i in range(n):
        r = (i + 0.5) * dr; dq = h(r) * (Tg - Ts) * 2 * math.pi * r * dr
        P += dq; Tg -= dq / mcp
    return P / EB.P_deliv(Ts, shape)
FC = min(cool_factor("uniform"), cool_factor("goldstein"))
GAS["gas cools across the floor (energy balance)"] = (1.0, 1.0 / FC)
G_HI = math.prod(v[1] for v in GAS.values()); G_OK = 2.0

def wilk_T(rock, v):
    """Steady spalling-surface T [K] at drilling speed v [m/h]: piecewise-linear in the means; below the slowest point
    the surface is at onset (returns None); above the fastest, held (Barre: Wilkinson p.57 says ~melting by SOD 22-26)."""
    pts = [(V, sum(T) / len(T)) for V, T in W3[rock]]
    if v < pts[0][0]: return None
    for (v1, t1), (v2, t2) in zip(pts, pts[1:]):
        if v <= v2: return 273.15 + t1 + (t2 - t1) * (v - v1) / (v2 - v1)
    return 273.15 + pts[-1][1]

def g_req_hoeser(T, V, A, shape):
    return EB.P_deliv(T, shape) * EB.T_EXP / sum(EB.E_terms(T, V, A))

def t_close(V, A, shape, G):
    f = lambda T: EB.P_deliv(T, shape) * EB.T_EXP / G - sum(EB.E_terms(T, V, A))
    lo, hi = T0 + 1, EB.T_F - 1e-6
    if f(lo) < 0: return T0
    for _ in range(200):
        m = 0.5 * (lo + hi); lo, hi = (m, hi) if f(m) > 0 else (lo, m)
    return 0.5 * (lo + hi)

def rng(xs): return min(xs), max(xs)
out = []; p = out.append
p("HEAT BUDGET -- rock-side hot surface + literature gas side vs the measured uptake gaps (offline; README.md)")
p("=" * 112)
p("A. HOESER: gas-side cut still needed, G_req(T_s) = delivered / E_max, if the flame-facing surface is at T_s")
p("   (P2-7 bound, generous to the rock: flakes leave at T_s, conduction from a face held at T_s, radiation at T_s)")
HC = {"0 MPa": (12.9e-6, (2.0e-3, 2.5e-3)), "47 MPa": (33.0e-6, (4.0e-3,))}
p(f"   {'load':7s} | " + " | ".join(f"T_s {T:6.0f} K" for T in (850, 1000, 1100) + T_MELT))
for L, (V, As) in HC.items():
    row = []
    for T in (850, 1000, 1100) + T_MELT:
        g = [g_req_hoeser(T, V, A, s) for A in As for s in ("uniform", "goldstein")]
        row.append("%4.2f-%4.2f" % rng(g))
    p(f"   {L:7s} | " + " | ".join(f"{r:>12s}" for r in row))
p(f"   (1,173 K = Hoeser's '>900 degC' melting; 1,323 K = Wilkinson's ~1,050 degC)")
p("")
p("B. WHERE EACH CASE SITS ON WILKINSON'S MEASURED SPEED -> STEADY SURFACE-T RELATION (Table 3, IR)")
CASES = {"Hoeser 0 MPa (axis)": (1.2, 2.0), "Hoeser 47 MPa (axis)": (3.6, 3.6), "Meier (Fig 8.8)": (1.28, 1.5),
         "Rauenzahn quarry (Barre)": (4.6, 7.3)}
TS = {}
for c, (v1, v2) in CASES.items():
    s = []
    for rock in ("Barre", "Westerly"):
        a, b = wilk_T(rock, v1), wilk_T(rock, v2)
        TS[(c, rock)] = (a, b)
        f = lambda x: "onset (<= ~500 degC)" if x is None else f"{x:6.0f} K ({x-273.15:4.0f} degC)"
        s.append(f"{rock}: {f(a)}" + ("" if v1 == v2 else f" .. {f(b)}"))
    p(f"   {c:26s} {v1:.2f}-{v2:.2f} m/h  " + " | ".join(s))
p("   ⚠ Wilkinson's slowest flame speed is 2.2 m/h; below it his surfaces sat at onset. Grimsel is unmeasured.")
p("")
p("C/E. PER CASE: total cut needed, rock-side share (hot surface from B), gas-side cut left, verdict")
p(f"   gas side available from the literature (product of factors): 1.00-{G_HI:.2f}; 'comfortable' <= {G_OK:.1f}")
p(f"   (gas cooling: delivered/uniform-T_f = {cool_factor('uniform'):.3f} uniform, {cool_factor('goldstein'):.3f} Goldstein; low end 1.0 if the ring recirculates hot gas)")
for k, v in GAS.items(): p(f"      {k:40s} {v[0]:.2f}-{v[1]:.2f}")
def verdict(g): return "closes comfortably" if g <= G_OK else ("only at the literature's extreme" if g <= G_HI else "DOES NOT CLOSE")
res = {}
# Hoeser: by T_s from B (floored at onset 850 K: generous)
for L, c in (("0 MPa", "Hoeser 0 MPa (axis)"), ("47 MPa", "Hoeser 47 MPa (axis)")):
    V, As = HC[L]
    tot = rng([g_req_hoeser(T_ON, V, A, s) for A in As for s in ("uniform", "goldstein")])
    for rock in ("Barre", "Westerly"):
        Ts = [max(T_ON, x if x is not None else T_ON) for x in TS[(c, rock)]]
        left = rng([g_req_hoeser(T, V, A, s) for T in Ts for A in As for s in ("uniform", "goldstein")])
        res[(c, rock)] = left
        p(f"   {c:26s} {rock:8s} total {tot[0]:4.2f}-{tot[1]:4.2f} | T_s {min(Ts):5.0f}-{max(Ts):5.0f} K | rock-side x{tot[0]/left[0]:4.2f}"
          f" | gas-side left {left[0]:4.2f}-{left[1]:4.2f} -> {verdict(left[0])} (low end) / {verdict(left[1])} (high end)")
# Meier: steady drilling. Absorbable = flakes at onset x (1 + wall loss <= 25 %, Rauenzahn) + floor radiation;
# delivered = Hoeser's calorimeter h at Meier's floor standoff 83-100 mm (her L table, D 7.9, 38 kW), Goldstein over
# the floor, CAPPED by the jet enthalpy above 850 K (13.6 g/s x 1250 J/kgK, 38 kW HHV, lambda 1.1).
RC = EB.RHO * EB.CP
hL = [(8, 7.1), (10, 5.7), (11, 5.0), (13, 3.7), (14, 3.2)]
def h_of_L(L):
    for (a, ha), (b, hb) in zip(hL, hL[1:]):
        if L <= b: return ha + (hb - ha) * (L - a) / (b - a)
cap = 13.6e-3 * 1250.0 * (EB.T_F - T_ON)
g_disc = EB.disc_mean(EB.goldstein, EB.R_DISC)
dl, ab = [], []
for v_mh in (1.28, 1.5):
    for R in (0.0425, 0.045):
        A = math.pi * R * R; rate = v_mh / 3600 * A
        for wall in (1.0, 1.25):
            ab.append(RC * (T_ON - T0) * rate * wall + EB.EPS * EB.SIG_SB * (T_ON**4 - T0**4) * A)
        for s in (0.083, 0.100):
            q0 = h_of_L(s / EB.D) * 1e3 * (EB.T_F - T_ON)
            dl.append(min(cap, q0 / g_disc * EB.disc_mean(EB.goldstein, R, 20000) * A))
mt = (min(dl) / max(ab), max(dl) / min(ab))
res[("Meier", "-")] = mt
p(f"   {'Meier (Fig 8.8)':26s} absorbable {min(ab)/1e3:4.2f}-{max(ab)/1e3:4.2f} kW, delivered {min(dl)/1e3:5.2f}-{max(dl)/1e3:5.2f} kW (jet cap {cap/1e3:4.1f})"
      f" -> total {mt[0]:4.2f}-{mt[1]:4.2f} | T_s onset (B) | rock-side x1.00 | gas-side left {mt[0]:4.2f}-{mt[1]:4.2f}"
      f" -> {verdict(mt[0])} / {verdict(mt[1])}")
# Rauenzahn: his gap 3-4x at his assumed T_s ~765 K (rise ~465 K, T0 300); torch 2,830 K; cold copper ~350 K.
for rock in ("Barre", "Westerly"):
    Ts = TS[("Rauenzahn quarry (Barre)", rock)]
    g = []
    for gap in (3.0, 4.0):
        for T in Ts:
            g.append(gap * (2830 - T) / (2830 - 765) * (765 - 300) / (T - 300))
    left = rng(g); res[("Rauenzahn", rock)] = left
    p(f"   {'Rauenzahn quarry':26s} {rock:8s} total 3.00-4.00 | T_s {min(Ts):5.0f}-{max(Ts):5.0f} K | rock-side x{3.0/left[0] if False else (3.0/min(g)):4.2f}"
      f" | gas-side left {left[0]:4.2f}-{left[1]:4.2f} -> {verdict(left[0])} / {verdict(left[1])}")
p(f"   {'Kant onset (before 1st spall)':26s} total 1 (measured h passes onset) | no flakes, no cavity -> every candidate gives 1x: consistent")
p("")
p("P3. HOESER 0 MPa: the flame-facing surface T it would need for the budget to close")
for G in (1.0, G_OK, G_HI):
    Tc = [t_close(12.9e-6, A, s, G) for A in (2.0e-3, 2.5e-3) for s in ("uniform", "goldstein")]
    p(f"   gas-side cut {G:4.2f}: T_s {min(Tc):6.0f}-{max(Tc):6.0f} K ({min(Tc)-273.15:5.0f}-{max(Tc)-273.15:5.0f} degC)")
p("=" * 112)
p("SCORING vs README")
h47 = res[("Hoeser 47 MPa (axis)", "Barre")]; h47w = res[("Hoeser 47 MPa (axis)", "Westerly")]
ra = res[("Rauenzahn", "Barre")]; raw = res[("Rauenzahn", "Westerly")]
p(f"   P1 fast cases close (G_req <= 2): Hoeser 47 MPa Barre {h47[0]:.2f}-{h47[1]:.2f}, Westerly {h47w[0]:.2f}-{h47w[1]:.2f};"
  f" Rauenzahn Barre {ra[0]:.2f}-{ra[1]:.2f}, Westerly {raw[0]:.2f}-{raw[1]:.2f}")
h0 = res[("Hoeser 0 MPa (axis)", "Barre")]
p(f"   P2 slow cases do not close (> {G_HI:.2f}): Hoeser 0 MPa {h0[0]:.2f}-{h0[1]:.2f}; Meier {mt[0]:.2f}-{mt[1]:.2f}")
txt = "\n".join(out); print(txt); open(os.path.join(HERE, "BUDGET.txt"), "w").write(txt + "\n")
