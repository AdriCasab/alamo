"""Radial floor-flux distributions: the present map against the options (2026-09-30).

Panel (a): flux delivered to the floor at a 100 mm axis standoff, per horizontal area.
Panel (b): the same, expressed as an EQUIVALENT flux at the present fixed 823 K flaking
temperature, so that options that change the rock response (C: T_s rising with flux,
D: spall-blowing reduction of h) can be drawn on the same axis as options that change the
map. Recession per column scales with this equivalent flux.

Numbers: map expression from studies/d2ac_feed/output/FH.log (Martin SRN local form);
gas-to-rock difference 957 K at 100 mm (D2ad C0 measured T_g 1684 K, T_s 823 K);
Rauenzahn 1986 Fig 6.16 spanwise distribution read at 0/1.5/3/6/9/12/15 nozzle radii;
Wilkinson 1989 Fig 1.8 slope for Barre granite, dT_s ~ q^(3/(m+3)) with m = 8, anchored at
823 K at Meier's 0.55 MW/m^2; Wilkinson §6.2.3 blowing factor 0.65.
"""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = 0.0075
def G(r, s):
    r = max(0.01875, r); s = min(0.09, max(0.015, s)); b = 0.1 * (s / D - 6.0)
    return (D / (2 * r)) * (r * r + 2 * b * D * r - 1.1 * b * D * D) / ((r + b * D) ** 2)
def w(r, a, bnd):
    x = min(bnd, max(a, r)); return 0.5 * (1 + math.cos(math.pi * (x - a) / (bnd - a)))
def h_map(r, s=0.10):
    M = 1449.0 * G(r, s) / 0.19763872491145218 * (0.09 / max(0.09, s))
    inner = M * w(r, 0.04, 0.056) + 142.0 * (1 - w(r, 0.04, 0.056))
    return inner * w(r, 0.06, 0.076) + 0.001 * (1 - w(r, 0.06, 0.076))

dT = 957.0
r_mm = np.linspace(0, 56, 561)
q_now = np.array([h_map(r * 1e-3) * dT / 1e6 for r in r_mm])          # MW/m^2
q_axis = q_now[0]
q_half = 0.5 * q_now                                                   # D2ad uniform scaling
q_cap50 = np.minimum(q_now, 0.50 * q_axis)                             # option A
q_cap45 = np.minimum(q_now, 0.45 * q_axis)

# Rauenzahn Fig 6.16 shape (fraction of axis vs r / R_noz), R_noz = 3.75 mm for Meier
rz_pts = np.array([0, 1.5, 3, 6, 9, 12, 15]) * 3.75
fz_pts = np.array([1.0, 0.72, 0.36, 0.28, 0.20, 0.15, 0.13])
f_rz = np.interp(r_mm, rz_pts, fz_pts)
q_rz_axis = q_axis * f_rz                                              # same axis level as now
# scaled so that the power inside 40 mm equals Meier's budget 2.8 kW
def power_inside(q, rmax=40.0):
    m = r_mm <= rmax
    return np.trapezoid(q[m] * 1e6 * 2 * np.pi * r_mm[m] * 1e-3, r_mm[m] * 1e-3)
q_rz_budget = q_rz_axis * (2800.0 / power_inside(q_rz_axis))
# Anderson & Stresino exponential, s0 = 40 mm, at the present axis level
q_as = q_axis * np.exp(-r_mm / 40.0)

# Option C + D: flux unchanged, but the rock needs dT_s(q) and h is cut by blowing.
# Equivalent flux at fixed 823 K: q_eq = 0.65 q * 530 / dT_s(q), dT_s = 530 (q/0.55)^(3/11)
def equiv_CD(q, blow=0.65, m=8):
    dTs = 530.0 * np.power(np.maximum(q, 1e-6) / 0.55, 3.0 / (m + 3))
    return blow * q * 530.0 / dTs
q_eq_CD = equiv_CD(q_now)
q_eq_C = equiv_CD(q_now, blow=1.0)

col = dict(now="#2b5c7a", half="#7a4a9c", cap="#3f8a63", rz="#b5461c", asx="#c8792a", cd="#1f7a8c")

fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), sharey=True)
for ax in axes:
    ax.axvspan(28, 40, color="0.75", alpha=0.35, lw=0, label="feet ring 28–40 mm")
    ax.axvline(43, color="0.4", ls=":", lw=1)
    ax.text(43.5, 0.03, "bore wall", fontsize=8, color="0.35")
    ax.axhspan(0.45, 0.87, color="#b5461c", alpha=0.10, lw=0)
    ax.text(0.8, 0.84, "Kant: flaking measured 0.45–0.87", fontsize=8, color="#8a3a1c", va="top")
    ax.axhspan(0.50, 0.60, color="0.3", alpha=0.18, lw=0)
    ax.text(0.8, 0.61, "Meier energy budget 0.50–0.60", fontsize=8, color="0.25", va="bottom")
    ax.axhspan(0.20, 0.25, color="#7a4a9c", alpha=0.15, lw=0)
    ax.text(0.8, 0.26, "model stops flaking (firing logs) 0.20–0.25", fontsize=8, color="#5a3a7c", va="bottom")
    ax.set_xlim(0, 56); ax.set_ylim(0, 1.3)
    ax.set_xlabel("radius from the hole axis, mm")
    ax.grid(alpha=0.25)

a = axes[0]
a.plot(r_mm, q_now, color=col["now"], lw=2.4, label="present map, full heat (D2ad C0)")
a.plot(r_mm, q_half, color=col["half"], lw=2.0, label="uniform ×0.5 (D2ad Q5): ring starves")
a.plot(r_mm, q_cap50, color=col["cap"], lw=2.4, label="A: cap at 0.50 of axis")
a.plot(r_mm, q_cap45, color=col["cap"], lw=1.6, ls="--", label="A: cap at 0.45 of axis")
a.plot(r_mm, q_rz_axis, color=col["rz"], lw=2.0, ls="-.", label="Rauenzahn CFD in-hole shape, same axis")
a.plot(r_mm, q_rz_budget, color=col["rz"], lw=1.4, ls=":", label="Rauenzahn shape scaled to 2.8 kW budget")
a.plot(r_mm, q_as, color=col["asx"], lw=1.4, ls=(0, (5, 2, 1, 2)), label="Anderson–Stresino exp(−r/40 mm), same axis")
a.set_ylabel("floor heat flux at 100 mm standoff, MW/m² (per horizontal area)")
a.set_title("(a) what the jet delivers: present map, scalings, literature shapes", fontsize=10.5, loc="left")
a.legend(fontsize=7.6, loc="upper right", framealpha=0.95)

b = axes[1]
b.plot(r_mm, q_now, color=col["now"], lw=2.4, label="present: fixed T_s = 823 K (equivalent = delivered)")
b.plot(r_mm, q_cap50, color=col["cap"], lw=2.4, label="A: cap 0.50 (map change)")
b.plot(r_mm, q_eq_C, color=col["cd"], lw=1.8, ls="--", label="C: T_s rising with flux (m = 8, Wilkinson) — equivalent")
b.plot(r_mm, q_eq_CD, color=col["cd"], lw=2.6, label="C + D: T_s(q) and spall blowing ×0.65 — equivalent")
b.set_title("(b) what drives recession: equivalent flux at a fixed 823 K flaking temperature", fontsize=10.5, loc="left")
b.legend(fontsize=7.6, loc="upper right", framealpha=0.95)

# annotate ring/axis ratios
def ratio(q):
    i40 = np.searchsorted(r_mm, 38.0); return q[i40] / q[0]
txt = ("ring(38 mm)/axis:  present %.2f · cap 0.50 %.2f · uniform ×0.5 %.2f · Rauenzahn %.2f\n"
       "                          C %.2f · C+D %.2f     (feed held when ≥ ~0.85–0.9)"
       % (ratio(q_now), ratio(q_cap50), ratio(q_half), ratio(q_rz_axis), ratio(q_eq_C), ratio(q_eq_CD)))
fig.text(0.5, -0.02, txt, ha="center", fontsize=8.5, family="monospace")
fig.suptitle("Floor heat-flux distributions: what the model has now and what each option changes  (2026-09-30)", fontsize=12)
fig.tight_layout()
out = "/Users/tzetze20/amr_tools/alamo/Claude_markdowns/figures/2026-09-30_flux_options.png"
import os; os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=170, bbox_inches="tight")
print(out)
print("axis: now %.2f  cap50 %.2f  C %.2f  C+D %.2f  rz_budget %.2f" % (q_now[0], q_cap50[0], q_eq_C[0], q_eq_CD[0], q_rz_budget[0]))
print("P<40: now %.0f  cap50 %.0f  cap45 %.0f  half %.0f  rz_axis %.0f  CD %.0f W" % tuple(power_inside(q) for q in (q_now, q_cap50, q_cap45, q_half, q_rz_axis, q_eq_CD)))
