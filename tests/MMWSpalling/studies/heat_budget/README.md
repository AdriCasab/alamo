# Heat budget: can a hot spalling surface + gas-side effects explain the uptake gap? (2026-10-26, user's request, OUTSIDE /plan)

Offline only: no model code, no runs, nothing tuned. Script `budget.py` -> `BUDGET.txt`.
Reuses P2-7's accepted energy bound (`p2_7_hoeser_heat/energy_bound.py`: E_terms, P_deliv) unchanged.
Literature: `validation/hoeser/literature_uptake_gap_notes.md`; Wilkinson & Tester 1993 Table 3 (°C, confirmed:
his energy balance reproduces his Q only in °C).

## Method
A. For each Höser load, the GAS-SIDE cut still needed if the surface the flame sees is at T_s:
   G_req(T_s) = delivered(T_s) / E_max(T_s)  (P2-7 bound; both shapes, both A_c). T_s ≤ melting
   (Höser ">900 °C" = 1,173 K; Wilkinson "~1,050 °C" = 1,323 K).
B. Where each case sits on Wilkinson's MEASURED drilling-speed -> steady-surface-T relation (Table 3, Barre and
   Westerly; IR T_s, V_dr from Rauenzahn's rate-vs-standoff curve). Below his slowest (2.2–2.5 m/h) the surface is
   at onset (~500 °C). Case speeds: Höser 0 MPa axis 1.2–2.0 m/h (no section exists: mean depth 10 mm x a 47 MPa-like
   axis/mean 1.6 = 16 mm in 30 s, low end = what-if s100 14 mm with onset); Höser 47 MPa axis 30 mm / 30 s = 3.6 m/h
   (Fig 5.10); Meier 1.28–1.5 m/h (Fig 8.8); Rauenzahn quarry 4.6–7.3 m/h (Table 2).
C. Rock-side share = G_req(onset) / G_req(T_s from B); the rest must be gas-side.
D. Gas-side available from the literature, as factors (low–high): cavity vortex 1.0–2.0 (Terekhov 2013; area offset),
   hot-wall gas properties 1.0–1.3, blowing 1.0–1.18, emissivity/catalysis 1.0–1.12, particles and roughness 1.0
   (they mostly ENHANCE). Product 1.0–3.4; "comfortable" ≤ 2.
E. Verdict per case: closes comfortably (G_req ≤ 2) / only at the literature's extreme (≤ 3.4) / does not close.

## Pre-registered predictions (written before running)
P1  FAST cases close: Höser 47 MPa (3.6 m/h, Barre-like T_s ~1,150–1,200 K) and Rauenzahn (4.6–7.3 m/h): G_req ≤ 2.
P2  SLOW cases do NOT close: Höser 0 MPa and Meier sit below Wilkinson's overheating speeds, so T_s ≈ onset and
    G_req stays 4.7–6.5 (Höser) / ~3–8 (Meier) > 3.4.
P3  Höser 0 MPa would close only if its flame-facing surface averaged ≥ ~1,150–1,250 K (at/above melting) with the
    gas side at 2–3.4 — a surface Wilkinson's speed relation says it should not have, and Höser reports no melt.
P4  Kant onset (before the first spall, no cavity): every candidate gives 1x — consistent, no discrimination.
⚠ Barre and Westerly differ strongly above 2.5 m/h; Grimsel is unmeasured, so B is bracketed by both.
⚠ Correction to my 10-26 claim "only rock-side fits the load trend": a spall temperature LOWERED by load also gives a
  smaller gap at 47 MPa (cheaper flakes + more delivered heat). The load trend does not discriminate by itself.

## Result (BUDGET.txt)
| # | measured | verdict |
|---|---|---|
| P1 | gas-side cut left: Höser 47 MPa Barre-like **1.01–1.14**, Westerly-like 2.11–2.39; Rauenzahn Barre 1.23–1.64, Westerly 1.73–2.55 | HOLDS for Barre-like rock; borderline for Westerly-like |
| P2 | Höser 0 MPa **4.99–6.40** vs gas side available ≤ 4.40; Meier 3.12–5.32 (its delivery is an upper bound at the jet's enthalpy cap) | HOLDS for Höser 0 MPa; Meier partly |
| P3 | Höser 0 MPa closes only with a flame-facing surface of **1,183–1,272 K** (gas side 2) or **894–985 K** (every gas effect at its maximum, 4.40) | HOLDS |
| P4 | trivially consistent | HOLDS |
⚠ Post-hoc: the "gas cools across the floor" factor (1.28, an energy balance on the wall jet; P2-7's delivery assumed T_f
everywhere over Ø40) was added AFTER the first run. It helps the gas side; no verdict changed.
⚠ Höser's thesis has IR only at ONSET (Ch. 4.2); the steady surface temperature of Grimsel during slow drilling is UNMEASURED.
**Conclusion:** energy conservation leaves two levers only — how hot the surface the flame sees is, and how well the gas
transfers heat. Fast drilling (Höser 47 MPa, Rauenzahn) closes with Wilkinson's measured hot surfaces plus a modest gas
side. Slow drilling (Höser 0 MPa, Meier) does not close with any measured surface temperature: Wilkinson's surfaces
below 2.2 m/h sat at onset (~500 °C), so stuck/overheated flakes have NO measured support in the regime of our main cases.
