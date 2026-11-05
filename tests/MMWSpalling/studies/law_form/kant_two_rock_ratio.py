#!/usr/bin/env python3
"""Two-rock threshold test on Kant & Rudolf von Rohr 2016 (IJRMMS 84), zero compute.

Question: does a thickness-type minimum-flux condition with ONE rock-independent length
(buckling span L_max, or equally a fixed reference depth d) reproduce the DIFFERENT
spall/no-spall standoffs of Gotthard (130 | 140 mm) and Bethel (150 | 160 mm) from the two
rocks' measured k and onset temperature alone?

Sources (all on disk):
  franck_material/MinimalBoundaryConditions.pdf  Table 1 (spall/no spall), Table 2 (k),
  Table 4 (h_fl, theta_f at 50/80/110/130/150 mm), Fig. 11 (theta_SP by SOD).
  Notes: tests/MMWSpalling/validation/hoeser/kant_2016_minimal_bc_notes.md
Assumptions, stated: h_fl and theta_f at 140 and 160 mm are linear in SOD between/after the
130 and 150 mm rows (Table 4 has no 140/160 row; Fig. 13 puts the no-spall cases at
h 600-700 W/m2K, theta_f 1100-1250 C, consistent with this). Onset at threshold: Gotthard
570 C at 130 mm, Bethel 540 C at 150 mm (Fig. 11 trend endpoints). T0 = 20 C.
"""
k   = {"Gotthard": 3.27, "Bethel": 2.62}            # W/m/K, Table 2
th  = {"Gotthard": 570.0, "Bethel": 540.0}          # C, onset at the last spalling SOD (Fig. 11)
T0  = 20.0
dT  = {r: th[r]-T0 for r in k}
# Table 4 rows: SOD -> (h_fl W/m2K, theta_f C)
tab = {130:(747,1250), 150:(633,1100)}
def interp(sod):
    (h1,f1),(h2,f2)=tab[130],tab[150]
    w=(sod-130)/20.0
    return h1+(h2-h1)*w, f1+(f2-f1)*w
rows={130:tab[130],140:interp(140),150:tab[150],160:interp(160)}
def q(sod,rock):            # flux into the rock surface at its onset temperature
    h,f=rows[sod]; return h*(f-th[rock])          # W/m2
brk={"Gotthard":(130,140),"Bethel":(150,160)}      # (last spall, first no-spall)
out=[]
P=lambda s: out.append(s)
P("Kant 2016 two-rock threshold test  (q = h_fl * (theta_f - theta_SP), W/m2)")
P("")
P("  SOD   h_fl   theta_f   q_Gotthard  q_Bethel   note")
for s in (130,140,150,160):
    h,f=rows[s]
    P(f"  {s:3d}   {h:5.0f}   {f:6.0f}    {q(s,'Gotthard')/1e3:7.0f}     {q(s,'Bethel')/1e3:6.0f}    "+("Table 4" if s in tab else "interpolated"))
P("")
thr={}
for r,(a,b) in brk.items():
    hi,lo=q(a,r),q(b,r)
    thr[r]=(lo,hi)
    P(f"  {r:9s} threshold flux bracket: {lo/1e3:5.0f} - {hi/1e3:5.0f} kW/m2   (spalls at {a} mm, not at {b} mm)")
P("")
# measured ratio bracket
rb=(thr['Bethel'][0]/thr['Gotthard'][1], thr['Bethel'][1]/thr['Gotthard'][0])
P(f"  MEASURED ratio q_min(Bethel)/q_min(Gotthard): {rb[0]:.2f} - {rb[1]:.2f}  (midpoint {0.5*(rb[0]+rb[1]):.2f})")
pred_kdT=(k['Bethel']*dT['Bethel'])/(k['Gotthard']*dT['Gotthard'])
pred_k=k['Bethel']/k['Gotthard']
P(f"  PREDICTED, one rock-independent length, q_min ∝ k*dT_s:  {pred_kdT:.2f}   (k alone: {pred_k:.2f})")
P(f"     -> inside the measured bracket: {rb[0] <= pred_kdT <= rb[1]}")
P("")
P("  Heated-layer depth at threshold, delta = k*dT_s/q  (the length both forms hold fixed):")
for r in k:
    lo,hi=thr[r]
    d_lo,d_hi=k[r]*dT[r]/hi*1e3, k[r]*dT[r]/lo*1e3
    P(f"    {r:9s} delta = {d_lo:4.1f} - {d_hi:4.1f} mm   -> buckling span C_L*c_f*delta (C_L 13, c_f 0.4) = {13*0.4*d_lo/10:4.1f} - {13*0.4*d_hi/10:4.1f} cm")
g=(k['Gotthard']*dT['Gotthard']/thr['Gotthard'][1]*1e3, k['Gotthard']*dT['Gotthard']/thr['Gotthard'][0]*1e3)
b=(k['Bethel']*dT['Bethel']/thr['Bethel'][1]*1e3, k['Bethel']*dT['Bethel']/thr['Bethel'][0]*1e3)
common=(max(g[0],b[0]),min(g[1],b[1]))
P(f"    common delta satisfying BOTH rocks: {common[0]:.1f} - {common[1]:.1f} mm  (non-empty: {common[0]<common[1]})")
P("")
P("  Alternative scalings, same bracket:")
P(f"    q_min ∝ UCS (Kant's own reading, 'Gotthard needs more'):  {152/214:.2f}")
P(f"    q_min ∝ k only (equal onset):                               {pred_k:.2f}")
P(f"    q_min ∝ k*dT_s (buckling OR fixed d, same length):           {pred_kdT:.2f}")
P("")
P("  VERDICT: the two-rock ratio is CONSISTENT with one rock-independent length (0.76 inside 0.55-0.85)")
P("  but NOT DISCRIMINATING: the fixed-d form predicts the identical ratio (both scale as k*dT_s), and a")
P("  UCS scaling lands inside the bracket too. Thresholds alone cannot separate the forms; the onset-vs-flux")
P("  SHAPE above threshold (flat vs climbing) remains the discriminator, and it already favours the step.")
txt="\n".join(out); print(txt)
open(__file__.replace("kant_two_rock_ratio.py","RATIO.txt"),"w").write(txt+"\n")
