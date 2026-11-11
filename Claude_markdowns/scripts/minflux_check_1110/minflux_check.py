"""Paper check 2026-11-10: the span the buckling condition needs at each measured low-flux spall
   L_req = C_L(T_s) * c_f * k * (T_s - T0) / Q     (flat face, mu = 1; tangent depth k dT/Q)
   C_L from the model's own law (studies/p2_15_buckle/reduce.py -> DS.C_L), c_f = 0.4."""
import importlib.machinery, importlib.util, os, sys
ROOT = "/Users/tzetze20/amr_tools/alamo"
P15 = os.path.join(ROOT, "tests", "MMWSpalling", "studies", "p2_15_buckle")
ld = importlib.machinery.SourceFileLoader("red15", os.path.join(P15, "reduce.py"))
R = importlib.util.module_from_spec(importlib.util.spec_from_loader("red15", ld)); ld.exec_module(R)
CL, CF, T0 = R.DS.C_L, R.CF, 300.0
print("c_f", CF, " C_L at 700/730/758/800/820/850 K:", [round(CL(T), 2) for T in (700, 730, 758.15, 800, 820, 850)])
# Wilkinson & Tester 1993 Table 5 (laser; Q from the T_s - sqrt(t) slope), T_s in deg C
barre = [(0.31,481),(0.35,495),(0.38,520),(0.43,495),(0.47,457),(0.53,497),(0.54,497),(0.55,423),(0.58,423),
         (0.59,630),(1.47,604),(1.55,870),(1.78,690),(1.92,685),(1.95,808),(2.28,686),(2.39,930),(2.41,890)]
west  = [(0.23,457),(0.40,524),(0.45,457),(0.52,423),(0.56,423),(1.07,442),(1.34,645),
         (1.48,598),(1.54,646),(1.64,645),(1.77,710),(1.95,767),(2.14,483),(2.50,551)]
def L(Q, TC, k):
    Ts = TC + 273.15
    return 1e3 * CL(Ts) * CF * k * (Ts - T0) / (Q * 1e6)
for name, data in (("Barre", barre), ("Westerly", west)):
    print("\n%s (laser, Table 5): Q [MW/m2], T_s [C], L_req [mm] at k 2.14 / 1.6" % name)
    for Q, TC in data:
        if Q < 0.62: print("   %.2f  %4d   %5.1f / %5.1f" % (Q, TC, L(Q, TC, 2.14), L(Q, TC, 1.6)))
    for k in (2.14, 1.6):
        need = max(L(Q, TC, k) for Q, TC in data)
        for span in (12.0, 13.2, 15.4, 17.8):
            nb = sum(L(Q, TC, k) > span for Q, TC in data)
            print("   k %.2f  span %.1f mm: blocked %d of %d" % (k, span, nb, len(data)), end="")
            print("   (all pass needs %.1f mm)" % need if span == 12.0 else "")
# Hu granite 2, his own properties: nu 0.19, beta 8e-6 constant, k 2.5, T0 280.15, T_s 773.15
import math
eps = 8e-6 * (773.15 - 280.15); cl = math.pi / math.sqrt(6 * 1.19 * eps)
print("\nHu granite 2 (his beta, nu): C_L %.1f;  L_req at 0.074 / 0.13 / 0.17 MW/m2: %s mm" % (
    cl, " / ".join("%.0f" % (1e3 * cl * CF * 2.5 * 493 / (q * 1e6)) for q in (0.074, 0.13, 0.17))))
print("   with the model's C_L(773 K) %.1f instead: %s mm" % (CL(773.15),
    " / ".join("%.0f" % (1e3 * CL(773.15) * CF * 2.5 * 493 / (q * 1e6)) for q in (0.074, 0.13, 0.17))))
