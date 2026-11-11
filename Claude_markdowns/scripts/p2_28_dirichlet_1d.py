import numpy as np
from math import erfc, exp, sqrt, pi
T0 = 280.15
def Ts(t): return 638.22 + 3.74*t
lam = lambda T: 2000.0/((T-273.15)+563.0)
cp  = lambda T: 1370.0 - 178000.0/((T-273.15)+271.0)
rho = 2650.0
def run(dz, factor, tend_list, L=0.1, dt=None, props="grimsel", k_c=2.5, rc_c=2640*827.2):
    n = int(round(L/dz)); T = np.full(n, T0)   # index 0 = top cell
    if props == "grimsel":
        kf = lambda T: lam(T); rcf = lambda T: rho*cp(T)
    else:
        kf = lambda T: np.full_like(T, k_c); rcf = lambda T: np.full_like(T, rc_c)
    amax = max(kf(np.array([T0]))[0]/rcf(np.array([T0]))[0], 1e-12)
    if dt is None: dt = 0.2*dz*dz/((2+factor)*amax)
    out = {}; t = 0.0; targets = sorted(tend_list)
    for tt in targets:
        while t < tt - 1e-12:
            h = min(dt, tt - t)
            k = kf(T); rc = rcf(T)
            kface = 2*k[:-1]*k[1:]/(k[:-1]+k[1:])
            flux = kface*(T[1:]-T[:-1])/dz**2          # into cell i from i+1
            div = np.zeros(n); div[:-1] += flux; div[1:] -= flux
            div[0] += factor*k[0]*(Ts(t)-T[0])/dz**2   # the prescribed top
            T = T + h*div/rc; t += h
        out[tt] = T.copy()
    return out
# ---- (1) the reviewer's table: Grimsel, 3.125 mm, factor 1 vs 2, + a fine reference
times = [37.0, 50.0, 66.0, 108.0]
c1 = run(3.125e-3, 1.0, times); c2 = run(3.125e-3, 2.0, times)
ref = run(0.05e-3, 2.0, times)
nref = int(round(3.125e-3/0.05e-3))
print("Grimsel, Hu granite-2 ramp, top cell (3.125 mm) average temperature [K]")
print("  t[s]   factor1   factor2   fine-ref cell avg   T_B1")
for t in times:
    print("  %5.0f   %7.1f   %7.1f   %7.1f            %7.1f" % (t, c1[t][0], c2[t][0], ref[t][:nref].mean(), Ts(t)))
# e-fold depth of the rise (fine reference)
print("e-fold depth of the rise (T - T0) = (Ts - T0)/e, fine reference:")
z = (np.arange(len(ref[37.0]))+0.5)*0.05e-3
for t in times:
    r = (ref[t]-T0)/(Ts(t)-T0); i = np.argmax(r < 1/np.e)
    print("  t %5.0f s: %.2f mm" % (t, z[i]*1e3))
# ---- (2) exact solution, constant properties (Hu granite k 2.5, rho cp 2640*827.2), cell averages
k_c, rc_c = 2.5, 2640*827.2; a = k_c/rc_c; A, B = 638.22, 3.74
def exact(zz, t):
    eta = zz/(2*sqrt(a*t))
    return T0 + (A-T0)*erfc(eta) + B*t*((1+2*eta*eta)*erfc(eta) - 2/sqrt(pi)*eta*exp(-eta*eta))
def exact_cellavg(dz, t, m=400):
    zs = (np.arange(m)+0.5)*dz/m
    return np.mean([exact(zz, t) for zz in zs])
print("\nexact solution check, constant properties, top-cell average: error in the RISE [%]")
print("  dz[mm]  t[s]  factor1   factor2")
errs = {}
for dz in (3.125e-3, 1.5625e-3, 0.78125e-3):
    tl = [10.0, 30.0, 60.0, 120.0]
    r1 = run(dz, 1.0, tl, props="const"); r2 = run(dz, 2.0, tl, props="const")
    for t in tl:
        ex = exact_cellavg(dz, t)
        e1 = (r1[t][0]-ex)/(ex-T0)*100; e2 = (r2[t][0]-ex)/(ex-T0)*100
        errs[(dz, t)] = (e1, e2)
        print("  %6.3f %5.0f  %7.2f   %7.2f" % (dz*1e3, t, e1, e2))
print("\nconvergence ratio coarse/fine of |err| (3.125 -> 1.5625 -> 0.78125 mm):")
for t in (30.0, 60.0, 120.0):
    f1 = [abs(errs[(d, t)][0]) for d in (3.125e-3, 1.5625e-3, 0.78125e-3)]
    f2 = [abs(errs[(d, t)][1]) for d in (3.125e-3, 1.5625e-3, 0.78125e-3)]
    print("  t %5.0f: factor1 %.2f, %.2f   factor2 %.2f, %.2f" % (t, f1[0]/f1[1], f1[1]/f1[2], f2[0]/f2[1], f2[1]/f2[2]))
# ---- (3) 3-D explicit limit for the top cell with the corrected term
aT0 = lam(T0)/(rho*cp(T0))
for dz in (3.125e-3,):
    print("\nalpha(280.15 K) = %.4e m2/s; 3-D top-cell limit dz^2/(7 alpha) = %.3f s (dz^2/(6 alpha) = %.3f s); Hu deck dt 0.5 s" % (aT0, dz*dz/(7*aT0), dz*dz/(6*aT0)))
