"""P2-9 rev 2 a4 predictions: per tilt, the per-face area error (P2-8 rule w = 1/|n|_1 on
every exposed face) and the line-rule removal excess (sum of 1/mu, mu = n.e* over the
exposed face most aligned with n), both with the normal from the solid indicator smoothed
by a Gaussian sigma = 1 cell truncated at +-2 cells.
PREDICTED STEADY NORMAL RECESSION-RATE ERROR = THE AREA ERROR. In steady state the heat in
equals the flake enthalpy out, and a line's flake carries enthalpy in proportion to the
volume it removes, so the volume removed per unit time is heat / (rho c dT_flake) whatever
the per-firing volume; the plane's mean normal speed is that volume / the TRUE area. The
removal excess (sum of 1/mu over the plane's area) is the volume of ONE simultaneous firing
of every surface cell -- it changes how often cells fire, not how fast the plane recedes,
and is reported only as that.
60 deg in x-z is the 30 deg case by the grid's cubic symmetry. Binary staircase, interior
window, 4 random plane offsets each."""
import numpy as np
from scipy.ndimage import gaussian_filter
rng = np.random.default_rng(7)
def one(theta, az, N=40, NZ=120, m=8, sig=1.0, trunc=2.0):
    t, a = np.radians(theta), np.radians(az)
    n = np.array([np.sin(t)*np.cos(a), np.sin(t)*np.sin(a), np.cos(t)])
    x = np.arange(N)+0.5; z = np.arange(NZ)+0.5
    X, Y, Z = np.meshgrid(x, x, z, indexing='ij')
    S = (n[0]*(X-N/2)+n[1]*(Y-N/2)+n[2]*Z < NZ/2+rng.uniform()).astype(float)
    Sp = np.pad(S, 4, mode='edge'); Sp[:, :, -4:] = 0; core = (slice(4, -4),)*3
    nbs = [((ax, s), np.roll(Sp, -s, axis=ax)[core] == 0) for ax in range(3) for s in (1, -1)]
    nf = sum(mm.astype(float) for _, mm in nbs)
    sel = np.zeros_like(S, bool); sel[m:N-m, m:N-m, 6:NZ-6] = True
    surf = (S == 1) & (nf > 0) & sel
    A = surf.sum()/n[2]                       # one surface cell per z-line on a plane
    F = gaussian_filter(Sp, sig, mode='nearest', truncate=trunc/sig)
    G = [np.gradient(F, axis=ax)[core] for ax in range(3)]
    gn = np.sqrt(sum(g*g for g in G)); ok = gn > 1e-12
    nn = [np.where(ok, -g/np.where(ok, gn, 1), 0) for g in G]
    l1 = sum(np.abs(c) for c in nn)
    area = np.sum((np.where(ok & (l1 > 0), nf/np.where(l1 > 0, l1, 1), nf))[surf])
    mu = np.full(S.shape, -2.0)
    for (ax, s), mm in nbs: mu = np.where(mm, np.maximum(mu, s*nn[ax]), mu)
    line = np.sum((1/np.maximum(mu, 1e-6))[surf])
    return area/A-1, line/A-1
print(__doc__)
print(f"{'tilt/az':>10s} {'area err = pred. rate err':>26s} {'one-firing volume excess':>26s}")
for th, az, lab in [(20,0,'20/0'),(30,0,'30/0 = 60'),(35,0,'35/0'),(45,0,'45/0'),
                    (35,45,'35/45 (TR)'),(np.degrees(np.arccos(1/np.sqrt(3))),45,'(1,1,1) TD')]:
    r = np.array([one(th, az) for _ in range(4)]).mean(axis=0)
    print(f"{lab:>10s} {100*r[0]:+25.2f}% {100*r[1]:+25.2f}%")
