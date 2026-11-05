import numpy as np
from scipy.ndimage import gaussian_filter
N=88; x=np.arange(N)+0.5-N/2; X,Y,Z=np.meshgrid(x,x,x,indexing='ij'); r=np.sqrt(X**2+Y**2+Z**2)
def check(R, label, normal='g1.0', exact=False):
    S=(r>R).astype(float)            # rock outside a spherical cavity
    Sp=np.pad(S,4,mode='edge'); core=(slice(4,-4),)*3
    nbs=[((ax,s), np.roll(Sp,-s,axis=ax)[core]==0) for ax in range(3) for s in (1,-1)]
    nf=sum(m.astype(float) for _,m in nbs); surf=(S==1)&(nf>0)
    if exact:   # the true normal of the sphere (pointing into the cavity = -r_hat)
        n=[-X/r,-Y/r,-Z/r]
    else:
        sig=float(normal[1:]); F=gaussian_filter(Sp,sig,mode='nearest',truncate=2.0/sig)
        G=[np.gradient(F,axis=ax)[core] for ax in range(3)]; gn=np.sqrt(sum(a*a for a in G))
        n=[-a/np.where(gn>0,gn,1) for a in G]
    mu=np.full(S.shape,-2.0)
    for (ax,s),m in nbs: mu=np.where(m,np.maximum(mu,s*n[ax]),mu)
    # removed volume per uniform normal flake h: each surface cell recedes h/mu along its line (cell area dx^2)
    vol=np.where(surf,1/np.maximum(mu,1e-6),0)
    # bin by angle of the cell direction: max component of r_hat ~ how 'diagonal'
    rh=np.stack([np.abs(X),np.abs(Y),np.abs(Z)])/r
    dia=rh.max(axis=0)                      # 1 on axes, 0.577 at the body diagonal
    bins=[(0.95,1.01),(0.85,0.95),(0.75,0.85),(0.65,0.75),(0.57,0.65)]
    tot=vol.sum()/(4*np.pi*R*R)-1
    # true area per bin: fraction of sphere area with dia in bin (Monte Carlo)
    rng=np.random.default_rng(0); v=rng.normal(size=(400000,3)); v/=np.linalg.norm(v,axis=1)[:,None]
    dv=np.abs(v).max(axis=1)
    s=f'{label:22s} total {100*tot:+6.2f}% |'
    for lo,hi in bins:
        sel=surf&(dia>=lo)&(dia<hi)
        A_bin=4*np.pi*R*R*np.mean((dv>=lo)&(dv<hi))
        s+=f' max-comp {lo:.2f}-{hi:.2f}: {100*(vol[sel].sum()/A_bin-1):+6.1f}%'
    print(s)
for R in (12,24,36):
    check(R,f'R{R} exact normal',exact=True)
    check(R,f'R{R} smoothed g1.0')
