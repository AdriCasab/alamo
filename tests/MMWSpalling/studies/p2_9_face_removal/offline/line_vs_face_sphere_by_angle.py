import numpy as np
from scipy.ndimage import gaussian_filter
N=88; x=np.arange(N)+0.5-N/2; X,Y,Z=np.meshgrid(x,x,x,indexing='ij'); r=np.sqrt(X**2+Y**2+Z**2)
rng=np.random.default_rng(0); v=rng.normal(size=(600000,3)); v/=np.linalg.norm(v,axis=1)[:,None]; dv=np.abs(v).max(axis=1)
bins=[(0.95,1.01),(0.85,0.95),(0.75,0.85),(0.65,0.75),(0.57,0.65)]
def run(R, rule, sig=1.0, cavity=True):
    S=(r>R).astype(float) if cavity else (r<R).astype(float)
    Sp=np.pad(S,4,mode='edge'); core=(slice(4,-4),)*3
    nbs=[((ax,s), np.roll(Sp,-s,axis=ax)[core]==0) for ax in range(3) for s in (1,-1)]
    nf=sum(m.astype(float) for _,m in nbs); surf=(S==1)&(nf>0)
    F=gaussian_filter(Sp,sig,mode='nearest',truncate=2.0/sig)
    G=[np.gradient(F,axis=ax)[core] for ax in range(3)]; gn=np.sqrt(sum(a*a for a in G))
    n=[-a/np.where(gn>0,gn,1) for a in G]
    if rule=='line':
        mu=np.full(S.shape,-2.0)
        for (ax,s),m in nbs: mu=np.where(m,np.maximum(mu,s*n[ax]),mu)
        a=np.where(surf,1/np.maximum(mu,1e-6),0)
    else:
        l1=sum(np.abs(c) for c in n); a=np.where(surf,nf/np.where(l1>0,l1,1),0)
    dia=(np.stack([np.abs(X),np.abs(Y),np.abs(Z)])/r).max(axis=0)
    A=4*np.pi*R*R; s=f'{"cavity" if cavity else "solid "} R{R:2d} {rule:5s} total {100*(a.sum()/A-1):+5.1f}% |'
    for lo,hi in bins:
        sel=surf&(dia>=lo)&(dia<hi); Ab=A*np.mean((dv>=lo)&(dv<hi))
        s+=f' {lo:.2f}-{hi:.2f} {100*(a[sel].sum()/Ab-1):+6.1f}%'
    print(s)
for R in (12,24,36):
    for rule in ('line','face'): run(R,rule)
for rule in ('line','face'): run(24,rule,cavity=False)
