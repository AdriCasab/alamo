import numpy as np
from scipy.ndimage import gaussian_filter
N=72
x=np.arange(N)+0.5-N/2
X,Y,Z=np.meshgrid(x,x,x,indexing='ij')
def estimate(S, A_true, label):
    Sp=np.pad(S,4,mode='edge'); core=(slice(4,-4),)*3
    nbs=[]
    for ax in range(3):
        for s in (1,-1):
            nbs.append(((ax,s), np.roll(Sp,-s,axis=ax)[core]==0))   # exposed face toward +s along ax
    nf=sum(m.astype(float) for _,m in nbs); surf=(S==1)&(nf>0)
    res={}
    # P2-8: central diff of indicator, w=1/|n|_1 on each exposed face
    g=[(np.roll(Sp,-1,axis=ax)-np.roll(Sp,1,axis=ax))[core] for ax in range(3)]
    gn=np.sqrt(sum(a*a for a in g)); l1=sum(np.abs(a) for a in g)
    w=np.where(gn>0,gn/np.where(l1>0,l1,1),1.0)
    res['P2-8 w']=np.sum((w*nf)[surf])
    for sig in (1.0,1.5):
        F=gaussian_filter(Sp,sig,mode='nearest',truncate=3.0/sig)
        G=[np.gradient(F,axis=ax)[core] for ax in range(3)]
        gn=np.sqrt(sum(a*a for a in G)); ok=gn>1e-12
        n=[np.where(ok,-a/np.where(ok,gn,1),0) for a in G]   # outward = toward void
        mu=np.full(S.shape,-2.0)
        for (ax,s),m in nbs:
            mu=np.where(m, np.maximum(mu, s*n[ax]), mu)   # n . e over EXPOSED faces
        bad=surf&(mu<=1e-6)
        a=np.where(mu>1e-6,1/np.maximum(mu,1e-6),nf)       # sheet/misaligned: count faces
        res[f'line s{sig}']=np.sum(a[surf])
        res[f'  misaligned s{sig}']=bad.sum()
    print(f'{label:34s} true {A_true:9.1f} | '+' | '.join(
        (f'{k} {100*(v/A_true-1):+6.2f}%' if 'mis' not in k else f'{k} {int(v)}') for k,v in res.items()))
R=20
r=np.sqrt(X**2+Y**2+Z**2)
estimate((r>R).astype(float), 4*np.pi*R**2, 'spherical cavity R20 (concave)')
estimate((r<R).astype(float), 4*np.pi*R**2, 'solid sphere R20 (convex)')
L=30; cube=(np.abs(X)<L/2)&(np.abs(Y)<L/2)&(np.abs(Z)<L/2)
estimate((~cube).astype(float), 6*L**2, 'cube cavity L30 (sharp concave)')
estimate(cube.astype(float), 6*L**2, 'solid cube L30 (sharp convex)')
H=30; Rc=15; cyl=(np.sqrt(X**2+Y**2)<Rc)&(np.abs(Z)<H/2)
estimate((~cyl).astype(float), 2*np.pi*Rc**2+2*np.pi*Rc*H, 'cylinder cavity R15 H30 (bore)')
