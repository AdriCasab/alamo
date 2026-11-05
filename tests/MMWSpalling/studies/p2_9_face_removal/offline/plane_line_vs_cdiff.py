import numpy as np
from scipy.ndimage import gaussian_filter
rng=np.random.default_rng(2)
N=48; NZ=64; m=8
def run(theta,az,ntr=5):
    t=np.radians(theta); a=np.radians(az)
    n=np.array([np.sin(t)*np.cos(a), np.sin(t)*np.sin(a), np.cos(t)])
    out={}
    for _ in range(ntr):
        c=0.5*NZ+rng.uniform(0,1)
        x=np.arange(N)+0.5; z=np.arange(NZ)+0.5
        X,Y,Z=np.meshgrid(x,x,z,indexing='ij')
        S=(n[0]*(X-N/2)+n[1]*(Y-N/2)+n[2]*Z < c).astype(float)
        Sp=np.pad(S,4,mode='edge'); Sp[:,:,-4:]=0
        core=(slice(4,-4),)*3
        nf=np.zeros_like(S)
        for ax in range(3):
            for s in (1,-1): nf+=(np.roll(Sp,-s,axis=ax)[core]==0)
        surf=(S==1)&(nf>0); inner=np.zeros_like(S,bool); inner[m:N-m,m:N-m,:]=True
        si=surf&inner; A_tr=si.sum()/n[2]
        def est(G,name):
            gn=np.sqrt(sum(g*g for g in G)); ok=gn>0
            l1=sum(np.abs(g) for g in G)
            w=np.where(ok,gn/np.where(ok,l1,1),1.0)
            mu=np.where(ok,np.maximum.reduce([np.abs(g) for g in G])/np.where(ok,gn,1),1.0)
            out.setdefault('W_'+name,[]).append(np.sum((w*nf)[si])/A_tr-1)        # P2-8 per-face w=1/|n|1
            out.setdefault('L_'+name,[]).append(np.sum((1/mu)[si])/A_tr-1)        # line rule dx^2/mu
            out.setdefault('rmsW_'+name,[]).append(np.sqrt(np.mean(((w*nf)[si]*n[2]-1)**2)))
        est([(np.roll(Sp,-1,axis=ax)-np.roll(Sp,1,axis=ax))[core] for ax in range(3)],'cdiff')
        for s in (1.0,1.5,2.0):
            F=gaussian_filter(Sp,s,mode='nearest',truncate=3.0/s)   # stencil +-3 cells
            est([np.gradient(F,axis=ax)[core] for ax in range(3)],f'g{s}')
    return {k:100*np.mean(v) for k,v in out.items()}
cases=[(5,0),(10,0),(20,0),(30,0),(40,0),(45,0),(15,30),(25,20),(35,45),(50,45)]
first=True
for th,az in cases:
    r=run(th,az)
    if first: print('tilt/az   '+' '.join(f'{k:>11s}' for k in r)); first=False
    print(f'{th:4.0f}/{az:3.0f}  '+' '.join(f'{v:+10.2f}%' for v in r.values()))
