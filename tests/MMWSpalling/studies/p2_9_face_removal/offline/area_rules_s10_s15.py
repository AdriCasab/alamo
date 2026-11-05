import numpy as np
from scipy.ndimage import gaussian_filter
def faces(S):
    Sp=np.pad(S,4,mode='edge'); core=(slice(4,-4),)*3
    nbs=[((ax,s), np.roll(Sp,-s,axis=ax)[core]==0) for ax in range(3) for s in (1,-1)]
    return Sp,core,nbs
def normals(Sp,core,kind):
    if kind=='cdiff':
        G=[(np.roll(Sp,-1,axis=ax)-np.roll(Sp,1,axis=ax))[core] for ax in range(3)]
    else:
        sig=float(kind[1:]); F=gaussian_filter(Sp,sig,mode='nearest',truncate=3.0/sig)
        G=[np.gradient(F,axis=ax)[core] for ax in range(3)]
    gn=np.sqrt(sum(a*a for a in G)); ok=gn>1e-12
    return [np.where(ok,-a/np.where(ok,gn,1),0) for a in G], ok
def est(S, sel=None):
    Sp,core,nbs=faces(S); out={}
    nf=sum(m.astype(float) for _,m in nbs); surf=(S==1)&(nf>0)
    if sel is not None: surf&=sel
    for kind in ('cdiff','g1.0','g1.5'):
        n,ok=normals(Sp,core,kind)
        l1=sum(np.abs(a) for a in n)
        A_w=0; A_dot=0; A_dots=0
        for (ax,s),m in nbs:
            mm=m&surf
            A_w+=np.sum(np.where(ok&(l1>0),1/np.where(l1>0,l1,1),1.0)[mm])      # P2-8 rule 1/|n|_1
            A_dot+=np.sum(np.where(ok,np.maximum(s*n[ax],0),1.0)[mm])           # projected, clamped
            A_dots+=np.sum(np.where(ok,s*n[ax],1.0)[mm])                          # projected, signed
        out[f'{kind}:1/|n|1']=A_w; out[f'{kind}:n.e+']=A_dot; out[f'{kind}:n.e']=A_dots
    return out
def show(label,S,A,sel=None):
    o=est(S,sel); print(f'{label:26s}'+' '.join(f'{k} {100*(v/A-1):+6.2f}%' for k,v in o.items()))
# planes: z-dominant tilts, interior window
N=48; NZ=64; m=8; rng=np.random.default_rng(3)
for th,az in [(10,0),(20,0),(35,0),(45,0),(35,45),(25,20)]:
    acc={}; 
    for _ in range(4):
        t,a=np.radians(th),np.radians(az); n=np.array([np.sin(t)*np.cos(a),np.sin(t)*np.sin(a),np.cos(t)])
        x=np.arange(N)+0.5; z=np.arange(NZ)+0.5; X,Y,Z=np.meshgrid(x,x,z,indexing='ij')
        S=(n[0]*(X-N/2)+n[1]*(Y-N/2)+n[2]*Z < NZ/2+rng.uniform()).astype(float)
        sel=np.zeros_like(S,bool); sel[m:N-m,m:N-m,:]=True; sel[:,:,NZ-4:]=False; sel[:,:,:4]=False
        Sp,core,nbs=faces(S); nf=sum(mm.astype(float) for _,mm in nbs); surf=(S==1)&(nf>0)&sel
        A=surf.sum()/n[2]
        for k,v in est(S,sel).items(): acc.setdefault(k,[]).append(v/A-1)
    print(f'plane {th:2d}/{az:2d}'+' '*14+' '.join(f'{k} {100*np.mean(v):+6.2f}%' for k,v in acc.items()))
N=72; x=np.arange(N)+0.5-N/2; X,Y,Z=np.meshgrid(x,x,x,indexing='ij'); r=np.sqrt(X**2+Y**2+Z**2)
show('sphere cavity R20',(r>20).astype(float),4*np.pi*400)
show('solid sphere R20',(r<20).astype(float),4*np.pi*400)
show('sphere cavity R10',(r>10).astype(float),4*np.pi*100)
cube=(np.abs(X)<15)&(np.abs(Y)<15)&(np.abs(Z)<15)
show('cube cavity L30',(~cube).astype(float),6*900)
cyl=(np.sqrt(X**2+Y**2)<15)&(np.abs(Z)<15)
show('cylinder cavity R15 H30',(~cyl).astype(float),2*np.pi*225+2*np.pi*15*30)
