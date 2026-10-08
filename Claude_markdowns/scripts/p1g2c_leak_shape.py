import os, sys, csv, math
import numpy as np
import yt
yt.set_log_level(50)
R="/Users/tzetze20/amr_tools/alamo/tests/MMWSpalling/studies"
LEGS={"CTRL2":R+"/p1g3_thickness/MEIER_CTRL2","CTRL4":R+"/p1g2b_subcell/CTRL4","FIXED":R+"/p1g2c_shiftfix/CTRL4_FIXED"}
DX=0.002; KTOP0=299
BANDS=((10,50),(50,90),(110,150),(150,200),(200,300),(300,400))
NX=60; NZ=300
ii,jj=np.meshgrid(np.arange(NX),np.arange(NX),indexing="ij")
RC=DX*np.hypot(ii+0.5,jj+0.5)*1e3          # mm
TH=np.degrees(np.arctan2(jj+0.5,ii+0.5))   # deg, 0 = +x axis
AZB=np.minimum((TH//15).astype(int),5)      # 6 bins
depth=(KTOP0-np.arange(NZ))*DX*1e3          # mm per k

def load(leg,step):
    ds=yt.load(f"{LEGS[leg]}/FH/{step:05d}cell")
    cg=ds.covering_grid(0,left_edge=ds.domain_left_edge,dims=ds.domain_dimensions)
    v=np.array(cg["boxlib","removed"])>0.5
    return v  # [i,j,k]

def floor_k(v):
    # deepest void per column; floor front = min over all columns of lowest void k
    anyv=v.any(axis=(0,1))
    ks=np.where(anyv)[0]
    return ks.min() if len(ks) else NZ

def layer_stats(v):
    # per k: N void, r_max, per-azimuth r_max
    N=v.sum(axis=(0,1))
    rmax=np.zeros(NZ); raz=np.zeros((NZ,6))
    for k in range(NZ):
        m=v[:,:,k]
        if not m.any(): continue
        rmax[k]=RC[m].max()
        for b in range(6):
            mb=m&(AZB==b)
            if mb.any(): raz[k,b]=RC[mb].max()
    return N,rmax,raz

def band_rows(N,rmax,raz,kfront):
    out=[]
    for a,c in BANDS:
        sel=(depth>=a)&(depth<c)&(depth<= (KTOP0-kfront)*DX*1e3-20)
        if not sel.any(): out.append(None); continue
        Deq=np.sqrt(16*N[sel]*DX*DX/math.pi).mean()*1e3
        Dmax=2*rmax[sel].max()
        Dmaxmean=2*rmax[sel].mean()
        Daz=2*raz[sel].mean()          # mean over layers and azimuth bins
        out.append((Deq,Dmax,Dmaxmean,Daz))
    return out

V={}; S={}
for leg in LEGS:
    V[leg]=load(leg,31250); S[leg]=layer_stats(V[leg])
kf={leg:floor_k(V[leg]) for leg in LEGS}
print("floor front k / depth mm:",{l:(int(kf[l]),round((KTOP0-kf[l])*DX*1e3,1)) for l in LEGS})
print("total void cells:",{l:int(V[l].sum()) for l in LEGS})

print("\n=== A. FINAL BORE PER BAND, four definitions (mm, full diameter) ===")
print("Dmax = 2x outermost void in band (the score's statistic); Dmm = 2x mean-over-layers of outermost; Daz = 2x mean over 6 azimuth bins & layers; Deq = area-equivalent")
hdr="%-9s"%"band"+"".join("%26s"%l for l in LEGS)
print(hdr); print("%-9s"%""+"".join("%26s"%"Deq  Dmax  Dmm  Daz" for l in LEGS))
BR={l:band_rows(*S[l],min(kf.values())) for l in LEGS}
for b,(a,c) in enumerate(BANDS):
    line="%3d-%-4d "%(a,c)
    for l in LEGS:
        r=BR[l][b]
        line+= "%26s"%("--" if r is None else "%5.1f %5.1f %5.1f %5.1f"%r)
    print(line)
print("\n   deltas (mm)        Deq    Dmax    Dmm    Daz")
for b,(a,c) in enumerate(BANDS):
    for (p,q) in (("CTRL2","CTRL4"),("CTRL4","FIXED"),("CTRL2","FIXED")):
        rp,rq=BR[p][b],BR[q][b]
        if rp is None or rq is None: continue
        d=[y-x for x,y in zip(rp,rq)]
        print("  %3d-%-4d %-6s->%-6s %6.2f %6.2f %6.2f %6.2f"%(a,c,p,q,*d))
    print()

print("=== B. OUTERMOST RADIUS BY AZIMUTH (mm), bands pooled 50-300 mm, mean over layers ===")
sel=(depth>=50)&(depth<300)
print("%-12s"%"azimuth"+"".join("%9s"%l for l in LEGS)+"   C4-C2  FX-C4  FX-C2")
for b in range(6):
    vals=[S[l][2][sel,b].mean() for l in LEGS]
    print("%3d-%3d deg  "%(15*b,15*b+15)+"".join("%9.2f"%x for x in vals)+"  %6.2f %6.2f %6.2f"%(vals[1]-vals[0],vals[2]-vals[1],vals[2]-vals[0]))
print("\n   per band, FX-C4 by azimuth bin:")
print("%-9s"%"band"+"".join("%9s"%("%d-%d"%(15*b,15*b+15)) for b in range(6)))
for a,c in BANDS[:5]:
    sel=(depth>=a)&(depth<c)
    print("%3d-%-4d "%(a,c)+"".join("%9.2f"%(S["FIXED"][2][sel,b].mean()-S["CTRL4"][2][sel,b].mean()) for b in range(6)))
print("   per band, C4-C2 by azimuth bin:")
for a,c in BANDS[:5]:
    sel=(depth>=a)&(depth<c)
    print("%3d-%-4d "%(a,c)+"".join("%9.2f"%(S["CTRL4"][2][sel,b].mean()-S["CTRL2"][2][sel,b].mean()) for b in range(6)))

print("\n=== C. CELL-BY-CELL DIFFERENCE CENSUS (void in one leg, solid in the other), depth 10-300 mm ===")
dsel=(depth>=10)&(depth<300)
for (p,q) in (("CTRL2","CTRL4"),("CTRL4","FIXED")):
    only_p=V[p]&~V[q]; only_q=V[q]&~V[p]
    print(f"  {p} vs {q}: void only in {p}: {int(only_p[:,:,dsel].sum())}   void only in {q}: {int(only_q[:,:,dsel].sum())}")
    print("    by azimuth bin  "+"  ".join("%d-%d"%(15*b,15*b+15) for b in range(6)))
    print("    only in %-6s "%p+"  ".join("%5d"%int((only_p[:,:,dsel]&(AZB==b)[:,:,None]).sum()) for b in range(6)))
    print("    only in %-6s "%q+"  ".join("%5d"%int((only_q[:,:,dsel]&(AZB==b)[:,:,None]).sum()) for b in range(6)))
    print("    by band        "+"  ".join("%d-%d"%bb for bb in BANDS[:5]))
    print("    only in %-6s "%p+"  ".join("%5d"%int(only_p[:,:,(depth>=a)&(depth<c)].sum()) for a,c in BANDS[:5]))
    print("    only in %-6s "%q+"  ".join("%5d"%int(only_q[:,:,(depth>=a)&(depth<c)].sum()) for a,c in BANDS[:5]))

print("\n=== D. x/y MIRROR SYMMETRY of the void field (voids with theta<45 vs >45, depth 10-300 mm; floor cells r<25 mm excluded) ===")
for l in LEGS:
    m=V[l][:,:,dsel]&(RC>25)[:,:,None]
    lo=int((m&(TH<45)[:,:,None]).sum()); hi=int((m&(TH>45)[:,:,None]).sum())
    asym=np.abs(V[l].astype(int)-np.transpose(V[l],(1,0,2)).astype(int))[:,:,dsel].sum()//2
    print("  %-6s theta<45: %6d  theta>45: %6d  ratio %.4f   cells not mirrored across the diagonal: %d"%(l,lo,hi,hi/max(lo,1),asym))

print("\n=== E. TIME EVOLUTION: Deq (mm) per band, every 50 s, bands >=20 mm above the shallowest front at that time ===")
steps=list(range(0,31251,3125))
for a,c in BANDS[:5]:
    print(f"  band {a}-{c} mm")
    print("    t[s]  "+"".join("%8s"%l for l in LEGS)+"   C4-C2   FX-C4")
    for st in steps:
        row=[]
        for l in LEGS:
            v=load(l,st); N=v.sum(axis=(0,1)); kfr=floor_k(v)
            sel=(depth>=a)&(depth<c)&(depth<=(KTOP0-kfr)*DX*1e3-20)
            row.append(np.sqrt(16*N[sel]*DX*DX/math.pi).mean()*1e3 if sel.any() else float("nan"))
        print("    %4d  "%(st*0.016)+"".join("%8.2f"%x for x in row)+"  %6.2f  %6.2f"%(row[1]-row[0],row[2]-row[1]))

print("\n=== F. EVENT LOG: lateral voids by azimuth bin (count) and median time-to-void after the floor passed, depth 50-300 mm ===")
def ev(leg):
    rows=[]; floor_t={}
    with open(LEGS[leg]+"/FH_removal_events.csv") as fh:
        for r in csv.DictReader(fh):
            reg=int(float(r["regime"])); t=float(r["time"]); k=int(float(r["k_top"]))
            i=int(r["col_i"]); j=int(r["col_j"])
            if reg==1:
                d=(KTOP0-k)*DX*1e3
                # first time ANY column reached depth d (floor pass)
                kk=k
                if kk not in floor_t: floor_t[kk]=t
            elif reg==6 and float(r["n_voided"] or 0)>=1:
                rows.append((t,i,j,k))
    out=[[] for _ in range(6)]
    for t,i,j,k in rows:
        d=(KTOP0-k)*DX*1e3
        if not (50<=d<300): continue
        th=math.degrees(math.atan2(j+0.5,i+0.5)); b=min(int(th//15),5)
        tf=floor_t.get(k,None)
        out[b].append((t,None if tf is None else t-tf))
    return out
EV={l:ev(l) for l in LEGS}
print("%-12s"%"azimuth"+"".join("%16s"%l for l in LEGS))
print("%-12s"%""+"".join("%16s"%"n   lag_p50[s]" for l in LEGS))
for b in range(6):
    line="%3d-%3d deg  "%(15*b,15*b+15)
    for l in LEGS:
        xs=EV[l][b]; lags=[x[1] for x in xs if x[1] is not None]
        line+="%7d %8.0f"%(len(xs),np.median(lags) if lags else float("nan"))
    print(line)
