# Tégument : enveloppe cutanée calculée autour de l'écorché et du squelette, et annexes (coussinets, griffes, vibrisses, pavillons, truffe)
from common import *
from common import parts
import pickle,json
from collections import defaultdict
from scipy import ndimage as ndi
OUT=defaultdict(list); LAY={}
def add(id,l,m): OUT[id].append(m); LAY[id]=l
def unit(v): v=np.asarray(v,float); return v/np.linalg.norm(v)
def ell(c,r,sub=2,R=None):
    m=trimesh.creation.icosphere(subdivisions=sub,radius=1.0); m.vertices*=np.asarray(r)
    if R is not None: m.vertices=m.vertices@R.T
    m.apply_translation(c); return m
lab2=json.load(open('data/bone_labels2.json')); G2=defaultdict(list)
for i,p in enumerate(parts): G2[lab2[str(i)]].append(p)
B2={k:trimesh.util.concatenate(v) for k,v in G2.items()}
# ---------- envelope ----------
ec=trimesh.load('sources/feline_ecorche_dec.glb').to_geometry()
P=[ec.sample(1500000)]+[p.sample(max(200,len(p.faces))) for p in parts]
# structures molles reconstruites : la peau doit les englober
EXC=('vibrisse','pavillon','truffe','coussinet','griffe','mamelle','peau','globe','cornee')
for f in ['net_branch.pkl','net_extra.pkl','net_more.pkl','net_trunk.pkl','net_tree.pkl','net_female.pkl']:
    d,_=pickle.load(open(f,'rb'))
    for k,ms in d.items():
        if any(e in k for e in EXC): continue
        for m in ms: P.append(m.sample(max(200,len(m.faces)//2)))
R=pickle.load(open('recon_thorax.pkl','rb'))
for k,v in R.items():
    m=v[1] if isinstance(v,tuple) else v
    P.append(m.sample(max(200,len(m.faces)//2)))
P=np.vstack(P)
h=0.008; lo=P.min(0)-0.06; hi=P.max(0)+0.06; shp=np.ceil((hi-lo)/h).astype(int)+1
G=np.zeros(shp,bool); ij=np.floor((P-lo)/h).astype(int); G[ij[:,0],ij[:,1],ij[:,2]]=True
print('grid',shp)
G=ndi.binary_closing(G,iterations=3)
G=ndi.binary_fill_holes(G)
for k in range(G.shape[2]): G[:,:,k]=ndi.binary_fill_holes(G[:,:,k])   # cavités ouvertes : pas de coque interne
for k in range(G.shape[0]): G[k]=ndi.binary_fill_holes(G[k])
G=ndi.binary_dilation(G,iterations=1)      # subcutis + skin thickness
# fentes palpébrales : la peau s'ouvre devant chaque globe oculaire
ii,jj,kk=np.indices(G.shape); XYZ=None
for s_ in (1,-1):
    c=np.array([s_*0.13,0.695,1.795]); r=0.056; gz=np.array([s_*0.5,0,1.0]); gz/=np.linalg.norm(gz)
    lo_i=np.floor((c-r-0.03-lo)/h).astype(int); hi_i=np.ceil((c+r+0.03-lo)/h).astype(int)
    sl=tuple(slice(a,b) for a,b in zip(lo_i,hi_i))
    P_=np.stack([ii[sl],jj[sl],kk[sl]],-1)*h+lo-c
    dist=np.linalg.norm(P_,axis=-1); front=(P_@gz)>0.25*r
    G[sl]&=~((dist<r+0.012)&front)
np.savez_compressed('skin_grid.npz',G=G,lo=lo,h=h)
skin=voxmesh(G,lo,h,smooth=15,sigma=1.4)
skin=max(skin.split(only_watertight=False),key=lambda c:len(c.faces))   # surface externe seule
print('skin faces',len(skin.faces)); add('peau','p',skin)
# ---------- pads, claws ----------
for s,S in [(1,'_g'),(-1,'_d')]:
    pads=[];claws=[]
    for ph,mc,src in [('m','mc','carpe_acc'),('p','mt',None)]:
        mids=[]
        for d in ['II','III','IV','V']:
            k2=f'phal2_{ph}{d}{S}'; k3=f'phal3_{ph}{d}{S}'; k1=f'phal1_{ph}{d}{S}'
            if k2 not in B2: continue
            p2=B2[k2].centroid; p3m=B2[k3]; p1=B2[k1].centroid
            low=min(np.asarray(B2[k2].vertices)[:,1].min(),np.asarray(p3m.vertices)[:,1].min())
            pads.append(ell([p2[0],low-0.006,p2[2]],[0.011,0.006,0.013]))
            mids.append([p1[0],np.asarray(B2[k1].vertices)[:,1].min()-0.008,p1[2]-0.01])
            v=np.asarray(p3m.vertices); c=v.mean(0); ax_=np.linalg.svd(v-c,full_matrices=False)[2][0]
            if ax_[2]<0: ax_=-ax_
            tip=c+ax_*0.03; base=c-ax_*0.005
            claws.append(tube([base,(base+tip)/2+np.array([0,0.004,0]),tip+np.array([0,-0.01,0])],0.0065,0.0008,seg=10))
        m=np.mean(mids,0); pads.append(ell(m+np.array([0,-0.002,-0.012]),[0.03,0.008,0.02]))  # metacarpal/metatarsal pad
        if src: c=B2[src+S].centroid; pads.append(ell(c+np.array([0,-0.03,-0.012]),[0.009,0.006,0.011]))  # carpal pad
    add('coussinets'+S,'p',trimesh.util.concatenate(pads)); add('griffes'+S,'p',trimesh.util.concatenate(claws))
# ---------- whiskers ----------
wk=[]
for s in (1,-1):
    for row,yy in enumerate([0.47,0.50,0.53,0.56]):
        for col in range(4):
            o=np.array([s*(0.075+0.008*col),yy,1.975-0.012*col])
            d=unit([s*(0.75+0.05*row),0.05-0.15*row,0.55]); L=0.30-0.02*row-0.015*col
            wk.append(tube([o,o+d*L*0.5+np.array([0,-0.01,0]),o+d*L+np.array([0,-0.035,0])],0.0013,0.0004,seg=5,step=0.02))
    for o,d,L in [((s*0.12,0.79,1.83),(s*0.5,0.7,0.5),0.14),((s*0.105,0.80,1.84),(s*0.35,0.8,0.5),0.12),((s*0.18,0.55,1.67),(s*0.9,0.1,-0.2),0.12)]:
        o=np.array(o); d=unit(d); wk.append(tube([o,o+d*L*0.5,o+d*L+np.array([0,-0.01,0])],0.0011,0.0004,seg=5,step=0.02))
add('vibrisses','p',trimesh.util.concatenate(wk))
# ---------- ear pinnae ----------
def pinna(s):
    base=np.array([s*0.15,0.80,1.60]); apex=np.array([s*0.21,1.00,1.59]); up=unit(apex-base)
    fwd=np.array([0,0,1.0]); lat=unit(np.cross(up,fwd)); fwd=unit(np.cross(lat,up))
    V=[];F=[]; nt,na=12,14
    for i in range(nt+1):
        t=i/nt; w=0.088*(1-t)**0.85+0.004; cup=0.03*(1-t)
        for j in range(na+1):
            a=-1+2*j/na; p=base+up*t*np.linalg.norm(apex-base)+lat*s*0*0 + fwd*(a*w) - lat*s*(cup*(1-a*a)) ; V.append(p)
    for i in range(nt):
        for j in range(na):
            a=i*(na+1)+j; b=a+1; c=a+na+1; d=c+1; F+=[[a,c,b],[b,c,d]]
    m=trimesh.Trimesh(np.array(V),np.array(F)); m2=m.copy(); m2.vertices=m2.vertices-m.vertex_normals*0.004; m2.invert()
    return trimesh.util.concatenate([m,m2])
add('pavillons','p',trimesh.util.concatenate([pinna(1),pinna(-1)]))
add('truffe','p',trimesh.util.concatenate([ell([0,0.50,1.962],[0.028,0.021,0.012]),ell([0.011,0.493,1.972],[0.006,0.005,0.004],sub=1),ell([-0.011,0.493,1.972],[0.006,0.005,0.004],sub=1)]))
pickle.dump((dict(OUT),LAY),open('net_skin.pkl','wb'))
for k,v in OUT.items(): print(k,sum(len(m.faces) for m in v))
