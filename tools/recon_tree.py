from common import *
import pickle
from collections import defaultdict
from scipy.spatial import cKDTree
rng=np.random.default_rng(7)
TH=pickle.load(open('recon_thorax.pkl','rb'))
OUT=defaultdict(list); LAY={}
def add(id,l,m): OUT[id].append(m); LAY[id]=l
def unit(v): v=np.asarray(v,float); return v/np.linalg.norm(v)
def colonize(lung, root, d0, N=900, di=0.15, dk=0.024, step=0.012, iters=260):
    lo,hi=lung.bounds; P=rng.uniform(lo,hi,(N*4,3)); P=P[lung.contains(P)][:N]
    nodes=[np.asarray(root,float)]; parent=[-1]
    # initial trunk segment
    for i in range(3): nodes.append(nodes[-1]+unit(d0)*step); parent.append(len(nodes)-2)
    alive=np.ones(len(P),bool)
    for it in range(iters):
        tr=cKDTree(np.array(nodes)); d,idx=tr.query(P[alive],distance_upper_bound=di)
        ap=P[alive]; grow=defaultdict(list)
        for j,(dd,ii) in enumerate(zip(d,idx)):
            if np.isfinite(dd): grow[ii].append(ap[j]-nodes[ii])
        if not grow: break
        for ii,vs in grow.items():
            dirv=unit(np.sum([unit(v) for v in vs],axis=0)+rng.normal(0,0.05,3))
            nodes.append(nodes[ii]+dirv*step); parent.append(ii)
        tr=cKDTree(np.array(nodes)); dd,_=tr.query(P); alive&=dd>dk
        if alive.sum()==0: break
    nodes=np.array(nodes); parent=np.array(parent)
    ch=defaultdict(list)
    for i,p in enumerate(parent):
        if p>=0: ch[p].append(i)
    r=np.zeros(len(nodes))
    for i in range(len(nodes)-1,-1,-1):
        r[i]=0.0012 if not ch[i] else (np.sum(r[ch[i]]**2.6))**(1/2.6)
    r=np.minimum(r,0.011)
    # chains
    chains=[]
    starts=[0]+[c for i in range(len(nodes)) if len(ch[i])>1 for c in ch[i]]
    for s in starts:
        seq=[parent[s]] if parent[s]>=0 else []
        seq.append(s); cur=s
        while len(ch[cur])==1: cur=ch[cur][0]; seq.append(cur)
        chains.append(seq)
    return nodes,r,chains
for sign,S in [(1,'_g'),(-1,'_d')]:
    lung=TH['poumon'+S][1]
    lung=trimesh.Trimesh(*fs.simplify(lung.vertices.astype(np.float32),lung.faces.astype(np.int32),target_reduction=0.85))
    hil=np.array([sign*0.10,0.41,0.68])
    nodes,r,chains=colonize(lung,hil,[sign*0.6,-0.2,-0.5])
    print(S,'nodes',len(nodes),'chains',len(chains))
    br=[];pa=[];pv=[]
    for c in chains:
        pts=nodes[c]
        if len(pts)<2: continue
        r0,r1=r[c[0]],r[c[-1]]
        br.append(tube(pts,r0*0.9,r1*0.9,seg=8 if r0>0.004 else 5,step=0.01))
        if r0>0.0022:
            t=unit(pts[-1]-pts[0]); o=unit(np.cross(t,[0,1,0])+1e-6)
            pa.append(tube(pts+o*(r0*0.9+0.002),r0*0.7,r1*0.7,seg=6,step=0.012))
            pv.append(tube(pts-o*(r0*0.9+0.002),r0*0.65,r1*0.65,seg=6,step=0.012))
    br.append(tube([[0,0.44,0.79],[sign*0.05,0.43,0.75],hil],0.012,0.011))
    pa.append(tube([[0.01,0.43,0.82],[sign*0.06,0.425,0.78],hil+np.array([0,0.012,0.01])],0.008,0.007))
    pv.append(tube([hil+np.array([0,-0.012,-0.01]),[sign*0.04,0.36,0.70],[0.0,0.34,0.72]],0.007,0.007))
    add('arbre_bronchique'+S,'o',trimesh.util.concatenate(br)); add('art_pulm'+S,'v',trimesh.util.concatenate(pa)); add('v_pulm'+S,'v',trimesh.util.concatenate(pv))
# ---------- coronaries on heart surface ----------
heart=TH['coeur'][1]
base=np.array([0.0,0.36,0.80]); apex=np.array([0.035,0.105,0.60]); ax=unit(apex-base); Lh=np.linalg.norm(apex-base)
u=unit(np.array([1,0,0])-ax*ax[0]); w=np.cross(ax,u)   # u ~ left, w ~ cranial-ish
if w[2]<0: w=-w
def surf(t,ang):
    c=base+ax*Lh*t; d=unit(np.cos(ang)*u+np.sin(ang)*w)
    loc,ri,_=heart.ray.intersects_location([c+d*0.3],[-d],multiple_hits=False)
    return (loc[0]+d*0.004) if len(loc) else c+d*0.1
cor=[]
# left coronary: paraconal interventricular (left-cranial face, base->apex)
cor.append(tube([surf(t,np.radians(55)) for t in np.linspace(0.18,0.97,12)],0.0035,0.0018))
# circumflex: coronary groove on left side going caudally/dorsally
cor.append(tube([surf(0.2,np.radians(a)) for a in np.linspace(60,-120,12)],0.0032,0.0018))
# subsinuosal interventricular (caudal face), from circumflex end
cor.append(tube([surf(t,np.radians(-125)) for t in np.linspace(0.2,0.85,9)],0.0026,0.0015))
# right coronary: groove on right side
cor.append(tube([surf(0.2,np.radians(a)) for a in np.linspace(80,250,12)],0.003,0.0018))
# diagonal branches
for t0,a0 in [(0.35,55),(0.55,55),(0.45,-60)]:
    cor.append(tube([surf(t0,np.radians(a0)),surf(t0+0.12,np.radians(a0-25)),surf(t0+0.22,np.radians(a0-45))],0.0018,0.0012))
add('coronaires','v',trimesh.util.concatenate(cor))
pickle.dump((dict(OUT),LAY),open('net_tree.pkl','wb'))
for k,v in OUT.items(): print(k,LAY[k],sum(len(m.faces) for m in v))
