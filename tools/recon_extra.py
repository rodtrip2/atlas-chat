from common import *
import pickle,json
from collections import defaultdict
from scipy.spatial import cKDTree
OUT=defaultdict(list); LAY={}
def add(id,l,m): OUT[id].append(m); LAY[id]=l
def unit(v): v=np.asarray(v,float); return v/np.linalg.norm(v)
O=organs(); LM=json.load(open('data/landmarks.json'))
def ell(c,r,R=None,sub=3):
    m=trimesh.creation.icosphere(subdivisions=sub,radius=1.0); m.vertices*=np.asarray(r)
    if R is not None: m.vertices=m.vertices@R.T
    m.apply_translation(c); return m
# ---------- urinary tract ----------
bv=np.asarray(O['vessie'].vertices); zmin=bv[:,2].min()
neck=bv[bv[:,2]<zmin+0.03]; neck=neck[neck[:,1]<np.percentile(neck[:,1],40)].mean(0)
print('col vésical',np.round(neck,3))
k=O['reins'].copy(); k.merge_vertices(); kp=sorted(k.split(only_watertight=False),key=lambda p:-len(p.faces))[:2]
for kk in kp:
    c=kk.centroid; S='_g' if c[0]>0 else '_d'; sg=np.sign(c[0])
    hil=np.array([c[0]-sg*0.075,c[1]-0.005,c[2]])
    add('uretere'+S,'o',tube([hil,[sg*0.055,0.515,c[2]-0.12],[sg*0.05,0.51,-0.32],[sg*0.04,0.50,-0.52],[sg*0.022,0.47,neck[2]+0.035],[sg*0.012,neck[1]+0.03,neck[2]+0.01]],0.0045,0.004))
ure=[neck,[0,0.43,-0.74],[0,0.425,-0.88],[0,0.43,-1.05],[0,0.44,-1.25],[0,0.43,-1.34]]
add('uretre','o',tube(ure,0.0065,0.005))
# ---------- male genital ----------
add('prostate','o',ell([0,0.43,-0.80],[0.022,0.018,0.02]))
add('bulbo_uretrales','o',trimesh.util.concatenate([ell([s*0.016,0.445,-1.29],[0.01,0.009,0.011],sub=2) for s in (1,-1)]))
add('penis','o',trimesh.util.concatenate([tube([[0,0.43,-1.33],[0,0.405,-1.39],[0,0.38,-1.44]],0.011,0.009),ell([0,0.375,-1.452],[0.008,0.008,0.012],sub=2)]))
for s,S in [(1,'_g'),(-1,'_d')]:
    add('testicule'+S,'o',ell([s*0.026,0.465,-1.425],[0.022,0.026,0.03]))
    add('deferent'+S,'o',tube([[s*0.03,0.49,-1.405],[s*0.05,0.42,-1.30],[s*0.075,0.36,-1.05],[s*0.07,0.36,-0.86],[s*0.05,0.46,-0.72],[s*0.012,0.445,-0.78]],0.0025))
# ---------- rectum ----------
iv=np.asarray(O['intestins'].vertices); cend=iv[np.argmin(iv[:,2])]
add('rectum','o',tube([cend,[-0.015,0.57,-0.72],[0,0.595,-0.95],[0,0.60,-1.20],[0,0.60,-1.36]],0.024,0.02))
# ---------- larynx and hyoid apparatus ----------
add('larynx','o',trimesh.util.concatenate([tube([[0,0.515,1.545],[0,0.505,1.515],[0,0.497,1.49]],0.036,0.03,seg=16),
    tube([[0,0.53,1.555],[0,0.565,1.575]],0.014,0.008)]))   # epiglottis
hy=[]
for s in (1,-1):
    hy.append(tube([[s*0.075,0.585,1.585],[s*0.065,0.53,1.585],[s*0.05,0.47,1.58],[s*0.03,0.44,1.585]],0.0045,0.004))   # stylo/epi/ceratohyoid chain
    hy.append(tube([[s*0.012,0.44,1.585],[s*0.025,0.47,1.55],[s*0.03,0.495,1.52]],0.0035))                              # thyrohyoid
hy.append(tube([[-0.032,0.44,1.585],[0,0.437,1.588],[0.032,0.44,1.585]],0.0045))                                     # basihyoid
add('hyoide','b',trimesh.util.concatenate(hy))
# ---------- eyes: largest sphere inscribed in each orbit ----------
sk=np.asarray(B['crane'].vertices); tr=cKDTree(sk)
for s,S in [(1,'_g'),(-1,'_d')]:
    c0=np.array(LM['cr_orbite'+S]); best=(0,c0)
    for dx in np.linspace(-0.05,0.06,12):
        for dy in np.linspace(-0.04,0.04,9):
            for dz in np.linspace(-0.04,0.06,11):
                c=c0+np.array([s*dx,dy,dz]); d=tr.query(c)[0]
                # must stay inside the skull's lateral extent
                if abs(c[0])>0.15 or not (1.74<=c[2]<=1.84) or not (0.67<=c[1]<=0.78): continue
                if d>best[0]: best=(d,c)
    r=0.056; c=np.array([s*0.13,0.695,1.795])
    gaze=unit([s*0.5,0.0,1.0])
    print(S,'œil centre',np.round(c,3),'rayon',round(r,3))
    R=np.eye(3)
    add('globe'+S,'o',ell(c,[r,r,r],sub=3))
    add('cornee'+S,'o',ell(c+gaze*r*0.62,[r*0.62]*3,sub=3))
    # optic nerve to posterior pole
    add('n_optique'+S,'n',tube([[s*0.01,0.70,1.74],[s*0.04,0.71,1.76],c-gaze*r*0.95],0.0028))
pickle.dump((dict(OUT),LAY),open('net_extra.pkl','wb'))
for k_,v in OUT.items(): print(k_,LAY[k_],sum(len(m.faces) for m in v))
