# Appareil génital femelle et mamelles (reconstruction schématique)
from common import *
import pickle,json
from collections import defaultdict
OUT=defaultdict(list); LAY={}
def add(id,l,m): OUT[id].append(m); LAY[id]=l
def ell(c,r,sub=2):
    m=trimesh.creation.icosphere(subdivisions=sub,radius=1.0); m.vertices*=np.asarray(r); m.apply_translation(c); return m
O=organs()
bv=np.asarray(O['vessie'].vertices); zmin=bv[:,2].min(); neck=bv[bv[:,2]<zmin+0.03]; neck=neck[neck[:,1]<np.percentile(neck[:,1],40)].mean(0)
for s,S in [(1,'_g'),(-1,'_d')]:
    ov=np.array([s*0.10,0.47,-0.19 if s>0 else -0.14])
    add('ovaire'+S,'o',ell(ov,[0.012,0.008,0.016]))
    add('trompe'+S,'o',tube([ov+np.array([s*0.012,0.008,0.01]),ov+np.array([s*0.018,-0.004,-0.004]),ov+np.array([s*0.008,-0.012,-0.018]),ov+np.array([-s*0.004,-0.01,-0.026])],0.0022))
    add('corne_uterine'+S,'o',tube([ov+np.array([-s*0.004,-0.012,-0.028]),[s*0.08,0.40,-0.32],[s*0.055,0.415,-0.48],[s*0.022,0.45,-0.63]],0.0075,0.0085))
add('uterus','o',tube([[0,0.455,-0.62],[0,0.465,-0.68],[0,0.47,-0.72]],0.012,0.011))
add('col_uterin','o',tube([[0,0.47,-0.72],[0,0.47,-0.75]],0.0095))
add('vagin','o',tube([[0,0.47,-0.75],[0,0.472,-0.92],[0,0.468,-1.10],[0,0.455,-1.25],[0,0.44,-1.33]],0.012,0.011))
add('vulve','o',trimesh.util.concatenate([ell([0.006,0.415,-1.40],[0.006,0.016,0.008]),ell([-0.006,0.415,-1.40],[0.006,0.016,0.008])]))
add('uretre_f','o',tube([neck,[0,0.43,-0.74],[0,0.428,-0.90],[0,0.435,-1.05],[0,0.445,-1.15]],0.0065,0.0055))
# mammary glands: 4 pairs on the ventral body wall (surface of the écorché)
EV=np.asarray(trimesh.load('sources/feline_ecorche_dec.glb').to_geometry().vertices)
for s,S in [(1,'_g'),(-1,'_d')]:
    ms=[]
    for z in [0.72,0.45,0.05,-0.40]:
        sel=EV[(np.abs(EV[:,2]-z)<0.02)&(np.abs(EV[:,0]-s*0.05)<0.015)]; y=sel[:,1].min()
        ms.append(ell([s*0.05,y-0.003,z],[0.028,0.006,0.035]))
        ms.append(tube([[s*0.05,y-0.007,z],[s*0.05,y-0.016,z]],0.0035,0.0025,seg=8))
    add('mamelles'+S,'o',trimesh.util.concatenate(ms))
pickle.dump((dict(OUT),LAY),open('net_female.pkl','wb'))
for k,v in OUT.items(): print(k,sum(len(m.faces) for m in v))
