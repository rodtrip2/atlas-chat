from common import *
import pickle
from collections import defaultdict
OUT=defaultdict(list); LAY={}
def add(id,layer,m): OUT[id].append(m); LAY[id]=layer
def unit(v): v=np.asarray(v,float); return v/np.linalg.norm(v)
TH=pickle.load(open('recon_thorax.pkl','rb'))
def aorta_pt(z): return np.array([0.02,float(body_ymin(z))-0.028,z])
# ================= A. intercostals =================
for sign,S in [(1,'_g'),(-1,'_d')]:
    for k in range(1,14):
        rib=centerline(B[f'cote{k}{S}'])
        if rib[0][1]<rib[-1][1]: rib=rib[::-1]
        pts=list(rib)
        if f'cart{k}{S}' in B:
            ca=centerline(B[f'cart{k}{S}'])
            if np.linalg.norm(ca[0]-rib[-1])>np.linalg.norm(ca[-1]-rib[-1]): ca=ca[::-1]
            pts+=list(ca[1:-2])
        pts=np.array(pts)[::2]
        # plane normal of rib, oriented caudally
        v=np.asarray(B[f'cote{k}{S}'].vertices); n=np.linalg.svd(v-v.mean(0),full_matrices=False)[2][2]
        if n[2]>0: n=-n
        def off(p,a,b):
            inw=np.array([0,0.33,p[2]])-p; inw[2]=0; inw=unit(inw)
            return p+n*a+inw*b
        z0=pts[0][2]
        art=[aorta_pt(z0),[sign*0.03,body_ymin(z0)+0.005,z0]]+[off(p,0.013,0.007) for p in pts]
        vei=[[-0.022,body_ymin(z0)-0.012,z0],[sign*0.03,body_ymin(z0)+0.012,z0]]+[off(p,0.009,0.005) for p in pts]
        ner=[[sign*0.012,canal_y(z0),z0],[sign*0.035,canal_y(z0)-0.02,z0-0.008]]+[off(p,0.017,0.006) for p in pts]
        add('art_intercost'+S,'v',tube(art,0.0028,0.0015)); add('v_intercost'+S,'v',tube(vei,0.0032,0.0018)); add('n_intercost'+S,'n',tube(ner,0.0022,0.0013))
# ================= B. internal thoracic + epigastric =================
SV=np.vstack([np.asarray(B[k].vertices) for k in ['manubrium','sternebre','xiphoide']])
def ster_y(z):
    sel=SV[np.abs(SV[:,2]-z)<0.03]; return sel[:,1].max() if len(sel) else 0.05
_ec=trimesh.load('sources/feline_ecorche_dec.glb').to_geometry(); EV=np.asarray(_ec.vertices)
def wall_y(z,x):
    sel=EV[(np.abs(EV[:,2]-z)<0.02)&(np.abs(EV[:,0]-x)<0.02)]; return sel[:,1].min()+0.024
for sign,S in [(1,'_g'),(-1,'_d')]:
    s_=lambda p: np.array([p[0]*sign,p[1],p[2]])
    zz=np.linspace(0.95,0.50,7)
    th=[s_([0.05,0.425,1.005]),s_([0.042,0.36,1.0]),s_([0.032,ster_y(0.97)+0.03,0.975])]+[s_([0.024,ster_y(z)+0.016,z]) for z in zz]
    # musculophrenic branch along the costal arch (cartilages 8 to 10)
    arch=[]
    for k in [8,9,10]:
        c=centerline(B[f'cart{k}{S}'])
        if abs(c[0][0])>abs(c[-1][0]): c=c[::-1]          # medial -> lateral
        arch+= list(c[len(c)//3::3])
    arch=[p+np.array([-sign*0.01,0.012,-0.006]) for p in arch]
    mp=[th[-1]]+arch
    # cranial epigastric inside the ventral abdominal wall, meets the caudal epigastric
    ze=[0.40,0.25,0.10,-0.05,-0.20,-0.35]
    ce=[th[-1]]+[s_([0.03+0.004*i,wall_y(z,0.035),z]) for i,z in enumerate(ze)]
    zc=[-0.35,-0.50,-0.65,-0.75]
    cd=[s_([0.054,wall_y(z,0.055),z]) for z in zc]+[s_([0.085,0.33,-0.83]),s_([0.11,0.355,-0.87])]
    add('art_thor_int'+S,'v',trimesh.util.concatenate([tube(th,0.004,0.0032),tube(mp,0.0024,0.0014),ball(th[-1],0.0034)]))
    add('art_epigastrique'+S,'v',tube(ce+cd[1:],0.0028,0.0028))
    vo=lambda pts,d=0.007: [np.asarray(p)+np.array([sign*d,0.002,0]) for p in pts]
    add('v_thor_int'+S,'v',trimesh.util.concatenate([tube(vo(th),0.0045,0.0036),tube(vo(mp,0.006),0.0026,0.0016)]))
    add('v_epigastrique'+S,'v',tube(vo(ce+cd[1:]),0.0032,0.0032))
# ================= C. azygos =================
az=[[-0.02,float(body_ymin(z))-0.012,z] for z in np.linspace(-0.05,0.86,10)]+[[-0.025,0.47,0.93],[-0.03,0.43,0.92]]
add('azygos','v',tube(az,0.0055,0.0045))
# ================= D. sympathetic chain =================
vseq=[f'T{i}' for i in range(1,14)]+[f'L{i}' for i in range(1,8)]
for sign,S in [(1,'_g'),(-1,'_d')]:
    zz=[B[v].centroid[2] for v in vseq]+[-0.80]
    pts=[[sign*0.03,float(body_ymin(z))+0.016,z] for z in zz]
    parts=[tube([[sign*0.045,0.50,1.30],[sign*0.04,0.49,1.15]]+pts,0.0025)]
    for i,p in enumerate(pts): parts.append(ball(p,0.0085 if i==0 else 0.0042))
    add('sympathique'+S,'n',trimesh.util.concatenate(parts))
# ================= E. spinal nerves (cervical, lumbar, caudal) =================
def jz(a,b): return (CAN[a][2]+CAN[b][2])/2
juncs=[(jz('C1','C2'),'c'),(jz('C2','C3'),'c'),(jz('C3','C4'),'c'),(jz('C4','C5'),'c'),(jz('C5','C6'),'c'),
       (jz('L1','L2'),'l'),(jz('L2','L3'),'l'),(jz('L3','L4'),'l'),(B['Cd1'].centroid[2],'s'),(B['Cd2'].centroid[2],'s')]
for sign,S in [(1,'_g'),(-1,'_d')]:
    for z,t in juncs:
        cy=canal_y(z) if t!='s' else 0.76
        if t=='l': p=[[sign*0.012,cy,z],[sign*0.05,cy-0.03,z-0.01],[sign*0.10,cy-0.10,z-0.04],[sign*0.15,cy-0.22,z-0.07]]
        else: p=[[sign*0.012,cy,z],[sign*0.045,cy-0.025,z-0.01],[sign*0.075,cy-0.07,z-0.03]]
        add('n_spinaux'+S,'n',tube(p,0.0026,0.0014))
        add('n_spinaux'+S,'n',tube([p[1],[sign*0.06,cy+0.04,z-0.02]],0.0018,0.0012))  # dorsal ramus
# ================= F. lumbar arteries =================
for sign,S in [(1,'_g'),(-1,'_d')]:
    for v in ['L1','L2','L3','L4','L5','L6']:
        z=B[v].centroid[2]; y0=float(body_ymin(z))-0.03
        add('art_lomb'+S,'v',tube([[-0.01,y0,z],[sign*0.035,body_ymin(z)+0.01,z],[sign*0.07,canal_y(z)+0.0,z-0.01],[sign*0.06,canal_y(z)+0.06,z-0.02]],0.0026,0.0015))
pickle.dump((dict(OUT),LAY),open('net_trunk.pkl','wb'))
for k,v in OUT.items(): print(k,LAY[k],sum(len(m.faces) for m in v))
