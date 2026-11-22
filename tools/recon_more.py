# Système lymphatique, tendons/ligaments/ménisques/disques, glandes et organes manquants
from common import *
from common import _G, parts
import pickle,json
from collections import defaultdict
from shapely.geometry import MultiPoint
import shapely
OUT=defaultdict(list); LAY={}
def add(id,l,m): OUT[id].append(m); LAY[id]=l
def unit(v): v=np.asarray(v,float); return v/np.linalg.norm(v)
def ell(c,r,sub=2,R=None):
    m=trimesh.creation.icosphere(subdivisions=sub,radius=1.0); m.vertices*=np.asarray(r)
    if R is not None: m.vertices=m.vertices@R.T
    m.apply_translation(c); return m
def ell_ax(c,r_long,r_small,d):
    d=unit(d); a=np.array([0,1,0]) if abs(d[1])<0.9 else np.array([1,0,0]); e1=unit(np.cross(d,a)); e2=np.cross(d,e1)
    R=np.c_[e1,e2,d]; return ell(c,[r_small,r_small*0.8,r_long],R=R)
def flatten(m,center_pts,axis,f):
    v=np.asarray(m.vertices); C=np.asarray(center_pts); from scipy.spatial import cKDTree
    i=cKDTree(C).query(v)[1]; off=v-C[i]; off[:,axis]*=f; m.vertices=C[i]+off; return m
LM=json.load(open('data/landmarks.json')); O=organs()
_lab2=json.load(open('data/bone_labels2.json')); _G2=defaultdict(list)
for _i,_p in enumerate(parts): _G2[_lab2[str(_i)]].append(_p)
B2={k:trimesh.util.concatenate(v) for k,v in _G2.items()}
P_=lambda a,b,t,off: a+(b-a)*t+off
# ---------- limb frames ----------
def frames(S,sign):
    def A(name):
        a,b=axis(B[name+S]); return (a,b) if a[1]>b[1] else (b,a)
    F={}; F['hp'],F['hd']=A('humerus'); F['rp'],F['rd']=A('radius'); F['fp'],F['fd']=A('femur'); F['tp'],F['td']=A('tibia'); F['bp'],F['bd']=A('fibula')
    def perps(p,d):
        dd=unit(d-p); n1=unit(np.array([0,-dd[2],dd[1]]))
        if n1[2]<0: n1=-n1
        return n1,-n1
    F['hc'],F['hk']=perps(F['hp'],F['hd']); F['rc'],F['rk']=perps(F['rp'],F['rd']); F['fc'],F['fk']=perps(F['fp'],F['fd']); F['tc'],F['tk']=perps(F['tp'],F['td'])
    if F['tk'][1]<0: F['tk'],F['tc']=-F['tk'],-F['tc']
    F['MED']=np.array([-sign,0,0]); F['LAT']=-F['MED']; return F
# =================== LYMPHATIC SYSTEM ===================
LN=lambda: None
venous=lambda s: np.array([s*0.05,0.42,1.06])
for s,S in [(1,'_g'),(-1,'_d')]:
    F=frames(S,s); sv=lambda p: np.array([p[0]*s,p[1],p[2]])
    nodes={'ln_mandibulaire':[ell(sv([0.11,0.46,1.66]),[0.011,0.008,0.015]),ell(sv([0.10,0.44,1.69]),[0.009,0.007,0.012])],
           'ln_retropharyngien':[ell(sv([0.04,0.57,1.50]),[0.008,0.012,0.022])],
           'ln_parotidien':[ell(sv([0.125,0.60,1.565]),[0.006,0.006,0.010])],
           'ln_cervical_sup':[ell(sv([0.16,0.45,1.18]),[0.010,0.012,0.020])],
           'ln_axillaire':[ell(sv([0.12,0.30,0.98]),[0.010,0.008,0.014])],
           'ln_iliaque_med':[ell(sv([0.03,0.55,-0.56]),[0.008,0.008,0.025])],
           'ln_inguinal_sup':[ell(sv([0.07,0.31,-0.82]),[0.010,0.006,0.014])]}
    cond=(np.array(LM['fem_cond_med'+S])+np.array(LM['fem_cond_lat'+S]))/2
    pop=cond+np.array([0,0.71,-0.71])*0.07; nodes['ln_popliteal']=[ell(pop,[0.012,0.014,0.018])]
    for k,ms in nodes.items():
        for m in ms: add(k+S,'l',m)
    # drainage vessels
    lv=[]
    lv.append(tube([sv([0.11,0.46,1.66]),sv([0.08,0.52,1.57]),sv([0.04,0.57,1.50])],0.0022))
    lv.append(tube([sv([0.04,0.56,1.49]),sv([0.032,0.49,1.30]),sv([0.03,0.47,1.15]),venous(s)],0.0026))   # tracheal duct
    lv.append(tube([sv([0.16,0.45,1.18]),sv([0.10,0.43,1.10]),venous(s)],0.0024))
    lv.append(tube([sv([0.12,0.30,0.98]),sv([0.08,0.38,1.04]),venous(s)],0.0024))
    mc=B['metacarpe'+S].centroid if 'metacarpe'+S in B else B['mcIII'+S].centroid
    lv.append(tube([mc+np.array([0,0.04,0.03]),P_(F['rp'],F['rd'],0.5,F['rc']*0.04),P_(F['hp'],F['hd'],0.92,F['MED']*0.03+F['hc']*0.02),P_(F['hp'],F['hd'],0.5,F['MED']*0.05),sv([0.12,0.30,0.98])],0.0018))
    mt=B2['mtIII'+S].centroid
    lv.append(tube([mt+np.array([0,0.02,-0.02]),P_(F['tp'],F['td'],0.6,F['tk']*0.045),pop],0.0018))
    lv.append(tube([pop,P_(F['fp'],F['fd'],0.6,F['MED']*0.025),sv([0.11,0.36,-0.88]),sv([0.06,0.48,-0.70]),sv([0.03,0.55,-0.56])],0.0022))
    lv.append(tube([sv([0.07,0.31,-0.82]),sv([0.08,0.40,-0.74]),sv([0.03,0.55,-0.56])],0.0022))
    lv.append(tube([sv([0.03,0.55,-0.58]),sv([0.012,0.585,-0.32]),[-0.008,0.60,-0.12]],0.0028))
    add('vx_lymph'+S,'l',trimesh.util.concatenate(lv))
# central lymphatics
cist=np.array([-0.012,0.60,-0.10]); add('citerne_chyle','l',ell(cist,[0.010,0.008,0.05]))
td=[cist+np.array([0,0,0.05])]+[[-0.012,float(body_ymin(z))-0.012,z] for z in np.linspace(0.10,0.72,6)]+[[0.008,0.53,0.82],[0.03,0.50,0.96],[0.05,0.45,1.04],venous(1)]
add('conduit_thoracique','l',tube(td,0.0035,0.003))
add('ln_tracheobronchiques','l',trimesh.util.concatenate([ell([0.035,0.45,0.79],[0.009,0.008,0.012]),ell([-0.035,0.45,0.79],[0.009,0.008,0.012]),ell([0,0.415,0.77],[0.008,0.007,0.010])]))
add('ln_mediastinaux','l',trimesh.util.concatenate([ell([0.004,0.42,0.96],[0.008,0.007,0.014]),ell([0.012,0.40,1.01],[0.007,0.006,0.011])]))
add('ln_sternal','l',ell([0.0,0.235,0.95],[0.007,0.006,0.010]))
add('ln_jejunaux','l',trimesh.util.concatenate([ell_ax([-0.01,0.40,-0.20],0.05,0.012,[0.2,0.3,1]),ell([0.02,0.385,-0.10],[0.01,0.009,0.02])]))
add('ln_coliques','l',ell([-0.03,0.45,-0.40],[0.008,0.008,0.012]))
add('vx_lymph_centraux','l',trimesh.util.concatenate([tube([[-0.01,0.40,-0.20],[-0.01,0.52,-0.15],cist],0.0028),
   tube([[0,0.44,0.79],[0.01,0.47,0.90],[0.03,0.50,0.96]],0.0022),tube([[0,0.235,0.95],[0.004,0.33,0.98],[0.004,0.42,0.96]],0.0018)]))
add('thymus','l',trimesh.util.concatenate([ell([0.012,0.30,0.93],[0.022,0.035,0.065]),ell([-0.012,0.31,0.95],[0.018,0.03,0.055])]))
add('tonsilles','l',trimesh.util.concatenate([ell([s*0.025,0.50,1.60],[0.006,0.004,0.010]) for s in (1,-1)]))
# =================== TENDONS, LIGAMENTS, MENISCI, DISCS ===================
def bc(name): return (B2[name] if name in B2 else B[name]).centroid
for s,S in [(1,'_g'),(-1,'_d')]:
    F=frames(S,s); MED,LAT=F['MED'],F['LAT']
    pat=bc('patella'+S); tt=np.array(LM['tib_tuberosite'+S]); ct=np.array(LM['calc_tuber'+S]); ol=np.array(LM['uln_olecrane'+S])
    add('t_calcaneen'+S,'t',tube([P_(F['tp'],F['td'],0.42,F['tk']*0.058),P_(F['tp'],F['td'],0.75,F['tk']*0.05),ct+np.array([0,0.01,0.005])],0.0075,0.0055))
    add('t_patellaire'+S,'t',tube([pat+np.array([0,-0.02,-0.004]),(pat+tt)/2+F['tc']*0.004,tt],0.0055,0.005))
    add('t_quadriceps'+S,'t',tube([P_(F['fp'],F['fd'],0.72,F['fc']*0.05),P_(F['fp'],F['fd'],0.9,F['fc']*0.045),pat+np.array([0,0.018,0])],0.008,0.006))
    add('t_triceps'+S,'t',tube([P_(F['hp'],F['hd'],0.62,F['hk']*0.05),P_(F['hp'],F['hd'],0.85,F['hk']*0.04),ol+np.array([0,0.006,-0.004])],0.008,0.0065))
    sg=np.array(LM['scap_glenoide'+S]); tm=np.array(LM['hum_tub_majeur'+S]); rh=np.array(LM['rad_tete'+S])
    add('t_biceps'+S,'t',tube([sg+F['hc']*0.02+np.array([0,0.01,0]),tm+MED*0.5*0.045+F['hc']*0.006,P_(F['hp'],F['hd'],0.5,F['hc']*0.036),P_(F['hp'],F['hd'],0.92,F['hc']*0.03+MED*0.3*0.045),P_(F['rp'],F['rd'],0.1,MED*0.35*0.045)],0.0042,0.0035))
    # digital tendons (each digit, real bones)
    for limb,mcp,ph,src_f,src_e in [('m','mc','m',bc('carpe_acc'+S)+np.array([0,-0.005,0.03]),bc('carpe_rad'+S)+np.array([0,0.02,0.01])),
                                    ('p','mt','p',np.array(LM['calc_tuber'+S])*0.3+bc('tarseIV'+S)*0.7+np.array([0,-0.02,0.02]),bc('tarse_cen'+S)+np.array([0,0.025,0.01]))]:
        fl=[];ex=[]
        for d in ['II','III','IV','V']:
            try: m1=bc(f'{mcp}{d}{S}'); p1=bc(f'phal1_{ph}{d}{S}'); p2=bc(f'phal2_{ph}{d}{S}'); p3=bc(f'phal3_{ph}{d}{S}')
            except KeyError: continue
            fl.append(tube([src_f,m1+np.array([0,-0.016,0]),p1+np.array([0,-0.013,0]),p2+np.array([0,-0.011,0]),p3+np.array([0,-0.008,-0.012])],0.0024,0.0016))
            ex.append(tube([src_e,m1+np.array([0,0.016,0]),p1+np.array([0,0.013,0]),p2+np.array([0,0.011,0]),p3+np.array([0,0.012,-0.012])],0.0019,0.0013))
        add(f't_flechisseurs_{limb}'+S,'t',trimesh.util.concatenate(fl)); add(f't_extenseurs_{limb}'+S,'t',trimesh.util.concatenate(ex))
    # stifle: cruciates, collaterals, menisci
    cm=np.array(LM['fem_cond_med'+S]); cl=np.array(LM['fem_cond_lat'+S]); J=(cm+cl)/2
    T=F['tp']; u=unit(F['td']-F['tp'])
    add('lig_croise_cranial'+S,'t',tube([J+LAT*0.012*0.0+np.array([s*0.012,0,0])+F['fk']*0.01,T+F['tc']*0.012+np.array([-s*0.006,0,0])],0.0034))
    add('lig_croise_caudal'+S,'t',tube([J+np.array([-s*0.012,0,0])+F['fc']*0.008,T+F['tk']*0.014],0.0034))
    add('lig_collat_genou'+S,'t',trimesh.util.concatenate([tube([cm+MED*0.012,P_(F['tp'],F['td'],0.09,MED*0.024)],0.003),tube([cl+LAT*0.012,F['bp']+LAT*0.006],0.003)]))
    e1=unit(np.cross(u,[1,0,0])); e2=np.cross(u,e1); men=[]
    for side_off in (MED,LAT):
        c0=T+side_off*0.014-u*0.006
        pts=[c0+0.012*(np.cos(a)*e1+np.sin(a)*e2) for a in np.linspace(0.5,2*np.pi-0.5,14)]
        men.append(tube(pts,0.0032,0.0032,seg=8))
    add('menisques'+S,'t',trimesh.util.concatenate(men))
    # elbow collaterals
    add('lig_collat_coude'+S,'t',trimesh.util.concatenate([tube([np.array(LM['hum_epi_med'+S]),rh+MED*0.85*0.045+F['rk']*0.004],0.0028),tube([np.array(LM['hum_epi_lat'+S]),rh+LAT*0.15*0.045],0.0028)]))
    add('lig_tete_femorale'+S,'t',tube([np.array(LM['fem_tete'+S])+MED*0.1*0.045,np.array(LM['cox_acetabulum'+S])+MED*0.25*0.045],0.0025))
# supraspinous ligament (cat: no nuchal ligament)
seq=[f'T{i}' for i in range(1,14)]+[f'L{i}' for i in range(1,8)]
tips=[]
for v in seq:
    vv=np.asarray(B[v].vertices); tips.append(vv[np.argmax(vv[:,1])])
add('lig_supra_epineux','t',tube(tips+[np.asarray(B['sacrum'].vertices)[np.argmax(np.asarray(B['sacrum'].vertices)[:,1])]],0.0028,step=0.006))
# intervertebral discs
def body(vname):
    vv=np.asarray(B[vname].vertices); c=CAN.get(vname)
    cy=c[1] if c is not None else canal_y(B[vname].centroid[2])
    sel=vv[(np.abs(vv[:,0])<0.028)&(vv[:,1]<cy-0.015)]
    if len(sel)<20: sel=vv[np.abs(vv[:,0])<0.03]
    return sel
order=[f'C{i}' for i in range(2,8)]+[f'T{i}' for i in range(1,14)]+[f'L{i}' for i in range(1,8)]+['sacrum']
discs=[]
for a,b in zip(order[:-1],order[1:]):
    A_=body(a); B_=body(b)
    za=A_[:,2].min(); zb=B_[:,2].max()          # caudal end of A, cranial end of B
    ya=np.median(A_[A_[:,2]<za+0.02][:,1]); yb=np.median(B_[B_[:,2]>zb-0.02][:,1])
    c=np.array([0,(ya+yb)/2,(za+zb)/2]); th=max(0.004,abs(za-zb)); r=max(0.012,(np.ptp(A_[A_[:,2]<za+0.02][:,1]))*0.45)
    d=unit([0,ya-yb,za-zb]) if abs(za-zb)>1e-4 else np.array([0,0,1.0])
    discs.append(tube([c-d*th/2,c+d*th/2],r,r,seg=16))
add('disques','t',trimesh.util.concatenate(discs))
# =================== GLANDS & OTHER ORGANS ===================
add('thyroide','o',trimesh.util.concatenate([ell([s*0.031,0.50,1.38],[0.006,0.008,0.025]) for s in (1,-1)]))
add('parathyroides','o',trimesh.util.concatenate([ell([s*0.037,0.505,1.40],[0.0032,0.0032,0.0032],sub=1) for s in (1,-1)]))
add('surrenales','o',trimesh.util.concatenate([ell([0.05,0.53,0.07],[0.006,0.006,0.012]),ell([-0.055,0.52,0.13],[0.006,0.006,0.012])]))
for s,S in [(1,'_g'),(-1,'_d')]:
    add('gl_mandibulaire'+S,'o',ell([s*0.10,0.505,1.59],[0.016,0.016,0.022]))
    add('gl_parotide'+S,'o',ell([s*0.132,0.58,1.585],[0.010,0.03,0.02]))
tg=[[0,0.47,1.58],[0,0.41,1.66],[0,0.35,1.76],[0,0.32,1.86],[0,0.31,1.90]]
t_=tube(tg,0.026,0.012,seg=16); t_=flatten(t_,spline(tg,0.005),1,0.5); add('langue','o',t_)
# diaphragm sheet from the computed diaphragm surface
D=pickle.load(open('diaph.pkl','rb')); zd,xs,ys=D['zd'],D['xs'],D['ys']
W=np.vstack([np.asarray(B2[k].vertices) for k in B2 if k.startswith(('cote','cart')) or k in('sternebre7','xiphoide')])
pts=[];idx={}
for i,x in enumerate(xs):
    for j,y in enumerate(ys):
        z=zd[i,j]; sl=W[np.abs(W[:,2]-z)<0.05]
        if len(sl)<20: continue
        idx[(i,j)]=None
grid_ok=np.zeros(zd.shape,bool)
for (i,j) in idx:
    z=zd[i,j]; sl=W[np.abs(W[:,2]-z)<0.05]; hull=MultiPoint(sl[:,:2]).convex_hull.buffer(-0.012)
    grid_ok[i,j]=shapely.contains_xy(hull,xs[i],ys[j]) and ys[j]<float(body_ymin(z))-0.005
V=[];Fc=[];vid=-np.ones(zd.shape,int)
for i in range(len(xs)):
    for j in range(len(ys)):
        if grid_ok[i,j]: vid[i,j]=len(V); V.append([xs[i],ys[j],zd[i,j]+0.004])
for i in range(len(xs)-1):
    for j in range(len(ys)-1):
        a,b,c,d=vid[i,j],vid[i+1,j],vid[i,j+1],vid[i+1,j+1]
        if min(a,b,c,d)>=0: Fc+= [[a,b,d],[a,d,c]]
dia=trimesh.Trimesh(np.array(V),np.array(Fc),process=True); trimesh.smoothing.filter_taubin(dia,iterations=10)
add('diaphragme','o',dia)
pickle.dump((dict(OUT),LAY),open('net_more.pkl','wb'))
for k,v in OUT.items(): print(k,LAY[k],sum(len(m.faces) for m in v))
