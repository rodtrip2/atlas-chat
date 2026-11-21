# Repères osseux calculés sur la géométrie réelle des os
import pickle,json,numpy as np,trimesh
from collections import defaultdict
parts=pickle.load(open('bones_parts.pkl','rb')); lab=json.load(open('data/bone_labels.json'))
G=defaultdict(list)
for i,p in enumerate(parts): G[lab[str(i)]].append(p)
Bm={k:trimesh.util.concatenate(v) for k,v in G.items()}
CAN={k:(np.array(v[0]) if v else None) for k,v in json.load(open('data/canal.json')).items()}
def unit(v): v=np.asarray(v,float); return v/np.linalg.norm(v)
def ax(m):
    v=np.asarray(m.vertices); c=v.mean(0); d=np.linalg.svd(v-c,full_matrices=False)[2][0]; t=(v-c)@d
    a,b=c+d*t.min(),c+d*t.max(); return (a,b) if a[1]>b[1] else (b,a)
LM={}
def put(name,p): LM[name]=[float(x) for x in p]
def pick(m,score,mask=None):
    v=np.asarray(m.vertices); s=score(v)
    if mask is not None: s=np.where(mask(v),s,-1e9)
    return v[np.argmax(s)]
for sign,S in [(1,'_g'),(-1,'_d')]:
    MED=np.array([-sign,0,0]); LAT=-MED
    def frame(name):
        m=Bm[name+S]; p,d=ax(m); u=unit(d-p); L=np.linalg.norm(d-p)
        c=unit(np.array([0,-u[2],u[1]])); c = c if c[2]>0 else -c
        t=lambda v:((v-p)@u)/L
        return m,p,d,u,c,t
    # humerus
    m,p,d,u,c,t=frame('humerus')
    put('hum_tete'+S,pick(m,lambda v:-(v@c)-0.5*t(v),lambda v:t(v)<0.15))
    put('hum_tub_majeur'+S,pick(m,lambda v:(v@c)-0.6*t(v)+0.3*(v@LAT),lambda v:t(v)<0.15))
    put('hum_epi_lat'+S,pick(m,lambda v:v@LAT,lambda v:t(v)>0.85))
    put('hum_epi_med'+S,pick(m,lambda v:v@MED,lambda v:t(v)>0.85))
    put('hum_foramen'+S,pick(m,lambda v:v@MED,lambda v:(t(v)>0.72)&(t(v)<0.82)))
    # radius / ulna
    m,p,d,u,c,t=frame('ulna')
    put('uln_olecrane'+S,pick(m,lambda v:-t(v)))
    put('uln_styloide'+S,pick(m,lambda v:t(v)+0.3*(v@LAT),lambda v:t(v)>0.9))
    m,p,d,u,c,t=frame('radius')
    put('rad_tete'+S,pick(m,lambda v:-t(v)))
    put('rad_styloide'+S,pick(m,lambda v:(v@MED)+0.5*t(v),lambda v:t(v)>0.9))
    # scapula
    m=Bm['scapula'+S]; v=np.asarray(m.vertices)
    put('scap_glenoide'+S,v[np.argmin(v[:,1])])
    put('scap_bord_dorsal'+S,v[np.argmax(v[:,1])])
    lat=v@LAT; hi=v[lat>np.percentile(lat,93)]
    put('scap_epine'+S,hi[np.argsort(hi[:,1])[len(hi)//2]])
    put('scap_acromion'+S,hi[np.argmin(hi[:,1])])
    # femur
    m,p,d,u,c,t=frame('femur')
    put('fem_tete'+S,pick(m,lambda v:(v@MED)-0.3*t(v),lambda v:t(v)<0.14))
    put('fem_gd_troch'+S,pick(m,lambda v:(v@LAT)-1.0*t(v),lambda v:t(v)<0.14))
    put('fem_pt_troch'+S,pick(m,lambda v:(v@MED)-0.8*(v@c),lambda v:(t(v)>0.08)&(t(v)<0.24)))
    put('fem_trochlee'+S,pick(m,lambda v:v@c,lambda v:t(v)>0.86))
    put('fem_cond_med'+S,pick(m,lambda v:-(v@c)+0.8*(v@MED),lambda v:t(v)>0.86))
    put('fem_cond_lat'+S,pick(m,lambda v:-(v@c)+0.8*(v@LAT),lambda v:t(v)>0.86))
    # tibia / fibula
    m,p,d,u,c,t=frame('tibia')
    put('tib_tuberosite'+S,pick(m,lambda v:(v@c)-0.3*t(v),lambda v:t(v)<0.18))
    put('tib_malleole'+S,pick(m,lambda v:(v@MED)+0.5*t(v),lambda v:t(v)>0.9))
    m,p,d,u,c,t=frame('fibula')
    put('fib_malleole'+S,pick(m,lambda v:(v@LAT)+0.5*t(v),lambda v:t(v)>0.9))
    # calcaneus
    m=Bm['calcaneus'+S]; v=np.asarray(m.vertices); put('calc_tuber'+S,v[np.argmax(-v[:,2]+0.5*v[:,1])])
    # pelvis (single mesh)
    m=Bm['coxal']; v=np.asarray(m.vertices); sd=v[:,0]*sign>0.01
    vv=v[sd]; put('cox_ilium'+S,vv[np.argmax(vv[:,2]+0.5*vv[:,1])]); put('cox_ischion'+S,vv[np.argmax(-vv[:,2]-0.3*vv[:,1])])
    fh=np.array(LM['fem_tete'+S]); put('cox_acetabulum'+S,vv[np.argmin(np.linalg.norm(vv-fh,axis=1))])
    sym=v[np.abs(v[:,0])<0.015]; symp=sym[np.argmin(sym[:,1])]
    put('cox_obture'+S,(np.array(LM['cox_acetabulum'+S])*0.45+symp*0.25+np.array(LM['cox_ischion'+S])*0.30))
    # skull
    m=Bm['crane']; v=np.asarray(m.vertices); s_=v[:,0]*sign>0.005
    vv=v[s_]; put('cr_arcade'+S,vv[np.argmax(np.abs(vv[:,0]))])
    reg=vv[(vv[:,2]>1.53)&(vv[:,2]<1.66)&(np.abs(vv[:,0])>0.04)&(np.abs(vv[:,0])<0.15)]; put('cr_bulle'+S,reg[np.argmin(reg[:,1])])
    reg=vv[(np.abs(vv[:,0])>0.02)&(np.abs(vv[:,0])<0.07)&(vv[:,1]>0.62)&(vv[:,1]<0.76)]; put('cr_condyle'+S,reg[np.argmin(reg[:,2])])
    reg=vv[(vv[:,1]>0.74)&(vv[:,2]>1.74)&(vv[:,2]<1.85)]; po=reg[np.argmax(np.abs(reg[:,0]))]; put('cr_postorb'+S,po)
    ar=np.array(LM['cr_arcade'+S]); put('cr_orbite'+S,(po+ar)/2+np.array([-sign*0.03,0.02,0.02]))
    tgt=np.array([sign*0.07,0.55,1.89]); put('cr_infraorb'+S,vv[np.argmin(np.linalg.norm(vv-tgt,axis=1))])
    tgt=np.array([sign*0.145,0.62,1.60]); put('cr_meat'+S,vv[np.argmin(np.linalg.norm(vv-tgt,axis=1))])
    # mandible
    m=Bm['mandibule']; v=np.asarray(m.vertices); vv=v[v[:,0]*sign>0.005]
    put('md_coronoide'+S,vv[np.argmax(vv[:,1])])
    reg=vv[(vv[:,1]>0.53)&(vv[:,1]<0.63)]; put('md_condyle'+S,reg[np.argmin(reg[:,2])])
    reg=vv[vv[:,1]<0.5]; put('md_angulaire'+S,reg[np.argmin(reg[:,2])])
    tgt=np.array([sign*0.06,0.29,1.85]); put('md_mental'+S,vv[np.argmin(np.linalg.norm(vv-tgt,axis=1))])
    # vertebrae
    for vb,nm in [('C1','atlas_aile'),('L4','L4_transverse'),('sacrum','sacrum_aile')]:
        v=np.asarray(Bm[vb].vertices); vv=v[v[:,0]*sign>0]; put(nm+S,vv[np.argmax(np.abs(vv[:,0]))])
v=np.asarray(Bm['crane'].vertices)
mid=v[np.abs(v[:,0])<0.012]; put('cr_crete_sag',mid[np.argmax(mid[:,1])])
hi=mid[mid[:,1]>0.75]; put('cr_protub_occ',hi[np.argmin(hi[:,2])])
put('cr_foramen_magnum',[0,float(CAN['C1'][1]),float(v[(np.abs(v[:,0])<0.03)&(np.abs(v[:,1]-CAN['C1'][1])<0.03)][:,2].min())+0.008])
v=np.asarray(Bm['mandibule'].vertices); mid=v[np.abs(v[:,0])<0.015]; put('md_symphyse',mid[np.argmax(mid[:,2])])
v=np.asarray(Bm['C2'].vertices); put('axis_dent',v[np.argmax(v[:,2])])
for vb in ['T1','T11']:
    v=np.asarray(Bm[vb].vertices); put(vb+'_epineux',v[np.argmax(v[:,1])])
m=Bm['coxal']; v=np.asarray(m.vertices); sym=v[np.abs(v[:,0])<0.015]; put('cox_symphyse',sym[np.argmin(sym[:,1])])
json.dump(LM,open('data/landmarks.json','w'),indent=0)
print(len(LM),'repères'); 
for k in ['fem_tete_g','fem_gd_troch_g','hum_tete_g','uln_olecrane_g','scap_acromion_g','cr_orbite_g','md_coronoide_g','cr_foramen_magnum','cox_obture_g']: print(k,np.round(LM[k],3))
