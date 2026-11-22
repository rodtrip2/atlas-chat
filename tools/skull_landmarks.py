import numpy as np, json
def skull_lms(crane_v, mand_v, canal_c1):
    LM={}
    v=crane_v
    for sign,S in [(1,'_g'),(-1,'_d')]:
        vv=v[v[:,0]*sign>0.005]
        LM['cr_arcade'+S]=vv[np.argmax(np.abs(vv[:,0]))]
        reg=vv[(vv[:,2]>1.53)&(vv[:,2]<1.66)&(np.abs(vv[:,0])>0.04)&(np.abs(vv[:,0])<0.15)]; LM['cr_bulle'+S]=reg[np.argmin(reg[:,1])]
        reg=vv[(np.abs(vv[:,0])>0.02)&(np.abs(vv[:,0])<0.07)&(vv[:,1]>0.62)&(vv[:,1]<0.76)]; LM['cr_condyle'+S]=reg[np.argmin(reg[:,2])]
        reg=vv[(vv[:,1]>0.74)&(vv[:,2]>1.74)&(vv[:,2]<1.85)]; po=reg[np.argmax(np.abs(reg[:,0]))]; LM['cr_postorb'+S]=po
        LM['cr_orbite'+S]=(po+LM['cr_arcade'+S])/2+np.array([-sign*0.03,0.02,0.02])
        t=np.array([sign*0.07,0.55,1.89]); LM['cr_infraorb'+S]=vv[np.argmin(np.linalg.norm(vv-t,axis=1))]
        t=np.array([sign*0.145,0.62,1.60]); LM['cr_meat'+S]=vv[np.argmin(np.linalg.norm(vv-t,axis=1))]
        mv=mand_v[mand_v[:,0]*sign>0.005]
        LM['md_coronoide'+S]=mv[np.argmax(mv[:,1])]
        reg=mv[(mv[:,1]>0.53)&(mv[:,1]<0.63)]; LM['md_condyle'+S]=reg[np.argmin(reg[:,2])]
        reg=mv[mv[:,1]<0.5]; LM['md_angulaire'+S]=reg[np.argmin(reg[:,2])]
        t=np.array([sign*0.06,0.29,1.85]); LM['md_mental'+S]=mv[np.argmin(np.linalg.norm(mv-t,axis=1))]
    mid=v[np.abs(v[:,0])<0.012]; LM['cr_crete_sag']=mid[np.argmax(mid[:,1])]
    hi=mid[mid[:,1]>0.75]; LM['cr_protub_occ']=hi[np.argmin(hi[:,2])]
    LM['cr_foramen_magnum']=np.array([0,canal_c1[1],v[(np.abs(v[:,0])<0.03)&(np.abs(v[:,1]-canal_c1[1])<0.03)][:,2].min()+0.008])
    mm=mand_v[np.abs(mand_v[:,0])<0.015]; LM['md_symphyse']=mm[np.argmax(mm[:,2])]
    return {k:[float(x) for x in p] for k,p in LM.items()}
