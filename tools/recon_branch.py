from common import *
from common import _G
import pickle
from collections import defaultdict
OUT=defaultdict(list); LAY={}
def add(id,l,m): OUT[id].append(m); LAY[id]=l
def unit(v): v=np.asarray(v,float); return v/np.linalg.norm(v)
def P(a,b,t,off): return a+(b-a)*t+off
for sign,S in [(1,'_g'),(-1,'_d')]:
    s=lambda p: np.array([p[0]*sign,p[1],p[2]])
    # ---------------- head ----------------
    car=s([0.045,0.53,1.47]); jug=s([0.08,0.50,1.47])
    add('art_faciale'+S,'v',tube([car,s([0.07,0.47,1.56]),s([0.09,0.42,1.64]),s([0.10,0.37,1.705]),s([0.14,0.42,1.78]),s([0.12,0.47,1.86]),s([0.08,0.49,1.93])],0.0035,0.002))
    add('art_linguale'+S,'v',tube([car,s([0.04,0.48,1.56]),s([0.03,0.42,1.64]),s([0.02,0.355,1.72]),s([0.012,0.305,1.82]),s([0.008,0.295,1.88])],0.0032,0.0018))
    add('art_maxillaire'+S,'v',tube([car,s([0.08,0.56,1.54]),s([0.12,0.57,1.60]),s([0.10,0.555,1.70]),s([0.08,0.56,1.80]),s([0.07,0.55,1.88])],0.0038,0.002))
    add('art_temporale'+S,'v',tube([s([0.08,0.56,1.54]),s([0.13,0.62,1.55]),s([0.16,0.72,1.58]),s([0.15,0.82,1.63])],0.0026,0.0016))
    add('v_linguofaciale'+S,'v',tube([jug,s([0.09,0.52,1.53]),s([0.11,0.49,1.59]),s([0.12,0.45,1.65]),s([0.115,0.40,1.70]),s([0.115,0.36,1.76]),s([0.14,0.42,1.84]),s([0.11,0.48,1.90])],0.0055,0.003))
    add('v_maxillaire'+S,'v',tube([jug,s([0.11,0.55,1.52]),s([0.14,0.62,1.56]),s([0.14,0.70,1.60])],0.005,0.003))
    fac0=s([0.10,0.60,1.55]); fac1=s([0.16,0.60,1.60])
    add('n_facial'+S,'n',trimesh.util.concatenate([tube([fac0,fac1],0.003),
        tube([fac1,s([0.205,0.55,1.70]),s([0.175,0.52,1.80]),s([0.12,0.515,1.89])],0.0022,0.0013),
        tube([fac1,s([0.178,0.475,1.70]),s([0.155,0.425,1.80]),s([0.12,0.41,1.88])],0.0022,0.0013),
        tube([fac1,s([0.18,0.68,1.64]),s([0.195,0.74,1.72])],0.0018,0.0012)]))
    io=[s([0.08,0.58,1.74]),s([0.07,0.55,1.86]),s([0.065,0.52,1.92])]
    add('n_infraorbitaire'+S,'n',trimesh.util.concatenate([tube(io,0.003,0.0024)]+[tube([io[-1],s(e)],0.0016,0.0009) for e in
        [[0.10,0.47,1.97],[0.105,0.52,1.98],[0.095,0.57,1.97],[0.08,0.44,1.95],[0.085,0.60,1.955]]]))
    add('n_optique'+S,'n',tube([s([0.008,0.70,1.74]),s([0.05,0.72,1.78]),s([0.09,0.73,1.81])],0.003,0.0028))
    add('n_alveolaire'+S,'n',tube([s([0.12,0.57,1.60]),s([0.14,0.53,1.635]),s([0.143,0.505,1.655]),s([0.122,0.463,1.68]),s([0.103,0.423,1.71]),s([0.08,0.385,1.74]),s([0.052,0.305,1.77]),s([0.035,0.272,1.80]),s([0.045,0.28,1.84])],0.0024,0.0015))
    add('n_hypoglosse'+S,'n',tube([s([0.03,0.62,1.53]),s([0.05,0.49,1.58]),s([0.035,0.40,1.68]),s([0.015,0.34,1.80])],0.0024,0.0016))
    add('n_lingual'+S,'n',tube([s([0.11,0.55,1.62]),s([0.075,0.45,1.70]),s([0.04,0.37,1.76]),s([0.02,0.32,1.84])],0.002,0.0013))
    # ---------------- limbs ----------------
    def A(name):
        a,b=axis(B[name+S]); return (a,b) if a[1]>b[1] else (b,a)
    hp,hd=A('humerus'); rp,rd=A('radius'); up,ud=A('ulna'); fp,fd=A('femur'); tp,td=A('tibia')
    mc=B['metacarpe'+S].centroid; mt=B['metatarse'+S].centroid
    med=np.array([-0.045*sign,0,0]); lat=-med
    def perps(p,d):
        dd=unit(d-p); n1=unit(np.array([0,-dd[2],dd[1]]))
        if n1[2]<0: n1=-n1
        return n1,-n1
    hc,hk=perps(hp,hd); rc,rk=perps(rp,rd); fc,fk=perps(fp,fd); tc,tk=perps(tp,td)
    if tk[1]<0: tk,tc=-tk,-tc
    # forelimb collaterals
    add('art_brach_prof'+S,'v',tube([P(hp,hd,0.3,med),P(hp,hd,0.36,hk*0.06),P(hp,hd,0.48,hk*0.085+lat*0.5)],0.0032,0.0018))
    add('art_coll_ulnaire'+S,'v',tube([P(hp,hd,0.68,med),P(hp,hd,0.9,med*0.8+hk*0.03),P(up,ud,0.06,med*0.4+rk*0.04),P(up,ud,0.25,rk*0.035)],0.0024,0.0014))
    add('art_brach_sup'+S,'v',tube([P(hp,hd,0.75,med),P(hp,hd,0.97,hc*0.035),P(rp,rd,0.4,rc*0.03+lat*0.3),P(rp,rd,0.8,rc*0.028+lat*0.2)],0.0024,0.0014))
    bv=[s([0.03,0.465,0.97])*[1,1,1],s([0.07,0.375,1.04]),P(hp,hd,0.08,med*1.3+hk*0.028),P(hp,hd,0.4,med+hk*0.009),P(hp,hd,0.8,med+hk*0.008),P(rp,rd,0.1,med*0.9+rk*0.024),P(rp,rd,0.5,med*0.8+rk*0.028)]
    add('v_brachiale'+S,'v',tube(bv,0.0068,0.0035))
    add('v_mediane_cub'+S,'v',tube([P(rp,rd,0.05,rc*0.035),P(hp,hd,0.98,med*0.6+hc*0.01),P(hp,hd,0.8,med+hk*0.008)],0.003))
    # hindlimb collaterals
    add('art_fem_prof'+S,'v',tube([s([0.12,0.36,-0.88]),P(fp,fd,0.2,med*0.6+fk*0.04),P(fp,fd,0.3,med*0.2+fk*0.08)],0.0034,0.002))
    add('art_saphene'+S,'v',tube([P(fp,fd,0.6,med),P(fp,fd,0.9,med*1.25),P(tp,td,0.3,med*1.1+tk*0.01),P(tp,td,0.8,med*0.85+tk*0.016),mt+np.array([-0.01*sign,-0.01,-0.03])],0.0028,0.0016))
    add('art_tib_cran'+S,'v',tube([P(tp,td,0.12,tk*0.03),P(tp,td,0.25,lat*0.5+tc*0.01),P(tp,td,0.6,tc*0.02+lat*0.4),P(tp,td,0.95,tc*0.02+lat*0.2),mt+np.array([0,0.02,0.015])],0.0028,0.0016))
    fv=[s([0.055,0.335,-0.835]),s([0.125,0.365,-0.885]),P(fp,fd,0.3,med*1.15+fc*0.008),P(fp,fd,0.6,med*1.15),P(fp,fd,0.92,med*0.4+fk*0.048),P(tp,td,0.15,tk*0.045),P(tp,td,0.5,med*0.4+tk*0.038)]
    add('v_femorale'+S,'v',tube(fv,0.0068,0.0035))
    add('v_saphene_lat'+S,'v',tube([mt+np.array([0.02*sign,0.02,-0.03]),P(tp,td,0.9,lat*0.9+tk*0.02),P(tp,td,0.4,lat*0.85+tk*0.035),P(tp,td,0.05,lat*0.4+tk*0.05),P(fp,fd,0.95,fk*0.055)],0.0035,0.0042))
    # digital branches, one per digit, using the real phalanges
    for limb,base_a,base_n,ph in [('mt',mc+np.array([0,-0.02,-0.01]),mc+np.array([0,-0.025,0]),'m'),('mp',mt+np.array([0,0.01,-0.02]),mt+np.array([0,-0.01,-0.02]),'p')]:
        p3=_G[f'phal3_{ph}{S}']; p1=_G[f'phal1_{ph}{S}']; p2=_G[f'phal2_{ph}{S}']
        c1=np.array([m.centroid for m in p1]); c2=np.array([m.centroid for m in p2])
        arts=[];ners=[]
        for m3 in p3:
            c3=m3.centroid; a1=c1[np.argmin(np.linalg.norm(c1-c3,axis=1))]; a2=c2[np.argmin(np.linalg.norm(c2-c3,axis=1))]
            pal=np.array([0,-0.012,0])
            for k,(lst,b0,dx,r) in enumerate([(arts,base_a,0.004,0.0016),(ners,base_n,-0.004,0.0013)]):
                o=np.array([dx*sign,0,0])
                lst.append(tube([b0,(b0+a1)/2+pal+o,a1+pal+o,a2+pal*0.8+o,c3+pal*0.5+o],r,r*0.6))
        add(f'art_digit_{limb}'+S,'v',trimesh.util.concatenate(arts)); add(f'n_digit_{limb}'+S,'n',trimesh.util.concatenate(ners))
pickle.dump((dict(OUT),LAY),open('net_branch.pkl','wb'))
print(len(OUT),'ids', sum(len(m.faces) for v in OUT.values() for m in v),'faces')
