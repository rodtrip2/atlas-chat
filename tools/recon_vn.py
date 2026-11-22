from common import *
import pickle
out={}
def ax_(k): a,b=axis(B[k]); return a,b
def P(a,b,t,off): return a+(b-a)*t+off
def unit(v): v=np.asarray(v,float); return v/np.linalg.norm(v)
def limbs(sign):
    s=lambda p: np.array([p[0]*sign,p[1],p[2]])
    S='_g' if sign>0 else '_d'
    # axes (computed on left bones then mirrored by using side-specific bones)
    def A(name):
        a,b=axis(B[name+S]); return a,b
    hp,hd=A('humerus'); 
    if hp[1]<hd[1]: hp,hd=hd,hp          # proximal = higher
    rp,rd=A('radius'); 
    if rp[1]<rd[1]: rp,rd=rd,rp
    up,ud=A('ulna');
    if up[1]<ud[1]: up,ud=ud,up
    fp,fd=A('femur');
    if fp[1]<fd[1]: fp,fd=fd,fp
    tp,td=A('tibia');
    if tp[1]<td[1]: tp,td=td,tp
    mc=B['metacarpe'+S].centroid; mt=B['metatarse'+S].centroid
    med=np.array([-0.045*sign,0,0]); lat=-med
    def perps(p,d):
        dd=unit(d-p); n1=np.array([0,-dd[2],dd[1]]); n1=unit(n1)
        if n1[2]<0: n1=-n1     # n1 points cranially
        return n1,-n1
    hc,hk=perps(hp,hd); rc,rk=perps(rp,rd); fc,fk=perps(fp,fd); tc,tk=perps(tp,td)
    if tk[1]<0: tk=-tk; tc=-tc
    R={}
    # ---- forelimb artery
    R['art_mt'+S]=('v',tube([s([0.03,0.47,0.97]),s([0.07,0.38,1.04]),P(hp,hd,0.08,med*1.3+hk*0.02),P(hp,hd,0.4,med+hc*0.01),P(hp,hd,0.8,med),
        P(rp,rd,0.1,med*0.9+rk*0.015),P(rp,rd,0.5,med*0.8+rk*0.02),P(rp,rd,0.9,med*0.6+rk*0.02),mc+np.array([0,-0.02,-0.01])],0.0065,0.003))
    # ---- forelimb nerves
    R['radial'+S]=('n',tube([s([0.10,0.45,1.06]),s([0.15,0.36,1.08]),P(hp,hd,0.25,med*0.5+hk*0.045),P(hp,hd,0.55,lat*0.6+hk*0.04),P(hp,hd,0.88,lat),
        P(rp,rd,0.2,lat*0.7+rc*0.02),P(rp,rd,0.6,lat*0.6+rc*0.02),P(rp,rd,0.95,lat*0.3+rc*0.02),mc+np.array([0,0.015,0.01])],0.0045,0.0025))
    R['median'+S]=('n',tube([s([0.15,0.36,1.08]),P(hp,hd,0.35,med*1.1+hc*0.008),P(hp,hd,0.8,med*1.1+hc*0.005),
        P(rp,rd,0.15,med*0.95+rk*0.008),P(rp,rd,0.6,med*0.85+rk*0.012),P(rp,rd,0.95,med*0.6+rk*0.012),mc+np.array([0,-0.025,0])],0.0035,0.002))
    R['ulnaire'+S]=('n',tube([s([0.15,0.36,1.08]),P(hp,hd,0.4,med*0.8+hk*0.035),P(hp,hd,0.85,med*0.7+hk*0.045),
        P(up,ud,0.1,med*0.5+rk*0.03),P(up,ud,0.5,rk*0.03),P(up,ud,0.9,lat*0.2+rk*0.02),mc+np.array([0.02*sign,-0.02,-0.02])],0.0032,0.002))
    # brachial plexus roots C6..T1
    zr=[(CAN['C6'][2]+CAN['C7'][2])/2,(CAN['C7'][2]+CAN['T1'][2])/2,(CAN['T1'][2]+CAN['T2'][2])/2,CAN['C6'][2]+0.02]
    roots=[tube([[0.012*sign,canal_y(z),z],[0.04*sign,canal_y(z)-0.02,z-0.01],s([0.10,0.45,1.06])],0.004,0.004) for z in zr]
    R['plexus'+S]=('n',trimesh.util.concatenate(roots+[tube([s([0.10,0.45,1.06]),s([0.15,0.36,1.08])],0.006)]))
    # ---- hindlimb artery
    R['art_mp'+S]=('v',tube([s([0.05,0.33,-0.83]),s([0.12,0.36,-0.88]),P(fp,fd,0.3,med+fc*0.01),P(fp,fd,0.6,med*1.0),P(fp,fd,0.92,med*0.3+fk*0.04),
        P(tp,td,0.15,tk*0.035),P(tp,td,0.5,med*0.3+tk*0.03),P(tp,td,0.9,med*0.3+tk*0.02),mt+np.array([0,0.01,-0.02])],0.0062,0.003))
    # ---- sciatic, tibial, fibular, femoral nerves
    zs7=[(CAN['L6'][2]+CAN['L7'][2])/2,(CAN['L7'][2]-0.10),-0.82]
    sroots=[tube([[0.012*sign,canal_y(z),z],[0.045*sign,canal_y(z)-0.02,z-0.02],s([0.09,0.66,-0.88])],0.004) for z in zs7]
    R['sciatique'+S]=('n',trimesh.util.concatenate(sroots+[tube([s([0.09,0.66,-0.88]),s([0.15,0.62,-1.04]),fp+np.array([0.0*sign,0.02,-0.08]),
        P(fp,fd,0.3,fk*0.06),P(fp,fd,0.6,fk*0.055),P(fp,fd,0.85,fk*0.05),P(tp,td,0.1,tk*0.045),P(tp,td,0.5,tk*0.04),P(tp,td,0.9,tk*0.025+med*0.2),mt+np.array([0,-0.01,-0.02])],0.0065,0.003)]))
    R['fibulaire'+S]=('n',tube([P(fp,fd,0.85,fk*0.05),P(fp,fd,1.0,lat*0.9+fk*0.03),P(tp,td,0.15,lat*1.1),P(tp,td,0.55,lat*0.9+tc*0.015),P(tp,td,0.95,lat*0.5+tc*0.02),mt+np.array([0,0.02,0.01])],0.0035,0.002))
    zf=[(CAN['L4'][2]+CAN['L5'][2])/2,(CAN['L5'][2]+CAN['L6'][2])/2]
    froots=[tube([[0.012*sign,canal_y(z),z],[0.05*sign,canal_y(z)-0.04,z-0.03],s([0.08,0.56,-0.66])],0.0035) for z in zf]
    R['femoral'+S]=('n',trimesh.util.concatenate(froots+[tube([s([0.08,0.56,-0.66]),s([0.12,0.42,-0.84]),P(fp,fd,0.35,med*0.7+fc*0.04),P(fp,fd,0.6,fc*0.05)],0.0042,0.0025)]))
    # ---- superficial veins used for venous access
    R['cephalique'+S]=('v',tube([mc+np.array([0,0.04,0.03]),P(rp,rd,0.85,rc*0.03+lat*0.2),P(rp,rd,0.45,rc*0.035+lat*0.2),P(rp,rd,0.05,rc*0.035),P(hp,hd,0.75,hc*0.04+lat*0.4),P(hp,hd,0.35,hc*0.045+lat*0.5),s([0.10,0.42,1.13])],0.004,0.005))
    R['saphene'+S]=('v',tube([mt+np.array([-0.0*sign,0.03,0.0]),P(tp,td,0.85,med*0.9),P(tp,td,0.5,med*1.0),P(tp,td,0.1,med*0.9),P(fp,fd,0.8,med*1.2),P(fp,fd,0.45,med*1.25),s([0.06,0.33,-0.86])],0.0038,0.0045))
    # carotid, jugular, vagus, phrenic
    R['carotide'+S]=('v',tube([[-0.005,0.47,1.03],s([0.035,0.47,1.08]),s([0.04,0.48,1.30]),s([0.045,0.53,1.47])],0.006,0.005))
    R['jugulaire'+S]=('v',tube([[-0.02,0.44,1.03],s([0.05,0.42,1.07]),s([0.078,0.42,1.16]),s([0.088,0.44,1.32]),s([0.08,0.50,1.47])],0.0085,0.007))
    R['vague'+S]=('n',tube([s([0.04,0.58,1.53]),s([0.045,0.50,1.30]),s([0.038,0.475,1.08]),s([0.03,0.47,0.92]),s([0.022,0.48,0.72]),s([0.018,0.48,0.56]),[0.024+0.012*sign,0.51,0.44]],0.0028))
    R['phrenique'+S]=('n',tube([s([0.035,0.56,1.20]),s([0.06,0.43,1.05]),s([0.12,0.33,0.82]),s([0.125,0.30,0.62]),s([0.11,0.30,0.42])],0.0022))
    return R
out.update(limbs(1)); out.update(limbs(-1))
# ---- right subclavian origin piece and brachiocephalic trunk
out['tronc_bc']=('v',tube([[0.006,0.495,0.865],[0.0,0.48,0.96],[-0.005,0.47,1.03]],0.009))
out['art_mt_d']=('v',trimesh.util.concatenate([out['art_mt_d'][1],tube([[-0.005,0.47,1.03],[-0.03,0.47,0.99],[-0.03,0.47,0.97]],0.0065)]))
out['art_mt_g']=('v',trimesh.util.concatenate([out['art_mt_g'][1],tube([[0.018,0.50,0.84],[0.025,0.49,0.92],[0.03,0.47,0.97]],0.0065)]))
# ---- aorta
desc=[[0.02,float(body_ymin(z))-0.028,z] for z in np.linspace(0.70,0.42,6)]
out['aorte']=('v',tube([[0.0,0.40,0.80],[0.0,0.46,0.86],[0.012,0.51,0.84],[0.022,0.53,0.76]]+desc+[[0.004,0.575,0.36]],0.017,0.015))
out['tronc_pulm']=('v',trimesh.util.concatenate([tube([[0.02,0.33,0.86],[0.015,0.40,0.86],[0.01,0.43,0.82]],0.013),
    tube([[0.01,0.43,0.82],[0.06,0.425,0.78],[0.10,0.41,0.72]],0.009,0.006),tube([[0.01,0.43,0.82],[-0.06,0.425,0.78],[-0.10,0.41,0.72]],0.009,0.006)]))
out['vcr']=('v',tube([[-0.035,0.33,0.80],[-0.03,0.42,0.92],[-0.02,0.44,1.03]],0.014))
out['vcc_thor']=('v',tube([[-0.02,0.42,0.40],[-0.03,0.36,0.55],[-0.035,0.32,0.68]],0.016))
# ---- spinal cord
cp=[[0,0.722,1.545]]+[list(CAN[k]) for k in [f'C{i}' for i in range(1,8)]+[f'T{i}' for i in range(1,14)]+[f'L{i}' for i in range(1,8)]]+[[0,0.76,-0.75],[0,0.77,-0.86]]
out['moelle']=('n',tube(cp,0.012,0.005,step=0.008))
pickle.dump(out,open('recon_vn.pkl','wb'))
for k,(l,m) in sorted(out.items()): print(k,l,len(m.faces))
