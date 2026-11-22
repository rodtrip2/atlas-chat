# Recalage du crâne scanné (Édimbourg) sur le crâne + mandibule + dents du squelette westerly
import pickle,json,numpy as np,trimesh,fast_simplification as fs
from scipy.spatial import cKDTree
parts=pickle.load(open('bones_parts.pkl','rb')); lab=json.load(open('data/bone_labels2.json'))
tgt=trimesh.util.concatenate([parts[int(i)] for i,v in lab.items() if v in('crane','mandibule') or v.startswith('dent_')])
T=tgt.sample(40000); tt=cKDTree(T)
_m=trimesh.load('sources/domestic_cat_felis_catus.glb').to_geometry(); _m.merge_vertices()
main=max(_m.split(only_watertight=False),key=lambda p:len(p.faces)); pickle.dump(main,open('edi_full.pkl','wb'))
Sv=main.sample(40000)
def umeyama(A,B):
    ma,mb=A.mean(0),B.mean(0); A0,B0=A-ma,B-mb
    U,S,Vt=np.linalg.svd(B0.T@A0/len(A)); D=np.eye(3); D[2,2]=np.sign(np.linalg.det(U@Vt))
    R=U@D@Vt; s=np.trace(np.diag(S)@D)/A0.var(0).sum(); t=mb-s*R@ma; return s,R,t
best=None
for flipx in [1,-1]:
    F=np.diag([flipx,1,1])
    # initial similarity: centroids + RMS scale, keep axes
    A=Sv@F.T
    s0=np.sqrt(((T-T.mean(0))**2).sum(1).mean()/((A-A.mean(0))**2).sum(1).mean())
    s,R,t=s0,np.eye(3),T.mean(0)-s0*A.mean(0)
    for it in range(60):
        X=(s*(R@A.T)).T+t; d,idx=tt.query(X)
        keep=d<np.percentile(d,80)
        s,R,t=umeyama(A[keep],T[idx[keep]])
    X=(s*(R@A.T)).T+t; d,_=tt.query(X); err=np.median(d)
    print('flipx',flipx,'scale',round(s,4),'median dist',round(err,4),'p80',round(np.percentile(d,80),4))
    if best is None or err<best[0]: best=(err,s,R,t,F)
err,s,R,t,F=best
M=np.eye(4); M[:3,:3]=s*R@F; M[:3,3]=t
np.save('data/edi_to_atlas.npy',M); print('saved transform; median residual',round(err,4))
