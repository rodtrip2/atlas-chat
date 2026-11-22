# Module « Crâne détaillé » : scanner d'Édimbourg recalé, os identifiés par transfert depuis l'atlas
import pickle,json,numpy as np,trimesh,fast_simplification as fs,os
from scipy.spatial import cKDTree
from skull_split import classify
parts=pickle.load(open('bones_parts.pkl','rb')); lab=json.load(open('data/bone_labels2.json'))
M=np.load('data/edi_to_atlas.npy')
main=pickle.load(open('edi_full.pkl','rb')).copy(); main.apply_transform(M)
# reference labels: cranium / mandible / each tooth / C1 / C2
ref_pts=[];ref_lab=[]
for i,v in lab.items():
    if v in('crane','mandibule','C1','C2') or v.startswith('dent_'):
        p=parts[int(i)]; n=max(200,int(len(p.faces)*0.6)); pts=p.sample(n); ref_pts.append(pts); ref_lab+= [v]*len(pts)
ref_pts=np.vstack(ref_pts); ref_lab=np.array(ref_lab,dtype=object); tr=cKDTree(ref_pts)
C=main.triangles_center; d,idx=tr.query(C); L=ref_lab[idx]
# cranium faces -> bone regions by the same rules as the atlas skull
cr=L=='crane'; sub=classify(C[cr]); L=L.astype(object); L[cr]=np.array(['cr_'+s for s in sub],dtype=object)
# pieces far from any atlas bone (e.g. extra fragments) go to nearest label anyway; report
print('faces',len(C),'labels',len(set(L)))
out={}
for k in sorted(set(L)):
    fi=np.nonzero(L==k)[0]
    if len(fi)<50: continue
    out[k]=main.submesh([fi],append=True,repair=False)
pickle.dump(out,open('skull_detail_parts.pkl','wb'))
for k,m in out.items(): print(k,len(m.faces))
