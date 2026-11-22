# Reporte sur le crâne de l'atlas les repères crâniens mesurés sur le scanner recalé
import json,pickle,numpy as np
from scipy.spatial import cKDTree
LM=json.load(open('data/landmarks.json')); SL=json.load(open('data/skull_detail_landmarks.json'))
parts=pickle.load(open('bones_parts.pkl','rb')); lab=json.load(open('data/bone_labels2.json'))
cr=np.vstack([parts[int(i)].vertices for i,v in lab.items() if v=='crane']); md=np.vstack([parts[int(i)].vertices for i,v in lab.items() if v=='mandibule'])
tc,tm=cKDTree(cr),cKDTree(md); moved=[]
for k,p in SL.items():
    p=np.array(p)
    if k.startswith('cr_orbite'): new=p
    elif k.startswith('md_'): new=md[tm.query(p)[1]]
    else: new=cr[tc.query(p)[1]]
    if k in LM: moved.append((k,round(float(np.linalg.norm(np.array(LM[k])-new)),3)))
    LM[k]=[float(x) for x in new]
json.dump(LM,open('data/landmarks.json','w'),indent=0)
print(sorted(moved,key=lambda x:-x[1])[:10])
