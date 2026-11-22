import pickle,json,numpy as np,trimesh,fast_simplification as fs,os,sys
from skull_landmarks import skull_lms
HD=os.environ.get('Q')=='hd'
P=pickle.load(open('skull_detail_parts.pkl','rb')); CAN=json.load(open('data/canal.json'))
def simp(m,r,minf=60):
    f=len(m.faces); tgt=max(minf,int(f*r))
    if tgt>=f: return m
    v,fa=fs.simplify(m.vertices.astype(np.float32),m.faces.astype(np.int32),target_reduction=1-tgt/f); return trimesh.Trimesh(v,fa,process=False)
sc=trimesh.Scene(); tot=0
for k,m in P.items():
    r=(0.75 if HD else 0.22) if not k.startswith('dent_') else (1.0 if HD else 0.5)
    q=simp(m,r); tot+=len(q.faces); sc.add_geometry(q,node_name='x:'+k,geom_name='x:'+k)
crv=np.vstack([np.asarray(m.vertices) for k,m in P.items() if k.startswith('cr_')]); mdv=np.asarray(P['mandibule'].vertices)
LM=skull_lms(crv,mdv,CAN['C1'][0])
json.dump(LM,open('data/skull_detail_landmarks.json','w'))
for k,p in LM.items():
    s=trimesh.creation.icosphere(subdivisions=2,radius=0.006); s.apply_translation(p); sc.add_geometry(s,node_name='y:'+k,geom_name='y:'+k)
out='skull_hd_raw.glb' if HD else 'skull_raw.glb'; sc.export(out); print(out,'faces',tot,'landmarks',len(LM))
