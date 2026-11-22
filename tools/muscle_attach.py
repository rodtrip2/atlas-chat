import pickle,json,numpy as np,trimesh
from scipy.spatial import cKDTree
from collections import Counter
parts=pickle.load(open('bones_parts.pkl','rb')); lab=json.load(open('data/bone_labels2.json'))
def coarse(l):
    import re
    l2=re.sub(r'_[dg]$','',l)
    if re.match(r'^C\d',l2): return 'cervicales'
    if re.match(r'^T\d',l2): return 'thoraciques'
    if re.match(r'^L\d',l2): return 'lombaires'
    if l2.startswith('Cd'): return 'caudales'
    if l2.startswith(('cote','cart')): return 'côtes'
    if l2.startswith(('sternebre','manubrium','xiphoide')): return 'sternum'
    if l2.startswith(('carpe','mc','phal1_m','phal2_m','phal3_m')): return 'main'
    if l2.startswith(('tarse','mt','phal1_p','phal2_p','phal3_p','talus')): return 'pied'
    if l2.startswith('dent'): return 'crane'
    return l2
pts=[];labs=[]
for i,p in enumerate(parts):
    l=coarse(lab[str(i)]); n=max(30,int(len(p.faces)*0.3)); s=p.sample(n); pts.append(s); labs+=[l]*n
pts=np.vstack(pts); labs=np.array(labs,dtype=object); tr=cKDTree(pts)
sc=trimesh.load('sources/feline_ecorche_dec.glb'); mm=sc.to_geometry(); mm.merge_vertices()
mp=[p for p in mm.split(only_watertight=False) if len(p.faces)>=40]
rows=[]
for k,p in enumerate(mp):
    v=np.asarray(p.vertices); c=v.mean(0); d=np.linalg.svd(v-c,full_matrices=False)[2][0]; t=(v-c)@d
    ends=[v[t<np.percentile(t,6)].mean(0), v[t>np.percentile(t,94)].mean(0)]
    att=[]
    for e in ends:
        dd,ii=tr.query(e,k=40); cnt=Counter(labs[ii]); att.append((cnt.most_common(1)[0][0],round(float(dd.min()),3)))
    side='g' if c[0]>0.01 else ('d' if c[0]<-0.01 else 'm')
    rows.append(dict(k=k,f=len(p.faces),c=c.round(3).tolist(),L=round(float(np.ptp(t)),3),a=att,side=side,area=round(float(p.area),4)))
json.dump(rows,open('data/muscle_rows.json','w'))
left=[r for r in rows if r['side']!='d']
print(len(rows),'pieces;',len(left),'left/mid')
pairs=Counter(tuple(sorted([r['a'][0][0],r['a'][1][0]])) for r in left)
for pr,n in pairs.most_common(70): print(n,pr)
