# Affinage : petits os nommés un par un (carpe, tarse, métacarpiens/tarsiens, phalanges par doigt, dents, sternèbres)
import pickle,json,numpy as np
from collections import defaultdict
parts=pickle.load(open('bones_parts.pkl','rb')); lab=json.load(open('data/bone_labels.json'))
C=np.array([p.centroid for p in parts]); E=np.array([p.extents for p in parts]); F=np.array([len(p.faces) for p in parts])
G=defaultdict(list)
for i in range(len(parts)): G[lab[str(i)]].append(i)
new=dict(lab)
def big(ids): return [i for i in ids if F[i]>=60]
def mirror(i):
    t=C[i]*[-1,1,1]; d=np.linalg.norm(C-t,axis=1); d[i]=9; return int(np.argmin(d))
def setboth(i,id_g):   # i is a LEFT (+x) piece; set left and its mirror
    new[str(i)]=id_g+'_g'; new[str(mirror(i))]=id_g+'_d'
# --- left side analysis (x>0 = left), medial = smaller x
def by_x(ids): return sorted(big(ids),key=lambda i:C[i][0])
# forelimb digits II..V
for row,key in [('phal1_m','phal1_m'),('phal2_m','phal2_m'),('phal3_m','phal3_m'),('phal1_p','phal1_p'),('phal2_p','phal2_p'),('phal3_p','phal3_p')]:
    for d,i in zip(['II','III','IV','V'],by_x(G[key+'_g'])): setboth(i,f'{row}{d}')
mc=by_x(G['metacarpe_g'])
small=[i for i in mc if E[i][1]<0.08]; long_=[i for i in mc if E[i][1]>=0.08]
for d,i in zip(['II','III','IV','V'],sorted(long_,key=lambda i:C[i][0])): setboth(i,f'mc{d}')
if small:
    s=sorted(small,key=lambda i:C[i][1],reverse=True)   # higher = metacarpal I, lower = proximal phalanx I
    setboth(s[0],'mcI')
    for i in s[1:]: setboth(i,'phal1_mI')
for d,i in zip(['II','III','IV','V'],by_x(G['metatarse_g'])): setboth(i,f'mt{d}')
# carpus: accessory = most caudal (min z) ; proximal row = higher y
ca=big(G['carpe_g']); acc=min(ca,key=lambda i:C[i][2]); setboth(acc,'carpe_acc'); ca.remove(acc)
radial=max(ca,key=lambda i:F[i]); setboth(radial,'carpe_rad'); ca.remove(radial)
prox=[i for i in ca if C[i][1]>-0.84 and C[i][0]>0.2]   # lateral proximal = ulnar carpal
uln=max(ca,key=lambda i:C[i][0]+ (C[i][1]+0.84)*2); setboth(uln,'carpe_uln'); ca.remove(uln)
for d,i in zip(['I','II','III','IV'],sorted(ca,key=lambda i:C[i][0])): setboth(i,f'carpe{d}')
# tarsus
ta=big(G['tarse_g']); t4=max(ta,key=lambda i:F[i]*(C[i][0]>0.22)); setboth(t4,'tarseIV'); ta.remove(t4)
cen=max(ta,key=lambda i:C[i][1]); setboth(cen,'tarse_cen'); ta.remove(cen)
for d,i in zip(['II','III'],sorted(ta,key=lambda i:C[i][0])): setboth(i,f'tarse{d}')
# sternebrae
for k,i in enumerate(sorted(G['sternebre'],key=lambda i:-C[i][2])): new[str(i)]=f'sternebre{k+2}'
new_m=[i for i in range(len(parts)) if lab[str(i)]=='manubrium']; new_x=[i for i in range(len(parts)) if lab[str(i)]=='xiphoide']
# teeth: jaw of origin = nearest bone (cranium vs mandible)
import trimesh
cr=trimesh.util.concatenate([parts[i] for i in G['crane']]); md=trimesh.util.concatenate([parts[i] for i in G['mandibule']])
from scipy.spatial import cKDTree
tc=cKDTree(np.asarray(cr.vertices)); tm=cKDTree(np.asarray(md.vertices))
teeth=[i for i in G['dent']]
quad=defaultdict(list)
for i in teeth:
    v=np.asarray(parts[i].vertices); dc=tc.query(v)[0].min(); dm=tm.query(v)[0].min()
    jaw='sup' if dc<dm else 'inf'; side='g' if C[i][0]>0.001 else ('d' if C[i][0]<-0.001 else ('g' if C[i][0]>=0 else 'd'))
    quad[(jaw,side)].append(i)
for (jaw,side),ids in quad.items():
    ids=sorted(ids,key=lambda i:-C[i][2])      # rostral -> caudal
    inc=[i for i in ids if abs(C[i][0])<0.04 and E[i].max()<0.03]
    rest=[i for i in ids if i not in inc]
    for k,i in enumerate(sorted(inc,key=lambda i:abs(C[i][0]))): new[str(i)]=f'dent_I{k+1}_{jaw}_{side}'
    can=max(rest,key=lambda i:E[i][2]*(C[i][2]>1.85)) if rest else None
    if can is not None: new[str(can)]=f'dent_C_{jaw}_{side}'; rest.remove(can)
    cheek=sorted(rest,key=lambda i:C[i][2])     # caudal -> rostral
    names=['M1','P4','P3','P2'] if jaw=='sup' else ['M1','P4','P3']
    for k,i in enumerate(cheek): new[str(i)]=f'dent_{names[k] if k<len(names) else "P"}_{jaw}_{side}'
    print(jaw,side,'incisives',len(inc),'canine',can is not None,'jugales',len(cheek))
json.dump(new,open('data/bone_labels2.json','w'),indent=0)
from collections import Counter
print(len(set(new.values())),'structures')
print(sorted(v for v in set(new.values()) if v.startswith(('dent','carpe','tarse','mc','mt','phal1_mI'))))
