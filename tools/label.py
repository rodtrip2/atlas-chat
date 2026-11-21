import pickle,numpy as np,json
parts=pickle.load(open('bones_parts.pkl','rb'))
C=np.array([p.centroid for p in parts]); E=np.array([p.extents for p in parts]); F=np.array([len(p.faces) for p in parts])
lab={}
def S(i,id): lab[i]=id
# midline
S(21,'crane'); S(67,'mandibule')
for k,i in enumerate([3,36,31,94,38,23,37]): S(i,f'C{k+1}')
for k,i in enumerate([40,83,33,6,2,14,88,12,44,10,17,81,78]): S(i,f'T{k+1}')
for k,i in enumerate([8,69,20,26,47,39,42]): S(i,f'L{k+1}')
S(16,'sacrum'); S(5,'coxal')
cau=[41,9,135,57,73,11,13,52,114,146,142,134,115,222,157,161,164]
for k,i in enumerate(cau): S(i,f'Cd{k+1}')
S(46,'manubrium'); [S(i,'sternebre') for i in [99,189,155,232,179,121]]; S(158,'xiphoide')
# right side explicit, then mirror
R={}
for k,i in enumerate([112,149,107,218,167,119,132,165,171,151,219,98,150]): R[i]=f'cote{k+1}'
for k,i in enumerate([214,145,168,185,56,110,138,7,66,180,90,61]): R[i]=f'cart{k+1}'
R[243]='cart13'
R.update({29:'scapula',74:'humerus',35:'ulna',76:'radius',25:'femur',101:'tibia',93:'fibula',131:'patella',58:'calcaneus',15:'talus'})
for i in [152,160,205,125]: R[i]='metacarpe'
for i in [108,96,143,183]: R[i]='phal1_m'
for i in [70,209,192,32]: R[i]='phal2_m'
for i in [200,136,163,120]: R[i]='phal3_m'
for i in [141,126,204,175]: R[i]='metatarse'
for i in [188,86,111,216]: R[i]='phal1_p'
for i in [190,201,215,230]: R[i]='phal2_p'
for i in [19,55,34,166]: R[i]='phal3_p'
def mirror(i):
    t=C[i]*[-1,1,1]; d=np.linalg.norm(C-t,axis=1); d[i]=9; j=int(np.argmin(d)); return j,d[j]
for i,id in R.items():
    S(i,id+'_g'); j,dd=mirror(i)
    if dd>0.02: print('mirror far',i,id,dd)
    S(j,id+'_d')
# region rules for the rest
for i in range(len(parts)):
    if i in lab: continue
    x,y,z=C[i]; side='_g' if x>0 else '_d'
    if F[i]<60: continue
    if z>1.7 and y<0.5: S(i,'dent'); continue
    if 0.75<z<1.0 and -0.9<y<-0.8 and abs(x)>0.05: S(i,'carpe'+side); continue
    if 0.85<z<1.0 and -0.96<y<-0.88: S(i,'metacarpe'+side); continue  # MC1 / dewclaw
    if -1.25<z<-1.1 and -0.72<y<-0.55: S(i,'tarse'+side); continue
    if z>1.3 and 0.2<y<0.6 and abs(x)<0.2: S(i,'hyoide'); continue
# tiny fragments -> nearest labeled
labeled=[i for i in lab]
for i in range(len(parts)):
    if i in lab: continue
    d=np.linalg.norm(C[labeled]-C[i],axis=1); j=labeled[int(np.argmin(d))]
    print('unlab',i,F[i],C[i].round(3),'->',lab[j], round(float(d.min()),3)) if F[i]>=60 else None
    lab[i]=lab[j]
json.dump({str(k):v for k,v in lab.items()},open('data/bone_labels.json','w'),indent=0)
from collections import Counter
print(Counter(v.rstrip('_dg') if v[-2:] in('_d','_g') else v for v in lab.values()))
