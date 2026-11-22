# Identification automatique proposée des muscles de l'écorché, d'après leurs attaches osseuses et leur position
import json,numpy as np
rows=json.load(open('data/muscle_rows.json')); LM=json.load(open('data/landmarks.json'))
def name(r):
    a=sorted([r['a'][0][0],r['a'][1][0]]); A=set(a); c=np.array(r['c']); x,y,z=c; s=1 if x>=0 else -1; S='_g' if s>0 else '_d'; ax=abs(x)
    P=lambda *k: set(k)==A or (len(k)==1 and A=={k[0]})
    # ---- head
    if A=={'crane','mandibule'}:
        if y>0.66: return 'temporal'
        if y<0.42 and z<1.7: return 'digastrique'
        return 'masseter'
    if A=={'mandibule'}: return 'masseter'
    if A=={'crane'}:
        if y>0.72 and z<1.78: return 'temporal'
        return 'faciaux'
    if A=={'crane','sternum'}: return 'sternocephalique'
    if A=={'crane','humerus'}: return 'brachiocephalique'
    if A=={'cervicales','crane'}: return 'splenius' if y>0.62 else 'longus_capitis'
    # ---- neck / shoulder girdle
    if A in ({'cervicales','scapula'},{'scapula','thoraciques'}): return 'rhomboide' if ax<0.1 else 'trapeze'
    if A in ({'cervicales'},{'cervicales','thoraciques'}): return 'cou_profonds'
    if A=={'humerus','scapula'} or A=={'scapula'}:
        sp=np.array(LM['scap_epine'+S])
        if ax < abs(sp[0])-0.03: return 'subscapulaire'
        if y<0.45: return 'deltoide'
        return 'supra_epineux' if z>sp[2] else 'infra_epineux'
    if A in ({'humerus','sternum'},{'sternum'},{'côtes','humerus'}): return 'pectoraux'
    if A=={'humerus','thoraciques'}: return 'grand_dorsal'
    # ---- forelimb
    if A in ({'humerus','ulna'},{'scapula','ulna'}): return 'triceps'
    if A=={'humerus'}: return 'triceps' if z<LM['hum_tete'+S][2]-0.02 else 'brachial'
    if A=={'humerus','radius'}: return 'biceps_brachial' if z>0.80 else 'brachial'
    if A & {'radius','ulna'} or A=={'humerus','main'} or A=={'main','ulna'}:
        return 'avant_bras_ext' if z>0.83 else 'avant_bras_flech'
    if A=={'main'}: return 'main_intrinseques'
    # ---- trunk
    if A in ({'thoraciques'},{'lombaires'},{'lombaires','thoraciques'},{'lombaires','sacrum'},{'sacrum','thoraciques'}):
        return 'epaxiaux'
    if A=={'côtes','thoraciques'}: return 'epaxiaux' if y>0.6 else 'dentele'
    if A=={'côtes'}: return 'dentele' if z>0.45 else 'obliques'
    if A=={'côtes','sternum'}: return 'pectoraux' if z>0.7 else 'droit_abdomen'
    if A in ({'coxal','côtes'},{'coxal','sternum'},{'lombaires','sternum'},{'coxal','lombaires'}): return 'obliques' if y>0.15 else 'droit_abdomen'
    if A in ({'caudales'},{'caudales','sacrum'},{'caudales','lombaires'}): return 'queue'
    # ---- hindlimb
    if A in ({'coxal','femur'},{'femur','sacrum'},{'coxal','sacrum'},{'coxal'}):
        return 'fessiers' if y>0.45 else 'adducteurs'
    if A in ({'coxal','tibia'},{'coxal','fibula'},{'caudales','femur'}):
        if ax>0.2: return 'biceps_femoral'
        if z>-0.85: return 'sartorius'
        return 'semitendineux'
    if A in ({'femur'},{'patella'}): return 'quadriceps'
    if A in ({'calcaneus','femur'},{'femur','pied'},{'calcaneus','fibula'}): return 'gastrocnemien'
    if A in ({'pied','tibia'},{'fibula','pied'},{'fibula'},{'calcaneus','pied'}):
        return 'jambe_cran' if z>-1.0 and y<-0.15 and r['c'][2]>-1.05 else 'jambe_caud'
    if A=={'pied'}: return 'pied_intrinseques'
    return None
out={}; from collections import Counter
for r in rows:
    n=name(r); side='' if abs(r['c'][0])<0.01 else ('_g' if r['c'][0]>0 else '_d')
    if n is None: n='non_identifie'
    out[r['k']]=n+side
json.dump(out,open('data/muscle_names.json','w'))
cnt=Counter(v.replace('_g','').replace('_d','') for v in out.values()); print(cnt.most_common())
