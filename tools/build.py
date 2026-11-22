import trimesh,numpy as np,json,pickle,fast_simplification as fs,os,sys
HD=os.environ.get('Q')=='hd'
def R(std,hd): return hd if HD else std
from collections import defaultdict
def simp(m,ratio,minf=150):
    f=len(m.faces); tgt=max(minf,int(f*ratio))
    if tgt>=f: return m
    v,fa=fs.simplify(m.vertices.astype(np.float32),m.faces.astype(np.int32),target_reduction=1-tgt/f)
    return trimesh.Trimesh(v,fa,process=False)
scene=trimesh.Scene(); meta={}
# ---- bones
parts=pickle.load(open('bones_parts.pkl','rb')); lab=json.load(open('data/bone_labels2.json'))
from skull_split import split as skull_split
groups=defaultdict(list)
for i,p in enumerate(parts): groups[lab[str(i)]].append(p)
nb=0
for id,ps in groups.items():
    m=trimesh.util.concatenate(ps)
    subs = {('cr_'+k):v for k,v in skull_split(m).items()} if id=='crane' else {id:m}
    for sid,sm in subs.items():
        sm=simp(sm,R(0.4,1.0),100); nb+=len(sm.faces)
        scene.add_geometry(sm,node_name='b:'+sid,geom_name='b:'+sid)
# ---- bony landmarks (markers)
for lid,pos in json.load(open('data/landmarks.json')).items():
    sph=trimesh.creation.icosphere(subdivisions=2,radius=0.0095); sph.apply_translation(pos)
    scene.add_geometry(sph,node_name='k:'+lid,geom_name='k:'+lid)
print('bones faces',nb,'groups',len(groups))
# ---- muscles
sc=trimesh.load('sources/feline_ecorche_dec.glb'); mm=sc.to_geometry(); mm.merge_vertices()
mparts=[p for p in mm.split(only_watertight=False) if len(p.faces)>=40]
MNAMES=json.load(open('data/muscle_names.json'))
def region(c):
    x,y,z=c; ax=abs(x)
    if z< -1.3 and ax<0.12: return 'queue'
    if z>1.45: return 'tete'
    if y<0.3 and z>0.55: return 'membre_thor'
    if y<0.45 and z< -0.45: return 'membre_pelv'
    if z< -0.55 and ax>0.12: return 'membre_pelv'
    if y>=0.3 and 0.7<z<1.2 and ax>0.13: return 'epaule'
    if z>1.05: return 'cou'
    return 'tronc'
nm=0; cnt=defaultdict(int)
for k,p in enumerate(mparts):
    mn=MNAMES[str(k)]; cnt[mn]+=1
    q=simp(p,R(0.28,0.9),80); nm+=len(q.faces)
    scene.add_geometry(q,node_name=f'm:{mn}:{k}',geom_name=f'm{k}')
print('muscle faces',nm,dict(cnt))
# ---- organs (Pixelbeaker), registered
sc=trimesh.load('/mnt/user-data/uploads/normal_feline_anatomy_-_digestive_system__1_.glb')
MAT={'default_08_Stomach_UpperGI_1':'estomac','default_10_Bladder_1':'vessie','default_01c_Blood_Volume_ab_dec':'vaisseaux_abdo',
     'default_05b_Liver_1':'foie','default_06b_Intestines_dec':'intestins','default_07b_Pancreas_dec':'pancreas','default_09c_Kidneys_1':'reins'}
RATIO={'foie':R(0.03,0.12),'intestins':R(0.08,0.3),'estomac':R(0.08,0.3),'vessie':R(0.06,0.25),'vaisseaux_abdo':R(0.12,0.45),'pancreas':R(0.12,0.4),'reins':R(0.06,0.25)}
s=0.0052; kx=0.85
def reg(v):
    X,Y,Z=v[:,0],v[:,1],v[:,2]
    zw=0.46-(X-165.5)*s; yw=0.49+(Y-216.7)*s; xw=(Z-5.9)*s*kx
    return np.c_[xw,yw,zw]
og=defaultdict(list)
for nn in sc.graph.nodes_geometry:
    T,gn=sc.graph[nn]; geo=sc.geometry[gn]; id=MAT.get(geo.visual.material.name)
    if id: og[id].append(geo.copy().apply_transform(T))
no=0
for id,ms in og.items():
    m=trimesh.util.concatenate(ms); m.merge_vertices()
    m=simp(m,RATIO[id],300); m.vertices=reg(np.asarray(m.vertices)); no+=len(m.faces)
    L='v' if id=='vaisseaux_abdo' else 'o'
    scene.add_geometry(m,node_name=L+':'+id,geom_name=L+':'+id)
    print(id,len(m.faces),np.round(m.bounds,2).tolist())
print('organ faces',no)
# ---- reconstructed structures
for pk in ['recon_thorax.pkl','recon_vn.pkl']:
    for id,(L,m) in pickle.load(open(pk,'rb')).items():
        if id in ('coeur','poumon_g','poumon_d'): continue
        r = R(0.12,0.35) if id.startswith(('poumon','lobe_')) else (R(0.5,1.0) if L in('v','n') and id!='encephale' else R(0.6,1.0))
        m=simp(m,r,200); scene.add_geometry(m,node_name=L+':'+id,geom_name=L+':'+id); no+=len(m.faces)
for pk,ratio in [('net_trunk.pkl',R(0.6,1.0)),('net_tree.pkl',R(0.55,1.0)),('net_branch.pkl',R(0.75,1.0)),('net_extra.pkl',R(0.8,1.0)),('net_more.pkl',R(0.8,1.0)),('net_female.pkl',R(0.8,1.0)),('net_skin.pkl',R(0.8,1.0)),('net_section.pkl',1.0)]:
    D,LAY=pickle.load(open(pk,'rb'))
    for id,ms in D.items():
        if pk=='net_branch.pkl' and id.startswith('n_optique'): continue
        m=trimesh.util.concatenate(ms); m=simp(m,R(0.06,0.2) if id=='peau' else ratio,150); L=LAY[id]
        scene.add_geometry(m,node_name=L+':'+id,geom_name=L+':'+id); no+=len(m.faces)
print('total with recon',no)
# ---- garde-fou : aucune structure molle ne doit traverser la peau finale
from scipy.spatial import cKDTree
SK=trimesh.load('cat_raw.glb').geometry['p:peau'] if HD else scene.geometry['p:peau']; sp,sf=trimesh.sample.sample_surface(SK,900000,seed=1); sn=SK.face_normals[sf]; ST=cKDTree(sp)
MARG=0.006 if HD else 0.004; moved=[]
for gn,g in scene.geometry.items():
    L=gn.split(':')[0]
    if L not in 'vnlot' or gn.startswith('o:mamelles'): continue
    v=np.asarray(g.vertices).copy(); tot=0
    for _ in range(4):
        d,i=ST.query(v); sd=((v-sp[i])*sn[i]).sum(1)
        cand=np.where((sd>-MARG)&(d<0.06))[0]
        if not len(cand): break
        out=~np.concatenate([SK.contains(v[cand[j:j+300]]) for j in range(0,len(cand),300)])   # vrai test intérieur/extérieur (rayons)
        near=(d[cand]<MARG*1.5)&(sd[cand]>-MARG)          # juste sous la surface : normale fiable
        sel=cand[out|near]
        if not len(sel): break
        v[sel]=sp[i[sel]]-sn[i[sel]]*MARG; tot+=len(sel)
    if tot: g.vertices=v; moved.append((gn,tot))
print('pushed inside skin:',moved)
scene.export('cat_raw_hd.glb' if HD else 'cat_raw.glb')
