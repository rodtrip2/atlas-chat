import pickle,json,numpy as np,trimesh,fast_simplification as fs
from collections import defaultdict
parts=pickle.load(open('bones_parts.pkl','rb')); lab=json.load(open('data/bone_labels.json'))
_G=defaultdict(list)
for i,p in enumerate(parts): _G[lab[str(i)]].append(p)
B={k:trimesh.util.concatenate(v) for k,v in _G.items()}
CAN={k:(np.array(v[0]) if v else None) for k,v in json.load(open('data/canal.json')).items()}
s=0.0052; kx=0.85
def reg(v):
    X,Y,Z=v[:,0],v[:,1],v[:,2]
    return np.c_[(Z-5.9)*s*kx, 0.49+(Y-216.7)*s, 0.46-(X-165.5)*s]
def organs():
    import os
    if os.path.exists('organs_reg.pkl'): return pickle.load(open('organs_reg.pkl','rb'))
    MAT={'default_08_Stomach_UpperGI_1':'estomac','default_10_Bladder_1':'vessie','default_01c_Blood_Volume_ab_dec':'vaisseaux_abdo',
     'default_05b_Liver_1':'foie','default_06b_Intestines_dec':'intestins','default_07b_Pancreas_dec':'pancreas','default_09c_Kidneys_1':'reins'}
    sc=trimesh.load('sources/normal_feline_anatomy_digestive_system.glb'); og=defaultdict(list)
    for nn in sc.graph.nodes_geometry:
        T,gn=sc.graph[nn]; geo=sc.geometry[gn]; id=MAT.get(geo.visual.material.name)
        if id: og[id].append(geo.copy().apply_transform(T))
    out={}
    for id,ms in og.items():
        m=trimesh.util.concatenate(ms); m.merge_vertices(); m.vertices=reg(np.asarray(m.vertices)); out[id]=m
    pickle.dump(out,open('organs_reg.pkl','wb')); return out
def axis(m):
    v=np.asarray(m.vertices); c=v.mean(0); u,sv,vt=np.linalg.svd(v-c,full_matrices=False); d=vt[0]
    t=(v-c)@d; return c+d*t.min(), c+d*t.max()

from scipy.interpolate import splprep, splev
def spline(pts, step=0.006):
    P=np.asarray(pts,float)
    if len(P)<3: P=np.vstack([P[0],(P[0]+P[1])/2,P[1]])
    k=min(3,len(P)-1); tck,u=splprep(P.T,s=0,k=k)
    L=np.sum(np.linalg.norm(np.diff(P,axis=0),axis=1)); n=max(8,int(L/step))
    return np.array(splev(np.linspace(0,1,n),tck)).T
def tube(pts, r0, r1=None, seg=None, step=None):
    r1=r0 if r1 is None else r1; rm=max(r0,r1)
    if seg is None: seg=16 if rm>=0.01 else 12 if rm>=0.005 else 8 if rm>=0.0025 else 6
    if step is None: step=0.005 if rm>=0.005 else 0.007
    C=spline(pts,step); n=len(C)
    R=np.linspace(r0,r1,n)
    T=np.gradient(C,axis=0); T/=np.linalg.norm(T,axis=1)[:,None]
    a=np.array([0,1,0]) if abs(T[0][1])<0.9 else np.array([1,0,0])
    N=np.cross(T[0],a); N/=np.linalg.norm(N); Ns=[N]
    for i in range(1,n):
        N=Ns[-1]-T[i]*np.dot(Ns[-1],T[i]); N/=np.linalg.norm(N); Ns.append(N)
    Ns=np.array(Ns); Bs=np.cross(T,Ns)
    ang=np.linspace(0,2*np.pi,seg,endpoint=False)
    V=(C[:,None,:]+R[:,None,None]*(np.cos(ang)[None,:,None]*Ns[:,None,:]+np.sin(ang)[None,:,None]*Bs[:,None,:])).reshape(-1,3)
    F=[]
    for i in range(n-1):
        for j in range(seg):
            a_=i*seg+j; b_=i*seg+(j+1)%seg; c_=(i+1)*seg+j; d_=(i+1)*seg+(j+1)%seg
            F+= [[a_,c_,b_],[b_,c_,d_]]
    V=np.vstack([V,C[0],C[-1]]); s0=len(V)-2; s1=len(V)-1
    for j in range(seg):
        F.append([s0,(j+1)%seg,j]); F.append([s1,(n-1)*seg+j,(n-1)*seg+(j+1)%seg])
    return trimesh.Trimesh(V,np.array(F),process=False)
def body_ymin(z):
    seq=[f'C{i}' for i in range(1,8)]+[f'T{i}' for i in range(1,14)]+[f'L{i}' for i in range(1,8)]
    zs=np.array([B[v].centroid[2] for v in seq]); ys=np.array([B[v].bounds[0][1] for v in seq])
    o=np.argsort(zs); return np.interp(z,zs[o],ys[o])
def canal_y(z):
    seq=[k for k,v in CAN.items() if v is not None]
    zs=np.array([CAN[k][2] for k in seq]); ys=np.array([CAN[k][1] for k in seq]); o=np.argsort(zs)
    return np.interp(z,zs[o],ys[o])
def voxmesh(mask, origin, h, smooth=8, sigma=1.0):
    from skimage.measure import marching_cubes
    from scipy.ndimage import gaussian_filter
    f=gaussian_filter(mask.astype(np.float32),sigma)
    v,fa,_,_=marching_cubes(np.pad(f,1),0.5)
    v=(v-1)*h+origin
    m=trimesh.Trimesh(v,fa[:,::-1],process=True)
    trimesh.smoothing.filter_taubin(m,iterations=smooth)
    return m

def ball(c,r,sub=1):
    m=trimesh.creation.icosphere(subdivisions=sub,radius=r); m.apply_translation(c); return m
def centerline(m, nb=24):
    """ordered centerline of an elongated curved bone (rib): PCA plane + angular binning around fitted circle"""
    v=np.asarray(m.vertices); c=v.mean(0); u,sv,vt=np.linalg.svd(v-c,full_matrices=False)
    e1,e2=vt[0],vt[1]; q=np.c_[(v-c)@e1,(v-c)@e2]
    A=np.c_[2*q,np.ones(len(q))]; bb=(q**2).sum(1); sol=np.linalg.lstsq(A,bb,rcond=None)[0]; cc=sol[:2]
    ang=np.unwrap(np.sort(np.arctan2(q[:,1]-cc[1],q[:,0]-cc[0])))
    a=np.arctan2(q[:,1]-cc[1],q[:,0]-cc[0])
    # choose cut at largest angular gap
    sa=np.sort(a); gaps=np.diff(np.r_[sa,sa[0]+2*np.pi]); cut=sa[(np.argmax(gaps)+1)%len(sa)]
    a=(a-cut)%(2*np.pi)
    edges=np.linspace(a.min(),a.max(),nb+1); pts=[]
    for i in range(nb):
        sel=(a>=edges[i])&(a<=edges[i+1])
        if sel.sum()>3: pts.append(v[sel].mean(0))
    return np.array(pts)
