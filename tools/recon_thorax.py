from common import *
from shapely.geometry import MultiPoint, Point
from shapely import vectorized
import pickle
O=organs(); out={}
h=0.007
xs=np.arange(-0.34,0.34,h); ys=np.arange(-0.05,0.70,h); zs=np.arange(0.05,1.14,h)
X,Y,Z=np.meshgrid(xs,ys,zs,indexing='ij')
# --- thoracic wall points
wall=[]
for k,m in B.items():
    if k.startswith(('cote','cart')) or k in('manubrium','sternebre','xiphoide'): wall.append(np.asarray(m.vertices))
W=np.vstack(wall)
cav=np.zeros(X.shape,bool)
for kz,z in enumerate(zs):
    sl=W[np.abs(W[:,2]-z)<0.05]
    if len(sl)<30: continue
    hull=MultiPoint(sl[:,:2]).convex_hull.buffer(-0.016)
    if hull.is_empty: continue
    cav[:,:,kz]=__import__('shapely').contains_xy(hull,X[:,:,kz],Y[:,:,kz])
# --- heart (parametric egg)
base=np.array([0.0,0.36,0.80]); apex=np.array([0.035,0.105,0.60]); ax=apex-base; Lh=np.linalg.norm(ax); ax/=Lh
P=np.stack([X,Y,Z],-1)-base; t=(P@ax)/Lh; rad=np.linalg.norm(P-(P@ax)[...,None]*ax,axis=-1)
R=0.112
prof=np.where(t<0.32, R*np.sqrt(np.clip(1-((0.32-t)/0.36)**2,0,1)), R*np.sqrt(np.clip(1-((t-0.32)/0.68)**2,0,1))*(1-0.25*np.clip(t-0.32,0,1)))
heart=(rad<prof)&(t>-0.05)&(t<1.0)
out['coeur']=('o',voxmesh(heart,np.array([xs[0],ys[0],zs[0]]),h,smooth=10,sigma=1.2))
_u=np.array([1.0,0,0])-ax*ax[0]; _u/=np.linalg.norm(_u); _w=np.cross(ax,_u); _w = _w if _w[2]>0 else -_w
radial=P-(P@ax)[...,None]*ax
cr=(_w-_u); cr/=np.linalg.norm(cr)          # cranial-right direction
score=radial@cr
atr=t<0.30
ch={'atrium_d':heart&atr&(score>0),'atrium_g':heart&atr&(score<=0),
    'ventricule_d':heart&~atr&(score>0.012)&(t<0.86),'ventricule_g':heart&~atr&((score<=0.012)|(t>=0.86))}
from scipy.ndimage import binary_erosion
for k_,m_ in ch.items():
    others=np.zeros_like(heart)
    for k2,m2 in ch.items():
        if k2!=k_: others|=m2
    from scipy.ndimage import binary_dilation as _bd
    m_=m_&~_bd(others,iterations=1)
    out[k_]=('o',voxmesh(m_,np.array([xs[0],ys[0],zs[0]]),h,smooth=8,sigma=0.9))
# --- diaphragm surface from liver+stomach cranial surface
V=np.vstack([np.asarray(O['foie'].vertices),np.asarray(O['estomac'].vertices)])
gx=np.round((V[:,0]-xs[0])/h).astype(int); gy=np.round((V[:,1]-ys[0])/h).astype(int)
hm=np.full((len(xs),len(ys)),-9.0)
np.maximum.at(hm,(gx.clip(0,len(xs)-1),gy.clip(0,len(ys)-1)),V[:,2])
ok=hm>-9; px,py=np.nonzero(ok); xv=xs[px]; yv=ys[py]; zv=hm[ok]
A=np.c_[np.ones_like(xv),xv,yv,xv**2,yv**2,xv*yv]; coef=np.linalg.lstsq(A,zv,rcond=None)[0]
XX,YY=np.meshgrid(xs,ys,indexing='ij'); zfit=np.c_[np.ones(XX.size),XX.ravel(),YY.ravel(),XX.ravel()**2,YY.ravel()**2,(XX*YY).ravel()]@coef
zfit=zfit.reshape(XX.shape)
from scipy.ndimage import grey_dilation, gaussian_filter
hmd=np.where(ok,hm,-9); hmd=grey_dilation(hmd,size=3)
zd=np.maximum(zfit,np.where(hmd>-9,hmd,-9))+0.012
zd=zd-np.clip(YY-0.45,0,None)*0.7
zd=np.clip(gaussian_filter(zd,1.5),0.16,0.70)
print('diaphragm z: ventral mid',zd[np.argmin(abs(xs)),np.argmin(abs(ys-0.15))].round(3),'dorsal mid',zd[np.argmin(abs(xs)),np.argmin(abs(ys-0.5))].round(3),'dorsal lat',zd[np.argmin(abs(xs-0.2)),np.argmin(abs(ys-0.5))].round(3))
below=Z<zd[:,:,None]
# --- exclusions
bym=body_ymin(Z)
vert=(np.abs(X)<0.05)&(Y>bym-0.012)
med=np.abs(X)<0.022
dmed=(np.abs(X)<0.04)&(Y>0.40)
from scipy.ndimage import binary_dilation
hd=binary_dilation(heart,iterations=2)
lung=cav&~below&~vert&~med&~dmed&~hd&(Z<1.10)
from scipy.ndimage import label
for side,cond,nm in [('g',X>0,'poumon_g'),('d',X<0,'poumon_d')]:
    m=lung&cond; lb,n=label(m); 
    if n>1:
        sizes=np.bincount(lb.ravel()); sizes[0]=0; m=lb==sizes.argmax()
    out[nm]=('o',voxmesh(m,np.array([xs[0],ys[0],zs[0]]),h,smooth=12,sigma=1.3))
    print(nm,len(out[nm][1].faces),np.round(out[nm][1].bounds,3).tolist())
    zf1=0.80+(Y-0.10)*0.22; zf2=0.56+(Y-0.05)*0.48
    if side=='d':
        lobes={'lobe_cranial_d':m&(Z>=zf1),'lobe_moyen_d':m&(Z<zf1)&(Z>=zf2),
               'lobe_accessoire':m&(Z<zf2)&(Z<0.62)&(X>-0.075)&(Y<0.40),
               'lobe_caudal_d':m&(Z<zf2)&~((Z<0.62)&(X>-0.075)&(Y<0.40))}
    else:
        lobes={'lobe_cranial_cr_g':m&(Z>=zf1),'lobe_cranial_cd_g':m&(Z<zf1)&(Z>=zf2),'lobe_caudal_g':m&(Z<zf2)}
    from scipy.ndimage import binary_dilation as _bd
    for k_,m_ in lobes.items():
        oth=np.zeros_like(m)
        for k2,m2 in lobes.items():
            if k2!=k_: oth|=m2
        mm_=m_&~_bd(oth,iterations=1)
        out[k_]=('o',voxmesh(mm_,np.array([xs[0],ys[0],zs[0]]),h,smooth=10,sigma=1.1)); print(' ',k_,len(out[k_][1].faces))
pickle.dump({'zd':zd,'xs':xs,'ys':ys},open('diaph.pkl','wb'))
# --- trachea, bronchi, esophagus
carina=np.array([0,0.44,0.79])
out['trachee']=('o',tube([[0,0.50,1.49],[0,0.475,1.32],[0,0.45,1.13],[0,0.435,0.95],carina],0.028,0.026))
out['bronches']=('o',trimesh.util.concatenate([tube([carina,[-0.05,0.43,0.75],[-0.10,0.41,0.68]],0.018,0.012),
                                              tube([carina,[0.05,0.43,0.75],[0.10,0.41,0.68]],0.018,0.012)]))
cardia=np.array([0.024,0.514,0.411])
out['oesophage']=('o',tube([[0,0.585,1.47],[0.02,0.545,1.30],[0.018,0.52,1.12],[0.012,0.505,0.95],[0.008,0.50,0.75],[0.006,0.49,0.58],[0.015,0.505,0.47],cardia],0.017))
# --- brain: union of ellipsoids, then conform to skull interior by ray casting
h2=0.005; bx=np.arange(-0.14,0.14,h2); by=np.arange(0.52,0.90,h2); bz=np.arange(1.48,1.92,h2)
BX,BY,BZ=np.meshgrid(bx,by,bz,indexing='ij')
def ell(c,r): return ((BX-c[0])/r[0])**2+((BY-c[1])/r[1])**2+((BZ-c[2])/r[2])**2<1
brain=(ell([0.032,0.768,1.675],[0.052,0.066,0.094])|ell([-0.032,0.768,1.675],[0.052,0.066,0.094])|
       ell([0,0.745,1.598],[0.05,0.042,0.03])|ell([0,0.722,1.565],[0.022,0.022,0.04])|
       ell([0.012,0.74,1.752],[0.012,0.018,0.022])|ell([-0.012,0.74,1.752],[0.012,0.018,0.022]))
del X,Y,Z,P,t,rad,prof,bym,vert,med,dmed,below,cav,lung
bm=voxmesh(brain,np.array([bx[0],by[0],bz[0]]),h2,smooth=10,sigma=1.0)
import fast_simplification as fs
vv,ff=fs.simplify(bm.vertices.astype(np.float32),bm.faces.astype(np.int32),target_reduction=0.85); bm=trimesh.Trimesh(vv,ff)
print('brain faces',len(bm.faces))
sk=B['crane'].copy(); vv,ff=fs.simplify(sk.vertices.astype(np.float32),sk.faces.astype(np.int32),target_reduction=0.8); sk=trimesh.Trimesh(vv,ff); ctr=np.array([0,0.76,1.67])
v=np.asarray(bm.vertices); d=v-ctr; dist=np.linalg.norm(d,axis=1); dirs=d/dist[:,None]
loc,ri,_=sk.ray.intersects_location(np.repeat(ctr[None],len(v),0),dirs,multiple_hits=False)
first=np.full(len(v),9.0); 
for p,i in zip(loc,ri): first[i]=min(first[i],np.linalg.norm(p-ctr))
# also conform to the real inner wall of the registered CT skull (Edinburgh)
ed=pickle.load(open('edi_small_reg.pkl','rb'))
loc2,ri2,_=ed.ray.intersects_location(np.repeat(ctr[None],len(v),0),dirs,multiple_hits=False)
first2=np.full(len(v),9.0)
for p_,i_ in zip(loc2,ri2): first2[i_]=min(first2[i_],np.linalg.norm(p_-ctr))
first=np.minimum(first,np.where(first2<9,first2,first))
lim=first-0.004; newd=np.where(dist>lim,lim,dist)
print('brain vertices clipped',(dist>lim).mean().round(3))
bm.vertices=ctr+dirs*newd[:,None]; trimesh.smoothing.filter_taubin(bm,iterations=5)
out['encephale']=('n',bm)
pickle.dump(out,open('recon_thorax.pkl','wb'))
for k,(l,m) in out.items(): print(k,l,len(m.faces))
