import numpy as np, trimesh
def classify(c):
    x,y,z=c[:,0],c[:,1],c[:,2]; ax=np.abs(x)
    side=np.where(x>=0,'_g','_d')
    lab=np.full(len(c),'base',dtype=object)
    # order matters: later rules override earlier
    lab[(z<1.705)&(y>=0.60)]='temporal'                       # lateral braincase default
    lab[(z<1.705)&(y>=0.70)&(ax<0.075)]='parietal'
    lab[(z>=1.705)&(z<1.845)&(y>=0.66)]='frontal'
    lab[(z>=1.705)&(z<1.845)&(y>=0.66)&(ax>0.14)]='frontal'   # postorbital process stays frontal
    lab[(z>=1.705)&(y<0.66)&(ax<0.09)]='base'
    lab[(z>=1.76)&(ax>=0.06)&(y<0.66)]='maxillaire'
    lab[(z>=1.845)&(ax>=0.035)]='maxillaire'
    lab[(z>=1.845)&(ax<0.035)&(y>=0.55)]='nasal'
    lab[(z>=1.845)&(ax<0.035)&(y<0.55)]='maxillaire'
    lab[(z>=1.93)]='incisif'
    lab[(ax>=0.15)&(z>=1.62)&(z<1.80)&(y<0.70)]='zygomatique'
    lab[(z<1.575)]='occipital'
    lab[(z<1.64)&(y<0.62)&(ax>=0.04)&(ax<0.15)]='bulle'
    lab[(z>=1.575)&(z<1.705)&(y<0.60)&(ax<0.04)]='base'
    lab[(lab=='frontal')&(z<1.72)&(ax<0.075)&(y>0.70)]='parietal'
    out=np.array([l+s if l not in('base','occipital') else l for l,s in zip(lab,side)],dtype=object)
    return out
def split(mesh):
    c=mesh.triangles_center; L=classify(c); res={}
    for k in np.unique(L):
        idx=np.nonzero(L==k)[0]; sub=mesh.submesh([idx],append=True,repair=False); res[k]=sub
    return res
