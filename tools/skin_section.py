# Module « Coupe de peau » : bloc schématique des couches et annexes cutanées
import numpy as np,trimesh,pickle
from common import tube, ball
from collections import defaultdict
OUT=defaultdict(list); LAY={}
O=np.array([0,0,3.4])
def add(id,m): m=m.copy(); m.apply_translation(O); OUT[id].append(m); LAY[id]='z'
def slab(y0,y1,inset=0.0):
    b=trimesh.creation.box(extents=[0.4-inset,abs(y1-y0),0.4-inset]); b.apply_translation([0,(y0+y1)/2,0]); return b
def slab_cut(y0,y1):
    t=abs(y1-y0); yc=(y0+y1)/2; bs=[]
    for (cx,cz,ex,ez) in [(-0.1,0,0.2,0.4),(0.1,-0.1,0.2,0.2)]:
        b=trimesh.creation.box(extents=[ex,t,ez]); b.apply_translation([cx,yc,cz]); bs.append(b)
    return trimesh.util.concatenate(bs)
def ell(c,r,sub=2):
    m=trimesh.creation.icosphere(subdivisions=sub,radius=1.0); m.vertices*=np.asarray(r); m.apply_translation(c); return m
add('sk_cornee',slab_cut(0.300,0.2965)); add('sk_epi_vivant',slab_cut(0.2965,0.2885)); add('sk_basale',slab_cut(0.2885,0.2855))
add('sk_derme_pap',slab_cut(0.2855,0.262)); add('sk_derme_ret',slab(0.262,0.18)); add('sk_hypoderme',slab(0.18,0.07)); add('sk_peaucier',slab(0.07,0.05))
rng=np.random.default_rng(3)
for i in range(70):
    c=[rng.uniform(-0.18,0.18),rng.uniform(0.085,0.165),rng.uniform(-0.18,0.18)]; add('sk_lobules',ell(c,[0.018,0.014,0.018],sub=1))
for gx in (-0.12,0,0.12):
    for gz in (-0.12,0,0.12):
        f=np.array([gx,0,gz]); pore=f+[0,0.300,0]
        add('sk_foll_primaire',tube([pore,f+[0,0.25,0],f+[0,0.205,0]],0.0065,0.0075))
        add('sk_bulbe',ell(f+[0,0.200,0],[0.011,0.012,0.011]))
        add('sk_poil_jarre',tube([f+[0,0.205,0],pore,pore+[0.03,0.08,0.01],pore+[0.075,0.16,0.02]],0.0028,0.0012,seg=8))
        for k in range(5):
            a=2*np.pi*k/5; o=np.array([np.cos(a)*0.018,0,np.sin(a)*0.018])
            add('sk_foll_secondaire',tube([f+o+[0,0.24,0],f+o*0.5+[0,0.285,0],pore+[0,-0.006,0]],0.0028,0.002,seg=8))
            add('sk_poil_bourre',tube([f+o+[0,0.24,0],pore+o*0.2,pore+o*0.4+[0.025,0.06,0.0],pore+o*0.6+[0.045,0.10,0.0]],0.0012,0.0006,seg=5))
        sg=f+[-0.012,0.262,0.006]
        add('sk_sebacee',trimesh.util.concatenate([ell(sg+o,[0.006,0.005,0.006],sub=1) for o in ([0,0,0],[-0.006,0.004,0.004],[-0.005,-0.004,-0.003],[-0.009,0,-0.001])]))
        add('sk_arrecteur',tube([f+[-0.007,0.235,0],f+[-0.025,0.258,0],f+[-0.04,0.283,0]],0.0018,0.0014,seg=6))
        gc=f+[0.03,0.165,0.02]; hel=[gc+[0.008*np.cos(t),0.0025*t/np.pi*0.6-0.006,0.008*np.sin(t)] for t in np.linspace(0,6*np.pi,60)]
        add('sk_sudoripare',tube(hel+[gc+[0.002,0.02,0],f+[0.01,0.25,0],f+[0.006,0.275,0]],0.0028,0.0018,seg=6,step=0.003))
# vascular plexuses
for y,rr in [(0.182,0.004),(0.268,0.0022)]:
    for z in np.linspace(-0.18,0.18,5):
        add('sk_arterioles',tube([[-0.2,y,z+0.006],[0,y+0.003,z+0.004],[0.2,y,z+0.006]],rr))
        add('sk_veinules',tube([[-0.2,y-0.004,z-0.008],[0,y-0.002,z-0.006],[0.2,y-0.004,z-0.008]],rr*1.15))
for x in np.linspace(-0.16,0.16,6):
    for z in np.linspace(-0.16,0.16,6):
        add('sk_arterioles',tube([[x,0.182,z],[x+0.004,0.225,z],[x,0.268,z]],0.0016))
        add('sk_capillaires',tube([[x+0.01,0.268,z],[x+0.01,0.283,z],[x+0.016,0.283,z],[x+0.016,0.268,z]],0.0009,seg=5,step=0.003))
# nerves and tylotrich pad
for z in (-0.06,0.06):
    add('sk_nerfs',tube([[-0.2,0.20,z],[0.2,0.20,z]],0.0022))
    for x in (-0.12,0,0.12):
        add('sk_nerfs',tube([[x,0.20,z],[x+0.004,0.225,z+(0.06 if z<0 else -0.06)*0.5],[x+0.004,0.245,(-0.12 if z<0 else 0.12)*0+z*2]],0.0012))
add('sk_tylotriche',ell([0.06,0.301,0.06],[0.012,0.004,0.012]))
add('sk_tylotriche',tube([[0.06,0.2,0.06],[0.06,0.28,0.06],[0.06,0.296,0.06]],0.0014))
pickle.dump((dict(OUT),LAY),open('net_section.pkl','wb'))
for k,v in OUT.items(): print(k,len(v),sum(len(m.faces) for m in v))
