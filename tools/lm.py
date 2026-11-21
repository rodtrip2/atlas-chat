import pickle,json,numpy as np,trimesh
from collections import defaultdict
parts=pickle.load(open('bones_parts.pkl','rb')); lab=json.load(open('data/bone_labels.json'))
G=defaultdict(list)
for i,p in enumerate(parts): G[lab[str(i)]].append(p)
B={k:trimesh.util.concatenate(v) for k,v in G.items()}
def canal(m):
    z=m.centroid[2]; best=None
    for dz in [0,-0.01,0.01,-0.02,0.02]:
        sec=m.section(plane_origin=[0,0,z+dz],plane_normal=[0,0,1])
        if sec is None: continue
        p2,T=sec.to_2D()
        for poly in p2.polygons_full:
            for h in poly.interiors:
                from shapely.geometry import Polygon
                hp=Polygon(h)
                if hp.area>1e-4 and (best is None or hp.area>best[0]):
                    c=np.array([hp.centroid.x,hp.centroid.y,0,1]); w=T@c
                    best=(hp.area,w[:3],np.sqrt(hp.area/np.pi))
        if best: break
    return best
if __name__=='__main__':
    seq=[f'C{i}' for i in range(1,8)]+[f'T{i}' for i in range(1,14)]+[f'L{i}' for i in range(1,8)]+['sacrum']+[f'Cd{i}' for i in range(1,4)]
    out={}
    for v in seq:
        r=canal(B[v]); m=B[v]
        print(v, np.round(m.centroid,3), 'ymin',round(m.bounds[0][1],3),'ymax',round(m.bounds[1][1],3), '-> canal', None if r is None else (np.round(r[1],3), round(r[2],3)))
        out[v]=None if r is None else [list(map(float,r[1])),float(r[2])]
    json.dump(out,open('data/canal.json','w'))
