# Étape 1 : sépare le squelette westerly en os individuels (composantes connexes)
import trimesh, pickle
sc = trimesh.load('sources/feline_skeleton_dec.glb')
m = sc.to_geometry(); m.merge_vertices()
parts = list(m.split(only_watertight=False))
pickle.dump(parts, open('bones_parts.pkl', 'wb'))
print(len(parts), 'pièces')
