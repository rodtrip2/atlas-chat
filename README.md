# Atlas 3D du chat

Atlas anatomique 3D interactif et gratuit du chat domestique (*Felis catus*), conçu pour les étudiants vétérinaires.

**➜ Consulter l'atlas : https://rodtrip2.github.io/atlas-chat/**

> ⚠️ **Version de travail.** Les noms des os et les régions musculaires ont été attribués automatiquement d'après leur position. Plusieurs structures (thorax, vaisseaux, nerfs) sont des reconstructions schématiques. Tout doit être relu par un enseignant en anatomie vétérinaire avant un usage pédagogique.

## Fonctionnalités

- 5 couches : squelette, muscles, organes, vaisseaux, nerfs. Le squelette, les muscles et les organes peuvent être opaques, transparents ou masqués.
- Environ 280 structures cliquables, chacune avec son nom en français et en latin et une courte description.
- Le bouton « ✦ Réseaux » affiche les vaisseaux et les nerfs à travers le squelette et les organes rendus transparents.
- Recherche par nom (par exemple « T5 », « fémur », « saphène »).
- Deux résolutions : **Standard** (environ 1 million de triangles, pour téléphone) et **HD** (environ 2,7 millions, pour ordinateur), à choisir avec le bouton « ◐ ».

## Origine des données

| Couche | Source | Nature |
|---|---|---|
| Squelette (144 structures) | « Feline Skeleton Dec », westerly | Modèle réel, os séparés et nommés automatiquement |
| Muscles (7 régions) | « Feline Ecorche Dec », westerly | Modèle réel, aligné sur le squelette, regroupé par régions |
| Organes abdominaux | « Normal Feline Anatomy – Digestive System », Pixelbeaker | Issu d'un scanner, **autre chat**, redimensionné et recalé approximativement |
| Cœur, poumons, trachée, œsophage, encéphale | Reconstruction | Ajustés sur les repères du squelette |
| Vaisseaux et nerfs (hors abdomen) | Reconstruction | Ajustés sur les repères du squelette |

### Méthode de reconstruction

- **Moelle épinière :** centre du canal vertébral relevé vertèbre par vertèbre (coupe de chaque vertèbre).
- **Poumons :** volume réel de la cage thoracique (côtes, cartilages, sternum), limité caudalement par la surface crâniale du foie, qui sert de diaphragme.
- **Encéphale :** dimensions mesurées de la cavité crânienne, puis forme ajustée à la paroi interne du crâne.
- **Vaisseaux et nerfs intercostaux :** tracés le long de la ligne centrale de chaque côte réelle.
- **Membres :** trajets calculés à partir de l'axe de chaque os (faces médiale, caudale, latérale) ; branches digitales tracées jusqu'aux phalanges réelles.
- **Arbre bronchique et vaisseaux pulmonaires :** générés par croissance (algorithme de colonisation de l'espace) dans le volume de chaque poumon.

Les scripts sont dans [`tools/`](tools/).

## Crédits et licences

Modèles 3D d'origine, sous licence [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/deed.fr) :

- « [Feline Skeleton Dec](https://sketchfab.com/3d-models/feline-skeleton-dec-a5d0bb8f55dc4f49b103cd20d65e0b17) », par [westerly](https://sketchfab.com/westerly)
- « [Feline Ecorche Dec](https://sketchfab.com/3d-models/feline-ecorche-dec-4c08824eb4914ff4b6e71b4dc2ba080e) », par [westerly](https://sketchfab.com/westerly)
- « [Normal Feline Anatomy – Digestive System](https://sketchfab.com/3d-models/normal-feline-anatomy-digestive-system-ab92a2a9d9024447b0ed4c080ad17106) », par [Pixelbeaker](https://sketchfab.com/pixelbeaker)

Modifications apportées : géométrie simplifiée, structures séparées et nommées, organes abdominaux redimensionnés et repositionnés, peau du modèle digestif retirée, structures thoraciques, vaisseaux et nerfs ajoutés par reconstruction.

- Les **modèles 3D dérivés** (`models/`) sont diffusés sous licence **CC BY 4.0**, avec les attributions ci-dessus.
- Le **code** (page web et scripts) est diffusé sous licence **MIT** (voir [`LICENSE`](LICENSE)).

## Contribuer

Une erreur de nomenclature, un trajet inexact ? Ouvrez un [ticket](https://github.com/rodtrip2/atlas-chat/issues) en indiquant la structure concernée (nom affiché dans sa fiche) et la correction proposée, avec une référence si possible. Les corrections de noms d'os se font dans [`tools/data/bone_labels.json`](tools/data/bone_labels.json) et dans les dictionnaires de `index.html`.

## Reconstruire les modèles

1. Télécharger les trois modèles Sketchfab ci-dessus au format glTF/GLB et les placer dans `tools/sources/` sous les noms `feline_skeleton_dec.glb`, `feline_ecorche_dec.glb` et `normal_feline_anatomy_digestive_system.glb`.
2. Installer les dépendances : `pip install trimesh numpy scipy scikit-image shapely rtree fast-simplification` et `npm i -g @gltf-transform/cli`.
3. Dans `tools/` : `python split_bones.py`, puis `python recon_thorax.py`, `python recon_vn.py`, `python recon_net.py`, `python recon_tree.py`, `python recon_branch.py`, puis `python build.py` (standard) et `Q=hd python build.py` (HD).
4. Compresser : `gltf-transform meshopt cat_raw.glb ../models/atlas-standard.glb` et `gltf-transform meshopt cat_raw_hd.glb ../models/atlas-hd.glb`.
