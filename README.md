# Atlas 3D du chat

Atlas anatomique 3D interactif et gratuit du chat domestique (*Felis catus*), conçu pour les étudiants vétérinaires.

**➜ Consulter l'atlas : https://rodtrip2.github.io/atlas-chat/**

> ⚠️ **Version de travail.** Les noms des os et les régions musculaires ont été attribués automatiquement d'après leur position. Plusieurs structures (thorax, vaisseaux, nerfs) sont des reconstructions schématiques. Tout doit être relu par un enseignant en anatomie vétérinaire avant un usage pédagogique.

## Fonctionnalités

- 9 couches : peau, squelette, muscles, organes, vaisseaux, nerfs, système lymphatique, tendons et ligaments, repères osseux. Le squelette, les muscles et les organes peuvent être opaques, transparents ou masqués.
- Environ 690 structures cliquables, chacune avec son nom en français et en latin et une fiche : articulations, insertions musculaires et intérêt clinique pour les os.
- Os détaillés : petits os du carpe et du tarse, métacarpiens et métatarsiens I à V, phalanges doigt par doigt, 28 dents (numérotation de Triadan), sternèbres, os du crâne (délimitation approximative, le crâne du modèle étant d'un seul bloc).
- 90 repères osseux cliquables, calculés sur la géométrie réelle des os (tête du fémur, trochanters, condyles, épine scapulaire, olécrâne, foramen magnum, processus coronoïde…).
- Le bouton « ✦ Réseaux » affiche les vaisseaux et les nerfs à travers le squelette et les organes rendus transparents.
- Recherche par nom (par exemple « T5 », « fémur », « saphène »).
- Module **« ◎ Crâne détaillé »** : scanner haute résolution d'un crâne de chat (Université d'Édimbourg), recalé sur l'atlas, avec os, mandibule et dents identifiés et 26 repères.
- Module **« ▦ Coupe de peau »** : bloc schématique de l'épiderme (couches), du derme, de l'hypoderme et des annexes (follicules pileux composés, glandes sébacées et sudoripares, muscles arrecteurs, plexus vasculaires, nerfs, coussinet tylotriche).
- Deux résolutions : **Standard** (environ 1 million de triangles, pour téléphone) et **HD** (environ 2,7 millions, pour ordinateur), à choisir avec le bouton « ◐ ».

## Origine des données

| Couche | Source | Nature |
|---|---|---|
| Squelette (environ 250 structures) | « Feline Skeleton Dec », westerly | Modèle réel, os séparés et nommés automatiquement ; crâne découpé en os de façon approximative |
| Repères osseux (90) | Calcul | Points extrêmes ou caractéristiques de la géométrie de chaque os |
| Muscles (environ 40 muscles ou groupes) | « Feline Ecorche Dec », westerly | Modèle réel, aligné sur le squelette ; identification automatique proposée d'après les attaches osseuses (à valider) |
| Crâne détaillé (module) | « Domestic Cat (Felis catus) », Université d'Édimbourg Open.Ed | Scanner d'un **autre chat**, recalé sur le crâne de l'atlas (ICP avec mise à l'échelle) ; os identifiés par transfert |
| Organes abdominaux | « Normal Feline Anatomy – Digestive System », Pixelbeaker | Issu d'un scanner, **autre chat**, redimensionné et recalé approximativement |
| Cœur (4 cavités), poumons (lobes), trachée, larynx, œsophage, encéphale, yeux | Reconstruction | Ajustés sur les repères du squelette |
| Voies urinaires basses, appareils génitaux mâle et femelle (bouton ♂/♀), mamelles, rectum, appareil hyoïde | Reconstruction | Ajustés sur les repères du squelette et du bassin |
| Système lymphatique (nœuds palpables et profonds, conduit thoracique, thymus…) | Reconstruction | Positions anatomiques usuelles, ajustées sur le squelette |
| Tendons, ligaments, ménisques, disques intervertébraux | Reconstruction | Ancrés sur les repères osseux et les phalanges ; disques mesurés sur les corps vertébraux |
| Peau (enveloppe calculée autour de l'écorché et du squelette), coussinets, griffes, vibrisses, pavillons, truffe | Calcul / reconstruction | Enveloppe : volume des muscles et des os, fermé puis légèrement dilaté |
| Thyroïde, parathyroïdes, surrénales, glandes salivaires, langue, diaphragme | Reconstruction | Ajustés sur le squelette, le foie et la cavité thoracique |
| Vaisseaux et nerfs (hors abdomen) | Reconstruction | Ajustés sur les repères du squelette |

### Méthode de reconstruction

- **Moelle épinière :** centre du canal vertébral relevé vertèbre par vertèbre (coupe de chaque vertèbre).
- **Poumons :** volume réel de la cage thoracique (côtes, cartilages, sternum), limité caudalement par la surface crâniale du foie, qui sert de diaphragme.
- **Encéphale :** forme ajustée à la paroi interne réelle de la cavité crânienne du scanner d'Édimbourg recalé.
- **Vaisseaux et nerfs intercostaux :** tracés le long de la ligne centrale de chaque côte réelle.
- **Membres :** trajets calculés à partir de l'axe de chaque os (faces médiale, caudale, latérale) ; branches digitales tracées jusqu'aux phalanges réelles.
- **Arbre bronchique et vaisseaux pulmonaires :** générés par croissance (algorithme de colonisation de l'espace) dans le volume de chaque poumon.

Les scripts sont dans [`tools/`](tools/).

### Référence de vérification

Les planches félines de *Hill's Atlas of Veterinary Clinical Anatomy* (Hill's Pet Nutrition, 2007) (cœur, thorax, denture, appareil urinaire, œil) ont servi à contrôler et compléter la reconstruction : lobes pulmonaires, cavités cardiaques, voies urinaires basses, larynx, appareil hyoïde, yeux, formule dentaire. Aucune illustration de cet atlas n'est reproduite dans ce projet.

## Crédits et licences

Modèles 3D d'origine, sous licence [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/deed.fr) :

- « [Feline Skeleton Dec](https://sketchfab.com/3d-models/feline-skeleton-dec-a5d0bb8f55dc4f49b103cd20d65e0b17) », par [westerly](https://sketchfab.com/westerly)
- « [Feline Ecorche Dec](https://sketchfab.com/3d-models/feline-ecorche-dec-4c08824eb4914ff4b6e71b4dc2ba080e) », par [westerly](https://sketchfab.com/westerly)
- « [Normal Feline Anatomy – Digestive System](https://sketchfab.com/3d-models/normal-feline-anatomy-digestive-system-ab92a2a9d9024447b0ed4c080ad17106) », par [Pixelbeaker](https://sketchfab.com/pixelbeaker)
- « [Domestic Cat (Felis catus)](https://sketchfab.com/3d-models/domestic-cat-felis-catus-7c99ca836d834c39872ecf5e9b5e2087) », par [The University of Edinburgh Open.Ed](https://sketchfab.com/openededinburgh) (scanner de Tobias Schwarz, Royal (Dick) School of Veterinary Studies)

Modifications apportées : géométrie simplifiée, structures séparées et nommées, organes abdominaux redimensionnés et repositionnés, peau du modèle digestif retirée, structures thoraciques, vaisseaux et nerfs ajoutés par reconstruction.

- Les **modèles 3D dérivés** (`models/`) sont diffusés sous licence **CC BY 4.0**, avec les attributions ci-dessus.
- Le **code** (page web et scripts) est diffusé sous licence **MIT** (voir [`LICENSE`](LICENSE)).

## Contribuer

Une erreur de nomenclature, un trajet inexact ? Ouvrez un [ticket](https://github.com/rodtrip2/atlas-chat/issues) en indiquant la structure concernée (nom affiché dans sa fiche) et la correction proposée, avec une référence si possible. Les corrections de noms d'os se font dans [`tools/data/bone_labels.json`](tools/data/bone_labels.json) et dans les dictionnaires de `index.html`.

## Reconstruire les modèles

1. Télécharger les trois modèles Sketchfab ci-dessus au format glTF/GLB et les placer dans `tools/sources/` sous les noms `feline_skeleton_dec.glb`, `feline_ecorche_dec.glb`, `normal_feline_anatomy_digestive_system.glb` et `domestic_cat_felis_catus.glb`.
2. Installer les dépendances : `pip install trimesh numpy scipy scikit-image shapely rtree fast-simplification` et `npm i -g @gltf-transform/cli`.
3. Dans `tools/` : `python split_bones.py`, `python refine_labels.py`, `python landmarks.py`, puis `python recon_thorax.py`, `python recon_vn.py`, `python recon_net.py`, `python recon_tree.py`, `python recon_branch.py`, `python recon_extra.py`, `python recon_more.py`, `python recon_female.py`, `python recon_skin.py`, `python skin_section.py`, `python muscle_attach.py`, `python muscle_names.py`, puis pour le crâne : `python register_skull.py`, `python build_skull.py`, `python export_skull.py` (et `Q=hd`), `python update_landmarks.py`, puis `python build.py` (standard) et `Q=hd python build.py` (HD).
4. Compresser : `gltf-transform meshopt cat_raw.glb ../models/atlas-standard.glb` et `gltf-transform meshopt cat_raw_hd.glb ../models/atlas-hd.glb`.
