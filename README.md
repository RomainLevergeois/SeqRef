# SeqRef

[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Description

Votre projet est un générateur de texte qui permet de générer une séquence FASTA d'un transcrit (à partir du nom RefSeq) avec la traductin alignée sur le cadre de lecture. Il utilise plusieurs dépendances pour gérer des fonctionnalités telles que l'accès à des bases de données biologiques comme NCBI et UCSC.

### Fonctionnalités principales :
- Requête de données génomiques à partir de **NCBI REST API Datasets v2**
- Requête de données via **UCSC REST API**

## Installation

### Prérequis

Avant d'installer ce projet, assurez-vous que les éléments suivants sont installés sur votre machine :
- Python 3.x
- Pip (gestionnaire de paquets Python)

### Étapes d'installation

1. Clonez le dépôt du projet :
```bash
git clone https://github.com/votre-utilisateur/votre-projet.git
```

2. Accédez au répertoire du projet :
```bash
cd SeqRef
```

3. Installez les dépendances nécessaires :
```bash
pip install -r requirements.txt
```

## Dépendances

Ce projet inclut les dépendances suivantes :
* autocomplete-lhc : Un outil d'autocomplétion développé par le Lister Hill National Center for Biomedical Communications (LHNCBC).
* NCBI REST API Datasets v2 : Une API permettant d'accéder aux métadonnées et séquences biologiques.
* UCSC REST API : Une API pour accéder aux données du UCSC Genome Browser.

Les détails complets des licences pour ces dépendances se trouvent dans le fichier LICENSE.

## Citations et Références

Si vous utilisez ce projet dans une publication académique ou un projet de recherche, merci de citer les dépendances comme suit :
- autocomplete-lhc : Lister Hill National Center for Biomedical Communications, National Library of Medicine, Bethesda, MD.
- NCBI REST API Datasets v2 :
O’Leary NA, et al. (2024). Sci Data. 11(1):732. doi: 10.1038/s41597-024-03571-y.
- UCSC REST API :
Lee CM, et al. (2020). Nucleic Acids Res. 48(D1).

## Licence

Ce projet est sous licence MIT. Voir le fichier LICENSE pour plus d'informations.
Auteurs

* Romain LEVERGEOIS – Contactez-moi

Contributions

Les contributions sont les bienvenues ! Si vous souhaitez contribuer à ce projet, merci de suivre ces étapes :
1. Forker le projet
2. Créer une nouvelle branche (git checkout -b nouvelle-fonctionnalité)
3. Commiter vos modifications (git commit -am 'Ajout d'une nouvelle fonctionnalité')
4. Pousser sur la branche (git push origin nouvelle-fonctionnalité)
5. Ouvrir une Pull Request

## Remerciements

Merci aux auteurs des bibliothèques et APIs utilisées dans ce projet.
