# 🌍 Qualité GPS dans les Enquêtes Agricoles : Pipeline de Validation Multi-Critères

Ce dépôt contient l'intégralité du code et des livrables du projet **"Analyse des Incohérences GPS dans les Données de l'EACI 2017 au Mali — Détection, Classification et Traitement des Anomalies Géographiques"**.

Ce travail a été réalisé par **ZONGO Abdoul Kader** dans le cadre de la formation des Ingénieurs Statisticiens Économistes (ISEP2) à l'**ENSAE Pierre Ndiaye** de Dakar, sous l'encadrement de M. **DIALLO Mouhamadou Hady** (Année académique 2025-2026).

---

## 🎯 Objectif du Projet

Le Mali est un pays à forte vocation agricole où la qualité des données GPS conditionne directement la fiabilité des analyses de superficie et de rendement. L'EACI 2017 (Enquête Agricole de Conjoncture Intégrée), conduite dans le cadre du programme **LSMS-ISA de la Banque mondiale**, fournit les coordonnées GPS de 953 exploitations agricoles géoréférencées.

L'objectif de ce projet est de construire un **pipeline de validation géographique systématique sous R**, capable de détecter, classifier et traiter les anomalies GPS en anticipant explicitement trois pièges méthodologiques majeurs :

- **Piège 1** : Confondre une erreur de projection avec une vraie anomalie GPS
- **Piège 2** : Exclure des doublons pourtant légitimes (jittering LSMS)
- **Piège 3** : Introduire un biais de sélection non documenté lors des exclusions

---

## 📊 Sources de Données Utilisées

| Source | Contenu | Rôle |
|---|---|---|
| `eaci_geovariables_2017.csv` | 953 exploitations + 36 variables géographiques (NDVI, TWI, élévation, précipitation...) | Base principale GPS |
| `eaci17_s01p1.csv` | Identifiants ménages : grappe, exploitation, **passage** | Traitement Piège 2 |
| `rnaturalearth` (package R) | Polygones des frontières administratives du Mali | Validation territoriale |

> **Source principale :** MLI_2017_EAC-I_v03, Banque mondiale — Programme LSMS-ISA

---

## ⚙️ Choix Méthodologiques Clés

Le pipeline intègre plusieurs techniques statistiques avancées pour garantir la rigueur de l'analyse :

- **Double projection spatiale (EPSG:3857 + UTM Zone 29N)** : Un point hors territoire en Web Mercator peut être une simple erreur de projection corrigible — la double vérification évite toute exclusion abusive *(Piège 1)*.

- **Clustering spatial DBSCAN (eps = 500 m en mètres)** : Contrairement à k-means, DBSCAN ne requiert pas de fixer le nombre de clusters et détecte les outliers naturellement. Le rayon de 500 m est calibré sur l'amplitude connue du jittering LSMS. Croisé avec la variable `passage`, il distingue jittering légitime, double visite et fraude probable *(Piège 2)*.

- **Quantile 95% non-paramétrique** : La distribution des distances au centroïde est asymétrique à droite — l'hypothèse de normalité n'est pas vérifiée. Le Q95% est robuste à toute forme de distribution, contrairement à la règle classique d̄ + 2σ.

- **Score composite additif (0 à 5)** : La convergence de plusieurs critères renforce la certitude de l'anomalie. Aucune observation n'est exclue sur la seule base d'un critère unique.

- **Vérification du biais de sélection** : Comparaison systématique des caractéristiques géographiques des exclus vs conservés sur cinq indicateurs *(Piège 3)*.

---

## 🏆 Résultats Clés

| Indicateur | Valeur |
|---|---|
| Total observations analysées | 953 |
| Observations de **bonne qualité** (Score 0) | **473 (49.6%)** |
| Observations à vérifier (Score 1-2) | 463 (48.6%) |
| Observations suspectes (Score > 2) | 17 (1.8%) |
| Doublons détectés — dont jittering légitime | **207 (dont 207, soit 100%)** |
| Outliers spatiaux (Q95%) | 48 (5%) |
| **Observations exclues définitivement** | **1 sur 953 (0.1%)** |
| Surface cultivée moyenne estimée | **0.56 ha** |
| Biais NDVI exclus vs conservés | **3% < seuil critique 10%** |

> **Résultat majeur :** Les 207 doublons détectés sont intégralement attribuables au jittering LSMS légitime. Une exclusion naïve aurait supprimé 21.7% de la base — un biais considérable évité grâce à la variable `passage`.

---

## 📂 Structure du Dépôt

```text
exposeZongoGPS2026ENSAE/
├── data/                        # Données brutes (EACI GPS, identifiants ménages)
├── scripts/                     # Scripts R modulaires
│   ├── 00_setup.R               # Installation et chargement des packages
│   ├── 01_validation_geo.R      # Étape 1 : Validation géographique (Piège 1)
│   ├── 02_taxonomie.R           # Étape 2 : Taxonomie des erreurs GPS
│   ├── 03_doublons_dbscan.R     # Étape 3 : Détection doublons DBSCAN (Piège 2)
│   ├── 04_outliers.R            # Étape 4 : Outliers spatiaux Q95%
│   ├── 05_score_decisions.R     # Étape 5 : Score composite et décisions
│   └── 06_sensibilite_viz.R     # Étape 6 : Analyse de sensibilité et visualisations
├── output/                      # Résultats générés
│   ├── figures/                 # Graphiques (carte, histogramme, sensibilité...)
│   └── base_nettoyee.csv        # Base GPS corrigée finale
├── beamer.pdf                   # Support de présentation (Slides Beamer)
├── rmd_WorkingFile.pdf          # Rapport méthodologique et analytique final
├── rmd_WorkingFile.Rmd          # Code source R Markdown du rapport
└── script_Rtraitement.R        # Script maître orchestrant l'ensemble du pipeline
```

---

## 🚀 Comment Reproduire l'Analyse ?

L'un des piliers de ce projet est la **reproductibilité totale**. Le pipeline est entièrement codé sous R, documenté et portable.

### 1. Cloner le dépôt

```bash
git clone https://github.com/[votre-username]/exposeZongoGPS2026ENSAE.git
cd exposeZongoGPS2026ENSAE
```

### 2. Installer les packages requis

```r
# Exécuter le script de setup
source("scripts/00_setup.R")
```

Packages principaux utilisés :

```r
install.packages(c("sf", "dplyr", "ggplot2", "dbscan",
                   "rnaturalearth", "rnaturalearthdata",
                   "units", "tidyr"))
```

### 3. Exécuter le pipeline complet

```r
# Lancer le script maître
source("script_Rtraitement.R")
```

> ⚠️ **Note importante :** Toutes les distances et opérations géométriques sont calculées en **EPSG:3857 (mètres)**. Les calculs en degrés décimaux (WGS84) sont métriquement incohérents et ne peuvent pas être utilisés pour DBSCAN.

---

## 🗺️ Pipeline en 6 Étapes

```
Étape 1 — Validation géographique    →  st_intersects() × 2 projections + TWI + NDVI   [Piège 1]
Étape 2 — Taxonomie des erreurs      →  Comparaison projections + NDVI + dist_road
Étape 3 — Doublons DBSCAN            →  eps=500m (EPSG:3857) + variable passage          [Piège 2]
Étape 4 — Outliers spatiaux          →  Quantile 95% des distances au centroïde
Étape 5 — Score et décisions         →  7 règles hiérarchisées et justifiées             [Piège 3]
Étape 6 — Analyse de sensibilité     →  Comparaison Q90%, Q95%, Q97.5%                   [Piège 3]
```

---

## 📝 Décharge

L'auteur se porte entièrement responsable des propos et positions tenus dans ce document. Ainsi, ceux-ci n'engagent ni l'École Nationale de la Statistique et de l'Analyse Économique (ENSAE Pierre Ndiaye), ni la Banque Mondiale, ni aucune autre institution mentionnée dans ce rapport.

---

*Rapport produit avec R Markdown — ENSAE Pierre Ndiaye, Dakar — 19 mai 2026*
