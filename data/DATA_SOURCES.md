# DATA_SOURCES.md
## Traçabilité des données — Projet Prédiction du sous-type moléculaire du cancer du sein

Ce fichier documente l'origine exacte de chaque fichier de données utilisé dans le projet. À compléter au fur et à mesure des téléchargements, pas a posteriori.

---

## 1. Expression génique (RNA-seq)

- **Fichier local** : `data/raw/tcga/TCGA-BRCA.star_tpm.tsv.gz`
- **Source** : UCSC Xena, hub GDC (`https://gdc.xenahubs.net`)
- **URL de téléchargement** : https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.star_tpm.tsv.gz
- **Cohorte** : GDC TCGA Breast Cancer (BRCA)
- **Date de téléchargement** : À COMPLÉTER (jj/mm/aaaa)
- **Nombre d'échantillons** : 1226
- **Version des données** : 05-20-2024 (version GDC affichée sur la page Xena)
- **Unité** : log2(TPM + 1)
- **Pipeline de génération** : STAR (voir https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)
- **Licence / conditions d'usage** : Données publiques GDC/TCGA, usage libre avec citation de la source (voir https://gdc.cancer.gov/about-data/publications/data-use)
- **Nombre de lignes/colonnes après chargement** : À COMPLÉTER

## 2. Table de correspondance ID Ensembl → nom de gène

- **Fichier local** : `data/raw/tcga/gencode.v36.annotation.gtf.gene.probemap`
- **Source** : UCSC Xena, hub GDC
- **URL de téléchargement** : https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap
- **Date de téléchargement** : À COMPLÉTER
- **Usage** : traduire les identifiants Ensembl (ENSG...) de la matrice d'expression en noms de gènes lisibles (ex. TP53, ESR1)
- **Licence / conditions d'usage** : GENCODE v36, domaine public / usage académique libre

## 3. Données cliniques (Phenotype)

- **Fichier local** : `data/raw/tcga/TCGA-BRCA.GDC_phenotype.tsv.gz`
- **Source** : UCSC Xena, hub GDC
- **URL de téléchargement** : À COMPLÉTER (récupérer le lien "download" de "Phenotype (n=1,255) GDC Hub" sur la page de la cohorte)
- **Date de téléchargement** : À COMPLÉTER
- **Nombre d'échantillons** : 1255
- **Contenu** : variables cliniques (âge, stade, statut ER/PR/HER2 si présent, etc.)
- **Licence / conditions d'usage** : Données publiques GDC/TCGA

## 4. Données de survie

- **Fichier local** : `data/raw/tcga/TCGA-BRCA.survival.tsv.gz`
- **Source** : UCSC Xena, hub GDC
- **URL de téléchargement** : À COMPLÉTER (lien "download" de "survival data (n=1,232) GDC Hub")
- **Date de téléchargement** : À COMPLÉTER
- **Nombre d'échantillons** : 1232
- **Contenu** : temps de suivi et statut vital (pour analyse de survie Kaplan-Meier / Cox)
- **Licence / conditions d'usage** : Données publiques GDC/TCGA

## 5. Labels de sous-type moléculaire PAM50

- **Fichier local** : `data/raw/cbioportal_tcga/data_clinical_patient.txt` (ou export .csv équivalent)
- **Source** : cBioPortal — étude "Breast Invasive Carcinoma (TCGA, PanCancer Atlas)"
- **URL** : https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
- **Date de téléchargement** : À COMPLÉTER
- **Colonne utilisée** : `SUBTYPE` (valeurs attendues : BRCA_LumA, BRCA_LumB, BRCA_Her2, BRCA_Basal, BRCA_Normal)
- **Licence / conditions d'usage** : cBioPortal, usage académique libre avec citation (voir https://www.cbioportal.org/faq)

---

## 6. Validation externe (à compléter en Phase 6 uniquement)

### Expression génique + clinique — METABRIC

- **Fichier local** : `data/raw/metabric/` (à créer au moment de l'utilisation)
- **Source** : cBioPortal — étude "Breast Cancer (METABRIC)"
- **URL** : https://www.cbioportal.org/study/summary?id=brca_metabric
- **Date de téléchargement** : À COMPLÉTER
- **Nombre d'échantillons** : ~2000
- **Plateforme technique** : microarray (à distinguer du RNA-seq de TCGA — attention lors de l'harmonisation, cf. phase 6 du guide)
- **Licence / conditions d'usage** : cBioPortal, usage académique libre avec citation

---

## Notes générales

- Toutes les données utilisées sont publiques et anonymisées, en accès libre pour usage académique/recherche.
- Aucune donnée d'identification patient n'est présente dans les fichiers (identifiants pseudonymisés type `TCGA-XX-XXXX`).
- En cas de mise à jour d'un fichier source (nouvelle version GDC, nouvelle release cBioPortal), noter ici l'ancienne et la nouvelle version pour garder une trace de reproductibilité.
- Toute transformation appliquée aux données brutes (filtrage, normalisation, imputation) doit être documentée séparément dans `report/gene_selection_justification.md` et dans le code source (`src/`), pas dans ce fichier — ce fichier ne concerne que la provenance des données brutes.