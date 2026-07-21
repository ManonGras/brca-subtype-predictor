"""
Phase 2 — Chargement, nettoyage et harmonisation des données TCGA-BRCA
========================================================================
Produit en sortie : data/processed/tcga_train.csv
(patients en lignes, gènes filtrés + variables cliniques en colonnes, + colonne 'subtype')

À exécuter depuis la racine du projet (breast-cancer-subtype/), avec :
    python src/prepare_dataset.py
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ----------------------------------------------------------------------------
# Chemins des fichiers (à ajuster si tes noms de fichiers diffèrent légèrement)
# ----------------------------------------------------------------------------
RAW_TCGA_DIR = Path("data/raw/tcga")
RAW_CBIO_DIR = Path("data/raw/cbioportal_tcga")
PROCESSED_DIR = Path("data/processed")
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

EXPRESSION_FILE = RAW_TCGA_DIR / "TCGA-BRCA.star_tpm.tsv"
PROBEMAP_FILE = RAW_TCGA_DIR / "gencode.v36.annotation.gtf.gene.probemap"
CLINICAL_FILE = RAW_TCGA_DIR / "TCGA-BRCA.clinical.tsv"
SURVIVAL_FILE = RAW_TCGA_DIR / "TCGA-BRCA.survival.tsv"
SUBTYPE_FILE = RAW_CBIO_DIR / "brca_pam50.tsv"  # export cBioPortal (colonnes: sample, PAM50)


# ----------------------------------------------------------------------------
# 2.1 — Chargement et inspection
# ----------------------------------------------------------------------------

def load_expression(path: Path) -> pd.DataFrame:
    """Charge la matrice d'expression et la transpose en (patients x gènes)."""
    print(f"Chargement de l'expression depuis {path} ...")
    # Fichier volumineux (~900 Mo) : on force les colonnes numériques en float32
    # pour réduire l'empreinte mémoire de moitié par rapport au float64 par défaut.
    expr = pd.read_csv(path, sep="\t", index_col=0)
    expr = expr.astype("float32")
    print(f"  -> forme brute (gènes x échantillons) : {expr.shape}")

    # Transposition : on veut patients en lignes, gènes en colonnes
    expr = expr.T
    expr.index.name = "sample_barcode"
    print(f"  -> forme après transposition (échantillons x gènes) : {expr.shape}")
    return expr


def load_probemap(path: Path) -> dict:
    """Charge la table de correspondance ID Ensembl -> nom de gène."""
    print(f"Chargement du probemap depuis {path} ...")
    probemap = pd.read_csv(path, sep="\t")
    # Le fichier probemap a généralement les colonnes : id, gene, chrom, chromStart, chromEnd, strand
    # On vérifie les noms de colonnes réels au cas où ils diffèrent :
    print(f"  -> colonnes disponibles : {list(probemap.columns)}")
    id_col = probemap.columns[0]
    gene_col = "gene" if "gene" in probemap.columns else probemap.columns[1]
    mapping = dict(zip(probemap[id_col], probemap[gene_col]))
    print(f"  -> {len(mapping)} correspondances ID -> gène chargées.")
    return mapping


def load_clinical(path: Path) -> pd.DataFrame:
    """Charge la matrice clinique (phenotype)."""
    print(f"Chargement des données cliniques depuis {path} ...")
    clinical = pd.read_csv(path, sep="\t")
    print(f"  -> forme : {clinical.shape}")
    print(f"  -> colonnes disponibles (premières 15) : {list(clinical.columns[:15])}")
    return clinical


def load_survival(path: Path) -> pd.DataFrame:
    """Charge les données de survie."""
    print(f"Chargement des données de survie depuis {path} ...")
    survival = pd.read_csv(path, sep="\t")
    print(f"  -> forme : {survival.shape}")
    print(f"  -> colonnes disponibles : {list(survival.columns)}")
    return survival


def load_subtype(path: Path) -> pd.DataFrame:
    """
    Charge le fichier PAM50 (colonnes attendues : 'sample' et 'PAM50').
    L'identifiant 'sample' est ici au format barcode complet (ex. TCGA-A7-A13F-01A),
    donc directement comparable aux barcodes de la matrice d'expression.
    """
    print(f"Chargement des labels PAM50 depuis {path} ...")
    subtype = pd.read_csv(path, sep="\t", comment="#")
    print(f"  -> colonnes disponibles : {list(subtype.columns)}")
    return subtype


def find_sample_id_column(df: pd.DataFrame, candidates=("sample", "sampleID", "submitter_id.samples", "_PATIENT")) -> str:
    """Cherche automatiquement la colonne d'identifiant échantillon/patient parmi les noms usuels."""
    # 1) Correspondance exacte avec la liste de candidats, dans l'ordre de priorité donné
    for c in candidates:
        if c in df.columns:
            return c
    # 2) Fallback : colonne contenant "sample" dans son nom (priorité sur "id" seul, trop générique)
    for c in df.columns:
        if "sample" in c.lower():
            return c
    # 3) Dernier recours : colonne se terminant par "id"
    for c in df.columns:
        if c.lower().endswith("id"):
            return c
    raise ValueError(f"Impossible de trouver une colonne d'identifiant parmi : {list(df.columns)}")


# ----------------------------------------------------------------------------
# 2.2 — Nettoyage
# ----------------------------------------------------------------------------

def barcode_to_patient_id(barcode: str) -> str:
    """Convertit un barcode échantillon (ex. TCGA-D8-A146-01A) en identifiant patient (TCGA-D8-A146)."""
    return "-".join(barcode.split("-")[:3])


def is_primary_tumor(barcode: str) -> bool:
    """
    Le code d'échantillon TCGA (2 chiffres après le 3e tiret) indique le type de tissu :
    01 = tumeur primaire, 11 = tissu normal, etc.
    On ne garde que les tumeurs primaires pour éviter les doublons normal/tumeur du même patient.
    """
    parts = barcode.split("-")
    if len(parts) < 4:
        return False
    sample_code = parts[3][:2]
    return sample_code == "01"


def merge_and_clean(expr, clinical, survival, subtype) -> pd.DataFrame:
    print("\n--- Fusion des jeux de données ---")

    # Filtrer l'expression pour ne garder que les échantillons de tumeur primaire
    expr = expr[[is_primary_tumor(b) for b in expr.index]]
    expr["patient_id"] = [barcode_to_patient_id(b) for b in expr.index]
    print(f"Après filtrage tumeur primaire : {expr.shape[0]} échantillons")

    # Identifier les colonnes d'ID dans chaque table clinique
    clinical_id_col = find_sample_id_column(clinical)
    survival_id_col = find_sample_id_column(survival)
    subtype_id_col = find_sample_id_column(
        subtype, candidates=("Patient ID", "PATIENT_ID", "_PATIENT", "bcr_patient_barcode", "sample")
    )

    clinical = clinical.rename(columns={clinical_id_col: "sample_barcode"})
    survival = survival.rename(columns={survival_id_col: "sample_barcode"})
    subtype = subtype.rename(columns={subtype_id_col: "patient_id"})

    # Harmonisation du nom de la colonne cible selon la source (export cBioPortal complet -> 'Subtype')
    for possible_name in ("Subtype", "PAM50"):
        if possible_name in subtype.columns and "SUBTYPE" not in subtype.columns:
            subtype = subtype.rename(columns={possible_name: "SUBTYPE"})

    if "SUBTYPE" not in subtype.columns:
        raise ValueError(
            f"Colonne 'SUBTYPE'/'Subtype'/'PAM50' introuvable dans le fichier PAM50. "
            f"Colonnes disponibles : {list(subtype.columns)}"
        )
    subtype = subtype[["patient_id", "SUBTYPE"]].dropna(subset=["SUBTYPE"])
    subtype = subtype.drop_duplicates(subset=["patient_id"])

    # Fusion expression <-> clinique (sur sample_barcode)
    merged = expr.reset_index().merge(clinical, on="sample_barcode", how="inner")
    print(f"Après fusion expression + clinique : {merged.shape[0]} échantillons")

    # Fusion avec survie (sur sample_barcode)
    if "sample_barcode" in survival.columns:
        merged = merged.merge(survival, on="sample_barcode", how="left")
        print(f"Après fusion avec survie : {merged.shape[0]} échantillons")

    # Fusion avec PAM50/SUBTYPE (sur patient_id — cet export cBioPortal est au niveau patient)
    merged = merged.merge(subtype, on="patient_id", how="inner")
    print(f"Après fusion avec SUBTYPE (PAM50) : {merged.shape[0]} échantillons")

    # Supprimer les échantillons sans label PAM50 (sécurité, déjà fait par le inner join ci-dessus)
    merged = merged.dropna(subset=["SUBTYPE"])

    # Supprimer les doublons patients éventuels (garder le premier échantillon rencontré)
    n_before = merged.shape[0]
    merged = merged.drop_duplicates(subset=["patient_id"], keep="first")
    print(f"Doublons patients supprimés : {n_before - merged.shape[0]}")

    return merged


def handle_missing_clinical(df: pd.DataFrame, clinical_cols: list) -> pd.DataFrame:
    """Impute les valeurs manquantes des variables cliniques : médiane si numérique, 'Unknown' sinon."""
    print("\n--- Gestion des valeurs manquantes cliniques ---")
    for col in clinical_cols:
        n_missing = df[col].isna().sum()
        if n_missing == 0:
            continue
        pct_missing = n_missing / len(df) * 100
        if pct_missing > 40:
            print(f"  '{col}' : {pct_missing:.1f}% manquant -> colonne exclue (trop de manquants)")
            df = df.drop(columns=[col])
            continue
        if pd.api.types.is_numeric_dtype(df[col]):
            median_val = df[col].median()
            df[col] = df[col].fillna(median_val)
            print(f"  '{col}' : {n_missing} valeurs imputées par la médiane ({median_val:.2f})")
        else:
            df[col] = df[col].fillna("Unknown")
            print(f"  '{col}' : {n_missing} valeurs remplacées par 'Unknown'")
    return df


# ----------------------------------------------------------------------------
# 2.3 — Filtrage des gènes
# ----------------------------------------------------------------------------

def filter_genes(expr_only: pd.DataFrame, variance_threshold=0.01, mean_expr_threshold=1.0) -> pd.DataFrame:
    """Supprime les gènes à variance quasi nulle et/ou à faible expression moyenne."""
    print("\n--- Filtrage des gènes ---")
    n_before = expr_only.shape[1]

    variances = expr_only.var(axis=0)
    means = expr_only.mean(axis=0)

    keep_mask = (variances > variance_threshold) & (means > mean_expr_threshold)
    expr_filtered = expr_only.loc[:, keep_mask]

    print(f"  Gènes avant filtrage : {n_before}")
    print(f"  Gènes après filtrage : {expr_filtered.shape[1]}")
    return expr_filtered


# ----------------------------------------------------------------------------
# Pipeline principal
# ----------------------------------------------------------------------------

def main():
    expr = load_expression(EXPRESSION_FILE)
    gene_mapping = load_probemap(PROBEMAP_FILE)
    expr = expr.rename(columns=gene_mapping)  # ID Ensembl -> nom de gène lisible

    clinical = load_clinical(CLINICAL_FILE)
    survival = load_survival(SURVIVAL_FILE)
    subtype = load_subtype(SUBTYPE_FILE)

    merged = merge_and_clean(expr, clinical, survival, subtype)

    # Séparer les colonnes d'expression (gènes) des colonnes cliniques/méta
    meta_cols = ["sample_barcode", "patient_id", "SUBTYPE"]
    known_clinical_cols = [c for c in clinical.columns if c != "sample_barcode" and c in merged.columns]
    known_survival_cols = [c for c in survival.columns if c != "sample_barcode" and c in merged.columns]
    clinical_cols = [c for c in known_clinical_cols + known_survival_cols if c not in meta_cols]

    gene_cols = [c for c in merged.columns if c not in meta_cols + clinical_cols]

    # Nettoyage des valeurs manquantes cliniques uniquement (l'expression ne devrait pas avoir de NaN)
    merged = handle_missing_clinical(merged, clinical_cols)
    # handle_missing_clinical peut avoir supprimé des colonnes (>40% manquant) : on met à jour la liste
    clinical_cols = [c for c in clinical_cols if c in merged.columns]

    # Filtrage des gènes
    expr_filtered = filter_genes(merged[gene_cols])

    # Reconstruction du tableau final
    final_df = pd.concat(
        [merged[["patient_id", "SUBTYPE"]], merged[clinical_cols], expr_filtered],
        axis=1,
    )
    final_df = final_df.rename(columns={"SUBTYPE": "subtype"})

    output_path = PROCESSED_DIR / "tcga_train.csv"
    final_df.to_csv(output_path, index=False)

    print(f"\n=== Terminé ===")
    print(f"Fichier final : {output_path}")
    print(f"Forme finale : {final_df.shape[0]} patients x {final_df.shape[1]} colonnes")
    print(f"Répartition des sous-types :\n{final_df['subtype'].value_counts()}")


if __name__ == "__main__":
    main()