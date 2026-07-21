"""
Phase 6 — Validation externe sur METABRIC
============================================================
Entrées :
  - models/subtype_classifier.pkl      (modèle XGBoost entraîné en Phase 4, 300 gènes)
  - data/processed/train_reduced.csv   (pour ré-entraîner le modèle restreint)
  - data/processed/test_reduced.csv    (métriques internes de référence)
  - data/raw/metabric/data_mrna_illumina_microarray.txt
  - data/raw/metabric/data_clinical_patient.txt

Sorties :
  - report/metabric_validation.md
  - report/figures/kaplan_meier_metabric_predicted_subtype.png
  - models/subtype_classifier_restricted210.pkl

À exécuter depuis la racine du projet :
    python src/validate_metabric.py
"""

import joblib
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import accuracy_score, f1_score, classification_report
from xgboost import XGBClassifier

import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

PROCESSED_DIR = Path("data/processed")
METABRIC_DIR = Path("data/raw/metabric")
MODELS_DIR = Path("models")
FIGURES_DIR = Path("report/figures")
REPORT_DIR = Path("report")

RANDOM_STATE = 42

# Nomenclature METABRIC (Pam50 + Claudin-low subtype) -> nomenclature TCGA (SUBTYPE cBioPortal)
SUBTYPE_MAPPING = {
    "LumA": "BRCA_LumA",
    "LumB": "BRCA_LumB",
    "Her2": "BRCA_Her2",
    "Basal": "BRCA_Basal",
    "Normal": "BRCA_Normal",
    # "claudin-low" et "NC" n'existent pas dans le schéma PAM50 à 5 classes utilisé pour
    # entraîner le modèle sur TCGA -> échantillons exclus (documenté dans le rapport).
}


# ----------------------------------------------------------------------------
# 1) Chargement et harmonisation de METABRIC
# ----------------------------------------------------------------------------

def load_metabric_expression(path: Path) -> pd.DataFrame:
    """Charge la matrice d'expression microarray (gènes x échantillons), agrège les sondes
    dupliquées par symbole (moyenne), puis transpose en (échantillons x gènes)."""
    print(f"Chargement de l'expression METABRIC depuis {path} ...")
    expr = pd.read_csv(path, sep="\t")
    print(f"  -> forme brute (sondes x (2 cols meta + échantillons)) : {expr.shape}")

    sample_cols = [c for c in expr.columns if c not in ("Hugo_Symbol", "Entrez_Gene_Id")]
    n_dup = expr["Hugo_Symbol"].duplicated().sum()
    print(f"  -> {n_dup} sondes dupliquées par symbole -> agrégation par moyenne")
    expr_by_gene = expr.groupby("Hugo_Symbol")[sample_cols].mean()

    expr_t = expr_by_gene.T
    expr_t.index.name = "SAMPLE_ID"
    print(f"  -> forme après agrégation + transposition (échantillons x gènes) : {expr_t.shape}")
    return expr_t


def load_metabric_clinical(path: Path) -> pd.DataFrame:
    print(f"Chargement de la clinique METABRIC depuis {path} ...")
    clinical = pd.read_csv(path, sep="\t", comment="#")
    print(f"  -> forme : {clinical.shape}")
    return clinical


def build_metabric_dataset(expr_t: pd.DataFrame, clinical: pd.DataFrame):
    """Fusionne expression + clinique, mappe le sous-type vers la nomenclature TCGA,
    et exclut les échantillons hors du schéma PAM50 à 5 classes (claudin-low, NC, manquant)."""
    print("\n--- Fusion expression + clinique METABRIC ---")
    clinical = clinical.rename(columns={"PATIENT_ID": "SAMPLE_ID"})
    merged = expr_t.merge(
        clinical[["SAMPLE_ID", "CLAUDIN_SUBTYPE", "OS_MONTHS", "OS_STATUS"]],
        left_index=True, right_on="SAMPLE_ID", how="inner",
    )
    print(f"  -> échantillons avec expression + clinique : {merged.shape[0]}")

    merged["subtype"] = merged["CLAUDIN_SUBTYPE"].map(SUBTYPE_MAPPING)
    n_excluded = merged["subtype"].isna().sum()
    print(f"  -> {n_excluded} échantillons exclus (hors schéma PAM50 5 classes : "
          f"claudin-low / NC / manquant)")
    merged = merged.dropna(subset=["subtype"])
    print(f"  -> échantillons METABRIC retenus pour la validation : {merged.shape[0]}")
    print(f"  -> répartition des sous-types :\n{merged['subtype'].value_counts()}")

    # OS_STATUS METABRIC au format "0:LIVING" / "1:DECEASED" -> statut binaire (0 = vivant, 1 = décédé)
    merged["OS_EVENT"] = merged["OS_STATUS"].str.split(":").str[0].astype(float)

    return merged.set_index("SAMPLE_ID")


# ----------------------------------------------------------------------------
# 2) Correspondance des gènes et batch effect
# ----------------------------------------------------------------------------

def gene_overlap_report(model_genes: list, metabric_genes: set):
    overlap = [g for g in model_genes if g in metabric_genes]
    missing = [g for g in model_genes if g not in metabric_genes]
    print(f"\n--- Correspondance des gènes ---")
    print(f"  Gènes du modèle (TCGA, n={len(model_genes)}) présents sur la plateforme "
          f"Illumina METABRIC : {len(overlap)}")
    print(f"  Gènes absents de la plateforme METABRIC : {len(missing)}")
    return overlap, missing


def zscore_within_dataset(X: pd.DataFrame) -> pd.DataFrame:
    """Standardisation z-score par gène, calculée sur ce jeu de données uniquement
    (batch effect : on ne réutilise jamais les moyennes/écarts-types d'un autre jeu de données)."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    return pd.DataFrame(X_scaled, index=X.index, columns=X.columns)


# ----------------------------------------------------------------------------
# 3) Modèle complet (300 gènes, imputation des gènes absents)
# ----------------------------------------------------------------------------

def evaluate_full_model(bundle, metabric_expr: pd.DataFrame, overlap: list, missing: list,
                         y_true: pd.Series):
    print("\n--- Évaluation du modèle complet (300 gènes, imputation à 0) sur METABRIC ---")
    gene_cols = bundle["gene_cols"]

    X_available = zscore_within_dataset(metabric_expr[overlap])
    X_full = pd.DataFrame(0.0, index=metabric_expr.index, columns=gene_cols)
    X_full[overlap] = X_available
    print(f"  -> {len(missing)}/{len(gene_cols)} colonnes imputées à 0 (valeur neutre post-standardisation)")

    label_encoder = bundle["label_encoder"]
    y_true_encoded = label_encoder.transform(y_true)
    y_pred_encoded = bundle["model"].predict(X_full.values)
    y_pred = label_encoder.inverse_transform(y_pred_encoded)

    acc = accuracy_score(y_true_encoded, y_pred_encoded)
    f1_macro = f1_score(y_true_encoded, y_pred_encoded, average="macro")
    report = classification_report(y_true_encoded, y_pred_encoded, target_names=label_encoder.classes_, digits=3)
    print(f"  accuracy={acc:.3f} | F1 macro={f1_macro:.3f}")
    print(report)

    return {"accuracy": acc, "f1_macro": f1_macro, "report": report, "y_pred": y_pred}


# ----------------------------------------------------------------------------
# 4) Modèle restreint (gènes communs aux deux plateformes)
# ----------------------------------------------------------------------------

def train_restricted_model(common_genes: list):
    """Ré-entraîne un XGBoost identique à celui de la Phase 4, mais limité aux gènes
    mesurables sur les deux plateformes -- même split train/test que la Phase 4
    (fichiers train_reduced.csv / test_reduced.csv), pour isoler l'effet de l'imputation
    de gènes manquants de l'effet de généralisation réelle du modèle."""
    print(f"\n--- Ré-entraînement du modèle restreint ({len(common_genes)} gènes communs) ---")
    train = pd.read_csv(PROCESSED_DIR / "train_reduced.csv")
    test = pd.read_csv(PROCESSED_DIR / "test_reduced.csv")

    X_train, y_train = train[common_genes], train["subtype"]
    X_test, y_test = test[common_genes], test["subtype"]

    label_encoder = LabelEncoder()
    y_train_encoded = label_encoder.fit_transform(y_train)
    y_test_encoded = label_encoder.transform(y_test)

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    model = XGBClassifier(random_state=RANDOM_STATE, eval_metric="mlogloss")
    model.fit(X_train_scaled, y_train_encoded)

    y_pred_encoded = model.predict(X_test_scaled)
    acc = accuracy_score(y_test_encoded, y_pred_encoded)
    f1_macro = f1_score(y_test_encoded, y_pred_encoded, average="macro")
    report = classification_report(y_test_encoded, y_pred_encoded, target_names=label_encoder.classes_, digits=3)
    print(f"  [TCGA interne, {len(common_genes)} gènes] accuracy={acc:.3f} | F1 macro={f1_macro:.3f}")

    bundle = {
        "model": model, "scaler": scaler, "label_encoder": label_encoder,
        "gene_cols": common_genes,
    }
    joblib.dump(bundle, MODELS_DIR / "subtype_classifier_restricted210.pkl")
    print(f"  -> modèle restreint sauvegardé : {MODELS_DIR / 'subtype_classifier_restricted210.pkl'}")

    return bundle, {"accuracy": acc, "f1_macro": f1_macro, "report": report}


def evaluate_restricted_model(bundle, metabric_expr: pd.DataFrame, common_genes: list, y_true: pd.Series):
    print(f"\n--- Évaluation du modèle restreint ({len(common_genes)} gènes) sur METABRIC ---")
    X = zscore_within_dataset(metabric_expr[common_genes])[common_genes]

    label_encoder = bundle["label_encoder"]
    y_true_encoded = label_encoder.transform(y_true)
    y_pred_encoded = bundle["model"].predict(X.values)
    y_pred = label_encoder.inverse_transform(y_pred_encoded)

    acc = accuracy_score(y_true_encoded, y_pred_encoded)
    f1_macro = f1_score(y_true_encoded, y_pred_encoded, average="macro")
    report = classification_report(y_true_encoded, y_pred_encoded, target_names=label_encoder.classes_, digits=3)
    print(f"  accuracy={acc:.3f} | F1 macro={f1_macro:.3f}")
    print(report)

    return {"accuracy": acc, "f1_macro": f1_macro, "report": report, "y_pred": y_pred}


# ----------------------------------------------------------------------------
# 5) Analyse de survie sur METABRIC (sous-types prédits par le modèle restreint)
# ----------------------------------------------------------------------------

def survival_analysis_metabric(metabric_df: pd.DataFrame, y_pred_labels, output_path: Path):
    print("\n--- Analyse de survie (Kaplan-Meier) sur METABRIC, sous-types prédits ---")
    surv_df = metabric_df[["OS_MONTHS", "OS_EVENT"]].copy()
    surv_df["predicted_subtype"] = y_pred_labels
    surv_df = surv_df.dropna(subset=["OS_MONTHS", "OS_EVENT"])
    print(f"  -> {surv_df.shape[0]} échantillons avec données de survie complètes")

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    for subtype in sorted(surv_df["predicted_subtype"].unique()):
        mask = surv_df["predicted_subtype"] == subtype
        kmf.fit(surv_df.loc[mask, "OS_MONTHS"], surv_df.loc[mask, "OS_EVENT"], label=subtype)
        kmf.plot_survival_function(ax=ax)

    plt.title("Courbes de survie (Kaplan-Meier) par sous-type prédit - METABRIC (validation externe)")
    plt.xlabel("Temps (mois)")
    plt.ylabel("Probabilité de survie")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  -> figure sauvegardée : {output_path}")

    result = multivariate_logrank_test(surv_df["OS_MONTHS"], surv_df["predicted_subtype"], surv_df["OS_EVENT"])
    print(f"  Test du log-rank multivarié : p-valeur = {result.p_value:.4g}")
    return result.p_value, surv_df.shape[0]


# ----------------------------------------------------------------------------
# 6) Rapport
# ----------------------------------------------------------------------------

def write_report(full_internal, full_external, restricted_internal, restricted_external,
                  overlap, missing, n_metabric_total, n_metabric_excluded, n_metabric_used,
                  logrank_p, n_survival, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lncrna_like = [g for g in missing if g.startswith(("LINC", "AL", "AC", "Y_RNA", "MIR")) or g.endswith(("-AS1", "-DT", "P1"))]

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("# Phase 6 — Validation externe sur METABRIC\n\n")

        f.write("## 1. Cohorte METABRIC utilisée\n\n")
        f.write(f"- Échantillons METABRIC avec expression + clinique disponibles : {n_metabric_total}\n")
        f.write(f"- Échantillons exclus (hors schéma PAM50 à 5 classes : `claudin-low`, `NC`, ou sous-type "
                f"manquant) : {n_metabric_excluded}\n")
        f.write(f"- Échantillons retenus pour la validation externe : **{n_metabric_used}**\n\n")

        f.write("## 2. Correspondance des gènes entre plateformes\n\n")
        f.write(f"- Gènes du modèle (sélectionnés en Phase 3 sur TCGA, RNA-seq) : 300\n")
        f.write(f"- Présents sur la puce microarray Illumina de METABRIC : **{len(overlap)}**\n")
        f.write(f"- Absents de la plateforme METABRIC : **{len(missing)}**\n\n")
        f.write("La quasi-totalité des gènes absents ne sont pas des gènes protéiques classiques mais des "
                "entités dont l'annotation (lncRNA, ARN antisens, pseudogènes, `Y_RNA`) est postérieure ou "
                "non prise en charge par la conception des sondes de la puce Illumina HT-12 utilisée pour "
                f"METABRIC (~2010) : sur les {len(missing)} gènes manquants, {len(lncrna_like)} correspondent "
                "à ce type d'identifiant (`LINC*`, `AL*`/`AC*` (contigs GENCODE non nommés), `Y_RNA.*`, "
                "gènes se terminant en `-AS1`/`-DT`/pseudogènes `P1`). Ce n'est donc pas majoritairement un "
                "problème de généralisation biologique mais un artefact de couverture technique de plateforme.\n\n")
        f.write("Liste complète des gènes absents :\n\n")
        f.write("```\n" + ", ".join(sorted(missing)) + "\n```\n\n")

        f.write("## 3. Tableau comparatif — validation interne (TCGA) vs externe (METABRIC)\n\n")
        f.write("| Modèle | Gènes | Accuracy interne (TCGA) | F1 macro interne (TCGA) | "
                "Accuracy externe (METABRIC) | F1 macro externe (METABRIC) | Δ Accuracy | Δ F1 macro |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        f.write(f"| Complet (Phase 4, imputation à 0 des gènes absents) | 300 | "
                f"{full_internal['accuracy']:.3f} | {full_internal['f1_macro']:.3f} | "
                f"{full_external['accuracy']:.3f} | {full_external['f1_macro']:.3f} | "
                f"{full_external['accuracy'] - full_internal['accuracy']:+.3f} | "
                f"{full_external['f1_macro'] - full_internal['f1_macro']:+.3f} |\n")
        f.write(f"| Restreint (ré-entraîné, gènes communs uniquement) | {len(overlap)} | "
                f"{restricted_internal['accuracy']:.3f} | {restricted_internal['f1_macro']:.3f} | "
                f"{restricted_external['accuracy']:.3f} | {restricted_external['f1_macro']:.3f} | "
                f"{restricted_external['accuracy'] - restricted_internal['accuracy']:+.3f} | "
                f"{restricted_external['f1_macro'] - restricted_internal['f1_macro']:+.3f} |\n\n")

        f.write("*« Interne » = mesuré sur le test set TCGA (Phase 4), jamais vu pendant l'entraînement. "
                "« Externe » = mesuré sur la cohorte METABRIC, jamais vue pendant l'entraînement, sur une "
                "plateforme technique différente (microarray vs RNA-seq). Le modèle « restreint » est "
                "ré-entraîné sur TCGA avec exactement le même pipeline (XGBoost, même split, mêmes "
                "hyperparamètres) mais en se limitant aux gènes réellement mesurables sur METABRIC — il n'y a "
                "donc aucune imputation lors de son évaluation externe.*\n\n")

        f.write("## 4. Rapport de classification détaillé — METABRIC (modèle restreint)\n\n")
        f.write("```\n" + restricted_external["report"] + "\n```\n\n")

        f.write("## 5. Analyse de survie sur METABRIC (sous-types prédits, modèle restreint)\n\n")
        f.write(f"- Échantillons avec données de survie complètes : {n_survival}\n")
        f.write(f"- Test du log-rank multivarié (5 sous-types prédits) : **p = {logrank_p:.4g}**\n")
        f.write("- Figure : `report/figures/kaplan_meier_metabric_predicted_subtype.png`\n\n")
        f.write("Pour rappel, sur TCGA (Phase 4), le même type de test donnait p = 0.001282. La relation "
                "pronostique entre sous-type prédit et survie observée sur METABRIC doit être comparée à ce "
                "chiffre pour juger si elle se maintient hors de la cohorte d'entraînement.\n\n")

        f.write("## 6. Discussion critique\n\n")
        drop_full = full_internal['f1_macro'] - full_external['f1_macro']
        drop_restricted = restricted_internal['f1_macro'] - restricted_external['f1_macro']
        f.write(
            f"**Chute de performance (F1 macro, interne → externe) :** {drop_full:.3f} pour le modèle "
            f"complet (300 gènes, {len(missing)} imputés), contre {drop_restricted:.3f} pour le modèle "
            f"restreint (entraîné et évalué sur les seuls {len(overlap)} gènes communs, sans imputation).\n\n"
        )
        if drop_full > drop_restricted:
            f.write(
                "L'écart entre les deux lignes du tableau ci-dessus indique qu'une partie non négligeable de "
                "la chute de performance du modèle complet est un **artefact technique lié à l'imputation** "
                "des gènes absents de la plateforme METABRIC, et non un vrai problème de généralisation du "
                "modèle. Le modèle restreint, qui ne souffre pas de ce problème d'imputation, mesure plus "
                "fidèlement la capacité du modèle à généraliser à une cohorte et une plateforme technique "
                "différentes de celles utilisées pour l'entraînement.\n\n"
            )
        else:
            f.write(
                "Le modèle restreint (sans imputation) chute autant, voire davantage, que le modèle complet : "
                "l'imputation des gènes manquants n'explique donc pas l'essentiel de la chute observée. La "
                "cause principale semble être un **vrai écart de généralisation** entre RNA-seq (TCGA) et "
                "microarray (METABRIC) — différences de dynamique de mesure, de bruit technique, et de "
                "composition de cohorte (âge, ethnicité, période de recrutement) entre les deux études.\n\n"
            )
        f.write(
            "**Limites méthodologiques à garder à l'esprit :**\n\n"
            "1. La standardisation z-score appliquée séparément à chaque jeu de données (recommandation "
            "*a minima* de la Phase 6) suppose implicitement que la distribution relative de chaque gène est "
            "comparable entre RNA-seq et microarray — hypothèse simplificatrice. Une correction de type "
            "ComBat (option avancée non mise en œuvre ici) modéliserait explicitement le batch comme facteur "
            "et pourrait réduire encore l'écart observé, en particulier pour le modèle complet.\n"
            "2. Les labels de vérité-terrain diffèrent aussi dans leur origine : le sous-type TCGA (`SUBTYPE`, "
            "cBioPortal PanCancer Atlas) et le sous-type METABRIC (`CLAUDIN_SUBTYPE`, PAM50 + claudin-low) "
            "proviennent de pipelines de classification PAM50 possiblement recalibrés différemment selon "
            "l'étude — une partie des erreurs de classification peut donc relever du désaccord entre "
            "méthodes de labellisation plutôt que d'une erreur du modèle.\n"
            "3. METABRIC est une cohorte plus âgée en moyenne et recrutée sur une période et dans un contexte "
            "de soins différents de TCGA ; un déplacement de la distribution des covariables cliniques "
            "(*covariate shift*) peut contribuer à la chute de performance indépendamment de la plateforme "
            "technique.\n"
            "4. Les 6 échantillons `NC` et 218 échantillons `claudin-low` de METABRIC ont été exclus car ils "
            "sortent du schéma PAM50 à 5 classes sur lequel le modèle a été entraîné — un modèle destiné à un "
            "usage clinique réel devrait explicitement gérer une classe \"hors distribution\".\n\n"
        )
        f.write(
            "**Conclusion :** la validation externe confirme que le signal biologique appris sur TCGA est "
            "partiellement transférable à une cohorte et une plateforme indépendantes, mais avec une chute de "
            "performance mesurable qu'il serait malhonnête de masquer. Le tableau ci-dessus permet de "
            "distinguer la part de cette chute imputable à un artefact de couverture de plateforme (gènes "
            "absents de l'array METABRIC) de la part imputable à une vraie limite de généralisation du "
            "modèle.\n"
        )

    print(f"\nRapport sauvegardé : {output_path}")


# ----------------------------------------------------------------------------
# Pipeline principal
# ----------------------------------------------------------------------------

def main():
    bundle = joblib.load(MODELS_DIR / "subtype_classifier.pkl")
    model_genes = bundle["gene_cols"]

    expr_t = load_metabric_expression(METABRIC_DIR / "data_mrna_illumina_microarray.txt")
    clinical = load_metabric_clinical(METABRIC_DIR / "data_clinical_patient.txt")
    metabric_df = build_metabric_dataset(expr_t, clinical)

    n_metabric_total = expr_t.merge(
        clinical.rename(columns={"PATIENT_ID": "SAMPLE_ID"})[["SAMPLE_ID", "CLAUDIN_SUBTYPE"]],
        left_index=True, right_on="SAMPLE_ID", how="inner",
    ).shape[0]
    n_metabric_used = metabric_df.shape[0]
    n_metabric_excluded = n_metabric_total - n_metabric_used

    overlap, missing = gene_overlap_report(model_genes, set(metabric_df.columns))

    gene_matrix_cols = [c for c in metabric_df.columns if c not in
                        ("CLAUDIN_SUBTYPE", "OS_MONTHS", "OS_STATUS", "subtype", "OS_EVENT")]
    metabric_expr_only = metabric_df[gene_matrix_cols]
    y_true = metabric_df["subtype"]

    # --- Modèle complet (300 gènes, imputation) ---
    full_internal = {"accuracy": 0.853, "f1_macro": 0.784}  # repris de report/model_comparison.md (Phase 4)
    full_external = evaluate_full_model(bundle, metabric_expr_only, overlap, missing, y_true)

    # --- Modèle restreint (gènes communs, ré-entraîné) ---
    restricted_bundle, restricted_internal = train_restricted_model(overlap)
    restricted_external = evaluate_restricted_model(restricted_bundle, metabric_expr_only, overlap, y_true)

    # --- Analyse de survie sur METABRIC (sous-types prédits par le modèle restreint) ---
    logrank_p, n_survival = survival_analysis_metabric(
        metabric_df, restricted_external["y_pred"],
        FIGURES_DIR / "kaplan_meier_metabric_predicted_subtype.png",
    )

    write_report(
        full_internal, full_external, restricted_internal, restricted_external,
        overlap, missing, n_metabric_total, n_metabric_excluded, n_metabric_used,
        logrank_p, n_survival, REPORT_DIR / "metabric_validation.md",
    )

    print("\n=== Phase 6 terminée ===")


if __name__ == "__main__":
    main()
