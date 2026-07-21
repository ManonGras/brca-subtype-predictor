"""
Phase 5 — Explicabilité biologique (SHAP)
============================================
Entrée : models/subtype_classifier.pkl (Phase 4), data/processed/train_reduced.csv
Sorties :
  - report/figures/shap_summary_global.png
  - report/figures/shap_summary_<classe>.png (un par sous-type)
  - report/interpretation_biologique.md

À exécuter depuis la racine du projet :
    python src/explainability.py
"""

import pandas as pd
import numpy as np
import joblib
import shap
import matplotlib.pyplot as plt
from pathlib import Path

PROCESSED_DIR = Path("data/processed")
MODELS_DIR = Path("models")
FIGURES_DIR = Path("report/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

N_TOP_GENES = 20
N_SAMPLES_FOR_SHAP = 300  # sous-échantillonnage pour accélérer le calcul SHAP si le train set est grand

# Liste officielle de la signature PAM50 (Parker et al., 2009, J Clin Oncol)
# Utilisée uniquement comme référence de comparaison, pas comme feature du modèle.
PAM50_GENES = {
    "ACTR3B", "ANLN", "BAG1", "BCL2", "BIRC5", "BLVRA", "CCNB1", "CCNE1", "CDC20", "CDC6",
    "CDH3", "CENPF", "CEP55", "CDCA1", "NUF2", "CXXC5", "EGFR", "ERBB2", "ESR1", "EXO1",
    "FGFR4", "FOXA1", "FOXC1", "GPR160", "GRB7", "KIF2C", "KRT14", "KRT17", "KRT5", "MAPT",
    "MDM2", "MELK", "MIA", "MKI67", "MLPH", "MMP11", "MYBL2", "MYC", "NAT1", "NDC80",
    "ORC6", "ORC6L", "PGR", "PHGDH", "PTTG1", "RRM2", "SFRP1", "SLC39A6", "TMEM45B", "TYMS",
    "UBE2C", "UBE2T",
}


def load_artifacts():
    bundle = joblib.load(MODELS_DIR / "subtype_classifier.pkl")
    model = bundle["model"]
    scaler = bundle["scaler"]
    label_encoder = bundle["label_encoder"]
    gene_cols = bundle["gene_cols"]
    print(f"Modèle chargé : {type(model).__name__} | {len(gene_cols)} gènes")
    return model, scaler, label_encoder, gene_cols


def compute_shap_values(model, scaler, gene_cols):
    train = pd.read_csv(PROCESSED_DIR / "train_reduced.csv")
    X = train[gene_cols]

    if len(X) > N_SAMPLES_FOR_SHAP:
        X_sample = X.sample(N_SAMPLES_FOR_SHAP, random_state=42)
    else:
        X_sample = X

    X_scaled = scaler.transform(X_sample)
    X_scaled_df = pd.DataFrame(X_scaled, columns=gene_cols)

    print(f"Calcul des valeurs SHAP sur {X_scaled_df.shape[0]} patients ...")
    explainer = shap.TreeExplainer(model)
    shap_values = explainer(X_scaled_df)
    print("  -> terminé.")
    return shap_values, X_scaled_df


def summarize_and_plot(shap_values, X_scaled_df, label_encoder):
    """
    Gère le format multi-classe de SHAP (values de forme n_samples x n_features x n_classes),
    produit un graphique global (importance moyenne toutes classes confondues) et un graphique
    par sous-type.
    """
    values = shap_values.values  # forme attendue : (n_samples, n_features, n_classes)
    n_classes = values.shape[-1] if values.ndim == 3 else 1
    class_names = label_encoder.classes_

    # Importance globale = moyenne des |SHAP| sur tous les échantillons ET toutes les classes
    if values.ndim == 3:
        global_importance = np.abs(values).mean(axis=(0, 2))
    else:
        global_importance = np.abs(values).mean(axis=0)

    importance_df = pd.DataFrame({
        "gene": X_scaled_df.columns,
        "mean_abs_shap": global_importance,
    }).sort_values("mean_abs_shap", ascending=False)

    # Graphique global (bar plot top N)
    top_global = importance_df.head(N_TOP_GENES)
    plt.figure(figsize=(8, 6))
    plt.barh(top_global["gene"][::-1], top_global["mean_abs_shap"][::-1], color="steelblue")
    plt.xlabel("Importance SHAP moyenne (|valeur SHAP|)")
    plt.title(f"Top {N_TOP_GENES} gènes les plus influents (toutes classes confondues)")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "shap_summary_global.png", dpi=150)
    plt.close()
    print(f"Graphique global sauvegardé : {FIGURES_DIR / 'shap_summary_global.png'}")

    # Graphiques par classe (summary plot SHAP natif)
    if values.ndim == 3:
        for i, class_name in enumerate(class_names):
            plt.figure()
            shap.summary_plot(
                values[:, :, i], X_scaled_df, show=False, max_display=15, plot_size=(8, 6)
            )
            plt.title(f"SHAP - sous-type {class_name}")
            plt.tight_layout()
            fname = FIGURES_DIR / f"shap_summary_{class_name}.png"
            plt.savefig(fname, dpi=150)
            plt.close()
            print(f"  Graphique par classe sauvegardé : {fname}")

    return importance_df


def cross_reference_pam50(importance_df: pd.DataFrame, n_top=N_TOP_GENES):
    top_genes = set(importance_df.head(n_top)["gene"])
    overlap = top_genes & PAM50_GENES
    print(f"\n--- Confrontation avec la signature PAM50 officielle ---")
    print(f"Parmi les {n_top} gènes les plus importants selon SHAP : {len(overlap)} appartiennent "
          f"à la signature PAM50 originale (Parker et al. 2009).")
    print(f"Gènes en commun : {sorted(overlap)}")
    return overlap


def write_report(importance_df: pd.DataFrame, overlap: set, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    top = importance_df.head(N_TOP_GENES).copy()
    top["dans_PAM50_officiel"] = top["gene"].apply(lambda g: "Oui" if g in PAM50_GENES else "Non")

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("# Interprétation biologique - Phase 5 (SHAP)\n\n")
        f.write(f"## Top {N_TOP_GENES} gènes les plus influents selon SHAP\n\n")
        f.write(top.to_markdown(index=False))
        f.write("\n\n")
        f.write(f"## Confrontation avec la signature PAM50 officielle\n\n")
        f.write(f"{len(overlap)} des {N_TOP_GENES} gènes les plus importants du modèle appartiennent "
                f"à la signature PAM50 originale (Parker et al., 2009) : {sorted(overlap)}.\n\n")
        f.write(
            "Ce recoupement, même partiel, constitue un indice que le modèle s'appuie sur des gènes "
            "biologiquement pertinents et déjà validés dans la littérature, plutôt que sur des "
            "corrélations arbitraires du jeu de données. Les gènes importants qui n'appartiennent "
            "pas à la liste PAM50 méritent une discussion prudente : ils peuvent refléter soit une "
            "information complémentaire réelle, soit un artefact du jeu de données à signaler comme "
            "limite.\n\n"
        )
        f.write(
            "**À compléter manuellement** : pour chaque gène important hors PAM50, rechercher sa "
            "fonction connue (GeneCards, Human Protein Atlas) et rédiger 2-3 phrases de contexte "
            "biologique, comme prévu dans le guide de réalisation.\n"
        )
    print(f"\nRapport sauvegardé : {output_path}")


def main():
    model, scaler, label_encoder, gene_cols = load_artifacts()
    shap_values, X_scaled_df = compute_shap_values(model, scaler, gene_cols)
    importance_df = summarize_and_plot(shap_values, X_scaled_df, label_encoder)
    overlap = cross_reference_pam50(importance_df)
    write_report(importance_df, overlap, Path("report/interpretation_biologique.md"))
    print("\n=== Phase 5 terminée ===")


if __name__ == "__main__":
    main()