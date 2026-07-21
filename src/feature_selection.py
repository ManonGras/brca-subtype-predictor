"""
Phase 3 — Sélection de features et réduction de dimension
============================================================
Entrée : data/processed/tcga_train.csv (sortie de la Phase 2)
Sorties :
  - data/processed/train_reduced.csv
  - data/processed/test_reduced.csv
  - report/gene_selection_justification.md
  - report/figures/pca_subtypes.png
  - report/figures/umap_subtypes.png

À exécuter depuis la racine du projet :
    python src/feature_selection.py
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import umap

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
PROCESSED_DIR = Path("data/processed")
FIGURES_DIR = Path("report/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

INPUT_FILE = PROCESSED_DIR / "tcga_train.csv"
ALPHA = 0.01           # seuil de significativité (p-valeur ajustée)
MAX_GENES_FINAL = 300  # plafond sur le nombre de gènes retenus, même si plus de gènes passent le seuil
TEST_SIZE = 0.2
RANDOM_STATE = 42

NON_GENE_COLS_HINTS = ("patient_id", "subtype")  # colonnes qu'on ne considère jamais comme des gènes


def identify_column_types(df: pd.DataFrame):
    """Sépare colonnes d'identification/label, colonnes cliniques, et colonnes de gènes."""
    meta_cols = [c for c in df.columns if c in NON_GENE_COLS_HINTS]
    # Heuristique par mots-clés pour repérer les colonnes cliniques les plus évidentes
    clinical_keywords = (
        "demographic", "diagnoses", "samples", "annotations", "OS.time", "OS", "_PATIENT"
    )
    clinical_cols = [c for c in df.columns if c not in meta_cols and any(k in c for k in clinical_keywords)]

    remaining = [c for c in df.columns if c not in meta_cols and c not in clinical_cols]

    # Parmi les colonnes restantes, seules celles réellement numériques peuvent être des gènes.
    # Toute colonne restante non-numérique (ex. 'id', 'disease_type', 'case_id', 'submitter_id',
    # 'primary_site') est en réalité une colonne clinique/administrative mal repérée par mots-clés,
    # donc reclassée ici plutôt que de fausser les tests statistiques sur les gènes.
    numeric_remaining = df[remaining].select_dtypes(include=[np.number]).columns.tolist()
    non_numeric_remaining = [c for c in remaining if c not in numeric_remaining]

    clinical_cols = clinical_cols + non_numeric_remaining
    gene_cols = numeric_remaining

    if non_numeric_remaining:
        print(f"  (reclassées en clinique car non-numériques : {non_numeric_remaining[:10]}"
              f"{'...' if len(non_numeric_remaining) > 10 else ''})")

    return meta_cols, clinical_cols, gene_cols


def differential_expression(X_train: pd.DataFrame, y_train: pd.Series, alpha=ALPHA, max_genes=MAX_GENES_FINAL):
    """
    Test de Kruskal-Wallis (ANOVA non-paramétrique) gène par gène, comparant les groupes de sous-types.
    Correction de tests multiples par Benjamini-Hochberg.
    IMPORTANT : ne s'exécute que sur le jeu d'entraînement pour éviter toute fuite de données.
    """
    print(f"\n--- Expression différentielle sur {X_train.shape[1]} gènes ({X_train.shape[0]} patients train) ---")

    # Pré-filtre : on écarte les gènes à variance nulle (non informatifs, et source de bugs
    # numériques dans certains tests statistiques sur des cas dégénérés).
    variances = X_train.var(axis=0)
    non_zero_var_genes = variances[variances > 0].index.tolist()
    n_dropped = X_train.shape[1] - len(non_zero_var_genes)
    if n_dropped > 0:
        print(f"  {n_dropped} gènes à variance nulle écartés avant le test statistique.")
    X_train = X_train[non_zero_var_genes]

    groups = [X_train.loc[y_train == cls] for cls in sorted(y_train.unique())]

    p_values = []
    gene_names = X_train.columns.tolist()

    for i, gene in enumerate(gene_names):
        if i % 5000 == 0:
            print(f"  ... {i}/{len(gene_names)} gènes testés")
        values_per_group = [g[gene].values for g in groups]
        try:
            stat, p = kruskal(*values_per_group)
            if not np.isfinite(p):
                p = 1.0
        except (ValueError, AttributeError, TypeError, ZeroDivisionError):
            # Cas dégénéré : valeurs identiques dans tous les groupes, variance nulle,
            # ou bug scipy sur certains cas limites -> gène traité comme non significatif
            p = 1.0
        p_values.append(p)

    p_values = np.array(p_values)
    reject, p_adj, _, _ = multipletests(p_values, alpha=alpha, method="fdr_bh")

    results = pd.DataFrame({
        "gene": gene_names,
        "p_value": p_values,
        "p_adj": p_adj,
        "significant": reject,
    }).sort_values("p_adj")

    n_significant = results["significant"].sum()
    print(f"  -> {n_significant} gènes significatifs (p_adj < {alpha})")

    selected = results[results["significant"]].head(max_genes)
    print(f"  -> {len(selected)} gènes retenus au final (plafond fixé à {max_genes})")

    return selected, results


def plot_pca(X_train_reduced: pd.DataFrame, y_train: pd.Series, output_path: Path):
    print("\n--- Génération de la PCA ---")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_train_reduced)

    pca = PCA(n_components=2, random_state=RANDOM_STATE)
    coords = pca.fit_transform(X_scaled)

    plot_df = pd.DataFrame(coords, columns=["PC1", "PC2"])
    plot_df["subtype"] = y_train.values

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=plot_df, x="PC1", y="PC2", hue="subtype", palette="Set2", alpha=0.7)
    plt.title(f"PCA des patients (gènes sélectionnés)\nVariance expliquée : PC1={pca.explained_variance_ratio_[0]:.1%}, PC2={pca.explained_variance_ratio_[1]:.1%}")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  -> figure sauvegardée : {output_path}")


def plot_umap(X_train_reduced: pd.DataFrame, y_train: pd.Series, output_path: Path):
    print("\n--- Génération de l'UMAP ---")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_train_reduced)

    reducer = umap.UMAP(random_state=RANDOM_STATE)
    coords = reducer.fit_transform(X_scaled)

    plot_df = pd.DataFrame(coords, columns=["UMAP1", "UMAP2"])
    plot_df["subtype"] = y_train.values

    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=plot_df, x="UMAP1", y="UMAP2", hue="subtype", palette="Set2", alpha=0.7)
    plt.title("UMAP des patients (gènes sélectionnés)")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  -> figure sauvegardée : {output_path}")


def write_justification(selected_genes: pd.DataFrame, all_results: pd.DataFrame, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("# Justification de la sélection de gènes\n\n")
        f.write(f"- Nombre total de gènes testés : {len(all_results)}\n")
        f.write(f"- Nombre de gènes significatifs (p_adj < {ALPHA}) : {all_results['significant'].sum()}\n")
        f.write(f"- Nombre de gènes retenus au final (plafond {MAX_GENES_FINAL}) : {len(selected_genes)}\n")
        f.write(f"- Méthode : test de Kruskal-Wallis par gène (comparaison des {5} sous-types PAM50), ")
        f.write("correction de tests multiples par Benjamini-Hochberg (FDR).\n\n")
        f.write("## Top 20 gènes les plus significatifs\n\n")
        f.write(selected_genes.head(20).to_markdown(index=False))
        f.write("\n")
    print(f"\nJustification sauvegardée : {output_path}")


def main():
    print(f"Chargement de {INPUT_FILE} ...")
    df = pd.read_csv(INPUT_FILE)
    print(f"  -> forme : {df.shape}")

    meta_cols, clinical_cols, gene_cols = identify_column_types(df)
    print(f"Colonnes méta : {len(meta_cols)} | cliniques : {len(clinical_cols)} | gènes : {len(gene_cols)}")

    X = df[gene_cols]
    y = df["subtype"]

    # Split train/test stratifié — la sélection de gènes ne doit voir QUE le train
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=TEST_SIZE, random_state=RANDOM_STATE, stratify=y
    )
    print(f"Train : {X_train.shape[0]} patients | Test : {X_test.shape[0]} patients")

    selected_genes, all_results = differential_expression(X_train, y_train)
    gene_list = selected_genes["gene"].tolist()

    X_train_reduced = X_train[gene_list]
    X_test_reduced = X_test[gene_list]

    # Sauvegarde des jeux réduits, avec les colonnes méta/cliniques ré-attachées
    train_out = pd.concat([
        df.loc[X_train.index, meta_cols + clinical_cols].reset_index(drop=True),
        X_train_reduced.reset_index(drop=True),
    ], axis=1)
    test_out = pd.concat([
        df.loc[X_test.index, meta_cols + clinical_cols].reset_index(drop=True),
        X_test_reduced.reset_index(drop=True),
    ], axis=1)

    train_out.to_csv(PROCESSED_DIR / "train_reduced.csv", index=False)
    test_out.to_csv(PROCESSED_DIR / "test_reduced.csv", index=False)
    print(f"\nFichiers sauvegardés : train_reduced.csv ({train_out.shape}), test_reduced.csv ({test_out.shape})")

    plot_pca(X_train_reduced, y_train, FIGURES_DIR / "pca_subtypes.png")
    plot_umap(X_train_reduced, y_train, FIGURES_DIR / "umap_subtypes.png")

    write_justification(selected_genes, all_results, Path("report/gene_selection_justification.md"))

    print("\n=== Phase 3 terminée ===")


if __name__ == "__main__":
    main()