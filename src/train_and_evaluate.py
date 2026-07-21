"""
Phase 4 — Modélisation : classification du sous-type + analyse de survie
==========================================================================
Entrées : data/processed/train_reduced.csv, data/processed/test_reduced.csv
Sorties :
  - models/subtype_classifier.pkl
  - report/model_comparison.md
  - report/figures/confusion_matrix.png
  - report/figures/kaplan_meier_by_predicted_subtype.png

À exécuter depuis la racine du projet :
    python src/train_and_evaluate.py
"""

import pandas as pd
import numpy as np
import joblib
from pathlib import Path

from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.metrics import classification_report, confusion_matrix, ConfusionMatrixDisplay, f1_score
from xgboost import XGBClassifier

import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

PROCESSED_DIR = Path("data/processed")
MODELS_DIR = Path("models")
FIGURES_DIR = Path("report/figures")
MODELS_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_STATE = 42
N_FOLDS = 5


def load_data():
    train = pd.read_csv(PROCESSED_DIR / "train_reduced.csv")
    test = pd.read_csv(PROCESSED_DIR / "test_reduced.csv")
    print(f"Train : {train.shape} | Test : {test.shape}")
    return train, test


def get_gene_columns(df: pd.DataFrame) -> list:
    """Les colonnes de gènes sont celles qui restent après avoir retiré les colonnes connues
    (identifiants, sous-type, variables cliniques). On les identifie comme les colonnes numériques
    qui ne font pas partie des mots-clés cliniques usuels."""
    non_gene_hints = ("patient_id", "subtype")
    clinical_keywords = ("demographic", "diagnoses", "samples", "annotations", "OS.time", "OS", "_PATIENT")
    candidates = [c for c in df.columns if c not in non_gene_hints and not any(k in c for k in clinical_keywords)]
    numeric_candidates = df[candidates].select_dtypes(include=[np.number]).columns.tolist()
    return numeric_candidates


def cross_validate_models(X_train, y_train_encoded):
    print("\n--- Validation croisée (5-fold stratifié) ---")
    cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)

    models = {
        "Baseline (classe majoritaire)": DummyClassifier(strategy="most_frequent"),
        "Régression logistique": LogisticRegression(max_iter=2000, random_state=RANDOM_STATE),
        "Random Forest": RandomForestClassifier(n_estimators=300, random_state=RANDOM_STATE),
        "XGBoost": XGBClassifier(random_state=RANDOM_STATE, eval_metric="mlogloss"),
    }

    results = {}
    for name, model in models.items():
        scores = cross_validate(
            model, X_train, y_train_encoded, cv=cv,
            scoring=["accuracy", "f1_macro"], n_jobs=-1
        )
        results[name] = {
            "accuracy_mean": scores["test_accuracy"].mean(),
            "accuracy_std": scores["test_accuracy"].std(),
            "f1_macro_mean": scores["test_f1_macro"].mean(),
            "f1_macro_std": scores["test_f1_macro"].std(),
        }
        print(f"  {name:35s} | accuracy={results[name]['accuracy_mean']:.3f} (+/-{results[name]['accuracy_std']:.3f})"
              f" | F1 macro={results[name]['f1_macro_mean']:.3f} (+/-{results[name]['f1_macro_std']:.3f})")

    return results, models


def evaluate_on_test(best_model, model_name, X_train, y_train_encoded, X_test, y_test_encoded, label_encoder):
    print(f"\n--- Évaluation finale sur le test set (modèle retenu : {model_name}) ---")
    best_model.fit(X_train, y_train_encoded)
    y_pred = best_model.predict(X_test)

    report = classification_report(
        y_test_encoded, y_pred, target_names=label_encoder.classes_, digits=3
    )
    print(report)

    cm = confusion_matrix(y_test_encoded, y_pred)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=label_encoder.classes_)
    fig, ax = plt.subplots(figsize=(7, 6))
    disp.plot(ax=ax, cmap="Blues", xticks_rotation=45)
    plt.title(f"Matrice de confusion - {model_name} (test set)")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "confusion_matrix.png", dpi=150)
    plt.close()
    print(f"Matrice de confusion sauvegardée : {FIGURES_DIR / 'confusion_matrix.png'}")

    return y_pred, report


def survival_analysis(test_df, y_pred_labels):
    """Courbes de Kaplan-Meier par sous-type PRÉDIT, + test du log-rank."""
    print("\n--- Analyse de survie (Kaplan-Meier par sous-type prédit) ---")

    if "OS.time" not in test_df.columns or "OS" not in test_df.columns:
        print("  Colonnes OS.time / OS introuvables dans le test set : analyse de survie ignorée.")
        return

    surv_df = test_df[["OS.time", "OS"]].copy()
    surv_df["predicted_subtype"] = y_pred_labels
    surv_df = surv_df.dropna(subset=["OS.time", "OS"])

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(8, 6))
    ax = plt.gca()

    for subtype in sorted(surv_df["predicted_subtype"].unique()):
        mask = surv_df["predicted_subtype"] == subtype
        kmf.fit(surv_df.loc[mask, "OS.time"], surv_df.loc[mask, "OS"], label=subtype)
        kmf.plot_survival_function(ax=ax)

    plt.title("Courbes de survie (Kaplan-Meier) par sous-type prédit - test set")
    plt.xlabel("Temps (jours)")
    plt.ylabel("Probabilité de survie")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "kaplan_meier_by_predicted_subtype.png", dpi=150)
    plt.close()
    print(f"Courbes de survie sauvegardées : {FIGURES_DIR / 'kaplan_meier_by_predicted_subtype.png'}")

    # Test du log-rank multivarié (comparaison des 5 courbes simultanément)
    result = multivariate_logrank_test(
        surv_df["OS.time"], surv_df["predicted_subtype"], surv_df["OS"]
    )
    print(f"Test du log-rank multivarié : p-valeur = {result.p_value:.4g}")
    return result.p_value


def write_report(cv_results, model_name, classification_rep, logrank_p, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("# Comparaison des modèles - Phase 4\n\n")
        f.write("## Résultats de validation croisée (5-fold, sur le train set uniquement)\n\n")
        f.write("| Modèle | Accuracy (moyenne ± écart-type) | F1 macro (moyenne ± écart-type) |\n")
        f.write("|---|---|---|\n")
        for name, res in cv_results.items():
            f.write(f"| {name} | {res['accuracy_mean']:.3f} ± {res['accuracy_std']:.3f} | "
                    f"{res['f1_macro_mean']:.3f} ± {res['f1_macro_std']:.3f} |\n")
        f.write(f"\n## Modèle retenu : {model_name}\n\n")
        f.write("## Rapport de classification sur le test set (jamais vu pendant l'entraînement)\n\n")
        f.write("```\n" + classification_rep + "\n```\n\n")
        if logrank_p is not None:
            f.write(f"## Analyse de survie\n\nTest du log-rank multivarié (sous-types prédits) : p = {logrank_p:.4g}\n")
    print(f"\nRapport sauvegardé : {output_path}")


def main():
    train, test = load_data()

    gene_cols = get_gene_columns(train)
    print(f"Nombre de gènes utilisés comme features : {len(gene_cols)}")

    X_train, y_train = train[gene_cols], train["subtype"]
    X_test, y_test = test[gene_cols], test["subtype"]

    # Encodage des labels (nécessaire pour XGBoost)
    label_encoder = LabelEncoder()
    y_train_encoded = label_encoder.fit_transform(y_train)
    y_test_encoded = label_encoder.transform(y_test)

    # Standardisation (utile pour la régression logistique, neutre pour les modèles à base d'arbres)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    cv_results, models = cross_validate_models(X_train_scaled, y_train_encoded)

    # Sélection du meilleur modèle (hors baseline) selon le F1 macro moyen en CV
    candidates = {k: v for k, v in cv_results.items() if "Baseline" not in k}
    best_model_name = max(candidates, key=lambda k: candidates[k]["f1_macro_mean"])
    best_model = models[best_model_name]
    print(f"\nMeilleur modèle retenu (F1 macro CV) : {best_model_name}")

    y_pred_encoded, classification_rep = evaluate_on_test(
        best_model, best_model_name, X_train_scaled, y_train_encoded, X_test_scaled, y_test_encoded, label_encoder
    )
    y_pred_labels = label_encoder.inverse_transform(y_pred_encoded)

    logrank_p = survival_analysis(test, y_pred_labels)

    # Sauvegarde du modèle + scaler + encoder (nécessaires pour réutiliser le modèle plus tard)
    joblib.dump({
        "model": best_model,
        "scaler": scaler,
        "label_encoder": label_encoder,
        "gene_cols": gene_cols,
    }, MODELS_DIR / "subtype_classifier.pkl")
    print(f"\nModèle sauvegardé : {MODELS_DIR / 'subtype_classifier.pkl'}")

    write_report(cv_results, best_model_name, classification_rep, logrank_p, Path("report/model_comparison.md"))

    print("\n=== Phase 4 terminée ===")


if __name__ == "__main__":
    main()