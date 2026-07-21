"""
Phase 7 — Démo interactive Streamlit
========================================
Utilise le modèle restreint (210 gènes, le plus robuste en généralisation externe, Phase 6).

Lancement depuis la racine du projet :
    streamlit run app/app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import joblib
import shap
import matplotlib.pyplot as plt
from pathlib import Path

MODELS_DIR = Path("models")
PROCESSED_DIR = Path("data/processed")

st.set_page_config(page_title="Prédiction du sous-type moléculaire - Démo", layout="wide")

# ----------------------------------------------------------------------------
# Avertissement éthique (toujours visible en haut de page)
# ----------------------------------------------------------------------------
st.warning(
    "⚠️ **Outil de démonstration à visée pédagogique.** Ce modèle a été entraîné sur des données "
    "publiques (TCGA) et validé sur une cohorte externe (METABRIC) dans le cadre d'un projet "
    "académique. Il ne constitue en aucun cas un dispositif de diagnostic médical, n'a pas été "
    "validé cliniquement, et ne doit jamais être utilisé pour une décision de soin réelle."
)

st.title("🧬 Prédiction du sous-type moléculaire du cancer du sein")
st.markdown(
    "Ce modèle (XGBoost, entraîné sur TCGA-BRCA, validé sur METABRIC) prédit le sous-type "
    "moléculaire PAM50 à partir d'un profil d'expression génique restreint à 210 gènes."
)


@st.cache_resource
def load_model():
    bundle = joblib.load(MODELS_DIR / "subtype_classifier_restricted210.pkl")
    return bundle


@st.cache_data
def load_example_patients():
    """Quelques patients du test set TCGA, utilisés comme exemples pré-chargés
    (on ne demande jamais à l'utilisateur de saisir manuellement des milliers de valeurs)."""
    test = pd.read_csv(PROCESSED_DIR / "test_reduced.csv")
    # On garde un échantillon varié : un patient par sous-type si possible
    examples = test.groupby("subtype", group_keys=False)[test.columns.tolist()].apply(
        lambda x: x.sample(1, random_state=42)
    )
    return examples.reset_index(drop=True)


bundle = load_model()
model = bundle["model"]
scaler = bundle["scaler"]
label_encoder = bundle["label_encoder"]
gene_cols = bundle["gene_cols"]

examples = load_example_patients()

# ----------------------------------------------------------------------------
# Choix du profil patient
# ----------------------------------------------------------------------------
st.header("1. Choisir un profil d'expression")

source = st.radio(
    "Source du profil",
    ["Exemple pré-chargé (cohorte TCGA, test set)", "Importer un fichier CSV"],
    horizontal=True,
)

profile = None
true_subtype = None

if source == "Exemple pré-chargé (cohorte TCGA, test set)":
    labels = [
        f"Patient {i+1} — sous-type réel : {row['subtype']}"
        for i, row in examples.iterrows()
    ]
    choice = st.selectbox("Sélectionner un patient exemple", labels)
    idx = labels.index(choice)
    profile = examples.loc[idx, gene_cols]
    true_subtype = examples.loc[idx, "subtype"]
else:
    uploaded = st.file_uploader(
        f"Fichier CSV avec une ligne et les {len(gene_cols)} colonnes de gènes attendues", type="csv"
    )
    if uploaded is not None:
        df_uploaded = pd.read_csv(uploaded)
        missing = [g for g in gene_cols if g not in df_uploaded.columns]
        if missing:
            st.error(f"{len(missing)} gènes attendus sont absents du fichier importé "
                      f"(ex. {missing[:5]}...). Impossible de continuer.")
        else:
            profile = df_uploaded.iloc[0][gene_cols]

# ----------------------------------------------------------------------------
# Prédiction + explication SHAP
# ----------------------------------------------------------------------------
if profile is not None:
    st.header("2. Prédiction du sous-type")

    X = profile.to_frame().T[gene_cols]
    X_scaled = scaler.transform(X)

    pred_encoded = model.predict(X_scaled)[0]
    pred_proba = model.predict_proba(X_scaled)[0]
    pred_label = label_encoder.inverse_transform([pred_encoded])[0]

    col1, col2 = st.columns([1, 2])
    with col1:
        st.metric("Sous-type prédit", pred_label)
        if true_subtype is not None:
            match = "✅ correspond" if true_subtype == pred_label else "❌ ne correspond pas"
            st.caption(f"Sous-type réel (référence) : {true_subtype} — {match}")

    with col2:
        proba_df = pd.DataFrame({
            "Sous-type": label_encoder.classes_,
            "Probabilité": pred_proba,
        }).sort_values("Probabilité", ascending=False)
        st.bar_chart(proba_df.set_index("Sous-type"))

    st.header("3. Gènes déterminants pour cette prédiction (SHAP)")
    st.markdown(
        "Les gènes ci-dessous ont le plus influencé la prédiction pour **ce patient précis** "
        "(pas une moyenne globale). Les valeurs positives poussent vers le sous-type prédit."
    )

    with st.spinner("Calcul de l'explication SHAP..."):
        explainer = shap.TreeExplainer(model)
        X_scaled_df = pd.DataFrame(X_scaled, columns=gene_cols)
        shap_values = explainer(X_scaled_df)

        values = shap_values.values
        if values.ndim == 3:
            local_shap = values[0, :, pred_encoded]
        else:
            local_shap = values[0]

        shap_df = pd.DataFrame({
            "gene": gene_cols,
            "shap_value": local_shap,
        })
        shap_df["abs_shap"] = shap_df["shap_value"].abs()
        top_local = shap_df.sort_values("abs_shap", ascending=False).head(10)

        fig, ax = plt.subplots(figsize=(7, 5))
        colors = ["#d62728" if v > 0 else "#1f77b4" for v in top_local["shap_value"][::-1]]
        ax.barh(top_local["gene"][::-1], top_local["shap_value"][::-1], color=colors)
        ax.set_xlabel("Valeur SHAP (impact sur la prédiction)")
        ax.set_title(f"Top 10 gènes déterminants — sous-type prédit : {pred_label}")
        plt.tight_layout()
        st.pyplot(fig)

    st.caption(
        "🔵 Bleu = pousse contre le sous-type prédit | 🔴 Rouge = pousse vers le sous-type prédit. "
        "Interprétation biologique détaillée disponible dans le rapport du projet "
        "(`report/interpretation_biologique.md`)."
    )
else:
    st.info("Sélectionne un patient exemple ou importe un fichier CSV pour lancer une prédiction.")

st.divider()
st.caption(
    "Projet académique — pipeline complet : TCGA-BRCA (entraînement) → sélection de gènes → "
    "XGBoost → validation externe METABRIC. Code source disponible sur demande."
)