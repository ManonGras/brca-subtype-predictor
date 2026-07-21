# Comparaison des modèles - Phase 4

## Résultats de validation croisée (5-fold, sur le train set uniquement)

| Modèle | Accuracy (moyenne ± écart-type) | F1 macro (moyenne ± écart-type) |
|---|---|---|
| Baseline (classe majoritaire) | 0.509 ± 0.001 | 0.135 ± 0.000 |
| Régression logistique | 0.867 ± 0.032 | 0.807 ± 0.032 |
| Random Forest | 0.865 ± 0.034 | 0.753 ± 0.050 |
| XGBoost | 0.899 ± 0.024 | 0.815 ± 0.033 |

## Modèle retenu : XGBoost

## Rapport de classification sur le test set (jamais vu pendant l'entraînement)

```
              precision    recall  f1-score   support

  BRCA_Basal      1.000     0.971     0.985        34
   BRCA_Her2      0.917     0.688     0.786        16
   BRCA_LumA      0.863     0.880     0.871       100
   BRCA_LumB      0.733     0.825     0.776        40
 BRCA_Normal      0.600     0.429     0.500         7

    accuracy                          0.853       197
   macro avg      0.823     0.758     0.784       197
weighted avg      0.855     0.853     0.852       197

```

## Analyse de survie

Test du log-rank multivarié (sous-types prédits) : p = 0.001282
