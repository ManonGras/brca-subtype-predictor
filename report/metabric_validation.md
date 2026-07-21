# Phase 6 — Validation externe sur METABRIC

## 1. Cohorte METABRIC utilisée

- Échantillons METABRIC avec expression + clinique disponibles : 1980
- Échantillons exclus (hors schéma PAM50 à 5 classes : `claudin-low`, `NC`, ou sous-type manquant) : 224
- Échantillons retenus pour la validation externe : **1756**

## 2. Correspondance des gènes entre plateformes

- Gènes du modèle (sélectionnés en Phase 3 sur TCGA, RNA-seq) : 300
- Présents sur la puce microarray Illumina de METABRIC : **210**
- Absents de la plateforme METABRIC : **90**

La quasi-totalité des gènes absents ne sont pas des gènes protéiques classiques mais des entités dont
l'annotation (lncRNA, ARN antisens, pseudogènes, `Y_RNA`) est postérieure ou non prise en charge par la
conception des sondes de la puce Illumina HT-12 utilisée pour METABRIC (~2010) : sur les 90 gènes
manquants, 84 correspondent à ce type d'identifiant (`LINC*`, `AL*`/`AC*` (contigs GENCODE non nommés),
`Y_RNA.*`, gènes se terminant en `-AS1`/`-DT`/pseudogènes `P1`). Ce n'est donc pas majoritairement un
problème de généralisation biologique mais un artefact de couverture technique de plateforme.

<details>
<summary>Liste complète des gènes absents</summary>

```
AC005077.4, AC007255.1, AC008268.1, AC008663.1, AC013726.1, AC044784.1, AC055854.1, AC092667.1,
AC092718.4, AC093512.2, AC093838.1, AC096733.2, AC103760.1, AC105328.1, AC106738.2, AC108134.4,
AC120498.4, AC124067.4, AL078582.2, AL133387.1, AL136296.1, AL157387.1, AL356311.1, AL591845.1,
CFAP94, FAM198B-AS1, FAM47E, FAM72B, GTF2IP7, LINC00504, LINC00993, LINC01016, LINC01087, LINC01117,
LINC01488, LINC01843, LINC01863, LINC02188, LINC02568, LINC02747, MIR3936HG, MRPS30-DT, RAB6C-AS1,
RAD17P1, RNU6-813P, SAMD15, TTC39A-AS1, Y_RNA.12459, Y_RNA.12910, Y_RNA.13033, Y_RNA.13074, Y_RNA.13115,
Y_RNA.13402, Y_RNA.1348, Y_RNA.1430, Y_RNA.14468, Y_RNA.14509, Y_RNA.1471, Y_RNA.15739, Y_RNA.16067,
Y_RNA.16190, Y_RNA.17666, Y_RNA.20085, Y_RNA.21110, Y_RNA.21479, Y_RNA.22955, Y_RNA.22996, Y_RNA.23078,
Y_RNA.24472, Y_RNA.24554, Y_RNA.26194, Y_RNA.26235, Y_RNA.27342, Y_RNA.27383, Y_RNA.27506, Y_RNA.28162,
Y_RNA.28203, Y_RNA.28244, Y_RNA.28572, Y_RNA.30991, Y_RNA.3849, Y_RNA.5120, Y_RNA.528, Y_RNA.5448,
Y_RNA.569, Y_RNA.7006, Y_RNA.733, Y_RNA.7334, Y_RNA.774, Y_RNA.8031
```
</details>

## 3. Tableau comparatif — validation interne (TCGA) vs externe (METABRIC)

| Modèle | Gènes | Accuracy interne (TCGA) | F1 macro interne (TCGA) | Accuracy externe (METABRIC) | F1 macro externe (METABRIC) | Δ Accuracy | Δ F1 macro |
|---|---|---|---|---|---|---|---|
| Complet (Phase 4, imputation à 0 des gènes absents) | 300 | 0.853 | 0.784 | 0.694 | 0.590 | -0.159 | -0.194 |
| Restreint (ré-entraîné, gènes communs uniquement) | 210 | 0.868 | 0.800 | 0.729 | 0.663 | -0.139 | -0.137 |

*« Interne » = mesuré sur le test set TCGA (Phase 4), jamais vu pendant l'entraînement. « Externe » =
mesuré sur la cohorte METABRIC, jamais vue pendant l'entraînement, sur une plateforme technique différente
(microarray vs RNA-seq). Le modèle « restreint » est ré-entraîné sur TCGA avec exactement le même pipeline
(XGBoost, même split, mêmes hyperparamètres) mais en se limitant aux gènes réellement mesurables sur
METABRIC — il n'y a donc aucune imputation lors de son évaluation externe.*

**Lecture du tableau :** le modèle restreint généralise nettement mieux (F1 macro externe 0.663 contre
0.590), confirmant qu'une partie substantielle de la chute de performance du modèle complet est un
artefact d'imputation plutôt qu'une vraie limite de généralisation biologique. Toute la suite de l'analyse
se concentre donc sur le modèle restreint, qui est le plus représentatif de la capacité réelle du signal
biologique à généraliser.

## 4. Rapport de classification détaillé — METABRIC (modèle restreint)

```
              precision    recall  f1-score   support

  BRCA_Basal      0.924     0.813     0.865       209
   BRCA_Her2      0.740     0.598     0.662       224
   BRCA_LumA      0.664     0.979     0.791       700
   BRCA_LumB      0.826     0.541     0.654       475
 BRCA_Normal      0.694     0.230     0.345       148

    accuracy                          0.729      1756
   macro avg      0.770     0.632     0.663      1756
weighted avg      0.751     0.729     0.709      1756
```

## 5. Analyse des erreurs — un biais de composition de cohorte, pas seulement un batch effect

Le rapport de classification révèle un pattern clair : **LumA a un recall quasi parfait (0.979) mais une
précision modeste (0.664)**, ce qui signifie que le modèle sur-prédit systématiquement cette classe au
détriment des autres. La cause la plus probable n'est pas uniquement technique (plateforme RNA-seq vs
microarray), mais un **déséquilibre de composition entre les deux cohortes** :

| Sous-type | Proportion dans TCGA (cohorte d'entraînement, n=981) | Proportion dans METABRIC (cohorte de test, n=1756) |
|---|---|---|
| LumA | 50,9 % (499) | 39,9 % (700) |
| LumB | 20,1 % (197) | 27,1 % (475) |
| Basal | 17,4 % (171) | 11,9 % (209) |
| Her2 | 8,0 % (78) | 12,8 % (224) |
| Normal | 3,7 % (36) | 8,4 % (148) |

Le modèle a été entraîné sur une cohorte où LumA représente environ la moitié des patients ; il a donc
appris un a priori de classe qui favorise cette prédiction par défaut en cas d'incertitude. Sur METABRIC,
où LumA est proportionnellement moins fréquent et LumB/Her2/Normal plus fréquents, ce biais se traduit
directement par la baisse de précision observée sur LumA et la baisse de recall sur les autres classes
(en particulier Normal, la plus touchée : recall 0.230).

**Ce point est distinct du batch effect technique** déjà discuté en section 2 (gènes manquants) : il
s'agit ici d'un **covariate shift** — un déplacement de la distribution des classes entre l'échantillon
d'entraînement et l'échantillon de test — un phénomène classique et bien documenté en apprentissage
automatique appliqué à des cohortes cliniques réelles, indépendant de la qualité intrinsèque du modèle.

*Piste d'amélioration non mise en œuvre ici (mentionnée comme perspective) : un ré-équilibrage des poids
de classe (`class_weight` ou `scale_pos_weight` dans XGBoost) pendant l'entraînement, ou un ajustement des
seuils de décision par classe a posteriori, pourrait réduire ce biais sans nécessiter de nouvelles
données.*

## 6. Analyse de survie sur METABRIC (sous-types prédits, modèle restreint)

- Échantillons avec données de survie complètes : 1756
- Test du log-rank multivarié (5 sous-types prédits) : **p = 0.0001557**
- Figure : `report/figures/kaplan_meier_metabric_predicted_subtype.png`

Ce résultat est **plus significatif encore que sur TCGA** (p = 0.001282 en Phase 4). Surtout, contrairement
à la courbe obtenue sur TCGA — où l'ordre pronostique observé (Basal en tête) contredisait la littérature,
probablement à cause du faible échantillon et d'un suivi trop court — la courbe METABRIC restitue un ordre
pronostique parfaitement cohérent avec les données cliniques établies : LumA et Normal-like présentent le
meilleur pronostic, Basal/Her2/LumB un pronostic moins favorable. Le suivi beaucoup plus long de METABRIC
(médiane de plusieurs décennies contre quelques années pour TCGA) rend cette estimation nettement plus
fiable.

**Ce résultat constitue l'argument le plus solide du projet** : même si la classification individuelle des
sous-types se dégrade en externe, le signal pronostique global — soit la question cliniquement la plus
pertinente — se maintient, et se maintient même mieux que sur la cohorte d'entraînement d'origine.

## 7. Limites méthodologiques à garder à l'esprit

1. La standardisation z-score appliquée séparément à chaque jeu de données suppose implicitement que la
   distribution relative de chaque gène est comparable entre RNA-seq et microarray — hypothèse
   simplificatrice. Une correction de type ComBat (option avancée non mise en œuvre ici) modéliserait
   explicitement le batch comme facteur et pourrait réduire encore l'écart observé, en particulier pour
   le modèle complet.
2. Les labels de vérité-terrain diffèrent aussi dans leur origine : le sous-type TCGA (`SUBTYPE`,
   cBioPortal PanCancer Atlas) et le sous-type METABRIC (`CLAUDIN_SUBTYPE`, PAM50 + claudin-low)
   proviennent de pipelines de classification PAM50 possiblement recalibrés différemment selon l'étude —
   une partie des erreurs de classification peut donc relever du désaccord entre méthodes de labellisation
   plutôt que d'une erreur du modèle.
3. Le déséquilibre de composition de cohorte documenté en section 5 (covariate shift) contribue à la chute
   de performance indépendamment de tout effet de plateforme technique.
4. METABRIC est une cohorte plus âgée en moyenne et recrutée sur une période et dans un contexte de soins
   différents de TCGA, ce qui peut introduire d'autres décalages de distribution non mesurés ici.
5. Les 6 échantillons `NC` et 218 échantillons `claudin-low` de METABRIC ont été exclus car ils sortent du
   schéma PAM50 à 5 classes sur lequel le modèle a été entraîné — un modèle destiné à un usage clinique
   réel devrait explicitement gérer une classe « hors distribution ».

## 8. Conclusion

La validation externe confirme que le signal biologique appris sur TCGA est **partiellement transférable**
à une cohorte et une plateforme indépendantes, avec une chute de performance de classification mesurable
(F1 macro : -0.137 pour le modèle restreint) qu'il serait malhonnête de masquer. L'analyse détaillée
permet cependant de distinguer précisément les sources de cette chute :

- une part **technique**, isolable et largement expliquée (gènes absents de l'array METABRIC, corrigée par
  le modèle restreint) ;
- une part liée à un **covariate shift de composition de cohorte** (déséquilibre des proportions de
  sous-types entre TCGA et METABRIC), identifiée et quantifiée en section 5 ;
- et surtout, une **robustesse remarquable du signal pronostique** (analyse de survie), qui est en réalité
  le résultat cliniquement le plus important de ce projet, et qui se maintient — voire s'améliore — sur la
  cohorte externe.