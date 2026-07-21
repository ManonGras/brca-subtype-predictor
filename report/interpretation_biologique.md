# Interprétation biologique - Phase 5 (SHAP)

## Top 20 gènes les plus influents selon SHAP

| gene       |   mean_abs_shap | dans_PAM50_officiel   |
|:-----------|----------------:|:----------------------|
| MLPH       |       0.311588  | Oui                   |
| SFRP1      |       0.299523  | Oui                   |
| FOXC1      |       0.28769   | Oui                   |
| GATA3      |       0.218003  | Non                   |
| PPP1R14C   |       0.198874  | Non                   |
| CENPA      |       0.194971  | Non                   |
| UBE2T      |       0.159541  | Oui                   |
| ESR1       |       0.158404  | Oui                   |
| CCDC170    |       0.138326  | Non                   |
| IL6ST      |       0.129552  | Non                   |
| XBP1       |       0.100387  | Non                   |
| NAT1       |       0.0960546 | Oui                   |
| AL356311.1 |       0.0847005 | Non                   |
| NDC80      |       0.0789337 | Oui                   |
| TNS2       |       0.0697862 | Non                   |
| F7         |       0.0605868 | Non                   |
| KRT16      |       0.0562864 | Non                   |
| CEP55      |       0.0551326 | Oui                   |
| HID1       |       0.0525381 | Non                   |
| EXO1       |       0.0487428 | Oui                   |

## Confrontation avec la signature PAM50 officielle

9 des 20 gènes les plus importants du modèle appartiennent à la signature PAM50 originale
(Parker et al., 2009) : `CEP55`, `ESR1`, `EXO1`, `FOXC1`, `MLPH`, `NAT1`, `NDC80`, `SFRP1`, `UBE2T`.

Ce recoupement, sur environ 20 000 gènes candidats possibles, constitue un indice fort que le
modèle s'appuie sur des gènes biologiquement pertinents et déjà validés dans la littérature,
plutôt que sur des corrélations arbitraires du jeu de données.

## Contexte biologique des gènes influents hors PAM50

**GATA3** — Facteur de transcription maître de la différenciation épithéliale luminale mammaire.
Bien qu'absent de la liste PAM50 originale de 2009, c'est aujourd'hui l'un des marqueurs
immunohistochimiques les plus utilisés en clinique pour confirmer une origine mammaire luminale,
et il est fréquemment muté dans les tumeurs luminales. Son rôle dominant dans le modèle est donc
parfaitement cohérent, simplement absent de la signature historique.

**CCDC170** — Particulièrement intéressant : ce gène est situé juste à côté du locus *ESR1* sur le
chromosome 6q25, et des fusions *ESR1-CCDC170* ont été décrites dans la littérature comme associées
à un phénotype plus agressif et à une résistance au tamoxifène dans les tumeurs luminales. Sa
proximité génomique avec *ESR1* explique probablement pourquoi il apparaît comme co-informatif dans
le modèle.

**XBP1** — Facteur de transcription de la réponse au stress du réticulum endoplasmique (UPR), connu
pour être régulé en aval de la signalisation des récepteurs aux œstrogènes et associé à la
résistance aux traitements endocriniens dans les tumeurs ER+.

**CENPA** — Protéine centromérique essentielle à la mitose, dans la même famille fonctionnelle que
plusieurs gènes de prolifération déjà présents dans PAM50 (*CDC20*, *MKI67*). Sa présence est
cohérente avec son rôle dans les gènes distinguant les tumeurs à forte/faible prolifération
(LumA vs LumB, notamment).

**IL6ST (gp130)** — Sous-unité du récepteur à l'IL-6, impliquée dans la voie de signalisation
JAK/STAT, dont l'activation a été rapportée dans certains contextes de résistance thérapeutique en
cancer du sein.

**KRT16** — Kératine associée à un phénotype épidermoïde/basal, cohérente avec son influence sur la
distinction Basal-like observée dans les graphiques SHAP par classe.

**PPP1R14C, TNS2, F7, HID1, AL356311.1** — Fonctions moins directement documentées dans la
littérature spécifique au cancer du sein. À mentionner avec prudence dans le rapport final :
possibles signaux biologiques réels non encore caractérisés dans ce contexte précis, ou possibles
artefacts liés à des facteurs techniques ou de composition tissulaire (`F7`, par exemple, code pour
un facteur de coagulation, plausiblement un signal de contamination vasculaire/stromale plutôt
qu'un signal tumoral intrinsèque — point à discuter comme limite du modèle).

## Cohérence par sous-type (résumé des graphiques SHAP par classe)

- **LumA** : dominé par `CENPA` (prolifération, effet négatif) et `GATA3`/`PGR` (effet positif) —
  cohérent avec un sous-type à faible prolifération et forte activité des récepteurs hormonaux.
- **LumB** : `SFRP1` et `ESR1` dominent, aux côtés de gènes de prolifération (`UBE2T`, `CENPA`,
  `BIRC5`) — cohérent avec la définition clinique de LumB comme "Luminal A mais plus proliférant".
- **Her2** : `ESR1` et `GREB1` (gène régulé par les œstrogènes) apparaissent en effet négatif —
  cohérent, Her2 étant généralement associé à des récepteurs hormonaux faibles ou négatifs.
- **Basal** : `MLPH` très fortement négatif — cohérent, `MLPH` étant un marqueur luminal classique,
  donc son absence est un signal attendu du phénotype basal.
- **Normal-like** : dominé par `GATA3` à l'opposé du profil Basal — cohérent avec le caractère
  "proche du tissu normal" attribué à cette catégorie ambiguë.

## Limite méthodologique à mentionner

La catégorie *Normal-like* (la plus petite de la cohorte, 36 patients à l'entraînement) est
cliniquement reconnue comme une catégorie ambiguë, reflétant parfois une faible pureté tumorale
(contamination par du tissu sain adjacent) plutôt qu'une entité biologique distincte à part entière.
Les performances plus faibles du modèle sur cette classe (recall 0.429 sur le test set) sont donc
attendues et documentées dans la littérature, pas un simple défaut du pipeline.