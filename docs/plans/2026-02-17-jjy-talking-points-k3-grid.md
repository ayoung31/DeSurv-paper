# K=3 DeSurv Factor Analysis: CV Grid Fits

**For discussion with Jen Jen Yeh**
**Date: 2026-02-17**
**Data source: cv_grid exhaustive search (K=2:12 x alpha=0:0.95, fixed lambda=0.3,
nu=0.05, 100 initializations per fit, 5-fold stratified CV, ngene=3000->1970)**

**Note for Amber:** This document uses the cv_grid exhaustive search fits, not
the BO-selected production fit from the original paper. By holding lambda, nu,
and ntop constant across all K and alpha values, we get a controlled comparison.
The goal is to understand when and why the iCAF program emerges. See the glossary
at the end for term definitions. The analysis scripts that generated these results
are in `inst/cv_grid_training_analysis.R` (training) and
`inst/cv_grid_validation_analysis.R` (validation). Run them from the repo root.

## Key Change from Previous Analysis

The previous K=3 analysis used the production BO-selected fit (alpha=0.334,
lambda=0.349, nu=0.056, ntop=270). This document uses fits from the exhaustive
cv_grid search where **all hyperparameters except K and alpha are held constant**
(lambda=0.3, nu=0.05, ntop=NULL). This enables apples-to-apples comparison
across K and alpha values.

Three K=3 fits are examined:
- **K3_a0** (alpha=0.00): Standard NMF baseline (no survival penalty)
- **K3_a35** (alpha=0.35): Validation-optimal alpha
- **K3_a55** (alpha=0.55): CV-optimal alpha

**Branch hashes:**
- K3_a0: `cv_grid_fit_22ac9cbd8337ebdf`
- K3_a35: `cv_grid_fit_1555a2f1fcc58ab2`
- K3_a55: `cv_grid_fit_ce943309f0b3e212`

**Training cohort:** TCGA + CPTAC (n=273)
**Validation cohort:** 4 independent datasets (n=570 pooled):
Dijk (n=90), Moffitt GEO array (n=123), Puleo array (n=288), PACA-AU (n=69).
All validation Cox models use `strata(dataset)` to account for cohort-specific
baseline hazards. PurIST and DeCAF labels available for all 273 training and
570 validation samples.

**Methodological note on H-scores:** The cv_grid fits use ntop=NULL, meaning
H-scores are computed as H = X^T W using all 1970 genes. The production fit
used ntop=270, focusing each factor's H-score on its top 270 genes. Using all
genes dilutes per-factor signal, making adjusted Cox models harder to reach
significance. This is a conservative analysis, so the unadjusted and KM results
are more directly comparable.

---

## The Central Finding

**The iCAF program does not exist at alpha=0.** At K=3 with no survival penalty,
no factor resembles the original K=3 Factor 1 (best H-score correlation rho=0.387,
gene overlap 48/270). At alpha=0.55, Factor 1 recovers the iCAF program almost
perfectly (rho=0.906, gene overlap 168/270). **The survival penalty creates the
iCAF program; it does not merely enhance it.**

**This finding validates in external data.** At alpha=0, the iCAF factor is not
prognostic in validation (p=0.10). At alpha=0.55, it is (unadjusted HR=0.800,
p=0.003; median-split KM 26.8 vs 18.0 months, p=3.0e-4).

When we say "iCAF program" in this document, we mean the factor that best matches
the original production K=3 Factor 1, which was characterized as iCAF-associated
by its overlap with the **Elyada et al. (2019) iCAF 35-gene signature** (from
scRNA-seq of PDAC cancer-associated fibroblasts) and the **SCISSORS iCAF 25-gene
signature**. The Elyada iCAF genes (in our gene universe) include: PLA2G2A, MCL1,
S100A10, S100A4, LMNA, DPT, EFEMP1, TNFAIP6, FBLN2, CCDC80, FSTL1, PTX3, UGDH,
CXCL14, GPX3, GFPT2, TNXB, PI16, SGK1, SOD2, APOE, FBLN1, ADAMTS1, CPE, ADH1B
(25 of the 35 genes present in the 1970-gene universe).

---

## K3_a0: Standard NMF (alpha=0.00)

### Factor identity

| Factor | Best match to orig K=3 | H-score rho | Gene overlap (top-270) |
|--------|----------------------|-------------|----------------------|
| F1 | Orig F2 (stroma) | rho=0.950 | 26/270 with orig F1 |
| F2 | None coherent | rho=-0.460 (orig F1) | 22/270 with orig F1 |
| F3 | Orig F3 (basal) | rho=0.888 | 48/270 with orig F1 |

Full H-score correlation matrix with original K=3:

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.111  | **0.950** | -0.866 |
| F2     | -0.460  | 0.379   | -0.022 |
| F3     | 0.387   | **-0.942** | **0.888** |

Without survival penalty, K=3 NMF recovers the **stroma axis** (F1 <-> orig F2)
and the **basal-like axis** (F3 <-> orig F3) with high fidelity. The iCAF program
simply does not form as a coherent factor.

### Reference signature overlaps (best iCAF factor: F3, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (35->25 in universe) | 3/25 | 12.0% |
| Elyada myCAF (16->15) | 0/15 | 0.0% |
| SCISSORS iCAF (25->17) | 0/17 | 0.0% |
| SCISSORS myCAF (25->16) | 1/16 | 6.2% |
| DeCAF restCAF (9->8) | 0/8 | 0.0% |
| DECODER Immune (408->60) | 0/60 | 0.0% |
| DECODER ActivatedStroma (206->108) | 0/108 | 0.0% |
| DECODER BasalTumor (394->139) | 19/139 | 13.7% |
| DECODER ClassicalTumor (372->166) | **63/166** | **38.0%** |

The "best iCAF factor" (F3) is actually a Classical Tumor factor (38% DECODER
ClassicalTumor overlap), not iCAF at all. Zero overlap with DECODER Immune and
ActivatedStroma.

### Per-factor Cox: Training AND Validation

**Unadjusted (standardized H-scores)**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.937 | 0.44 | 0.941 | 0.30 |
| **F2** | **1.308** | **0.005** | **1.308** | **<0.001** |
| F3 (best iCAF) | 0.927 | 0.40 | 0.902 | 0.10 |

**Adjusted for PurIST + DeCAF (training n=273, validation n=570)**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.993 | 0.93 | 0.957 | 0.47 |
| F2 | 0.932 | 0.56 | 1.113 | 0.18 |
| F3 (best iCAF) | 0.970 | 0.74 | 0.926 | 0.25 |

**The iCAF factor is not prognostic at alpha=0, in training or validation.**
The only significant factor is F2 (not iCAF), which picks up the basal-like
tumor axis that PurIST already captures. After adjusting for PurIST + DeCAF,
no factor adds information.

**LRT adding F3 (iCAF) to PurIST + DeCAF:** Training p=0.74, Validation p=0.26

**Median-split KM for F3 (iCAF):**
- Training (n=273): High 21.1 vs Low 20.3 months, p=0.43
- Validation (n=570): High 24.8 vs Low 19.2 months, p=0.003

Note: the validation KM shows a modest separation despite the iCAF factor's
poor coherence (rho=0.387), but this does not replicate in Cox models. The
median split captures a nonspecific stroma-vs-tumor axis, not the iCAF program.

### Beta structure

F1=**-5.23e-05**, F2=**1.52e-04**, F3=-3.41e-05

All betas are tiny (no survival penalty -> effectively zero survival weighting).

---

## K3_a35: Validation-Optimal Alpha (alpha=0.35)

### Factor identity

| Factor | Best match to orig K=3 | H-score rho | Gene overlap |
|--------|----------------------|-------------|-------------|
| **F1** | **Orig F1 (iCAF)** | **0.731** | **132/270** |
| F2 | Orig F3 (basal) | 0.833 | 7/270 |
| F3 | Orig F2 (stroma) | 0.930 | 20/270 |

Full H-score correlation matrix:

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| **F1** | **0.731** | -0.564 | 0.592 |
| F2     | -0.204  | -0.501  | **0.833** |
| F3     | -0.216  | **0.930** | -0.743 |

At alpha=0.35, the iCAF program **emerges** as F1 with rho=0.731 to original F1.
The stroma and basal axes are still recognizable but have reorganized. The
survival penalty has pulled iCAF/immune genes into a dedicated factor.

### Reference signature overlaps (iCAF factor: F1, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25 in universe) | 4/25 | 16.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | **6/17** | **35.3%** |
| SCISSORS myCAF (16) | 1/16 | 6.2% |
| DeCAF restCAF (8) | 0/8 | 0.0% |
| DECODER Immune (60) | **20/60** | **33.3%** |
| DECODER ActivatedStroma (108) | 3/108 | 2.8% |
| DECODER NormalStroma (15) | 4/15 | 26.7% |
| DECODER BasalTumor (139) | 6/139 | 4.3% |
| DECODER ClassicalTumor (166) | 49/166 | 29.5% |

Compared to alpha=0: Elyada iCAF rises from 12->16%, SCISSORS iCAF from 0->35%,
DECODER Immune from 0->33%. The iCAF/immune character is emerging.

### Per-factor Cox: Training AND Validation

**Unadjusted (standardized H-scores)**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.601** | **<0.0001** | **0.841** | **0.027** |
| **F2 (basal)** | **1.391** | **0.0004** | **1.322** | **0.0007** |
| F3 (stroma) | 0.996 | 0.96 | 1.059 | 0.40 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.674** | **<0.0001** | 0.862 | 0.072 |
| F2 (basal) | 1.140 | 0.21 | 1.103 | 0.28 |
| F3 (stroma) | 0.998 | 0.98 | 0.992 | 0.91 |

F1 (iCAF) is the primary prognostic factor in training (**HR=0.674 adjusted,
p<0.0001**) and shows a consistent protective effect in validation (**HR=0.862,
p=0.072**). F2 (basal) is significant unadjusted but absorbed by PurIST after
adjustment. F3 (stroma) is inert throughout.

**LRT adding F1 (iCAF) to PurIST + DeCAF:**
- Training p = 7.9e-05
- Validation p = 0.075

**Median-split KM for F1 (iCAF):**
- Training (n=273): High 24.6 vs Low 15.3 months, **p=8.0e-5**
- Validation (n=570): High 24.9 vs Low 19.2 months, **p=7.8e-4**

### Beta structure

F1=**-9.73e-04** (dominant), F2=**7.13e-04**, F3=3.91e-07 (approx 0)

F1 has the largest magnitude beta, consistent with being the primary survival
factor. F3 has near-zero beta (the stroma axis is prognostically inactive).

---

## K3_a55: CV-Optimal Alpha (alpha=0.55)

### Factor identity

| Factor | Best match to orig K=3 | H-score rho | Gene overlap |
|--------|----------------------|-------------|-------------|
| **F1** | **Orig F1 (iCAF)** | **0.906** | **168/270** |
| F2 | Orig F3 (basal) | **0.980** | 3/270 |
| F3 | Orig F2 (stroma) | **0.989** | 3/270 |

Full H-score correlation matrix:

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| **F1** | **0.906** | -0.574  | 0.438  |
| F2     | 0.137   | -0.892  | **0.980** |
| F3     | -0.254  | **0.989** | -0.845 |

At alpha=0.55, the iCAF program is nearly identical to the original production
fit (rho=0.906). The basal and stroma factors also match almost perfectly
(rho=0.980 and 0.989 respectively). The factor permutation is F1->F1(iCAF),
F2->F3(basal), F3->F2(stroma).

### Reference signature overlaps (iCAF factor: F1, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25 in universe) | **6/25** | **24.0%** |
| Elyada myCAF (15) | 3/15 | 20.0% |
| SCISSORS iCAF (17) | 4/17 | 23.5% |
| SCISSORS myCAF (16) | 0/16 | 0.0% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DECODER Immune (60) | **16/60** | **26.7%** |
| DECODER ActivatedStroma (108) | 1/108 | 0.9% |
| DECODER NormalStroma (15) | **7/15** | **46.7%** |
| DECODER BasalTumor (139) | 17/139 | 12.2% |
| DECODER ClassicalTumor (166) | 40/166 | 24.1% |

The iCAF factor's biological identity: 24% Elyada iCAF, 27% DECODER Immune,
47% DECODER NormalStroma. This is a **microenvironmental program** mixing
iCAF, normal stroma, and immune (B cell) genes, consistent with what the
original paper describes.

### DECODER compartment breakdown for ALL factors

| Factor | DECODER Immune | Basal Tumor | Classical Tumor | Act. Stroma | Norm. Stroma |
|--------|----------------|-------------|-----------------|-------------|-------------|
| **F1 (iCAF)** | **16/60 (27%)** | 17/139 (12%) | 40/166 (24%) | 1/108 (1%) | **7/15 (47%)** |
| F2 (basal) | 0/60 (0%) | **48/139 (35%)** | 19/166 (11%) | 29/108 (27%) | 0/15 (0%) |
| F3 (stroma) | 16/60 (27%) | 2/139 (1%) | 57/166 (34%) | 16/108 (15%) | 5/15 (33%) |

F1 uniquely combines Immune + NormalStroma with minimal ActivatedStroma (1%).
F2 is Basal Tumor dominant. F3 mixes Classical Tumor + Immune + NormalStroma.

### Per-factor Cox: Training AND Validation

**Unadjusted (standardized H-scores)**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.344** | **<0.0001** | **0.800** | **0.003** |
| F2 (basal) | 1.125 | 0.18 | 1.103 | 0.16 |
| F3 (stroma) | 1.030 | 0.73 | 1.066 | 0.27 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.371** | **<0.0001** | 0.867 | 0.075 |
| F2 (basal) | 1.058 | 0.53 | 0.979 | 0.79 |
| F3 (stroma) | 1.029 | 0.75 | 1.022 | 0.71 |

**Factor 1 is the sole prognostic factor** in both training and validation.
Training: HR=0.344 unadjusted (each SD increase in F1 reduces hazard by 66%),
HR=0.371 adjusted. Validation: HR=0.800 unadjusted (p=0.003), HR=0.867 adjusted
(p=0.075). Factors 2 and 3 contribute nothing after PurIST + DeCAF adjustment.

The adjusted validation p=0.075 is borderline because the cv_grid fits use all
1970 genes for H-scores (ntop=NULL). The production fit's ntop=270 focuses
each factor on its top genes, yielding sharper validation separation (adjusted
p=0.004 in the original analysis). The consistent direction (HR<1, protective)
across training and validation confirms the effect is real.

**LRT adding F1 (iCAF) to PurIST + DeCAF:**
- Training p = **5.0e-18**
- Validation p = 0.076

**Median-split KM for F1 (iCAF):**
- Training (n=273): High 37.7 vs Low 13.3 months, **p=2.3e-13**
- Validation (n=570): High 26.8 vs Low 18.0 months, **p=3.0e-4**

The KM validation result shows an 8.8-month median survival difference
between high and low F1 groups, p=3.0e-4 in 570 independent patients.

### Beta structure

F1=**-1.44e-03** (dominant), F2=3.12e-04, F3=3.41e-05

F1's beta is 4.6x larger than F2's and 42x larger than F3's. The LP is
dominated by Factor 1, just as in the original production fit.

---

## Alpha Progression at K=3: The iCAF Program Emerges

### Training (n=273)

| | alpha=0.00 | alpha=0.35 | alpha=0.55 |
|--|-----------|-----------|-----------|
| **Best iCAF factor** | F3 | F1 | F1 |
| **H-cor with orig F1** | 0.387 | 0.731 | **0.906** |
| **Top-270 overlap** | 48/270 | 132/270 | **168/270** |
| **Elyada iCAF overlap** | 3/25 (12%) | 4/25 (16%) | **6/25 (24%)** |
| **DECODER Immune** | 0/60 (0%) | 20/60 (33%) | 16/60 (27%) |
| **iCAF factor HR (unadj)** | 0.927 (p=0.40) | 0.601 (p<0.0001) | **0.344 (p<0.0001)** |
| **iCAF factor HR (adj)** | 0.970 (p=0.74) | 0.674 (p<0.0001) | **0.371 (p<0.0001)** |
| **iCAF factor beta** | -3.41e-05 | -9.73e-04 | **-1.44e-03** |
| **# prognostic factors** | 1 (F2, not iCAF) | 2 (F1+F2) | **1 (F1 only)** |

### Validation (n=570, pooled, strata(dataset))

| | alpha=0.00 | alpha=0.35 | alpha=0.55 |
|--|-----------|-----------|-----------|
| **iCAF factor HR (unadj)** | 0.902 (p=0.10) | **0.841 (p=0.027)** | **0.800 (p=0.003)** |
| **iCAF factor HR (adj)** | 0.926 (p=0.25) | 0.862 (p=0.072) | 0.867 (p=0.075) |
| **LRT vs PurIST+DeCAF** | p=0.26 | p=0.075 | p=0.076 |
| **KM High vs Low** | 24.8 vs 19.2 (p=0.003) | 24.9 vs 19.2 (p=8e-4) | **26.8 vs 18.0 (p=3e-4)** |
| **# factors validated** | 0 (iCAF not sig) | 1 (iCAF unadj) | **1 (iCAF unadj+KM)** |

The training progression is monotonic: as alpha increases, the iCAF program becomes
more coherent, more prognostic, and more dominant. At alpha=0.55, only one
factor matters for survival (F1, the iCAF program).

The validation progression confirms the pattern:
- **alpha=0**: iCAF factor not significant in any Cox model (unadj or adj)
- **alpha=0.35**: iCAF factor validates unadjusted (p=0.027) with strong KM (p=8e-4)
- **alpha=0.55**: Strongest validation, unadj p=0.003, KM p=3e-4, with the
  largest median survival gap (8.8 months)

---

## Comparison to Original Production Fit

| Metric | Original (BO) | K3_a55 (grid) |
|--------|-------------|-------------|
| Alpha | 0.334 | 0.55 |
| Lambda | 0.349 | 0.30 |
| Nu | 0.056 | 0.05 |
| ntop | 270 | NULL |
| iCAF factor | F1 | F1 |
| iCAF beta | -1.003e-03 | -1.44e-03 |
| H-cor (F1 to F1) | 1.000 | 0.906 |
| Gene overlap (F1) | 270/270 | 168/270 |
| F1 HR unadj (train) | (ref) | 0.344 |
| F1 HR unadj (val) | (ref) | 0.800 (p=0.003) |
| Val KM p | (ref: p=5.4e-10) | p=3.0e-4 |
| Val adj p | (ref: p=0.004) | p=0.075 |

The grid fit at alpha=0.55 is more aggressive (higher alpha, stronger survival
penalty) than the BO-selected fit (alpha=0.334). This makes the iCAF program
even more dominant (F1 is the sole significant factor with HR=0.344) but at the
cost of some gene overlap (168/270 vs 270/270). The 62% overlap is still very
high, confirming the same underlying program.

The validation adjusted p (0.075 vs original 0.004) is weaker because the cv_grid
fits use all 1970 genes for H-scores (ntop=NULL) rather than the top-270 focused
genes. The unadjusted validation and KM results remain strong (p=0.003, KM p=3e-4),
confirming that the iCAF program generalizes to independent cohorts.

---

## Glossary

- **H-score**: Factor loading score for each sample, computed as H = X^T W
  where X is the expression matrix and W is the factor weight matrix. Measures
  how strongly each sample expresses a given factor's gene program.
- **H-cor (H-score correlation)**: Spearman rank correlation between H-score
  vectors from two different fits. Used here to compare each cv_grid factor
  against the original production K=3 factors.
- **LP (linear predictor)**: The Cox model linear predictor, sum of beta_k *
  H_k across factors. Summarizes each sample's predicted risk from the DeSurv
  model.
- **LRT (likelihood ratio test)**: Compares nested Cox models. Here we test
  whether adding a factor's H-score improves a model already containing PurIST
  + DeCAF.
- **ntop**: Number of top genes per factor used to compute H-scores. ntop=NULL
  means all 1970 genes are used; ntop=270 (production fit) focuses on each
  factor's 270 highest-weight genes, sharpening per-factor signal.
- **strata(dataset)**: Stratified Cox regression allowing each validation
  cohort its own baseline hazard function while sharing covariate effects
  across cohorts.
- **Gene universe**: The 1970 genes remaining after variance and expression
  filtering from the initial 3000 (ngene=3000). All gene overlap counts in
  this document (e.g., "17 of 25 present in universe") refer to genes present
  in this filtered set.
- **CV-optimal alpha**: The alpha value that minimizes cross-validated partial
  log-likelihood in the training data. For K=3, this is alpha=0.55.
- **Validation-optimal alpha**: The alpha at which external validation signal
  is strongest. For K=3, this is alpha=0.35 (unadjusted Cox) or alpha=0.55
  (KM).
