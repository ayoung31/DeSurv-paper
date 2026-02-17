# K=2 DeSurv Factor Analysis: CV Grid Fits

**For discussion with Jen Jen Yeh**
**Date: 2026-02-17**
**Companion to: 2026-02-17-jjy-talking-points-k3-grid.md**
**Data source: cv_grid exhaustive search (fixed lambda=0.3, nu=0.05, 100 inits)**

**Note for Amber:** This document uses the cv_grid exhaustive search fits, not
the BO-selected production fit. By holding lambda, nu, and ntop constant, we get
a controlled comparison across K values. See the glossary at the end for term
definitions. The analysis scripts are in `inst/cv_grid_training_analysis.R` (training)
and `inst/cv_grid_validation_analysis.R` (validation). Run from the repo root.

## Key Change from Previous K=2 Analysis

The previous K=2 analysis used K=3's BO-selected hyperparameters (alpha=0.33,
lambda=0.35, nu=0.056) with only K changed to 2. This made it difficult to
attribute differences to K alone, since the hyperparameters were optimized for K=3.

This document uses cv_grid fits where **all hyperparameters except K and alpha
are identical** to the K=3 and K=5 grid fits (lambda=0.3, nu=0.05, ntop=NULL,
100 initializations). Three K=2 fits are examined:

- **K2_a0** (alpha=0.00): Standard NMF baseline
- **K2_a35** (alpha=0.35): Validation-optimal alpha
- **K2_a55** (alpha=0.55): Matched to K=3 CV-optimal alpha

**Branch hashes:**
- K2_a0: `cv_grid_fit_083dae3cd3bba1c5`
- K2_a35: `cv_grid_fit_e675ee353d14ef89`
- K2_a55: `cv_grid_fit_7ced5aaef3b594d7`

**Training cohort:** TCGA + CPTAC (n=273)
**Validation cohort:** 4 independent datasets (n=570 pooled):
Dijk (n=90), Moffitt GEO array (n=123), Puleo array (n=288), PACA-AU (n=69).
All validation Cox models use `strata(dataset)`.

**iCAF reference**: When we refer to the "iCAF program," we mean the factor
that best matches the original K=3 Factor 1. That factor was characterized
by its overlap with the **Elyada et al. (2019) iCAF 35-gene signature** from
scRNA-seq of PDAC CAFs (25 of 35 genes present in the 1970-gene universe).

---

## Summary: K=2 Partially Recovers the iCAF Program But Cannot Isolate It

At K=2, there are only two factors, not enough capacity to separate the iCAF
program from the tumor subtype axis. The result is that one factor absorbs
iCAF + Classical tumor genes while the other absorbs Basal-like tumor genes.
The iCAF signal is present but diluted.

**The K=2 iCAF factor does not validate independently.** While the
training signal is strong (HR=0.634, p<0.0001 at alpha=0.55), the unadjusted
validation is borderline (HR=0.876, p=0.088) and the adjusted validation is not
significant (HR=0.891, p=0.16). This is because K=2's iCAF factor is entangled
with the PurIST axis. Once PurIST is included as a covariate, the residual
iCAF signal is too dilute to reach significance.

| | alpha=0.00 | alpha=0.35 | alpha=0.55 |
|--|-----------|-----------|-----------|
| **Best iCAF factor** | F2 | F2 | F2 |
| **H-cor with orig F1** | 0.108 | 0.518 | **0.664** |
| **Top-270 overlap** | 42/270 | 61/270 | **92/270** |
| **Elyada iCAF overlap** | 6/25 (24%) | **7/25 (28%)** | 3/25 (12%) |
| **DECODER Immune** | 4/60 (7%) | 6/60 (10%) | **10/60 (17%)** |
| **Train HR (unadj)** | 0.907 (p=0.20) | **0.728 (p<0.0001)** | **0.634 (p<0.0001)** |
| **Val HR (unadj)** | 1.018 (p=0.81) | 0.918 (p=0.28) | 0.876 (p=0.088) |
| **Val KM p** | 0.13 | 0.051 | **0.012** |
| **# prognostic factors** | 0 | 2 (both) | 2 (both) |

---

## K2_a0: Standard NMF (alpha=0.00)

### Factor identity

H-score correlation with original K=3:

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.056  | 0.305   | 0.124  |
| F2     | 0.108   | **0.589** | -0.283 |

Neither factor matches original F1. F2 best correlates with original F2
(stroma, rho=0.589). The iCAF program is absent at K=2 with no survival penalty.

### Reference signature overlaps (F2, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 6/25 | 24.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | 4/17 | 23.5% |
| DECODER Immune (60) | 4/60 | 6.7% |
| DECODER ActivatedStroma (108) | 18/108 | 16.7% |
| DECODER BasalTumor (139) | 12/139 | 8.6% |
| DECODER ClassicalTumor (166) | 17/166 | 10.2% |

F2 is a diffuse mix with no dominant biological identity. Moderate Elyada iCAF
overlap (24%) but this is deceiving because the Elyada iCAF genes are broadly
stroma-associated and overlap with many general stroma factors.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.030 | 0.73 | 1.131 | 0.14 |
| F2 (best iCAF) | 0.907 | 0.20 | 1.018 | 0.81 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.928 | 0.41 | 1.004 | 0.96 |
| F2 (best iCAF) | 0.900 | 0.21 | 0.953 | 0.53 |

**Neither factor is prognostic** in training or validation. Standard NMF at K=2
captures no survival signal.

**LRT adding F2 (iCAF) to PurIST + DeCAF:** Training p=0.22, Validation p=0.54

**Median-split KM for F2 (iCAF):**
- Training (n=273): High 23.2 vs Low 17.7 months, p=0.15
- Validation (n=570): High 22.5 vs Low 21.0 months, p=0.13

### Beta structure

F1=2.68e-04, F2=-2.79e-04. Both near zero (no survival penalty).

---

## K2_a35: Validation-Optimal Alpha (alpha=0.35)

### Factor identity

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.435  | **0.548** | -0.100 |
| **F2** | **0.518** | 0.258  | -0.068 |

F2 now correlates meaningfully with original F1 (rho=0.518). F1 aligns with
original F2 (stroma, rho=0.548). The survival penalty is pulling iCAF-related
genes toward F2.

### Reference signature overlaps (F2, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | **7/25** | **28.0%** |
| Elyada myCAF (15) | 3/15 | 20.0% |
| SCISSORS iCAF (17) | 4/17 | 23.5% |
| DECODER Immune (60) | 6/60 | 10.0% |
| DECODER ActivatedStroma (108) | 13/108 | 12.0% |
| DECODER BasalTumor (139) | 12/139 | 8.6% |
| DECODER ClassicalTumor (166) | 31/166 | 18.7% |

Elyada iCAF overlap peaks at 28%, the highest across all K=2 fits. But this
factor also contains 19% Classical Tumor genes, reflecting the K=2 limitation:
with only 2 factors, iCAF and Classical tumor signals are mixed.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1** | **1.330** | **0.005** | **1.222** | **0.009** |
| **F2 (iCAF)** | **0.728** | **<0.0001** | 0.918 | 0.28 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.124 | 0.28 | 1.052 | 0.53 |
| **F2 (iCAF)** | **0.760** | **0.001** | 0.906 | 0.24 |

F2 (iCAF-like) is strongly prognostic in training (unadj HR=0.728, adj HR=0.760)
but **does not replicate in validation** (unadj p=0.28, adj p=0.24). This is
the K=2 problem: the iCAF signal is entangled with the PurIST axis, and once
the larger validation samples provide more precise PurIST estimates, the
residual F2 signal is absorbed.

**LRT adding F2 (iCAF) to PurIST + DeCAF:** Training p=0.002, Validation p=0.24

**Median-split KM for F2 (iCAF):**
- Training (n=273): High 27.0 vs Low 15.3 months, **p=6.9e-6**
- Validation (n=570): High 23.7 vs Low 20.1 months, p=0.051 (borderline)

### Beta structure

F1=**8.54e-04**, F2=**-1.10e-03** (dominant)

F2 (protective, iCAF-like) has the larger beta, consistent with the original
K=3 pattern where the iCAF factor dominates.

---

## K2_a55: Matched to K=3 CV-Optimal Alpha (alpha=0.55)

### Factor identity

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.461  | **0.520** | -0.065 |
| **F2** | **0.664** | 0.207  | -0.102 |

F2 reaches its highest correlation with original F1 (rho=0.664) at this alpha.
But notably F1 also anti-correlates with original F1 (rho=-0.461), suggesting
the two K=2 factors capture opposite ends of the *same* iCAF axis.

### Reference signature overlaps (F2, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 3/25 | 12.0% |
| Elyada myCAF (15) | 4/15 | 26.7% |
| SCISSORS iCAF (17) | 3/17 | 17.6% |
| DECODER Immune (60) | **10/60** | **16.7%** |
| DECODER ActivatedStroma (108) | 8/108 | 7.4% |
| DECODER NormalStroma (15) | **5/15** | **33.3%** |
| DECODER BasalTumor (139) | 7/139 | 5.0% |
| DECODER ClassicalTumor (166) | 28/166 | 16.9% |

Interesting shift: Elyada iCAF overlap drops to 12% while DECODER Immune rises
to 17% and NormalStroma to 33%. The factor is becoming more immune/stroma-heavy
but less specifically iCAF. This suggests that at high alpha with K=2, the
survival penalty pushes toward immune/stroma genes broadly rather than the
specific iCAF program.

### DECODER compartment breakdown for both factors

| Factor | Immune | Basal Tumor | Classical Tumor | Act. Stroma | Norm. Stroma |
|--------|--------|-------------|-----------------|-------------|-------------|
| F1 | 1/60 (2%) | **26/139 (19%)** | 24/166 (14%) | 19/108 (18%) | 0/15 (0%) |
| **F2** | **10/60 (17%)** | 7/139 (5%) | 28/166 (17%) | 8/108 (7%) | **5/15 (33%)** |

F2 is Immune + NormalStroma + Classical. F1 is Basal + ActivatedStroma. At K=2,
the model creates a single tumor-stroma contrast axis rather than isolating
the iCAF program from the tumor subtype axis.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1** | **1.377** | **0.002** | **1.237** | **0.005** |
| **F2 (iCAF)** | **0.634** | **<0.0001** | 0.876 | 0.088 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.164 | 0.17 | 1.058 | 0.49 |
| **F2 (iCAF)** | **0.672** | **<0.0001** | 0.891 | 0.16 |

The K=2 pattern at alpha=0.55: extremely strong training signal (F2 adj HR=0.672,
p<0.0001) but **weak validation** (F2 adj HR=0.891, p=0.16). The unadjusted
validation is also borderline (p=0.088). Compare to K=3 at alpha=0.55 where
unadjusted validation reaches p=0.003.

**LRT adding F2 (iCAF) to PurIST + DeCAF:**
- Training p = **1.8e-5**
- Validation p = 0.16

**Median-split KM for F2 (iCAF):**
- Training (n=273): High 33.4 vs Low 15.1 months, **p=1.8e-8**
- Validation (n=570): High 24.1 vs Low 20.1 months, **p=0.012**

The KM validation p=0.012 shows there IS an iCAF signal in K=2's Factor 2, but
the 4.0-month gap (24.1 vs 20.1) is much smaller than K=3's 8.8-month gap
(26.8 vs 18.0). The iCAF program is diluted by mixing with the Classical tumor axis.

### Beta structure

F1=**7.00e-04**, F2=**-1.13e-03** (dominant)

---

## The K=2 Limitation: Why Two Factors Cannot Isolate iCAF

At K=2, the model must capture all biological variation in just two factors.
The result:

1. **One factor mixes iCAF + Classical tumor + Normal stroma + Immune**
   (the "good prognosis" end)
2. **One factor mixes Basal-like tumor + Activated stroma**
   (the "bad prognosis" end)

This is essentially the **PurIST/DeCAF axis** (Basal/proCAF vs Classical/restCAF)
with iCAF genes layered on top of the Classical end. The iCAF program cannot be
separated from the Classical tumor program because there aren't enough factors.

### Validation evidence of mixing

The validation data makes this limitation concrete. At K=2 alpha=0.55:

| Metric | K=2 a=0.55 | K=3 a=0.55 |
|--------|-----------|-----------|
| iCAF factor HR (train unadj) | 0.634 | **0.344** |
| iCAF factor HR (val unadj) | 0.876 (p=0.088) | **0.800 (p=0.003)** |
| iCAF factor HR (val adj) | 0.891 (p=0.16) | 0.867 (p=0.075) |
| Val KM gap | 4.0 months | **8.8 months** |
| Val KM p | 0.012 | **3.0e-4** |
| Non-iCAF factors val significant? | **Yes (F1 p=0.005)** | No |
| Val LRT vs PurIST+DeCAF | 0.16 | 0.076 |
| iCAF beta dominance | 1.6x F1 | **4.6x F2, 42x F3** |

At K=2, both factors contribute to survival. At K=3, only the iCAF factor
matters. This means K=3's iCAF program is a purer survival signal; it is not
entangled with the tumor subtype axis. The validation data confirms this:
K=3's iCAF factor reaches significance independently, while K=2's does not.

---

## Glossary

- **H-score**: Factor loading score for each sample, computed as H = X^T W.
  Measures how strongly each sample expresses a given factor's gene program.
- **H-cor**: Spearman rank correlation between H-score vectors from two fits.
  Used to compare cv_grid factors against original production K=3 factors.
- **LP**: Cox model linear predictor, sum of beta_k * H_k across factors.
- **LRT**: Likelihood ratio test comparing nested Cox models. Tests whether
  adding a factor improves a model already containing PurIST + DeCAF.
- **ntop**: Number of top genes per factor for H-scores. ntop=NULL = all 1970
  genes; ntop=270 (production fit) = top 270 genes per factor.
- **strata(dataset)**: Stratified Cox allowing each validation cohort its own
  baseline hazard while sharing covariate effects.
- **Gene universe**: The 1970 genes remaining after filtering from the initial
  3000. All overlap counts (e.g., "25 in universe") refer to this set.
- **CV-optimal alpha**: Alpha minimizing cross-validated partial log-likelihood.
- **Validation-optimal alpha**: Alpha with strongest external validation signal.
