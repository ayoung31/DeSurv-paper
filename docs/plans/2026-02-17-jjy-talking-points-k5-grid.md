# K=5 DeSurv Factor Analysis: CV Grid Fits

**For discussion with Jen Jen Yeh**
**Date: 2026-02-17**
**Companion to: 2026-02-17-jjy-talking-points-k3-grid.md,
2026-02-17-jjy-talking-points-k2-grid.md**
**Data source: cv_grid exhaustive search (fixed lambda=0.3, nu=0.05, 100 inits)**

**Note for Amber:** This document uses the cv_grid exhaustive search fits, not
the BO-selected production fit. By holding lambda, nu, and ntop constant, we get
a controlled comparison across K values. See the glossary at the end for term
definitions. The analysis scripts are in `inst/cv_grid_training_analysis.R` (training)
and `inst/cv_grid_validation_analysis.R` (validation). Run from the repo root.

## Key Change from Previous K=5 Analysis

The previous K=5 analysis used the production HPC fit with BO-selected
hyperparameters (alpha=0.362, lambda=0.314, **nu=0.005**). The critical
parameter was nu=0.005, 11x smaller than K=3's nu=0.056, meaning the BO
had essentially turned off survival weighting at K=5, making it converge
toward standard NMF. That made it difficult to attribute differences to K
alone.

This document uses cv_grid fits where **nu=0.05 is fixed across all K values**.
The BO's conclusion that "K=5 needs less survival weighting" is set aside;
we ask instead: at a fixed survival penalty, what happens to the iCAF program
at K=5? Three fits are examined:

- **K5_a0** (alpha=0.00): Standard NMF baseline
- **K5_a25** (alpha=0.25): Validation-optimal alpha for K=5
- **K5_a55** (alpha=0.55): Matched to K=3 CV-optimal alpha

**Branch hashes:**
- K5_a0: `cv_grid_fit_128a06d24f80875e`
- K5_a25: `cv_grid_fit_5c725504aa99734e`
- K5_a55: `cv_grid_fit_ed2d68b41575585f`

**Training cohort:** TCGA + CPTAC (n=273)
**Validation cohort:** 4 independent datasets (n=570 pooled):
Dijk (n=90), Moffitt GEO array (n=123), Puleo array (n=288), PACA-AU (n=69).
All validation Cox models use `strata(dataset)`.

**iCAF reference**: The "iCAF program" refers to the factor best matching the
original K=3 Factor 1, characterized by overlap with the **Elyada et al. (2019)
iCAF 35-gene signature** (25 of 35 present in the 1970-gene universe).

---

## Summary: K=5 Recovers iCAF But Fragments It Among Competing Factors

At K=5, with survival penalty, the iCAF program emerges (rho=0.737 with original
F1 at alpha=0.55) but is accompanied by additional prognostic factors. Unlike
K=3 where the iCAF factor is the sole survival signal, K=5 spreads prognostic
weight across 2-3 factors.

**The validation reveals an interesting pattern.** K5_a25's iCAF factor has the
strongest unadjusted validation (HR=0.797, p=0.0005) and the strongest KM separation
(25.6 vs 18.0 months, p=4.4e-6) of any fit. But it does not survive PurIST+DeCAF
adjustment (p=0.069), and its H-correlation with the original iCAF program is
only 0.449. It is a "purer" immune factor that shares survival signal with
PurIST/DeCAF rather than an independent iCAF program.

| | alpha=0.00 | alpha=0.25 | alpha=0.55 |
|--|-----------|-----------|-----------|
| **Best iCAF factor** | F1 | F4 | F4 |
| **H-cor with orig F1** | 0.313 | 0.449 | **0.737** |
| **Top-270 overlap** | 47/270 | 121/270 | **146/270** |
| **Elyada iCAF overlap** | 3/25 (12%) | 6/25 (24%) | 6/25 (24%) |
| **DECODER Immune** | 1/60 (2%) | **27/60 (45%)** | 15/60 (25%) |
| **Train HR (unadj)** | 0.969 (p=0.72) | **0.652 (p<0.0001)** | **0.529 (p<0.0001)** |
| **Val HR (unadj)** | 0.971 (p=0.65) | **0.797 (p=0.0005)** | 0.884 (p=0.13) |
| **Val KM p** | 0.59 | **4.4e-6** | **0.005** |
| **# prognostic factors** | 1 (not iCAF) | 2 | **3** |

---

## K5_a0: Standard NMF (alpha=0.00)

### Factor identity

H-score correlation with original K=3:

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | 0.313   | **-0.961** | **0.925** |
| F2     | -0.490  | 0.463   | -0.141 |
| F3     | 0.112   | 0.274   | -0.243 |
| F4     | -0.133  | **0.964** | **-0.884** |
| F5     | 0.292   | -0.471  | 0.554 |

The dominant correlations are with original F2 and F3 (the stroma/basal axis),
not F1 (iCAF). F1 and F4 are anti-correlated versions of the stroma axis
(F1 <-> -origF2, F4 <-> +origF2). No factor strongly matches the iCAF program.

### Reference signature overlaps (F1 as best iCAF match, top-270)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 3/25 | 12.0% |
| SCISSORS iCAF (17) | 0/17 | 0.0% |
| DECODER Immune (60) | 1/60 | 1.7% |
| DECODER BasalTumor (139) | 19/139 | 13.7% |
| DECODER ClassicalTumor (166) | 35/166 | 21.1% |

F1 is dominated by Classical Tumor genes (21%), not iCAF. Where do the Immune
genes go? F4 has 33/60 DECODER Immune genes (55%), but F4 anti-correlates
with original F1 (rho=-0.133). At K=5 with alpha=0, the immune program forms
its own factor but is not linked to the iCAF program.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 (best iCAF) | 0.969 | 0.72 | 0.971 | 0.65 |
| **F2** | **1.320** | **0.004** | **1.301** | **<0.0001** |
| F3 | 0.912 | 0.29 | 0.961 | 0.45 |
| F4 | 0.954 | 0.58 | 0.940 | 0.32 |
| F5 | 0.898 | 0.21 | 0.901 | 0.13 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 (best iCAF) | 0.972 | 0.75 | 0.962 | 0.57 |
| F2 | 0.952 | 0.68 | 1.138 | 0.068 |
| F3 | 1.048 | 0.61 | 1.022 | 0.68 |
| F4 | 0.986 | 0.87 | 0.941 | 0.34 |
| F5 | 0.886 | 0.20 | 0.839 | 0.022 |

**The iCAF factor is not prognostic at alpha=0.** Only F2 is significant, and
F2 captures the basal-like tumor axis, not iCAF. In validation, F2 replicates
strongly (p<0.0001), confirming that standard NMF captures tumor subtype but
not the iCAF program.

**LRT adding F1 (iCAF) to PurIST + DeCAF:** Training p=0.75, Validation p=0.57

**Median-split KM for F1 (iCAF):**
- Training (n=273): High 20.5 vs Low 23.4 months, p=0.79 (reversed!)
- Validation (n=570): High 22.5 vs Low 21.4 months, p=0.59

### Beta structure

F1=0, F2=1.34e-04, F3=-1.37e-05, F4=-4.28e-05, F5=-1.13e-04

All tiny (no survival penalty -> near-zero betas).

---

## K5_a25: Validation-Optimal Alpha (alpha=0.25)

### Factor identity

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.246  | **0.992** | -0.857 |
| F2     | -0.345  | 0.281   | 0.051  |
| F3     | 0.225   | -0.907  | **0.936** |
| **F4** | **0.449** | 0.484  | -0.512 |
| F5     | 0.078   | -0.213  | 0.538  |

F4 is the best iCAF match (rho=0.449). F1 perfectly matches original F2
(rho=0.992) and F3 matches original F3 (rho=0.936). The stroma/basal axis is
well recovered; the iCAF program is partially forming in F4.

### Reference signature overlaps (F4, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 6/25 | 24.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | **7/17** | **41.2%** |
| SCISSORS myCAF (16) | 1/16 | 6.2% |
| DeCAF restCAF (8) | 2/8 | 25.0% |
| DECODER Immune (60) | **27/60** | **45.0%** |
| DECODER ActivatedStroma (108) | 5/108 | 4.6% |
| DECODER NormalStroma (15) | 5/15 | 33.3% |
| DECODER BasalTumor (139) | 0/139 | 0.0% |
| DECODER ClassicalTumor (166) | 16/166 | 9.6% |

F4 is remarkably pure: **45% DECODER Immune, 0% Basal Tumor**. The SCISSORS
iCAF overlap (41%) and DeCAF restCAF overlap (25%) are the highest across all
fits. This is arguably the purest immune/iCAF factor in the entire grid.

### Where do the iCAF-adjacent genes go?

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 6/25 | 23/60 | 2/139 | Stroma/Immune |
| F2 | 3/25 | 0/60 | 23/139 | Basal Tumor |
| F3 | 0/25 | 0/60 | 43/139 | Basal Tumor (dominant) |
| **F4** | **6/25** | **27/60** | **0/139** | **iCAF/Immune (purest)** |
| F5 | 7/25 | 4/60 | 13/139 | Mixed |

But note: F1 also has 23/60 DECODER Immune genes and 6/25 Elyada iCAF genes.
The immune/iCAF signal is split between F1 and F4. Both capture parts of it.
At K=5, the model has enough capacity to separate Immune from iCAF, but this
fragmentation means no single factor captures the full iCAF program as K=3 does.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.011 | 0.90 | 1.055 | 0.38 |
| **F2** | **1.336** | **0.002** | **1.224** | **0.003** |
| F3 | 1.009 | 0.92 | 1.022 | 0.73 |
| **F4 (iCAF)** | **0.652** | **<0.0001** | **0.797** | **0.0005** |
| F5 | 1.009 | 0.91 | 1.173 | 0.052 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.003 | 0.97 | 1.008 | 0.89 |
| **F2** | **1.256** | **0.017** | 1.093 | 0.21 |
| F3 | 0.955 | 0.61 | 0.941 | 0.38 |
| **F4 (iCAF)** | **0.767** | **0.007** | 0.879 | 0.069 |
| F5 | 0.860 | 0.10 | 1.037 | 0.67 |

F4 (iCAF) validates unadjusted (HR=0.797, p=0.0005) but is borderline adjusted
(HR=0.879, p=0.069). Two factors are independently prognostic in training
(F4 + F2), and the iCAF factor's training signal (HR=0.652, p<0.0001) matches
K=3's pattern. In validation, F4 is the strongest single-factor prognostic signal
in the entire grid (unadj p=0.0005).

**LRT adding F4 (iCAF) to PurIST + DeCAF:** Training p=0.008, Validation p=0.072

**Median-split KM for F4 (iCAF):**
- Training (n=273): High 24.8 vs Low 15.7 months, **p=6.9e-5**
- Validation (n=570): High 25.6 vs Low 18.0 months, **p=4.4e-6**

The validation KM is the strongest in the entire grid, with a 7.6-month survival
gap, p=4.4e-6 in 570 patients. But this factor only has rho=0.449 with the
original iCAF program. It is a purer immune factor, not the full iCAF program.

### Beta structure

F1=6.02e-05, F2=**5.18e-04**, F3=-6.72e-05, F4=**-6.77e-04** (dominant), F5=-9.62e-05

F4 (iCAF) has the dominant beta, but F2 (Basal) also has substantial weight.

---

## K5_a55: Matched to K=3 CV-Optimal Alpha (alpha=0.55)

### Factor identity

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.081  | -0.646  | **0.920** |
| F2     | -0.223  | **0.886** | -0.614 |
| F3     | -0.054  | -0.256  | 0.630  |
| **F4** | **0.737** | -0.652 | 0.654  |
| F5     | -0.089  | **0.920** | -0.813 |

F4 matches original F1 well (rho=0.737). F1 matches original F3 (rho=0.920).
F2 and F5 both match original F2 (rho=0.886 and 0.920). The stroma axis has
split into two factors at K=5.

### Reference signature overlaps (F4, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 6/25 | 24.0% |
| Elyada myCAF (15) | 0/15 | 0.0% |
| SCISSORS iCAF (17) | **7/17** | **41.2%** |
| SCISSORS myCAF (16) | 2/16 | 12.5% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | **2/4** | **50.0%** |
| DECODER Immune (60) | 15/60 | 25.0% |
| DECODER ActivatedStroma (108) | 4/108 | 3.7% |
| DECODER NormalStroma (15) | **6/15** | **40.0%** |
| DECODER BasalTumor (139) | 15/139 | 10.8% |
| DECODER ClassicalTumor (166) | 42/166 | 25.3% |

F4 maintains SCISSORS iCAF overlap (41%) and gains NormalStroma (40%). But
unlike K5_a25 where F4 had 45% DECODER Immune overlap with zero Basal Tumor,
here F4 has only 25% Immune and picks up 11% Basal Tumor and 25% Classical
Tumor. The factor is becoming less pure; it is absorbing more tumor genes as
alpha increases.

### DECODER compartment breakdown for ALL factors

| Factor | Immune | Basal | Classical | Act.Stroma | Norm.Stroma | Identity |
|--------|--------|-------|-----------|------------|-------------|----------|
| F1 | 0/60 | **39/139** | 12/166 | 22/108 | 0/15 | Basal/Stroma |
| F2 | 5/60 | 12/139 | **36/166** | 10/108 | 1/15 | Classical Tumor |
| F3 | 5/60 | **31/139** | 13/166 | 10/108 | 0/15 | Mixed/Basal |
| **F4** | **15/60** | 15/139 | **42/166** | 4/108 | **6/15** | **iCAF/Mixed** |
| F5 | **24/60** | 1/139 | 18/166 | 2/108 | 3/15 | Immune/Classical |

At K=5 alpha=0.55, the Immune genes are split: F4 gets 15/60, F5 gets 24/60.
Compare to K=3 alpha=0.55 where F1 (iCAF) gets 16/60 and F3 gets 16/60, a
similar split, but at K=3 the iCAF factor combines Immune + NormalStroma into
a single coherent program, while K=5's F4 dilutes this with Classical Tumor
genes (42/166 = 25%).

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1** | **1.318** | **0.003** | **1.166** | **0.039** |
| F2 | 1.008 | 0.92 | 1.119 | 0.13 |
| **F3** | **1.218** | **0.029** | 1.166 | 0.054 |
| **F4 (iCAF)** | **0.529** | **<0.0001** | 0.884 | 0.13 |
| F5 | 0.920 | 0.32 | 0.955 | 0.50 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1** | 1.178 | 0.091 | 1.001 | 0.99 |
| F2 | 0.941 | 0.51 | 1.014 | 0.85 |
| F3 | 1.047 | 0.64 | 1.014 | 0.87 |
| **F4 (iCAF)** | **0.562** | **<0.0001** | 0.888 | 0.17 |
| F5 | 0.986 | 0.88 | 0.962 | 0.58 |

**Three training factors are prognostic but only F4 survives PurIST+DeCAF
adjustment in training.** In validation, the iCAF factor DOES NOT reach
significance, neither unadjusted (p=0.13) nor adjusted (p=0.17). This is a
critical failure: K=5 at high alpha overtrains the survival signal. The
extremely strong training HR (0.529) does not replicate (val HR=0.884).

**LRT adding F4 (iCAF) to PurIST + DeCAF:** Training p=4.3e-9, Validation p=0.17

**Median-split KM for F4 (iCAF):**
- Training (n=273): High 24.6 vs Low 14.2 months, **p=3.0e-7**
- Validation (n=570): High 25.3 vs Low 19.2 months, **p=0.005**

The KM still validates (6.1-month gap, p=0.005) despite the Cox model failure,
suggesting the iCAF signal exists but is spread across multiple factors that
individually don't reach significance.

### Beta structure

F1=7.25e-04, F2=-1.12e-04, F3=4.15e-04, F4=**-1.53e-03** (dominant), F5=-4.42e-05

F4 (iCAF) dominates the LP, but F1 and F3 also have substantial betas.

---

## The K=5 Limitation: Fragmentation and Validation Failure

### The fragmentation problem, confirmed by validation

At K=5, the model has enough capacity to separate biological programs into
finer components. But this creates two problems visible in validation:

1. **Immune/iCAF splitting**: The immune component that is part of the K=3 iCAF
   program gets split across multiple K=5 factors (F4 gets 15/60 Immune genes,
   F5 gets 24/60, F1 gets 0). At K=3, F1 captures the combined iCAF+Immune
   signal coherently.

2. **Validation instability across alpha**: K5_a25's iCAF factor validates
   strongly (unadj p=0.0005), but K5_a55's does not (unadj p=0.13). At K=3,
   both alpha values validate (p=0.027 and p=0.003). The K=5 validation signal
   is fragile. It depends on the specific alpha, whereas K=3 is robust.

3. **Multiple survival axes hurt replication**: At K=5 alpha=0.55, 3 of 5
   factors are independently prognostic in training, meaning the survival signal
   is distributed. None of the individual factors reliably replicate in validation.

### Comparison across K at alpha=0.55 (with validation)

| Metric | K=2 | K=3 | K=5 |
|--------|-----|-----|-----|
| iCAF H-cor with orig F1 | 0.664 | **0.906** | 0.737 |
| iCAF gene overlap | 92/270 | **168/270** | 146/270 |
| iCAF Train HR (unadj) | 0.634 | **0.344** | 0.529 |
| iCAF Val HR (unadj) | 0.876 (p=0.088) | **0.800 (p=0.003)** | 0.884 (p=0.13) |
| iCAF Val HR (adj) | 0.891 (p=0.16) | 0.867 (p=0.075) | 0.888 (p=0.17) |
| Val KM gap | 4.0 mo (p=0.012) | **8.8 mo (p=3e-4)** | 6.1 mo (p=0.005) |
| # prognostic factors | 2 | **1** | 3 |
| iCAF beta dominance | 1.6x | **4.6x** | 2.2x |

K=3 wins on validation: strongest unadjusted validation (p=0.003), largest KM
gap (8.8 months), and the only K where the iCAF factor is the sole prognostic
signal. K=5's stronger training signal (HR=0.529) doesn't replicate, suggesting
overfitting at higher K with strong survival penalty.

### The "purest immune factor" at K=5_a25

Worth noting: K5_a25's F4 (45% DECODER Immune, 0% Basal Tumor) is arguably
the purest immune/iCAF factor in the entire grid, and it has the best unadjusted
validation KM (p=4.4e-6). But it only captures rho=0.449 of the original K=3 F1
program, because the original F1 was not purely immune; it was a mixed
iCAF + immune + normal stroma program. K=5 separates what K=3 combined, and
the original paper's claim is about the combined program, not just the immune
component.

The K5_a25 result does suggest that the immune component of Factor 1 is the
primary driver of the survival signal, and that this signal exists independently
of PurIST/DeCAF (since the unadjusted validation is highly significant). But
the adjusted validation borderline (p=0.069) indicates substantial overlap with
existing classifiers. The immune signal alone is not fully independent of PurIST.

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
