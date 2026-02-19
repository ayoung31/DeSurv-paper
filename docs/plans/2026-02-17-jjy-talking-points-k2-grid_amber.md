# K=2 DeSurv Factor Analysis: CV Grid Fits (Expanded Grid)

**For discussion with Jen Jen Yeh**
**Date: 2026-02-18 (amber revision)**
**Companion to: 2026-02-17-k-sensitivity-synthesis_amber.md**
**Data source: cv_grid exhaustive search (fixed lambda=0.3, nu=0.05, 100 inits)**

**Note for Amber:** This is the expanded-grid version. It covers 7 K=2 fits
(alpha=0, 0.25, 0.35, 0.55, 0.75, 0.85, 0.95) instead of the original 3.
All fits are loaded by parameter matching from the targets store (no hardcoded
branch hashes). Analysis scripts: `inst/cv_grid_training_analysis_amber.R` and
`inst/cv_grid_validation_analysis_amber.R`.

## Key Change from Previous K=2 Analysis

The previous analysis examined 3 alpha values (0, 0.35, 0.55). This version
extends to 7 alpha values (0, 0.25, 0.35, 0.55, 0.75, 0.85, 0.95) to test:

1. Does iCAF coherence (H-cor) keep improving beyond alpha=0.55?
2. Does validation signal eventually emerge at high alpha?
3. Is the K=2 mixing problem (iCAF entangled with PurIST) inherent to the rank
   or alpha-dependent?

**Training cohort:** TCGA + CPTAC (n=273)
**Validation cohort:** 4 independent datasets (n=570 pooled):
Dijk (n=90), Moffitt GEO array (n=123), Puleo array (n=288), PACA-AU (n=69).
All validation Cox models use `strata(dataset)`.

**iCAF reference**: The "iCAF program" = the factor best matching original K=3
Factor 1 (Elyada iCAF + DECODER Immune + Normal Stroma).

---

## Summary Table: K=2 Across 7 Alpha Values

| | a=0.00 | a=0.25 | a=0.35 | a=0.55 | a=0.75 | a=0.85 | a=0.95 |
|--|--------|--------|--------|--------|--------|--------|--------|
| **Best iCAF factor** | F2 | F1 | F2 | F2 | F2 | F1 | F2 |
| **H-cor with orig F1** | 0.108 | 0.409 | 0.518 | 0.664 | 0.674 | 0.671 | 0.728 |
| **Top-270 overlap** | 42/270 | 52/270 | 61/270 | 92/270 | 85/270 | 97/270 | 100/270 |
| **Elyada iCAF overlap** | 6/25 (24.0%) | 3/25 (12.0%) | 7/25 (28.0%) | 3/25 (12.0%) | 5/25 (20.0%) | 4/25 (16.0%) | 5/25 (20.0%) |
| **DECODER Immune** | 4/60 (6.7%) | 5/60 (8.3%) | 6/60 (10.0%) | 10/60 (16.7%) | 5/60 (8.3%) | 2/60 (3.3%) | 6/60 (10.0%) |
| **Train HR (unadj)** | 0.907 | 0.786 | 0.728 | 0.634 | 0.619 | 0.607 | 0.544 |
| **Train HR (adj)** | 0.900 | 0.804 | 0.760 | 0.672 | 0.638 | 0.615 | 0.550 |
| **Val HR (unadj)** | 1.018 | 0.951 | 0.918 | 0.876 | 0.878 | 0.881 | 0.886 |
| **Val HR (adj)** | 0.953 | 0.920 | 0.906 | 0.891 | 0.883 | 0.883 | 0.892 |
| **Val p (unadj)** | 8.14e-01 | 5.35e-01 | 2.83e-01 | 8.75e-02 | 8.13e-02 | 1.19e-01 | 1.06e-01 |
| **Val p (adj)** | 5.35e-01 | 3.22e-01 | 2.36e-01 | 1.57e-01 | 1.12e-01 | 1.47e-01 | 1.45e-01 |
| **Train LRT p** | 2.22e-01 | 1.02e-02 | 1.66e-03 | 1.84e-05 | 7.87e-07 | 8.51e-08 | 1.81e-10 |
| **Val LRT p** | 5.37e-01 | 3.25e-01 | 2.40e-01 | 1.61e-01 | 1.16e-01 | 1.51e-01 | 1.49e-01 |
| **Val KM p** | 1.33e-01 | 1.27e-01 | 5.10e-02 | 1.17e-02 | 4.27e-03 | 2.40e-03 | 1.74e-02 |
| **Val KM gap (mo)** | 1.5 | 3.4 | 3.6 | 4.0 | 4.3 | 4.9 | 4.1 |
| **# prognostic factors** | 0 | 1 | 2 | 2 | 2 | 2 | 2 |

The training HR monotonically decreases from 0.907 to 0.544 as alpha rises,
yet the validation HR plateaus around 0.88 for alpha >= 0.55. The adjusted
validation p-value never reaches significance (best: 0.112 at alpha=0.75).
The LRT tells the same story: training p-values plummet to 1.8e-10 while
validation LRT p-values hover around 0.12-0.54, never significant. The KM
median-split p-values do reach significance at alpha >= 0.55, but the survival
gap is modest (4-5 months) compared to K=3's 8.8-month gap.

---

## Per-Alpha Analysis

### K2_a000: Standard NMF (alpha=0.00)

#### Factor identity

Best iCAF match: **F2** (rho=0.108). Neither factor meaningfully matches
original F1.

**H-score correlation with original K=3:**

| | Orig_F1 | Orig_F2 | Orig_F3 |
|--|---------|---------|---------|
| F1 | -0.056 | 0.305 | 0.124 |
| F2 | 0.108 | **0.589** | -0.283 |

Neither factor matches original F1. F2 best correlates with original F2
(stroma, rho=0.589). The iCAF program is absent at K=2 with no survival penalty.

#### Reference signature overlaps (F2, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 6/25 | 24.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | 4/17 | 23.5% |
| SCISSORS myCAF (16) | 3/16 | 18.8% |
| DeCAF restCAF (8) | 0/8 | 0.0% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 4/60 | 6.7% |
| DECODER ActStroma (108) | 18/108 | 16.7% |
| DECODER NormStroma (15) | 1/15 | 6.7% |
| DECODER BasalTumor (139) | 12/139 | 8.6% |
| DECODER ClassTumor (166) | 17/166 | 10.2% |

F2 is a diffuse mix with no dominant biological identity. The 24% Elyada iCAF
overlap is misleading because those genes overlap broadly with stroma programs.

#### DECODER compartment breakdown (both factors, top-270)

| Factor | Immune | BasalTumor | ClassTumor | ActStroma | NormStroma |
|--------|--------|------------|------------|-----------|------------|
| F1 | 1/60 (2%) | 15/139 (11%) | 17/166 (10%) | 18/108 (17%) | 1/15 (7%) |
| **F2** | 4/60 (7%) | 12/139 (9%) | 17/166 (10%) | 18/108 (17%) | 1/15 (7%) |

Both factors are diffuse across all compartments. No clear biological identity.

#### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.030 | 0.73 | 1.131 | 0.14 |
| F2 (iCAF) | 0.907 | 0.20 | 1.018 | 0.81 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.928 | 0.41 | 1.004 | 0.96 |
| F2 (iCAF) | 0.900 | 0.21 | 0.953 | 0.53 |

**Neither factor is prognostic** in training or validation. Standard NMF at K=2
captures no survival signal.

**LRT adding F2 (iCAF) to PurIST + DeCAF:** Training p=0.22, Validation p=0.54

**Median-split KM for F2 (iCAF):**
- Training (n=273): High 23.2 vs Low 17.7 months, p=0.15
- Validation (n=570): High 22.5 vs Low 21.0 months, p=0.13

#### Beta structure

F1=2.68e-04, F2=-2.79e-04. Both near zero (no survival penalty).

---

### K2_a025: alpha=0.25

#### Factor identity

Best iCAF match: **F1** (rho=0.409). The survival penalty begins pulling
iCAF-related genes into F1, which also partially correlates with original F2.

**H-score correlation with original K=3:**

| | Orig_F1 | Orig_F2 | Orig_F3 |
|--|---------|---------|---------|
| **F1** | **0.409** | 0.314 | -0.069 |
| F2 | -0.283 | **0.554** | -0.122 |

F1 starts to pick up original F1 signal (rho=0.409) while F2 aligns with
original F2 (stroma, rho=0.554). Notably, the iCAF factor is F1 here (not F2
as in most other alphas), reflecting the instability of factor ordering at
moderate alpha.

#### Reference signature overlaps (F1, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 3/25 | 12.0% |
| Elyada myCAF (15) | 1/15 | 6.7% |
| SCISSORS iCAF (17) | 4/17 | 23.5% |
| SCISSORS myCAF (16) | 0/16 | 0.0% |
| DeCAF restCAF (8) | 2/8 | 25.0% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 5/60 | 8.3% |
| DECODER ActStroma (108) | 8/108 | 7.4% |
| DECODER NormStroma (15) | 1/15 | 6.7% |
| DECODER BasalTumor (139) | 10/139 | 7.2% |
| DECODER ClassTumor (166) | 21/166 | 12.7% |

Low Elyada iCAF overlap (12%). The factor is picking up a mild Classical Tumor
(12.7%) and SCISSORS iCAF (23.5%) signal but remains poorly defined.

#### DECODER compartment breakdown (both factors, top-270)

| Factor | Immune | BasalTumor | ClassTumor | ActStroma | NormStroma |
|--------|--------|------------|------------|-----------|------------|
| **F1** | 5/60 (8%) | 10/139 (7%) | 21/166 (13%) | 8/108 (7%) | 1/15 (7%) |
| F2 | 3/60 (5%) | 25/139 (18%) | -- | -- | -- |

F1 is mildly Classical Tumor-enriched. F2 absorbs Basal Tumor genes (18%).
The compartment separation is beginning but incomplete.

#### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.786** | **7.18e-04** | 0.951 | 0.54 |
| F2 | 1.164 | 0.10 | **1.172** | **0.040** |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.804** | **0.007** | 0.920 | 0.32 |
| F2 | 1.025 | 0.80 | 1.026 | 0.75 |

F1 (iCAF) is prognostic in training (unadj HR=0.786, adj HR=0.804) but not
in validation (unadj p=0.54, adj p=0.32). Interestingly, F2 reaches unadjusted
validation significance (p=0.040) but loses it after adjustment (p=0.75),
suggesting its signal is collinear with PurIST/DeCAF.

**LRT adding F1 (iCAF) to PurIST + DeCAF:** Training p=0.010, Validation p=0.33

**Median-split KM for F1 (iCAF):**
- Training (n=273): High 25.4 vs Low 15.3 months, **p=5.84e-05**
- Validation (n=570): High 23.6 vs Low 20.2 months, p=0.13

#### Beta structure

F1=**-9.52e-04** (dominant, protective), F2=6.90e-04

F1 (iCAF) has the larger absolute beta, consistent with the protective iCAF
pattern. But the magnitude is still modest compared to higher-alpha fits.

---

### K2_a035: alpha=0.35

#### Factor identity

Best iCAF match: **F2** (rho=0.518). The factors now clearly separate: F1
aligns with original F2 (stroma) and F2 begins capturing the iCAF program.

**H-score correlation with original K=3:**

| | Orig_F1 | Orig_F2 | Orig_F3 |
|--|---------|---------|---------|
| F1 | -0.435 | **0.548** | -0.100 |
| **F2** | **0.518** | 0.258 | -0.068 |

F2 now correlates meaningfully with original F1 (rho=0.518). F1 aligns with
original F2 (stroma, rho=0.548). The survival penalty is pulling iCAF-related
genes toward F2.

#### Reference signature overlaps (F2, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | **7/25** | **28.0%** |
| Elyada myCAF (15) | 3/15 | 20.0% |
| SCISSORS iCAF (17) | 4/17 | 23.5% |
| SCISSORS myCAF (16) | 2/16 | 12.5% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 6/60 | 10.0% |
| DECODER ActStroma (108) | 13/108 | 12.0% |
| DECODER NormStroma (15) | 1/15 | 6.7% |
| DECODER BasalTumor (139) | 12/139 | 8.6% |
| DECODER ClassTumor (166) | 31/166 | 18.7% |

Elyada iCAF overlap peaks at **28%**, the highest across all K=2 fits. But this
factor also contains 18.7% Classical Tumor genes, reflecting the K=2 limitation:
with only 2 factors, iCAF and Classical tumor signals are mixed.

#### DECODER compartment breakdown (both factors, top-270)

| Factor | Immune | BasalTumor | ClassTumor | ActStroma | NormStroma |
|--------|--------|------------|------------|-----------|------------|
| F1 | 3/60 (5%) | 19/139 (14%) | -- | 13/108 (12%) | -- |
| **F2** | 6/60 (10%) | 12/139 (9%) | 31/166 (19%) | 13/108 (12%) | 1/15 (7%) |

F2 is enriched for Classical Tumor (19%) alongside iCAF genes. F1 absorbs
Basal Tumor (14%). The model creates a Basal-vs-Classical+iCAF axis.

#### Per-factor Cox: Training AND Validation

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
PurIST is included as a covariate, the residual F2 signal is absorbed.

**LRT adding F2 (iCAF) to PurIST + DeCAF:** Training p=0.002, Validation p=0.24

**Median-split KM for F2 (iCAF):**
- Training (n=273): High 27.0 vs Low 15.3 months, **p=6.87e-06**
- Validation (n=570): High 23.7 vs Low 20.1 months, p=0.051 (borderline)

The validation KM is right at the significance threshold, but the 3.6-month gap
is modest.

#### Beta structure

F1=**8.54e-04**, F2=**-1.10e-03** (dominant)

F2 (protective, iCAF-like) has the larger beta, consistent with the original
K=3 pattern where the iCAF factor dominates.

---

### K2_a055: alpha=0.55

#### Factor identity

Best iCAF match: **F2** (rho=0.664). This is the alpha matched to K=3's
CV-optimal value.

**H-score correlation with original K=3:**

| | Orig_F1 | Orig_F2 | Orig_F3 |
|--|---------|---------|---------|
| F1 | -0.461 | **0.520** | -0.065 |
| **F2** | **0.664** | 0.207 | -0.102 |

F2 reaches high correlation with original F1 (rho=0.664). But F1 also
anti-correlates with original F1 (rho=-0.461), suggesting the two K=2 factors
capture opposite ends of the *same* iCAF axis.

#### Reference signature overlaps (F2, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 3/25 | 12.0% |
| Elyada myCAF (15) | 4/15 | 26.7% |
| SCISSORS iCAF (17) | 3/17 | 17.6% |
| SCISSORS myCAF (16) | 1/16 | 6.2% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | **10/60** | **16.7%** |
| DECODER ActStroma (108) | 8/108 | 7.4% |
| DECODER NormStroma (15) | **5/15** | **33.3%** |
| DECODER BasalTumor (139) | 7/139 | 5.0% |
| DECODER ClassTumor (166) | 28/166 | 16.9% |

Interesting shift: Elyada iCAF overlap drops to 12% while DECODER Immune rises
to 16.7% and NormalStroma to 33.3%. The factor is becoming more immune/stroma-heavy
but less specifically iCAF. At high alpha with K=2, the survival penalty pushes
toward immune/stroma genes broadly rather than the specific iCAF program.

#### DECODER compartment breakdown (both factors, top-270)

| Factor | Immune | BasalTumor | ClassTumor | ActStroma | NormStroma |
|--------|--------|------------|------------|-----------|------------|
| F1 | 1/60 (2%) | **26/139 (19%)** | 24/166 (14%) | 19/108 (18%) | 0/15 (0%) |
| **F2** | **10/60 (17%)** | 7/139 (5%) | 28/166 (17%) | 8/108 (7%) | **5/15 (33%)** |

F2 is Immune + NormalStroma + Classical. F1 is Basal + ActivatedStroma. At K=2,
the model creates a single tumor-stroma contrast axis rather than isolating
the iCAF program from the tumor subtype axis.

#### Per-factor Cox: Training AND Validation

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
- Training p = **1.84e-05**
- Validation p = 0.16

**Median-split KM for F2 (iCAF):**
- Training (n=273): High 33.4 vs Low 15.1 months, **p=1.80e-08**
- Validation (n=570): High 24.1 vs Low 20.1 months, **p=0.012**

The KM validation p=0.012 shows there IS an iCAF signal in K=2's Factor 2, but
the 4.0-month gap (24.1 vs 20.1) is much smaller than K=3's 8.8-month gap
(26.8 vs 18.0). The iCAF program is diluted by mixing with the Classical tumor axis.

#### Beta structure

F1=**7.00e-04**, F2=**-1.13e-03** (dominant)

---

### K2_a075: alpha=0.75

#### Factor identity

Best iCAF match: **F2** (rho=0.674). Very similar structure to alpha=0.55,
with slightly higher H-cor.

**H-score correlation with original K=3:**

| | Orig_F1 | Orig_F2 | Orig_F3 |
|--|---------|---------|---------|
| F1 | -0.482 | **0.542** | -0.090 |
| **F2** | **0.674** | 0.192 | -0.075 |

The correlation structure is nearly identical to alpha=0.55. F2 aligns with
original F1 (rho=0.674, up slightly from 0.664). F1 continues to anti-correlate
with original F1 (rho=-0.482).

#### Reference signature overlaps (F2, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 5/25 | 20.0% |
| Elyada myCAF (15) | 4/15 | 26.7% |
| SCISSORS iCAF (17) | 2/17 | 11.8% |
| SCISSORS myCAF (16) | 3/16 | 18.8% |
| DeCAF restCAF (8) | 0/8 | 0.0% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 5/60 | 8.3% |
| DECODER ActStroma (108) | 7/108 | 6.5% |
| DECODER NormStroma (15) | **6/15** | **40.0%** |
| DECODER BasalTumor (139) | 9/139 | 6.5% |
| DECODER ClassTumor (166) | 26/166 | 15.7% |

DECODER NormStroma peaks at 40.0%. Immune drops back to 8.3% (from 16.7% at
alpha=0.55). The factor is increasingly dominated by NormalStroma + myCAF genes
rather than immune or iCAF-specific genes.

#### DECODER compartment breakdown (both factors, top-270)

| Factor | Immune | BasalTumor | ClassTumor | ActStroma | NormStroma |
|--------|--------|------------|------------|-----------|------------|
| F1 | 2/60 (3%) | **31/139 (22%)** | -- | -- | -- |
| **F2** | 5/60 (8%) | 9/139 (6%) | 26/166 (16%) | 7/108 (6%) | **6/15 (40%)** |

F1 is strongly Basal Tumor (22%). F2 continues the NormalStroma + Classical
Tumor pattern. The axis is increasingly polarized.

#### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1** | **1.457** | **3.44e-04** | **1.254** | **0.005** |
| **F2 (iCAF)** | **0.619** | **<0.0001** | 0.878 | 0.081 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.256 | 0.047 | 1.074 | 0.40 |
| **F2 (iCAF)** | **0.638** | **<0.0001** | 0.883 | 0.11 |

Training signal continues to strengthen (F2 adj HR=0.638) while validation
remains non-significant (adj p=0.11). F1 reaches adjusted training significance
for the first time (p=0.047), reflecting the increasing polarization of the
two-factor axis.

**LRT adding F2 (iCAF) to PurIST + DeCAF:**
- Training p = **7.87e-07**
- Validation p = 0.12

**Median-split KM for F2 (iCAF):**
- Training (n=273): High 33.4 vs Low 14.2 months, **p=1.17e-09**
- Validation (n=570): High 24.1 vs Low 19.8 months, **p=0.004**

The KM validation strengthens (p=0.004, gap=4.3 months), the best KM result
so far. But the Cox model still cannot disentangle the iCAF factor from
PurIST/DeCAF in validation.

#### Beta structure

F1=**7.79e-04**, F2=**-1.26e-03** (dominant)

---

### K2_a085: alpha=0.85

#### Factor identity

Best iCAF match: **F1** (rho=0.671). Factor ordering flips again at this alpha.

**H-score correlation with original K=3:**

| | Orig_F1 | Orig_F2 | Orig_F3 |
|--|---------|---------|---------|
| **F1** | **0.671** | 0.232 | -0.120 |
| F2 | -0.501 | **0.508** | -0.048 |

F1 is the iCAF factor (rho=0.671) and F2 is the stroma/Basal factor
(rho=-0.501 with original F1, 0.508 with original F2). The flip in factor
ordering (iCAF is F1 here vs F2 at alpha=0.55/0.75) does not change the
underlying biological structure.

#### Reference signature overlaps (F1, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 4/25 | 16.0% |
| Elyada myCAF (15) | **5/15** | **33.3%** |
| SCISSORS iCAF (17) | **5/17** | **29.4%** |
| SCISSORS myCAF (16) | 0/16 | 0.0% |
| DeCAF restCAF (8) | 0/8 | 0.0% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 2/60 | 3.3% |
| DECODER ActStroma (108) | 6/108 | 5.6% |
| DECODER NormStroma (15) | **6/15** | **40.0%** |
| DECODER BasalTumor (139) | 8/139 | 5.8% |
| DECODER ClassTumor (166) | 25/166 | 15.1% |

High overlap with myCAF (33.3%), SCISSORS iCAF (29.4%), and NormalStroma (40.0%).
Immune drops to 3.3%. The factor at this alpha is NormalStroma + CAF (both
iCAF and myCAF) with minimal immune component.

#### DECODER compartment breakdown (both factors, top-270)

| Factor | Immune | BasalTumor | ClassTumor | ActStroma | NormStroma |
|--------|--------|------------|------------|-----------|------------|
| **F1** | 2/60 (3%) | 8/139 (6%) | 25/166 (15%) | 6/108 (6%) | **6/15 (40%)** |
| F2 | 1/60 (2%) | **39/139 (28%)** | -- | -- | -- |

F2 is very strongly Basal Tumor (28%). F1 is NormalStroma + Classical Tumor.
The polarization into Basal-vs-Classical/Stroma is now extreme.

#### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.607** | **<0.0001** | 0.881 | 0.12 |
| **F2** | **1.525** | **<0.0001** | **1.236** | **0.005** |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.615** | **<0.0001** | 0.883 | 0.15 |
| F2 | 1.318 | 0.018 | 1.064 | 0.43 |

F1 (iCAF) has the strongest training signal yet (adj HR=0.615, p<0.0001) but
validation remains non-significant (adj p=0.15). F2 (Basal-like) validates
unadjusted (p=0.005) but not adjusted (p=0.43), confirming it captures the
PurIST axis.

**LRT adding F1 (iCAF) to PurIST + DeCAF:**
- Training p = **8.51e-08**
- Validation p = 0.15

**Median-split KM for F1 (iCAF):**
- Training (n=273): High 33.4 vs Low 15.1 months, **p=1.27e-09**
- Validation (n=570): High 24.1 vs Low 19.2 months, **p=0.002**

The KM validation reaches its best p-value (0.002) with a 4.9-month gap.
But the Cox-based adjusted test still fails.

#### Beta structure

F1=**-1.28e-03** (dominant, protective), F2=**7.96e-04**

The iCAF factor (F1) has the largest absolute beta across all K=2 fits so far.

---

### K2_a095: alpha=0.95

#### Factor identity

Best iCAF match: **F2** (rho=0.728). This is the highest H-cor with original
F1 across all K=2 fits.

**H-score correlation with original K=3:**

| | Orig_F1 | Orig_F2 | Orig_F3 |
|--|---------|---------|---------|
| F1 | -0.479 | **0.521** | -0.063 |
| **F2** | **0.728** | 0.192 | -0.113 |

F2 reaches the peak correlation with original F1 (rho=0.728), the strongest
iCAF match across all K=2 fits. Yet this high H-cor still does not translate
to validation significance.

#### Reference signature overlaps (F2, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 5/25 | 20.0% |
| Elyada myCAF (15) | 3/15 | 20.0% |
| SCISSORS iCAF (17) | **6/17** | **35.3%** |
| SCISSORS myCAF (16) | 4/16 | 25.0% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | 1/4 | 25.0% |
| DECODER Immune (60) | 6/60 | 10.0% |
| DECODER ActStroma (108) | 12/108 | 11.1% |
| DECODER NormStroma (15) | **8/15** | **53.3%** |
| DECODER BasalTumor (139) | 5/139 | 3.6% |
| DECODER ClassTumor (166) | 28/166 | 16.9% |

DECODER NormStroma reaches **53.3%**, the highest for any K=2 fit. SCISSORS
iCAF also peaks at 35.3%. The factor is increasingly dominated by
stroma/CAF-related genes. BasalTumor drops to its lowest (3.6%), showing
maximal separation from the tumor axis in gene space -- even though the H-score
axis cannot be separated from PurIST.

#### DECODER compartment breakdown (both factors, top-270)

| Factor | Immune | BasalTumor | ClassTumor | ActStroma | NormStroma |
|--------|--------|------------|------------|-----------|------------|
| F1 | 0/60 (0%) | **31/139 (22%)** | -- | 12/108 (11%) | -- |
| **F2** | 6/60 (10%) | 5/139 (4%) | 28/166 (17%) | 12/108 (11%) | **8/15 (53%)** |

F2 is half NormalStroma genes (53%) plus Classical Tumor (17%) and Immune (10%).
F1 is purely Basal Tumor (22%) + ActivatedStroma (11%). Maximum compartment
separation, but still only 2 axes.

#### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1** | **1.520** | **<0.0001** | **1.222** | **0.012** |
| **F2 (iCAF)** | **0.544** | **<0.0001** | 0.886 | 0.11 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.332 | 0.014 | 1.051 | 0.55 |
| **F2 (iCAF)** | **0.550** | **<0.0001** | 0.892 | 0.14 |

Training reaches the most extreme HR across all K=2 fits (F2 adj HR=0.550),
yet validation remains stubbornly non-significant (adj p=0.14). This is the
clearest evidence that K=2's training-validation gap is **structural**, not
fixable by increasing alpha.

**LRT adding F2 (iCAF) to PurIST + DeCAF:**
- Training p = **1.81e-10**
- Validation p = 0.15

**Median-split KM for F2 (iCAF):**
- Training (n=273): High 33.5 vs Low 13.3 months, **p=1.94e-12**
- Validation (n=570): High 23.9 vs Low 19.8 months, **p=0.017**

The training KM gap reaches 20.2 months (33.5 vs 13.3), a massive separation
that far exceeds anything seen at K=3. But the validation gap is only 4.1 months.
This 5:1 ratio of training-to-validation gap size is the hallmark of
overfitting to the training data's PurIST axis.

#### Beta structure

F1=**7.72e-04**, F2=**-1.45e-03** (dominant)

The iCAF beta magnitude increases monotonically with alpha: -2.79e-04 (a=0),
-9.52e-04 (a=0.25), -1.10e-03 (a=0.35), -1.13e-03 (a=0.55), -1.26e-03 (a=0.75),
-1.28e-03 (a=0.85), -1.45e-03 (a=0.95).

---

## The K=2 Limitation: Why Two Factors Cannot Isolate iCAF

At K=2, the model must capture all biological variation in just two factors.
The result across all 7 alpha values is consistent:

1. **One factor mixes iCAF + Classical tumor + Normal stroma + Immune**
   (the "good prognosis" end)
2. **One factor mixes Basal-like tumor + Activated stroma**
   (the "bad prognosis" end)

This is essentially the **PurIST/DeCAF axis** with iCAF genes layered on top.
The iCAF program cannot be separated from the Classical tumor program because
there aren't enough factors.

### The expanded grid answers all three questions

**Q1: Does iCAF coherence keep improving beyond alpha=0.55?**
Yes. H-cor rises from 0.664 (a=0.55) to 0.728 (a=0.95). Gene overlap rises
from 92/270 to 100/270. NormalStroma overlap rises from 33% to 53%. Higher
alpha pushes more iCAF/stroma genes into the protective factor.

**Q2: Does validation signal eventually emerge at high alpha?**
No. The adjusted validation p-value never reaches significance (best: 0.11 at
a=0.75). The unadjusted validation HR improves marginally from 0.876 (a=0.55)
to 0.878-0.886 (a=0.75-0.95) but is not trending toward significance. The KM
median-split does reach significance (p=0.002 at a=0.85), but this coarser
test cannot distinguish iCAF from PurIST.

**Q3: Is the mixing problem inherent to K=2 or alpha-dependent?**
Inherent to K=2. The DECODER compartment breakdown shows that even at a=0.95
where the iCAF factor has 53% NormalStroma overlap and only 4% BasalTumor, the
factor still contains 17% Classical Tumor genes. The iCAF and Classical tumor
programs share the same factor because K=2 forces them together.

### Training vs validation divergence summary

| Alpha | Train adj HR | Val adj HR | Train LRT p | Val LRT p |
|-------|-------------|-----------|-------------|----------|
| 0.00 | 0.900 | 0.953 | 0.22 | 0.54 |
| 0.25 | 0.804 | 0.920 | 0.010 | 0.33 |
| 0.35 | 0.760 | 0.906 | 0.002 | 0.24 |
| 0.55 | 0.672 | 0.891 | 1.84e-05 | 0.16 |
| 0.75 | 0.638 | 0.883 | 7.87e-07 | 0.12 |
| 0.85 | 0.615 | 0.883 | 8.51e-08 | 0.15 |
| 0.95 | 0.550 | 0.892 | 1.81e-10 | 0.15 |

The training HR decreases monotonically (0.900 to 0.550) while the validation
HR barely moves (0.953 to 0.883). The training LRT p-value drops 11 orders of
magnitude while the validation LRT p-value improves only from 0.54 to ~0.15.
This scissors pattern is the signature of a model fitting to confounded
training signal (iCAF + PurIST) that does not generalize.

---

## Glossary

- **H-score**: Factor loading score, H = X^T W. How strongly each sample expresses a factor.
- **H-cor**: Spearman correlation between H-score vectors from two fits.
- **LP**: Cox model linear predictor, sum of beta_k * H_k across factors.
- **LRT**: Likelihood ratio test comparing nested Cox models. Tests whether
  adding a factor improves a model already containing PurIST + DeCAF.
- **ntop**: NULL = all genes; 270 (production) = top 270 per factor.
- **strata(dataset)**: Stratified Cox allowing per-cohort baseline hazards.
- **Gene universe**: 1970 genes after variance/expression filtering from 3000.
- **CV-optimal alpha**: Alpha minimizing cross-validated partial log-likelihood.
- **Validation-optimal alpha**: Alpha with strongest external validation signal.
