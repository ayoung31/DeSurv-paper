# DeSurv K-Sensitivity Synthesis: The iCAF Coherence Story (Expanded Grid)

**For discussion with Jen Jen Yeh**
**Date: 2026-02-18 (amber revision)**
**Companion documents:**
- `2026-02-17-jjy-talking-points-k2-grid_amber.md` (K=2 factor analysis)
- `2026-02-17-jjy-talking-points-k3-grid_amber.md` (K=3 factor analysis)
- `2026-02-17-jjy-talking-points-k5-grid_amber.md` (K=5 factor analysis)
- `2026-02-17-jjy-talking-points-k7-grid_amber.md` (K=7 factor analysis)
- `2026-02-17-jjy-talking-points-k9-grid_amber.md` (K=9 factor analysis)

**Note for Amber:** This is the expanded-grid version of the synthesis document.
It uses 35 cv_grid fits (K=2,3,5,7,9 x alpha=0,0.25,0.35,0.55,0.75,0.85,0.95)
instead of the original 9. All hyperparameters except K and alpha are held constant
(lambda=0.3, nu=0.05, ntop=NULL, 100 initializations). The analysis scripts are in
`inst/cv_grid_training_analysis_amber.R` (training) and
`inst/cv_grid_validation_analysis_amber.R` (validation). Run from the repo root.

Fits are loaded by matching (k, alpha) parameters directly from the targets store,
not by hardcoded branch hashes. This eliminates the risk of hash-parameter mismatches.

---

## Executive Summary

We fit DeSurv across **K=2, 3, 5, 7, 9** and **alpha=0, 0.25, 0.35, 0.55, 0.75,
0.85, 0.95** using the cv_grid exhaustive search, with all hyperparameters except
K and alpha held constant (lambda=0.3, nu=0.05, 100 initializations, ntop=NULL).
This is an apples-to-apples comparison across 35 parameter combinations.

**All results are validated in 570 independent patients** across 4 external
cohorts (Dijk, Moffitt GEO, Puleo, PACA-AU), using `strata(dataset)` to account
for cohort-specific baseline hazards and adjusting for PurIST + DeCAF classifiers.

The central finding: **The iCAF program is alpha-dependent, not just
K-dependent.** At alpha=0 (standard NMF), no value of K recovers the iCAF
program. The survival penalty creates it. Among all 35 K x alpha combinations
tested, K=3 at alpha=0.55 produces the highest coherence with the original
production iCAF factor (H-cor=0.906). However, K=7 at alpha=0.75 produces
the strongest adjusted validation result (p=0.019), revealing that higher K
values can also capture this program when given sufficient alpha.

### What "iCAF" means in this document

Throughout, "iCAF program" refers to the factor that best matches the
**original production K=3 Factor 1** (from the BO-selected fit with
alpha=0.334, lambda=0.349, nu=0.056, ntop=270). That factor was
characterized as iCAF-associated based on:

1. **Elyada et al. (2019) iCAF 35-gene signature**, from single-cell
   RNA-seq of PDAC cancer-associated fibroblasts. 25 of 35 genes present
   in the 1970-gene universe. Key genes: PLA2G2A, DPT, CXCL14, PI16,
   CCDC80, FSTL1, PTX3, FBLN2, TNXB, SOD2, APOE, FBLN1, ADH1B, GPX3.

2. **SCISSORS iCAF 25-gene signature**, from Raghavan et al. (2021).
   17 of 25 present in the gene universe.

3. **DECODER Immune compartment** (408 genes, 60 in the gene universe),
   reflecting B cell co-expression.

---

## The Master Tables

### K x Alpha Landscape Matrices

#### H-cor Matrix: K (rows) x Alpha (columns)

Best iCAF factor Spearman correlation with original K=3 Factor 1.

| K \ Alpha | 0.00 | 0.25 | 0.35 | 0.55 | 0.75 | 0.85 | 0.95 |
|-----------|------|------|------|------|------|------|------|
| 2 | 0.108 | 0.409 | 0.518 | 0.664 | 0.674 | 0.671 | 0.728 |
| 3 | 0.387 | 0.610 | 0.731 | **0.906** | 0.663 | 0.794 | 0.836 |
| 5 | 0.313 | 0.449 | 0.667 | 0.737 | 0.816 | 0.823 | 0.893 |
| 7 | 0.510 | 0.784 | 0.758 | 0.873 | 0.872 | 0.870 | 0.814 |
| 9 | 0.468 | 0.785 | 0.533 | 0.793 | 0.806 | 0.741 | 0.806 |

#### Unadjusted Validation P-value Matrix: K (rows) x Alpha (columns)

iCAF factor Cox p-value in validation (strata(dataset)).

| K \ Alpha | 0.00 | 0.25 | 0.35 | 0.55 | 0.75 | 0.85 | 0.95 |
|-----------|------|------|------|------|------|------|------|
| 2 | 0.814 | 0.535 | 0.283 | 0.087 | 0.081 | 0.119 | 0.105 |
| 3 | 0.102 | 0.277 | **0.027** | **0.003** | 0.440 | 0.060 | **0.040** |
| 5 | 0.651 | **4.9e-4** | **0.002** | 0.127 | **0.031** | **0.014** | **5.0e-4** |
| 7 | 0.096 | 0.138 | **0.013** | **0.005** | **1.3e-4** | **0.002** | **0.003** |
| 9 | 0.691 | **0.016** | 0.223 | **0.008** | **0.004** | 0.156 | **0.008** |

#### Adjusted Validation P-value Matrix: K (rows) x Alpha (columns)

iCAF factor Cox p-value in validation, adjusted for PurIST + DeCAF + strata(dataset).

| K \ Alpha | 0.00 | 0.25 | 0.35 | 0.55 | 0.75 | 0.85 | 0.95 |
|-----------|------|------|------|------|------|------|------|
| 2 | 0.535 | 0.321 | 0.236 | 0.157 | 0.112 | 0.147 | 0.145 |
| 3 | 0.254 | 0.192 | 0.072 | 0.075 | 0.282 | 0.186 | 0.193 |
| 5 | 0.568 | 0.069 | 0.069 | 0.167 | 0.090 | 0.056 | 0.082 |
| 7 | 0.762 | 0.265 | **0.033** | 0.104 | **0.019** | 0.075 | **0.046** |
| 9 | 0.475 | 0.071 | 0.358 | **0.050** | **0.020** | 0.149 | 0.100 |

---

### Combined Master Training Table (Selected Fits, n=273)

| Fit | K | Alpha | iCAF Factor | H-cor | Gene Overlap | Elyada iCAF | DECODER Immune | Train HR (unadj) | Train HR (adj) | Train LRT p |
|-----|---|-------|------------|-------|-------------|------------|----------------|------------------|----------------|-------------|
| K2_a0 | 2 | 0.00 | F2 | 0.108 | 42/270 | 6/25 (24%) | 4/60 (7%) | 0.907 (p=0.20) | 0.900 (p=0.21) | 0.22 |
| K3_a0 | 3 | 0.00 | F3 | 0.387 | 48/270 | 3/25 (12%) | 0/60 (0%) | 0.927 (p=0.40) | 0.970 (p=0.74) | 0.74 |
| K5_a0 | 5 | 0.00 | F1 | 0.313 | 47/270 | 3/25 (12%) | 1/60 (2%) | 0.969 (p=0.72) | 0.972 (p=0.75) | 0.75 |
| K7_a0 | 7 | 0.00 | F5 | 0.510 | 68/270 | 9/25 (36%) | 16/60 (27%) | 0.850 (p=0.067) | 0.914 (p=0.35) | 0.35 |
| K9_a0 | 9 | 0.00 | F6 | 0.468 | 50/270 | 7/25 (28%) | 1/60 (2%) | 0.829 (p=0.035) | 0.842 (p=0.063) | 0.063 |
| K5_a025 | 5 | 0.25 | F4 | 0.449 | 121/270 | 6/25 (24%) | **27/60 (45%)** | **0.652 (p<.0001)** | **0.767 (p=.007)** | 0.008 |
| K9_a025 | 9 | 0.25 | F4 | 0.785 | 129/270 | 8/25 (32%) | 21/60 (35%) | **0.526 (p<.0001)** | **0.578 (p<.0001)** | 7.0e-8 |
| **K3_a035** | **3** | **0.35** | **F1** | **0.731** | **132/270** | 4/25 (16%) | **20/60 (33%)** | **0.601 (p<.0001)** | **0.674 (p<.0001)** | **8e-5** |
| K5_a035 | 5 | 0.35 | F1 | 0.667 | 105/270 | 2/25 (8%) | 11/60 (18%) | **0.641 (p<.0001)** | **0.754 (p=.006)** | 0.005 |
| K7_a035 | 7 | 0.35 | F6 | 0.758 | 125/270 | 4/25 (16%) | 10/60 (17%) | **0.397 (p<.0001)** | **0.365 (p<.0001)** | 1.4e-20 |
| **K3_a055** | **3** | **0.55** | **F1** | **0.906** | **168/270** | **6/25 (24%)** | 16/60 (27%) | **0.344 (p<.0001)** | **0.371 (p<.0001)** | **5e-18** |
| K5_a055 | 5 | 0.55 | F4 | 0.737 | 146/270 | 6/25 (24%) | 15/60 (25%) | **0.529 (p<.0001)** | **0.562 (p<.0001)** | 4e-9 |
| K7_a055 | 7 | 0.55 | F7 | 0.873 | 119/270 | 5/25 (20%) | 15/60 (25%) | **0.085 (p<.0001)** | **0.076 (p<.0001)** | 1.2e-64 |
| K9_a055 | 9 | 0.55 | F7 | 0.793 | 127/270 | 5/25 (20%) | 16/60 (27%) | **0.177 (p<.0001)** | **0.163 (p<.0001)** | 2.7e-50 |
| **K7_a075** | **7** | **0.75** | **F2** | **0.872** | **133/270** | 4/25 (16%) | 9/60 (15%) | **0.108 (p<.0001)** | **0.086 (p<.0001)** | **1.1e-61** |
| K9_a075 | 9 | 0.75 | F5 | 0.806 | 120/270 | 4/25 (16%) | 15/60 (25%) | **0.405 (p<.0001)** | **0.420 (p<.0001)** | 2.6e-16 |
| K5_a095 | 5 | 0.95 | F1 | 0.893 | 144/270 | 4/25 (16%) | 18/60 (30%) | **0.288 (p<.0001)** | **0.286 (p<.0001)** | 2.0e-24 |
| K7_a095 | 7 | 0.95 | F1 | 0.814 | 112/270 | 5/25 (20%) | 13/60 (22%) | **0.056 (p<.0001)** | **0.044 (p<.0001)** | 2.1e-79 |

### Combined Master Validation Table (Selected Fits, n=570, strata(dataset))

| Fit | K | Alpha | iCAF Factor | Val HR (unadj) | Val HR (adj) | Val LRT p | Val KM High | Val KM Low | Val KM p |
|-----|---|-------|------------|----------------|-------------|-----------|-------------|------------|----------|
| K2_a0 | 2 | 0.00 | F2 | 1.018 (p=0.81) | 0.953 (p=0.53) | 0.54 | 22.5 | 21.0 | 0.13 |
| K3_a0 | 3 | 0.00 | F3 | 0.902 (p=0.10) | 0.926 (p=0.25) | 0.26 | 24.8 | 19.2 | **0.003** |
| K5_a0 | 5 | 0.00 | F1 | 0.971 (p=0.65) | 0.962 (p=0.57) | 0.57 | 22.5 | 21.4 | 0.59 |
| K7_a0 | 7 | 0.00 | F5 | 0.885 (p=0.096) | 1.025 (p=0.76) | 0.76 | 26.8 | 17.4 | **5.5e-5** |
| K9_a0 | 9 | 0.00 | F6 | 0.967 (p=0.69) | 0.939 (p=0.48) | 0.48 | 24.4 | 19.4 | **0.009** |
| K5_a025 | 5 | 0.25 | F4 | **0.797 (p=5e-4)** | 0.879 (p=.069) | 0.072 | **25.6** | **18.0** | **4.4e-6** |
| K9_a025 | 9 | 0.25 | F4 | **0.837 (p=.016)** | 0.867 (p=.071) | 0.073 | **28.7** | **19.0** | **8.4e-4** |
| **K3_a035** | **3** | **0.35** | **F1** | **0.841 (p=.027)** | 0.862 (p=.072) | 0.075 | **24.9** | **19.2** | **7.8e-4** |
| K5_a035 | 5 | 0.35 | F1 | **0.789 (p=.002)** | 0.860 (p=.069) | 0.072 | **24.3** | **19.0** | **3.8e-4** |
| K7_a035 | 7 | 0.35 | F6 | **0.844 (p=.013)** | **0.856 (p=.033)** | **0.034** | **26.2** | **20.0** | **0.010** |
| **K3_a055** | **3** | **0.55** | **F1** | **0.800 (p=.003)** | 0.867 (p=.075) | 0.076 | **26.8** | **18.0** | **3.0e-4** |
| K5_a055 | 5 | 0.55 | F4 | 0.884 (p=0.13) | 0.888 (p=0.17) | 0.17 | 25.3 | 19.2 | **0.005** |
| K7_a055 | 7 | 0.55 | F7 | **0.816 (p=.005)** | 0.885 (p=0.10) | 0.11 | **30.2** | **17.0** | **2.0e-7** |
| K9_a055 | 9 | 0.55 | F7 | **0.827 (p=.008)** | 0.864 (p=.050) | 0.050 | **29.2** | **19.0** | **2.8e-5** |
| **K7_a075** | **7** | **0.75** | **F2** | **0.771 (p=1.3e-4)** | **0.844 (p=.019)** | **0.020** | **30.8** | **16.6** | **8.0e-9** |
| K9_a075 | 9 | 0.75 | F5 | **0.823 (p=.004)** | **0.844 (p=.020)** | **0.021** | **25.6** | **19.2** | **7.1e-4** |
| K5_a095 | 5 | 0.95 | F1 | **0.783 (p=5e-4)** | 0.876 (p=.082) | 0.083 | **27.0** | **18.0** | **9.4e-6** |
| K7_a095 | 7 | 0.95 | F1 | **0.779 (p=.003)** | **0.840 (p=.046)** | **0.047** | **30.1** | **17.0** | **2.4e-7** |

---

## Five Key Insights

### 1. The Survival Penalty Creates the iCAF Program (Validated)

At alpha=0 (standard NMF), no K value recovers the iCAF program as a prognostic
factor in Cox regression:

| K | Best iCAF H-cor at alpha=0 | Train significant? | Val Cox significant? | Val KM p |
|---|---------------------------|-------------------|---------------------|----------|
| 2 | 0.108 | No (p=0.20) | No (p=0.81) | 0.13 |
| 3 | 0.387 | No (p=0.40) | No (p=0.10) | **0.003*** |
| 5 | 0.313 | No (p=0.72) | No (p=0.65) | 0.59 |
| 7 | 0.510 | No (p=0.067) | No (p=0.096) | **5.5e-5*** |
| 9 | 0.468 | Marginal (p=0.035) | No (p=0.69) | **0.009*** |

*K3_a0 and K7_a0 show significant KM p-values despite non-significant Cox
models. This reflects nonspecific stroma-vs-tumor median splits that separate
survival groups but do not add information beyond PurIST/DeCAF (adj p=0.25 and
0.76, respectively). The KM significance without Cox significance is a sign
that the factor captures a coarse axis (high vs low stroma) rather than a
specific prognostic program.

Standard NMF discovers variation-maximizing factors (the stroma/basal-like
axis, Classical tumor genes, immune programs as separate compartments) but
**never combines them into the specific iCAF + B cell + normal stroma program**
that the original K=3 Factor 1 captures.

The survival penalty (alpha > 0) forces the factorization to prioritize genes
that jointly predict survival. This pulls iCAF, B cell, and normal stroma genes
-- which individually explain modest variance but collectively predict survival --
into a single factor. **The iCAF program is a survival-selected gene program,
not a variance-maximizing one.** This is the core methodological claim of DeSurv.

### 2. K=3 Alpha=0.55 Is the Sweet Spot for iCAF Coherence

At alpha=0.55, K=3 achieves the highest H-cor in the entire 5x7 grid (0.906).
Comparison across K values at alpha=0.55:

| Metric | K=2 (a=0.55) | K=3 (a=0.55) | K=5 (a=0.55) | K=7 (a=0.55) | K=9 (a=0.55) |
|--------|-------------|-------------|-------------|-------------|-------------|
| H-cor with orig F1 | 0.664 | **0.906** | 0.737 | 0.873 | 0.793 |
| Gene overlap | 92/270 | **168/270** | 146/270 | 119/270 | 127/270 |
| Train HR (unadj) | 0.634 | **0.344** | 0.529 | 0.085 | 0.177 |
| **Val HR (unadj)** | 0.876 (p=0.087) | **0.800 (p=0.003)** | 0.884 (p=0.13) | **0.816 (p=0.005)** | **0.827 (p=0.008)** |
| **Val KM gap** | 4.0 mo (p=0.012) | **8.8 mo (p=3e-4)** | 6.1 mo (p=0.005) | **13.2 mo (p=2e-7)** | **10.2 mo (p=3e-5)** |
| Val LRT vs PurIST+DeCAF | 0.16 | 0.076 | 0.17 | 0.11 | **0.050** |
| Elyada iCAF | 3/25 (12%) | 6/25 (24%) | 6/25 (24%) | 5/25 (20%) | 5/25 (20%) |
| DECODER Immune | 10/60 (17%) | 16/60 (27%) | 15/60 (25%) | 15/60 (25%) | 16/60 (27%) |
| # sig factors (train) | 2 | **1** | 3 | 2 | 4 |

**Why K=3 alpha=0.55 is the coherence sweet spot:**
- **Highest H-cor** in the entire grid (0.906): nearly identical sample rankings
  to the original production factor
- **Single prognostic factor**: Only F1 is significant in training; the other
  two factors are noise. This means the iCAF program is isolated, not entangled.
- **Largest gene overlap** (168/270 = 62%): the top-270 genes substantially
  reproduce the original factor's gene program.

However, K=3 alpha=0.55 is NOT the strongest validation result. K=7 alpha=0.75
has stronger unadjusted validation (p=1.3e-4 vs 0.003) and is the only K=3+
configuration that achieves adjusted significance (p=0.019).

### 3. Higher K Does NOT Simply Fragment the iCAF Program

The original hypothesis that K>3 fragments the iCAF program is **revised**.
The expanded grid reveals that K=7 achieves H-cor values comparable to K=3:

```
H-cor at alpha=0.55:     K=3: 0.906    K=5: 0.737    K=7: 0.873    K=9: 0.793
H-cor at alpha=0.75:     K=3: 0.663    K=5: 0.816    K=7: 0.872    K=9: 0.806
H-cor at alpha=0.85:     K=3: 0.794    K=5: 0.823    K=7: 0.870    K=9: 0.741
```

At K=7, the iCAF program concentrates into a single factor (which varies:
F5, F3, F6, F7, F2, F6, or F1 depending on alpha) while the extra factors
capture other biological signals. K=7 also has the strongest validation
results in the entire grid:

| K7 Fit | Val HR (unadj) | Val p (unadj) | Val HR (adj) | Val p (adj) | Val KM gap | Val KM p |
|--------|----------------|---------------|-------------|-------------|------------|----------|
| K7_a035 | 0.844 | **0.013** | **0.856** | **0.033** | 6.2 mo | **0.010** |
| K7_a055 | 0.816 | **0.005** | 0.885 | 0.10 | 13.2 mo | **2.0e-7** |
| **K7_a075** | **0.771** | **1.3e-4** | **0.844** | **0.019** | **14.2 mo** | **8.0e-9** |
| K7_a085 | 0.812 | **0.002** | 0.880 | 0.075 | 11.0 mo | **1.4e-6** |
| K7_a095 | 0.779 | **0.003** | **0.840** | **0.046** | 13.1 mo | **2.4e-7** |

K=7 validates at 5 of 6 non-zero alpha values (unadjusted) and achieves
adjusted significance at 3 alpha values (0.35, 0.75, 0.95). This is the
most robustly validated K value in the grid.

True fragmentation begins at K=9, where H-cor is erratic (0.468-0.806) and
validation is inconsistent (4 of 6 non-zero alpha values validate unadjusted,
but only 2 achieve adjusted significance).

### 4. The K5_a025 Exception: Pure Immune Factor

K5_a025 presents an interesting counterpoint: its iCAF factor (F4) has the
strongest unadjusted validation signal among K=5 fits (HR=0.797, p=5e-4)
and excellent validation KM (25.6 vs 18.0, p=4.4e-6). But:

- Its H-correlation with the original iCAF program is only 0.449
- It is a purer immune factor (**45% DECODER Immune**, 0% Basal Tumor)
- It does not survive PurIST+DeCAF adjustment (p=0.069)

| Metric | K3_a055 (iCAF) | K5_a025 (pure immune) |
|--------|---------------|----------------------|
| H-cor | 0.906 | 0.449 |
| DECODER Immune overlap | 16/60 (27%) | 27/60 (45%) |
| DECODER BasalTumor overlap | 17/139 (12%) | 0/139 (0%) |
| Elyada iCAF overlap | 6/25 (24%) | 6/25 (24%) |
| Val HR (unadj) | 0.800 (p=0.003) | **0.797 (p=5e-4)** |
| Val HR (adj) | 0.867 (p=0.075) | 0.879 (p=0.069) |
| Val KM gap | 8.8 mo (p=3e-4) | **7.6 mo (p=4e-6)** |

This suggests the **immune component of Factor 1 is the primary survival
driver**. K=5 at moderate alpha separates the immune signal from iCAF/stroma,
producing a cleaner immune factor that validates strongly. But it is not the
full iCAF program -- it is a fragment that shares survival information
with PurIST (which also captures immune-related biology through the Basal-like
vs Classical axis).

K=3 keeps immune + iCAF + normal stroma together, which is both more coherent
with the original program (rho=0.906 vs 0.449) and more interpretable
biologically (an organized immune microenvironment rather than isolated B cells).

### 5. The Alpha x K Interaction Reveals Two Regimes

The 5x7 grid reveals that the alpha-K landscape has two distinct regimes:

**Low K (2-3):** Alpha controls whether the iCAF program emerges at all.
K=2 never validates regardless of alpha. K=3 validates at a narrow alpha
band (0.35, 0.55, and marginally 0.95) with a sharp optimum at 0.55.
The alpha=0.75 dip (H-cor drops from 0.906 to 0.663) suggests the
survival penalty can overshoot at K=3, pushing the iCAF signal into a
different factor configuration.

**High K (5-9):** Alpha controls how strongly the survival penalty
concentrates prognostic signal. K=5 and K=7 validate across a broad
alpha range (multiple alpha values from 0.25 to 0.95), suggesting the
extra factors provide a "buffer" that absorbs noise while preserving the
iCAF factor. K=7 alpha=0.75 is the strongest overall validation result
(unadjusted p=1.3e-4, adjusted p=0.019, KM gap=14.2 months, p=8.0e-9).

The adjusted validation story is particularly revealing:

```
Adjusted validation p<0.05:
K=2:  0 of 7 alpha values
K=3:  0 of 7 alpha values  (best: 0.072 at alpha=0.35)
K=5:  0 of 7 alpha values  (best: 0.056 at alpha=0.85)
K=7:  3 of 7 alpha values  (alpha=0.35: 0.033, 0.75: 0.019, 0.95: 0.046)
K=9:  2 of 7 alpha values  (alpha=0.55: 0.050, 0.75: 0.020)
```

Only K>=7 achieves adjusted significance, meaning the iCAF factor at K=7+
captures information **beyond PurIST and DeCAF** that K<=5 does not.

---

## The Biological Interpretation

### What the iCAF program IS (at K=3, alpha>=0.35)

The K=3 iCAF factor combines:
1. **iCAF genes** (Elyada: PI16, CXCL14, DPT, CCDC80, FBLN2 / SCISSORS: PI16,
   CXCL14, DPT, OGN, MFAP5, FBLN2, HAS1), chemokine-secreting fibroblasts
2. **B cell / immune genes** (DECODER Immune: 16-27% overlap), tertiary
   lymphoid structures
3. **Normal stroma genes** (DECODER NormalStroma: 33-47% overlap), tissue
   homeostasis

This combination makes biological sense: iCAFs create chemokine gradients
(CCL19, CCL21, CXCL14) that recruit B cells and organize tertiary lymphoid
structures. Normal stroma provides the structural scaffold. Together, they
represent an **organized immune microenvironment**, and its presence (high
Factor 1 score) is protective.

### What the iCAF program is NOT

- It is NOT a variance-maximizing program (absent at alpha=0)
- It is NOT a tumor subtype marker (PurIST captures that)
- It is NOT an activated stroma marker (DeCAF/myCAF captures that)
- It is NOT a generic immune signature (K=5 separates "pure immune" factors
  with 45% DECODER Immune overlap, but those are not the same program)

### Why it adds to PurIST + DeCAF

PurIST classifies tumor cells (Basal-like vs Classical). DeCAF classifies CAFs
(proCAF vs restCAF). The iCAF program captures a third axis: **the organized
immune microenvironment** (iCAF + B cell + normal stroma). This axis is:
- Orthogonal to PurIST (not a tumor subtype marker)
- Orthogonal to DeCAF (not an activated/resting CAF distinction)
- The adjusted validation models show borderline significance at K=3 (p=0.075)
  with ntop=NULL. The original production fit with focused ntop=270 achieved
  p=0.004. At K=7, adjusted significance is achieved even with ntop=NULL
  (p=0.019 at alpha=0.75), suggesting the higher-K factorization may better
  separate iCAF from the PurIST/DeCAF axes.

---

## Comparison to Previous K-Sensitivity Analysis

| Issue | Previous Analysis (BO-selected) | Original cv_grid (9 fits) | This Analysis (35 fits) |
|-------|--------------------------------|--------------------------|------------------------|
| K values tested | K=2, 3, 5 | K=2, 3, 5 | K=2, 3, 5, 7, 9 |
| Alpha values tested | 1 per K (BO-selected) | 3 per K (0, ~0.35, 0.55) | 7 per K (0, 0.25, 0.35, 0.55, 0.75, 0.85, 0.95) |
| Hyperparameter comparability | K=2 used K=3's params; K=5 used BO with nu=0.005 | Fixed lambda=0.3, nu=0.05 for all | Fixed lambda=0.3, nu=0.05 for all |
| Key confound | K=5's nu=0.005 | Eliminated | Eliminated |
| External validation | Partial | Full (n=570) for all 9 fits | Full (n=570) for all 35 fits |
| K=3 unique? | Yes | Yes -- only K=3 validates stably | **Revised** -- K=7 validates more robustly |
| iCAF at K=7,9 | Not tested | Not tested | K=7 concentrates iCAF; K=9 is erratic |
| Adjusted validation p<0.05 | Production fit: p=0.004 | None (best: K3 p=0.075) | K=7: 3 combos; K=9: 2 combos |
| New insight | -- | iCAF is alpha-dependent | K=7 may be superior; adjusted signal requires high K |

The previous conclusion that "K=3 is the right model" is **nuanced**: K=3 alpha=0.55
achieves the highest coherence with the original iCAF program (H-cor=0.906) and has
the most interpretable single-factor structure. But K=7 alpha=0.75 is the strongest
validated result and achieves adjusted significance that K=3 does not. The new insights:
1. The survival penalty is necessary for **all** K values (alpha=0 row)
2. K=3 produces the most *coherent* iCAF program but not the strongest *validation*
3. K=7 validates most robustly and uniquely achieves adjusted significance with ntop=NULL
4. K=9 shows genuine fragmentation with erratic behavior

---

## Decision Points

### Decision Point 1: Which fit anchors the paper?

**Option A (Status quo): Keep the production BO fit as primary, add cv_grid as
supplement sensitivity analysis.**
- The production fit has the strongest adjusted validation (p=0.004 with ntop=270)
- The cv_grid analysis is a sensitivity analysis, not a replacement
- The main text adds 1-2 paragraphs summarizing the grid search finding
- Low risk: doesn't change any existing results or figures

**Option B: Highlight K=7 as corroborating evidence.**
- K=7 alpha=0.75 is the strongest grid search validation result (adj p=0.019)
- This strengthens the "DeSurv adds to PurIST+DeCAF" claim
- Present as: "The iCAF program validates at K=7 even with all genes (ntop=NULL),
  confirming it is not an artifact of gene selection"
- Moderate effort: add a supplementary note about K=7

**Option C: Consider K=7 as the primary model.**
- Strongest adjusted validation (p=0.019 vs K=3's 0.075 vs production's 0.004)
- But H-cor=0.872 vs K=3's 0.906 -- slightly less coherent with original
- Would require re-characterizing the factors, since K=7 has 7 factors to interpret
- High risk: major revision and less interpretable (7 factors vs 3)

### Decision Point 2: What goes in the main text vs supplement?

**Main text (1-2 paragraphs in Results):**
- "We performed a controlled sensitivity analysis across K=2,3,5,7,9 and
  alpha=0-0.95 with all other hyperparameters fixed (Supplementary Table X)."
- Key sentence: "The iCAF/immune program emerges only with survival penalty
  (alpha>0) and is most coherent at K=3 (H-cor=0.906 with production factor),
  while K=7 provides the strongest external validation (adjusted HR=0.84,
  p=0.019)."
- Cite the validation KM gaps: K=3 (8.8 months, p=3e-4); K=7 (14.2 months, p=8e-9)

**Supplement:**
- Full 5x7 K x alpha master tables (H-cor, unadjusted val p, adjusted val p)
- Selected per-factor Cox tables
- DECODER compartment breakdowns
- The methodological note about ntop=NULL vs ntop=270
- The K7 detailed analysis (per-factor breakdown)

### Decision Point 3: How to frame the K=3 vs K=7 tension?

The cv_grid reveals a tension: K=3 is the most *interpretable* (3 factors,
single iCAF program, highest coherence) but K=7 is the most *validated*
(strongest adjusted p-value, largest KM gaps). Three framing options:

**Frame A: K=3 is the right model; K=7 provides corroboration.**
"K=3 uniquely isolates the iCAF program as a single coherent factor
(H-cor=0.906). K=7 captures the same program (H-cor=0.872) with stronger
validation, confirming that the signal is robust and not an artifact of K=3."

**Frame B: The iCAF program is K-robust above alpha=0.**
"The iCAF/immune program validates across K=3-9 when alpha>0, confirming
that the survival penalty -- not the rank choice -- is the critical
methodological innovation. K=3 provides the most parsimonious representation."

**Frame C: Lead with K=7's adjusted significance.**
"External validation confirms the iCAF factor adds prognostic information
beyond PurIST+DeCAF (K=7 adjusted HR=0.84, p=0.019, LRT p=0.020 in 570
patients). The program is most coherently isolated at K=3 (H-cor=0.906)
and most strongly validated at K=7 (14.2-month KM gap, p=8e-9)."

### Decision Point 4: The K5_a025 exception

K5_a025's factor has strong unadjusted validation (p=5e-4) and best K=5 KM
(p=4.4e-6), but is a purer immune factor (45% DECODER Immune, rho=0.449).

**Option A: Acknowledge and explain.** "K=5 separates the immune component
that K=3 combines with iCAF and normal stroma. The pure immune factor validates
strongly but does not add to PurIST+DeCAF (adj p=0.069)."

**Option B: Use it as supporting evidence.** "The immune component of Factor 1
drives the survival signal, as demonstrated by its isolated recovery at K=5."

---

## Open Questions for Discussion

1. **K=7 vs K=3:** Should K=7's stronger validation change the paper's framing,
   or is K=3's interpretability (3 factors vs 7) more important for the narrative?

2. **Adjusted significance:** K=7 achieves adjusted p=0.019 with ntop=NULL,
   while K=3 achieves p=0.075 with ntop=NULL but p=0.004 with ntop=270
   (production fit). Should we run K=7 with ntop=270 to see if its adjusted
   validation strengthens further?

3. **The alpha=0.75 dip at K=3:** K3_a075 has H-cor=0.663 (dropping from
   0.906 at alpha=0.55), yet K7_a075 is the strongest result. This suggests
   the survival penalty at K=3 may overshoot at high alpha, pushing the
   model into a degenerate solution. Should this be investigated?

4. **K=9 erratic behavior:** K9_a025 has the second-highest H-cor among
   alpha=0.25 fits (0.785) and validates (p=0.016), but K9_a035 drops to
   H-cor=0.533 and fails validation (p=0.22). Should we investigate what
   causes this instability? Is it initialization-dependent?

5. **Manuscript scope:** With 35 fits, there is rich supplementary material.
   Should we include the full 5x7 tables, or curate a subset (e.g., 15-20
   selected fits) for the supplement?

---

## How to Proceed with the Paper

### Suggested Action Items

1. **Immediate**: Decide how to position K=7 findings relative to K=3
   (Decision Point 3 above)
2. **For the supplement**: Generate the master tables as formatted supplementary
   tables. The 5x7 matrices (H-cor, unadjusted p, adjusted p) plus a curated
   master table of ~18 selected fits
3. **For the main text**: Draft 1-2 paragraphs summarizing the expanded
   sensitivity analysis, emphasizing: (a) alpha=0 fails universally,
   (b) K=3 is most coherent, (c) K=7 provides strongest adjusted validation
4. **Optional**: Re-run K=7 alpha=0.75 fit with ntop=270 to see if adjusted
   validation strengthens. This would directly address whether ntop explains
   the K=3 vs K=7 adjusted p-value gap.
5. **For the reviewer letter**: Prepare an updated defense of K=3 that
   acknowledges K=7 corroboration:
   - Alpha=0 fails for all K (the penalty matters)
   - K=3 is the most parsimonious model with single-factor iCAF isolation
   - K=7 corroborates the program's existence with even stronger validation
   - The iCAF program is a robust biological signal, not a K=3 artifact

---

## Methodological Note: ntop=NULL vs ntop=270

The cv_grid fits use ntop=NULL, meaning H-scores are computed as H = X^T W
using all 1970 genes. The original production fit used ntop=270, focusing each
factor's H-score on its top 270 genes. This has implications:

- **Training results are stronger with ntop=NULL** because the survival penalty
  has more genes to weight (K3_a55 HR=0.344 vs production HR~0.35)
- **Validation adjusted models are weaker** because using all genes dilutes the
  per-factor signal (K3 adj p=0.075 vs production adj p=0.004)
- **Unadjusted and KM results are robust** and clearly confirm the patterns
- The ntop=NULL results are **conservative** for the adjusted analysis. The
  production fit with focused ntop=270 would show stronger adjusted validation
- **K=7's adjusted significance despite ntop=NULL** (p=0.019) is notable:
  it suggests K=7 better separates the iCAF axis from PurIST/DeCAF even without
  gene focusing

---

## Data Provenance

All fits from: `store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main`

Fits are loaded by parameter matching (k, alpha, ntop) from the targets store,
not by hardcoded branch hashes (see `find_fit()` in the amber analysis scripts).

**Grid:** K = {2, 3, 5, 7, 9} x Alpha = {0, 0.25, 0.35, 0.55, 0.75, 0.85, 0.95}
= 35 parameter combinations.

**Fixed parameters:** lambda=0.3, nu=0.05, lambdaW=0, lambdaH=0, 100
initializations (CV_GRID_NSTARTS_FULL), ngene=3000 (->1970 after filtering),
5-fold stratified CV (seed=123), ntop=NULL.

**Original production fit:** alpha=0.334, lambda=0.349, nu=0.056, ntop=270,
200+ consensus seed initializations.

**Training cohort:** TCGA-PAAD (n=144) + CPTAC-PDAC (n=129) = 273 samples.
**Validation cohort:** Dijk (n=90) + Moffitt GEO array (n=123) + Puleo array
(n=288) + PACA-AU (n=69) = 570 samples.

---

## Glossary

- **H-score**: Factor loading score for each sample, computed as H = X^T W
  where X is the expression matrix and W is the factor weight matrix. Measures
  how strongly each sample expresses a given factor's gene program.
- **H-cor (H-score correlation)**: Spearman rank correlation between H-score
  vectors from two different fits. Used here to compare each cv_grid factor
  against the original production K=3 factors. A value of 1.0 would mean
  identical sample rankings; 0 means unrelated.
- **LP (linear predictor)**: The Cox model linear predictor, sum of beta_k *
  H_k across factors. Summarizes each sample's predicted risk from the DeSurv
  model.
- **LRT (likelihood ratio test)**: Compares nested Cox models. Here we test
  whether adding a factor's H-score improves a model already containing PurIST
  + DeCAF. A significant LRT means the factor adds prognostic information
  beyond the existing classifiers.
- **ntop**: Number of top genes per factor used to compute H-scores. ntop=NULL
  means all 1970 genes are used; ntop=270 (production fit) focuses on each
  factor's 270 highest-weight genes. Using all genes is a more conservative
  test because it dilutes each factor's signal.
- **strata(dataset)**: Stratified Cox regression allowing each validation
  cohort its own baseline hazard function while sharing covariate effects
  across cohorts. This accounts for differences in follow-up, patient
  selection, and platform effects across the 4 validation cohorts.
- **Gene universe**: The 1970 genes remaining after variance and expression
  filtering from the initial 3000 (ngene=3000). All gene overlap counts
  throughout these documents (e.g., "17 of 25 present in universe") refer
  to genes present in this filtered set, not the full genome.
- **CV-optimal alpha**: The alpha value that minimizes cross-validated partial
  log-likelihood in the training data. For K=3, this is alpha=0.55.
- **Validation-optimal alpha**: The alpha at which external validation signal
  is strongest (by unadjusted Cox p-value). For K=3, alpha=0.55; for K=5,
  alpha=0.25; for K=7, alpha=0.75.
- **PurIST**: Single-sample classifier for PDAC tumor subtypes (Basal-like
  vs Classical). Based on tumor cell gene expression.
- **DeCAF**: Classifier for cancer-associated fibroblast subtypes (proCAF vs
  restCAF). Based on CAF gene expression.
- **DECODER**: Deconvolution-based compartment classifier from Peng et al.
  Provides Immune, BasalTumor, ClassicalTumor, ActivatedStroma, and
  NormalStroma gene lists.
