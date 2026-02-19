# K=7 DeSurv Factor Analysis: CV Grid Fits (Expanded Grid)

**For discussion with Jen Jen Yeh**
**Date: 2026-02-18 (amber revision)**
**Companion to: 2026-02-17-k-sensitivity-synthesis_amber.md**
**Data source: cv_grid exhaustive search (fixed lambda=0.3, nu=0.05, 100 inits)**

**Note for Amber:** This is a new document -- K=7 was not examined in the original
3-K analysis. Seven fits are examined: alpha=0, 0.25, 0.35, 0.55, 0.75, 0.85, 0.95.
Fits loaded by parameter matching. Scripts: `inst/cv_grid_training_analysis_amber.R`
and `inst/cv_grid_validation_analysis_amber.R`.

## Motivation for K=7

K=5 showed that higher rank fragments the iCAF program across multiple factors.
K=7 tests whether this fragmentation continues. With 7 factors, the model has
even more capacity to split biological programs into fine-grained components.

**Expectations:**
- H-cor with original F1 should be lower than K=5 (more fragmentation)
- More factors should be prognostic in training (survival signal dispersed)
- Validation should be weaker or absent (individual factors too dilute)
- Possibly one very pure immune or iCAF sub-factor at moderate alpha

**Training cohort:** TCGA + CPTAC (n=273)
**Validation cohort:** 4 independent datasets (n=570 pooled):
Dijk (n=90), Moffitt GEO array (n=123), Puleo array (n=288), PACA-AU (n=69).
All validation Cox models use `strata(dataset)`.

**iCAF reference**: The "iCAF program" refers to the factor best matching the
original K=3 Factor 1, characterized by overlap with the **Elyada et al. (2019)
iCAF 35-gene signature** (25 of 35 present in the 1970-gene universe).

---

## Summary Table

| | a=0.00 | a=0.25 | a=0.35 | a=0.55 | a=0.75 | a=0.85 | a=0.95 |
|--|--------|--------|--------|--------|--------|--------|--------|
| **Best iCAF factor** | F5 | F3 | F6 | F7 | F2 | F6 | F1 |
| **H-cor with orig F1** | 0.510 | 0.784 | 0.758 | 0.873 | 0.872 | 0.870 | 0.814 |
| **Top-270 overlap** | 68/270 | 117/270 | 125/270 | 119/270 | 133/270 | 128/270 | 112/270 |
| **Elyada iCAF overlap** | 9/25 (36%) | 7/25 (28%) | 4/25 (16%) | 5/25 (20%) | 4/25 (16%) | 4/25 (16%) | 5/25 (20%) |
| **DECODER Immune** | 16/60 (26.7%) | 22/60 (36.7%) | 10/60 (16.7%) | 15/60 (25.0%) | 9/60 (15.0%) | 16/60 (26.7%) | 13/60 (21.7%) |
| **DECODER Basal** | 0/139 (0.0%) | 8/139 (5.8%) | 12/139 (8.6%) | 16/139 (11.5%) | 17/139 (12.2%) | 19/139 (13.7%) | 21/139 (15.1%) |
| **Train HR (unadj)** | 0.850 (p=0.067) | 0.557 (p<0.0001) | 0.397 (p<0.0001) | 0.085 (p<0.0001) | 0.108 (p<0.0001) | 0.074 (p<0.0001) | 0.056 (p<0.0001) |
| **Val HR (unadj)** | 0.885 (p=0.096) | 0.890 (p=0.138) | 0.844 (p=0.013) | 0.816 (p=0.005) | **0.771 (p=0.0001)** | 0.812 (p=0.002) | 0.779 (p=0.003) |
| **Val HR (adj)** | 1.025 (p=0.762) | 0.912 (p=0.265) | 0.856 (p=0.033) | 0.885 (p=0.104) | **0.844 (p=0.020)** | 0.880 (p=0.075) | 0.840 (p=0.046) |
| **Val LRT p** | 0.761 | 0.269 | **0.034** | 0.105 | **0.020** | 0.075 | **0.047** |
| **Val KM p** | 5.5e-5 | 9.7e-3 | 1.0e-2 | **2.0e-7** | **8.0e-9** | **1.4e-6** | **2.4e-7** |
| **Val KM gap** | 9.4 mo | 4.5 mo | 6.2 mo | **13.2 mo** | **14.2 mo** | **11.0 mo** | **13.1 mo** |
| **# train prognostic** | 2/7 | 3/7 | 2/7 | 3/7 | 3/7 | 2/7 | 3/7 |

**Key pattern:** H-cor plateaus around 0.87 for alpha >= 0.55, confirming K=7 can
recover the iCAF program at high alpha. The dominant iCAF factor migrates across
factor indices (F5 -> F3 -> F6 -> F7 -> F2 -> F6 -> F1) as alpha increases,
reflecting instability in which factor "captures" the program. Training HRs become
extreme (0.056-0.108) at high alpha, far exceeding K=3's 0.344, raising overfitting
concerns. Yet validation KM separation is remarkably strong (p < 1e-6 for alpha >= 0.55).
The adjusted Cox validation is significant only at alpha=0.35, 0.75, and 0.95.

---

## K7_a000: Standard NMF (alpha=0.00)

### H-score correlation with original K=3

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | 0.480   | **-0.826** | **0.773** |
| F2     | -0.475  | **0.777** | -0.493 |
| F3     | -0.020  | **-0.793** | **0.903** |
| F4     | -0.056  | 0.670   | -0.569 |
| **F5** | **0.510** | -0.741 | 0.632  |
| F6     | -0.144  | **0.930** | **-0.877** |
| F7     | 0.027   | 0.298   | -0.251 |

Best iCAF match is F5 (rho=0.510). No factor strongly captures the iCAF program.
The dominant correlations are with original F2 and F3 (stroma/basal axis). F6
perfectly tracks original F2 (rho=0.930) and F3 matches original F3 (rho=0.903).

### Reference signature overlaps (F5 as best iCAF match, top-270)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 9/25 | 36.0% |
| Elyada myCAF (15) | 3/15 | 20.0% |
| SCISSORS iCAF (17) | 5/17 | 29.4% |
| SCISSORS myCAF (16) | 0/16 | 0.0% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 16/60 | 26.7% |
| DECODER ActivatedStroma (108) | 14/108 | 13.0% |
| DECODER NormalStroma (15) | 6/15 | 40.0% |
| DECODER BasalTumor (139) | 0/139 | 0.0% |
| DECODER ClassicalTumor (166) | 3/166 | 1.8% |

F5 is the cleanest iCAF-like factor at alpha=0: strong Elyada iCAF (36%),
NormalStroma (40%), zero Basal Tumor. But its H-cor with original F1 is only 0.510.

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 2/25 | 0/60 | 4/139 | Mixed |
| F2 | 1/25 | 1/60 | 21/139 | Basal Tumor |
| F3 | 2/25 | 0/60 | 97/139 | **Basal Tumor (dominant)** |
| F4 | 0/25 | 1/60 | 3/139 | Uncharacterized |
| **F5** | **9/25** | **16/60** | **0/139** | **iCAF/Immune** |
| F6 | 8/25 | 38/60 | 0/139 | **Immune (purest)** |
| F7 | 3/25 | 2/60 | 0/139 | Mixed |

Notably, F6 has the most DECODER Immune genes (38/60 = 63%) but does not correlate
with original F1 (rho=-0.144). At K=7 alpha=0, the immune program forms its own
separate factor (F6) distinct from the iCAF program (F5).

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.832 | 0.046 | **0.800** | **0.0003** |
| **F2** | **1.232** | **0.026** | **1.275** | **0.0001** |
| F3 | 1.157 | 0.095 | **1.238** | **0.001** |
| F4 | 0.917 | 0.301 | **0.830** | **0.009** |
| F5 (iCAF) | 0.850 | 0.067 | 0.885 | 0.096 |
| F6 | 0.951 | 0.559 | 0.999 | 0.987 |
| F7 | 0.953 | 0.578 | 0.982 | 0.733 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.952 | 0.604 | **0.827** | **0.006** |
| F2 | 0.985 | 0.891 | 1.101 | 0.187 |
| F3 | 0.977 | 0.801 | 1.060 | 0.433 |
| F4 | 0.955 | 0.618 | **0.790** | **0.003** |
| F5 (iCAF) | 0.914 | 0.347 | 1.025 | 0.762 |
| F6 | 0.972 | 0.748 | 1.004 | 0.938 |
| F7 | 1.071 | 0.451 | 1.035 | 0.522 |

**The iCAF factor (F5) is not prognostic in training or validation.** Instead, F1
and F2 are the prognostic pair in training. In validation, 4 of 7 factors reach
significance unadjusted, and F1 and F4 survive PurIST+DeCAF adjustment -- but
neither is the iCAF factor. Standard NMF at K=7 finds prognostic structure but
it does not align with the iCAF program.

**LRT adding F5 (iCAF) to PurIST + DeCAF:** Training p=0.346, Validation p=0.761

**Median-split KM for F5 (iCAF):**
- Training (n=273): High 20.9 vs Low 20.3 months, p=0.321
- Validation (n=570): High 26.8 vs Low 17.4 months, **p=5.5e-5**

The validation KM is significant despite the Cox model failure. This suggests the
iCAF factor has a threshold effect (high vs low separation) even though its
continuous HR is not significant. But the training KM shows no separation (0.6-month
gap), confirming this is not a reliable iCAF signal.

### Beta structure

F1=-7.75e-05, F2=6.88e-05, F3=5.63e-05, F4=-9.75e-05, F5=-1.08e-04, F6=-4.04e-05, F7=1.77e-06

All betas are tiny (no survival penalty). No factor dominates the LP.

---

## K7_a025: alpha=0.25

### H-score correlation with original K=3

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.198  | **0.852** | -0.601 |
| F2     | -0.189  | **0.810** | -0.524 |
| **F3** | **0.784** | -0.402 | 0.435  |
| F4     | -0.029  | 0.407   | -0.067 |
| F5     | -0.236  | -0.394  | **0.757** |
| F6     | -0.067  | **0.793** | -0.522 |
| F7     | 0.434   | -0.192  | 0.387  |

Best iCAF match is F3 (rho=0.784). F1, F2, and F6 all track original F2
(rho=0.852, 0.810, 0.793), showing the stroma axis has fragmented across three
K=7 factors. F5 matches original F3 (rho=0.757).

### Reference signature overlaps (F3 as best iCAF match, top-270)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 7/25 | 28.0% |
| Elyada myCAF (15) | 1/15 | 6.7% |
| SCISSORS iCAF (17) | 6/17 | 35.3% |
| SCISSORS myCAF (16) | 0/16 | 0.0% |
| DeCAF restCAF (8) | 2/8 | 25.0% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 22/60 | 36.7% |
| DECODER ActivatedStroma (108) | 5/108 | 4.6% |
| DECODER NormalStroma (15) | 6/15 | 40.0% |
| DECODER BasalTumor (139) | 8/139 | 5.8% |
| DECODER ClassicalTumor (166) | 41/166 | 24.7% |

F3 has the highest Immune overlap at any K=7 alpha (36.7%), plus strong SCISSORS
iCAF (35.3%) and NormalStroma (40%). But it also absorbs substantial Classical
Tumor genes (24.7%), diluting the iCAF purity.

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 0/25 | 6/60 | 8/139 | Mixed |
| F2 | 2/25 | 16/60 | 24/139 | Immune + Basal mix |
| **F3** | **7/25** | **22/60** | **8/139** | **iCAF/Immune** |
| F4 | 3/25 | 1/60 | 16/139 | Basal |
| F5 | 1/25 | 3/60 | 68/139 | **Basal Tumor (dominant)** |
| F6 | 5/25 | 7/60 | 13/139 | Mixed |
| F7 | 3/25 | 1/60 | 13/139 | Mixed/Basal |

The Immune genes are concentrated in F3 (22/60) and F2 (16/60). The iCAF program
is partially captured by F3 but the immune component is leaking into F2.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.993 | 0.936 | 1.047 | 0.571 |
| F2 | 1.012 | 0.884 | 1.142 | 0.077 |
| **F3 (iCAF)** | **0.557** | **<0.0001** | 0.890 | 0.138 |
| F4 | 1.009 | 0.914 | 1.124 | 0.155 |
| **F5** | **1.387** | **0.0004** | **1.326** | **0.0003** |
| F6 | 0.957 | 0.588 | 1.020 | 0.781 |
| **F7** | **0.815** | **0.009** | 0.884 | 0.098 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.974 | 0.768 | 0.944 | 0.495 |
| F2 | 0.927 | 0.392 | 1.040 | 0.608 |
| **F3 (iCAF)** | **0.612** | **<0.0001** | 0.912 | 0.265 |
| F4 | 0.974 | 0.759 | 1.026 | 0.754 |
| F5 | 1.130 | 0.244 | 1.121 | 0.187 |
| F6 | 0.958 | 0.623 | 0.941 | 0.407 |
| F7 | 0.876 | 0.131 | 0.854 | 0.049 |

F3 (iCAF) has a strong training signal (HR=0.557, p<0.0001) that survives
PurIST+DeCAF adjustment (HR=0.612, p<0.0001). However, **validation fails
completely**: unadjusted p=0.138, adjusted p=0.265. Three factors are prognostic
in training (F3, F5, F7), but only F5 validates unadjusted (p=0.0003). The iCAF
factor does not generalize.

**LRT adding F3 (iCAF) to PurIST + DeCAF:** Training p=7.7e-7, Validation p=0.269

**Median-split KM for F3 (iCAF):**
- Training (n=273): High 25.4 vs Low 15.2 months, **p=9.4e-6**
- Validation (n=570): High 24.3 vs Low 19.8 months, **p=9.7e-3**

The validation KM (4.5-month gap, p=0.010) validates despite the Cox failure,
mirroring the K5_a25 pattern where the median-split captures signal that the
continuous HR misses.

### Beta structure

F1=-3.61e-07, F2=0.00e+00, F3=**-1.19e-03** (dominant), F4=1.25e-04, F5=**6.75e-04**, F6=2.16e-05, F7=-8.85e-05

F3 (iCAF) dominates the LP, but F5 also has substantial weight.

---

## K7_a035: alpha=0.35

### H-score correlation with original K=3

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | 0.270   | **-0.936** | **0.925** |
| F2     | -0.170  | -0.335  | 0.609  |
| F3     | -0.150  | **0.924** | **-0.810** |
| F4     | 0.208   | 0.231   | -0.251 |
| F5     | -0.421  | **0.865** | -0.652 |
| **F6** | **0.758** | -0.393 | 0.347  |
| F7     | -0.068  | -0.156  | 0.359  |

Best iCAF match is F6 (rho=0.758). F1 matches original F3 (rho=0.925), F3
matches original F2 (rho=0.924), and F5 also tracks original F2 (rho=0.865).
The stroma axis is splitting across F3 and F5.

### Reference signature overlaps (F6 as best iCAF match, top-270)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 4/25 | 16.0% |
| Elyada myCAF (15) | 1/15 | 6.7% |
| SCISSORS iCAF (17) | 1/17 | 5.9% |
| SCISSORS myCAF (16) | 0/16 | 0.0% |
| DeCAF restCAF (8) | 0/8 | 0.0% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 10/60 | 16.7% |
| DECODER ActivatedStroma (108) | 0/108 | 0.0% |
| DECODER NormalStroma (15) | 3/15 | 20.0% |
| DECODER BasalTumor (139) | 12/139 | 8.6% |
| DECODER ClassicalTumor (166) | 34/166 | 20.5% |

F6 has relatively low signature overlap across the board: only 16% iCAF, 5.9%
SCISSORS iCAF, 16.7% Immune. Despite having the highest H-cor with original F1
(0.758), the gene content is diluted with Classical Tumor genes (20.5%) and Basal
Tumor (8.6%).

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 0/25 | 0/60 | 26/139 | Basal Tumor |
| F2 | 3/25 | 0/60 | 18/139 | Basal Tumor |
| F3 | 1/25 | 13/60 | 3/139 | Immune |
| F4 | 5/25 | 10/60 | 2/139 | iCAF/Immune mix |
| F5 | 4/25 | 23/60 | 11/139 | **Immune (most)** |
| **F6** | **4/25** | **10/60** | **12/139** | **iCAF/Mixed** |
| F7 | 10/25 | 11/60 | 18/139 | Mixed (highest iCAF overlap!) |

Strikingly, F7 has the most Elyada iCAF genes (10/25 = 40%) but is not the best
H-cor match. The Immune program is split between F5 (23/60) and F3 (13/60). The
iCAF program is scattered across F4, F6, and F7.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.990 | 0.913 | 0.961 | 0.522 |
| **F2** | **1.529** | **0.0001** | 1.058 | 0.439 |
| F3 | 0.970 | 0.718 | 0.934 | 0.308 |
| F4 | 0.869 | 0.115 | 0.926 | 0.146 |
| F5 | 1.166 | 0.093 | **1.225** | **0.001** |
| **F6 (iCAF)** | **0.397** | **<0.0001** | **0.844** | **0.013** |
| F7 | 1.118 | 0.193 | **1.253** | **0.0002** |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.977 | 0.799 | 0.924 | 0.233 |
| **F2** | **1.537** | **0.0001** | 0.958 | 0.588 |
| F3 | 0.979 | 0.816 | 0.893 | 0.109 |
| F4 | 1.037 | 0.709 | 1.013 | 0.813 |
| F5 | 1.010 | 0.922 | 1.080 | 0.269 |
| **F6 (iCAF)** | **0.365** | **<0.0001** | **0.856** | **0.033** |
| F7 | 0.908 | 0.321 | 1.134 | 0.053 |

F6 (iCAF) is the standout: extremely strong training (HR=0.397, p<0.0001) that
**validates both unadjusted (p=0.013) and adjusted (p=0.033)**. This is one of
only three alpha values where the iCAF factor survives PurIST+DeCAF adjustment in
validation. F2 is strongly prognostic in training (HR=1.529) but does not validate
(p=0.439).

**LRT adding F6 (iCAF) to PurIST + DeCAF:** Training p=1.4e-20, **Validation p=0.034**

**Median-split KM for F6 (iCAF):**
- Training (n=273): High 33.4 vs Low 12.4 months, **p=6.7e-16**
- Validation (n=570): High 26.2 vs Low 20.0 months, **p=1.0e-2**

Training KM shows a massive 21-month gap. Validation KM confirms with a 6.2-month
gap. This is the lowest alpha at which K=7 achieves adjusted validation significance.

### Beta structure

F1=1.50e-05, F2=**9.17e-04**, F3=-2.67e-06, F4=-1.01e-04, F5=0.00e+00, F6=**-1.86e-03** (dominant), F7=1.54e-05

F6 (iCAF) dominates, with F2 as secondary.

---

## K7_a055: alpha=0.55

### H-score correlation with original K=3

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | 0.362   | **-0.926** | **0.890** |
| F2     | -0.445  | **0.810** | -0.513 |
| F3     | 0.051   | 0.665   | -0.601 |
| F4     | -0.082  | **0.897** | **-0.824** |
| F5     | 0.180   | -0.140  | 0.134  |
| F6     | -0.409  | -0.290  | 0.577  |
| **F7** | **0.873** | -0.095 | -0.093 |

Best iCAF match is F7 (rho=0.873). F7 is remarkable: it correlates strongly with
original F1 (0.873) but has near-zero correlation with original F2 (-0.095) and F3
(-0.093). This is the purest recovery of the iCAF program at K=7 -- isolated from
the stroma/basal axis. F1 captures original F3 (0.890), F2 and F4 capture original
F2 (0.810 and 0.897).

### Reference signature overlaps (F7 as best iCAF match, top-270)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 5/25 | 20.0% |
| Elyada myCAF (15) | 1/15 | 6.7% |
| SCISSORS iCAF (17) | 5/17 | 29.4% |
| SCISSORS myCAF (16) | 2/16 | 12.5% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | 1/4 | 25.0% |
| DECODER Immune (60) | 15/60 | 25.0% |
| DECODER ActivatedStroma (108) | 3/108 | 2.8% |
| DECODER NormalStroma (15) | 6/15 | 40.0% |
| DECODER BasalTumor (139) | 16/139 | 11.5% |
| DECODER ClassicalTumor (166) | 17/166 | 10.2% |

Moderate overlaps across the board. NormalStroma (40%) and SCISSORS iCAF (29.4%)
are the strongest markers. The 11.5% Basal Tumor contamination shows the factor is
absorbing some tumor genes.

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 0/25 | 0/60 | 11/139 | Basal |
| F2 | 1/25 | 0/60 | 16/139 | Basal Tumor |
| F3 | 13/25 | 18/60 | 0/139 | **iCAF/Immune (purest)** |
| F4 | 1/25 | 23/60 | 0/139 | **Immune** |
| F5 | 3/25 | 4/60 | 4/139 | Mixed |
| F6 | 2/25 | 1/60 | 75/139 | **Basal Tumor (dominant)** |
| **F7** | **5/25** | **15/60** | **16/139** | **iCAF/Mixed** |

The immune program is splitting: F3 gets 18/60 Immune + 13/25 Elyada iCAF (the
most iCAF-enriched factor!), F4 gets 23/60 Immune with zero Basal. F7 (the best
H-cor match) gets only 15/60 Immune and 5/25 Elyada iCAF. The program content is
more coherently captured by F3, but the H-score pattern is best captured by F7.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.920 | 0.348 | 0.885 | 0.043 |
| **F2** | **1.219** | **0.034** | **1.219** | **0.002** |
| F3 | 0.942 | 0.465 | 1.018 | 0.752 |
| F4 | 0.949 | 0.536 | 0.865 | 0.035 |
| F5 | 0.950 | 0.564 | 0.947 | 0.303 |
| **F6** | **1.517** | **<0.0001** | **1.390** | **<0.0001** |
| **F7 (iCAF)** | **0.085** | **<0.0001** | **0.816** | **0.005** |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.944 | 0.527 | 0.878 | 0.049 |
| F2 | 0.997 | 0.981 | 1.054 | 0.463 |
| F3 | 0.976 | 0.786 | 1.024 | 0.678 |
| F4 | 0.980 | 0.825 | 0.853 | 0.028 |
| F5 | 1.114 | 0.253 | 1.018 | 0.743 |
| F6 | 1.163 | 0.205 | **1.221** | **0.011** |
| **F7 (iCAF)** | **0.076** | **<0.0001** | 0.885 | 0.104 |

F7 (iCAF) has the most extreme training HR in the entire grid (0.085), but
**validation adjusted p=0.104 -- borderline failure**. Three factors are prognostic
in training (F2, F6, F7), and in validation F6 (Basal axis) is actually the only
factor to survive PurIST+DeCAF adjustment (p=0.011). The iCAF factor's training
signal is 12x stronger than K=3's (HR=0.085 vs 0.344), a hallmark of overfitting.

**LRT adding F7 (iCAF) to PurIST + DeCAF:** Training p=1.2e-64, Validation p=0.105

**Median-split KM for F7 (iCAF):**
- Training (n=273): High 44.4 vs Low 11.6 months, **p<1e-16**
- Validation (n=570): High 30.2 vs Low 17.0 months, **p=2.0e-7**

The validation KM is strikingly strong (13.2-month gap, p=2.0e-7), despite the
Cox model borderline. This pattern -- strong KM, weak Cox -- suggests the factor
captures a threshold-based survival difference but its continuous score does not
add beyond PurIST+DeCAF.

### Beta structure

F1=3.29e-07, F2=0.00e+00, F3=3.17e-05, F4=-8.78e-07, F5=0.00e+00, F6=9.57e-05, F7=**-2.93e-03** (dominant)

F7 completely dominates the LP. Multiple betas are exactly zero, meaning those
factors are penalized out of the survival model.

---

## K7_a075: alpha=0.75

### H-score correlation with original K=3

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | 0.207   | 0.621   | -0.631 |
| **F2** | **0.872** | -0.055 | -0.144 |
| F3     | -0.348  | **0.946** | **-0.765** |
| F4     | -0.397  | -0.171  | 0.502  |
| F5     | 0.244   | 0.591   | -0.555 |
| F6     | 0.148   | 0.203   | -0.192 |
| F7     | 0.253   | **-0.940** | **0.943** |

Best iCAF match is F2 (rho=0.872). Like K7_a055's F7, F2 is nearly orthogonal
to the stroma axis (orig F2 rho=-0.055, orig F3 rho=-0.144). F3 captures original
F2 (0.946), F7 captures original F3 (0.943). Clean separation of the three
original programs.

### Reference signature overlaps (F2 as best iCAF match, top-270)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 4/25 | 16.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | 4/17 | 23.5% |
| SCISSORS myCAF (16) | 1/16 | 6.2% |
| DeCAF restCAF (8) | 0/8 | 0.0% |
| DeCAF proCAF (4) | 1/4 | 25.0% |
| DECODER Immune (60) | 9/60 | 15.0% |
| DECODER ActivatedStroma (108) | 3/108 | 2.8% |
| DECODER NormalStroma (15) | 5/15 | 33.3% |
| DECODER BasalTumor (139) | 17/139 | 12.2% |
| DECODER ClassicalTumor (166) | 35/166 | 21.1% |

F2 has the lowest DECODER Immune overlap of any iCAF factor at alpha >= 0.25
(15.0%). It absorbs substantial Classical Tumor (21.1%) and Basal Tumor (12.2%).
Despite high H-cor with original F1, the gene content is increasingly diluted.

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 9/25 | 21/60 | 2/139 | **iCAF/Immune** |
| **F2** | **4/25** | **9/60** | **17/139** | **iCAF (H-cor match)/Mixed** |
| F3 | 1/25 | 11/60 | 3/139 | Immune |
| F4 | 3/25 | 0/60 | 67/139 | **Basal Tumor (dominant)** |
| F5 | 6/25 | 17/60 | 0/139 | Immune/iCAF |
| F6 | 1/25 | 2/60 | 4/139 | Uncharacterized |
| F7 | 1/25 | 0/60 | 28/139 | Basal Tumor |

The Immune genes are split three ways: F1 (21/60), F5 (17/60), F3 (11/60). The
Elyada iCAF genes go to F1 (9/25) and F5 (6/25). F2 (the best H-cor match) has
surprisingly low immune/iCAF gene content, suggesting it captures the iCAF
H-score pattern through a different gene set.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1** | **0.788** | **0.003** | 0.998 | 0.975 |
| **F2 (iCAF)** | **0.108** | **<0.0001** | **0.771** | **0.0001** |
| F3 | 1.078 | 0.396 | 1.120 | 0.062 |
| **F4** | **1.416** | **0.0002** | **1.340** | **<0.0001** |
| F5 | 0.921 | 0.312 | 0.888 | 0.042 |
| F6 | 0.890 | 0.188 | 0.927 | 0.160 |
| F7 | 1.025 | 0.777 | 0.994 | 0.918 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.819 | 0.023 | 1.047 | 0.540 |
| **F2 (iCAF)** | **0.086** | **<0.0001** | **0.844** | **0.020** |
| F3 | 1.013 | 0.890 | 1.011 | 0.864 |
| F4 | 1.040 | 0.753 | 1.118 | 0.196 |
| F5 | 1.097 | 0.321 | 0.924 | 0.195 |
| F6 | 1.037 | 0.699 | 0.996 | 0.948 |
| F7 | 1.009 | 0.920 | 0.952 | 0.455 |

**This is the best-validating K=7 fit.** F2 (iCAF) validates both unadjusted
(HR=0.771, p=0.0001) and adjusted (HR=0.844, p=0.020). It is the only factor
to survive PurIST+DeCAF adjustment in validation. Three training factors are
prognostic (F1, F2, F4), but only F2 carries through to validation adjusted.
The training HR (0.108) is still extreme vs validation (0.771), indicating
substantial overfitting, but the validated effect is real.

**LRT adding F2 (iCAF) to PurIST + DeCAF:** Training p=1.1e-61, **Validation p=0.020**

**Median-split KM for F2 (iCAF):**
- Training (n=273): High 44.4 vs Low 11.6 months, **p<1e-16**
- Validation (n=570): High 30.8 vs Low 16.6 months, **p=8.0e-9**

The strongest validation KM in the entire K=7 grid: 14.2-month gap, p=8.0e-9.
The training KM gap (32.8 months) is implausibly large, confirming overfitting,
but the validation effect is robust.

### Beta structure

F1=-2.15e-04, F2=**-2.79e-03** (dominant), F3=1.56e-05, F4=4.46e-05, F5=2.41e-04, F6=-2.60e-05, F7=1.23e-06

F2 dominates overwhelmingly. All other betas are at least 10x smaller.

---

## K7_a085: alpha=0.85

### H-score correlation with original K=3

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.280  | **0.990** | **-0.831** |
| F2     | 0.093   | 0.415   | -0.320 |
| F3     | -0.181  | **0.893** | **-0.712** |
| F4     | -0.014  | -0.365  | 0.618  |
| F5     | 0.683   | -0.229  | 0.269  |
| **F6** | **0.870** | -0.109 | -0.096 |
| F7     | 0.086   | **-0.900** | **0.973** |

Best iCAF match is F6 (rho=0.870). Like the a055 and a075 best matches, F6 is
nearly orthogonal to original F2 (-0.109) and F3 (-0.096). F1 almost perfectly
captures original F2 (0.990), F7 captures original F3 (0.973). Clean three-program
separation persists.

### Reference signature overlaps (F6 as best iCAF match, top-270)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 4/25 | 16.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | 5/17 | 29.4% |
| SCISSORS myCAF (16) | 1/16 | 6.2% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | 1/4 | 25.0% |
| DECODER Immune (60) | 16/60 | 26.7% |
| DECODER ActivatedStroma (108) | 3/108 | 2.8% |
| DECODER NormalStroma (15) | 7/15 | 46.7% |
| DECODER BasalTumor (139) | 19/139 | 13.7% |
| DECODER ClassicalTumor (166) | 32/166 | 19.3% |

F6 has the highest NormalStroma overlap of any K=7 iCAF factor (46.7%) and
moderate Immune (26.7%), SCISSORS iCAF (29.4%). Basal contamination continues
to rise (13.7%).

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 1/25 | 18/60 | 5/139 | Immune |
| F2 | 2/25 | 9/60 | 2/139 | Mixed |
| F3 | 4/25 | 6/60 | 6/139 | Mixed |
| F4 | 0/25 | 3/60 | 18/139 | Basal Tumor |
| F5 | 10/25 | 15/60 | 4/139 | **iCAF/Immune** |
| **F6** | **4/25** | **16/60** | **19/139** | **iCAF (H-cor match)/Mixed** |
| F7 | 0/25 | 0/60 | 59/139 | **Basal Tumor (dominant)** |

F5 has the most Elyada iCAF genes (10/25 = 40%) and good Immune overlap (15/60),
but its H-cor with original F1 is only 0.683. F6 (the H-cor match) splits its
content between Immune (16/60) and Basal (19/139). The iCAF program is again
fragmented between the factor that best matches the H-score pattern (F6) and the
factor with the most iCAF gene content (F5).

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.035 | 0.694 | 1.074 | 0.268 |
| F2 | 0.914 | 0.295 | 0.962 | 0.505 |
| F3 | 1.032 | 0.716 | 0.947 | 0.430 |
| F4 | 1.071 | 0.438 | **1.249** | **0.004** |
| **F5** | **0.673** | **<0.0001** | **0.842** | **0.016** |
| **F6 (iCAF)** | **0.074** | **<0.0001** | **0.812** | **0.002** |
| F7 | 1.151 | 0.109 | **1.176** | **0.014** |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.016 | 0.860 | 0.999 | 0.986 |
| F2 | 1.036 | 0.699 | 1.002 | 0.969 |
| F3 | 1.083 | 0.379 | 0.905 | 0.166 |
| F4 | 0.838 | 0.075 | 1.072 | 0.416 |
| F5 | 0.794 | 0.025 | 0.892 | 0.131 |
| **F6 (iCAF)** | **0.057** | **<0.0001** | 0.880 | 0.075 |
| F7 | 1.059 | 0.526 | 1.038 | 0.609 |

F6 (iCAF) has the most extreme training HR in the entire document (0.074) and
validates unadjusted (p=0.002), but **misses adjusted significance (p=0.075)**.
F5 is also prognostic in training (HR=0.673) but loses significance adjusted in
validation (p=0.131). The iCAF signal is split between F5 and F6, diluting each.

**LRT adding F6 (iCAF) to PurIST + DeCAF:** Training p=2.6e-70, Validation p=0.075

**Median-split KM for F6 (iCAF):**
- Training (n=273): High 44.4 vs Low 11.2 months, **p<1e-16**
- Validation (n=570): High 29.0 vs Low 18.0 months, **p=1.4e-6**

Strong validation KM (11.0-month gap) despite the adjusted Cox borderline,
following the same pattern as K7_a055.

### Beta structure

F1=1.68e-06, F2=0.00e+00, F3=9.02e-05, F4=0.00e+00, F5=0.00e+00, F6=**-3.36e-03** (dominant), F7=5.48e-05

F6 overwhelmingly dominates. Four of seven factors have zero or near-zero betas.

---

## K7_a095: alpha=0.95

### H-score correlation with original K=3

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| **F1** | **0.814** | -0.081 | -0.078 |
| F2     | -0.358  | **0.940** | **-0.746** |
| F3     | -0.051  | **0.825** | **-0.772** |
| F4     | 0.412   | -0.668  | 0.686  |
| F5     | -0.180  | -0.376  | 0.696  |
| F6     | 0.283   | **-0.932** | **0.938** |
| F7     | 0.089   | -0.424  | 0.701  |

Best iCAF match is F1 (rho=0.814). F1 is orthogonal to original F2 (-0.081) and
F3 (-0.078), the cleanest separation in the grid. F2 captures original F2 (0.940),
F6 captures original F3 (0.938). But the stroma axis is now split: F2 and F3 both
track original F2 (0.940 and 0.825).

### Reference signature overlaps (F1 as best iCAF match, top-270)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 5/25 | 20.0% |
| Elyada myCAF (15) | 1/15 | 6.7% |
| SCISSORS iCAF (17) | 5/17 | 29.4% |
| SCISSORS myCAF (16) | 1/16 | 6.2% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 13/60 | 21.7% |
| DECODER ActivatedStroma (108) | 7/108 | 6.5% |
| DECODER NormalStroma (15) | 5/15 | 33.3% |
| DECODER BasalTumor (139) | 21/139 | 15.1% |
| DECODER ClassicalTumor (166) | 23/166 | 13.9% |

F1 has the highest Basal Tumor contamination of any K=7 iCAF factor (15.1%).
At alpha=0.95 the survival penalty is so strong it pulls in any prognostic genes
regardless of compartment. The iCAF purity is lowest here.

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| **F1** | **5/25** | **13/60** | **21/139** | **iCAF (H-cor match)/Mixed** |
| F2 | 3/25 | 9/60 | 9/139 | Mixed |
| F3 | 6/25 | 19/60 | 0/139 | **Immune** |
| F4 | 6/25 | 10/60 | 18/139 | iCAF/Mixed |
| F5 | 0/25 | 5/60 | 38/139 | **Basal Tumor** |
| F6 | 3/25 | 0/60 | 31/139 | Basal Tumor |
| F7 | 3/25 | 1/60 | 15/139 | Basal |

F3 has the most DECODER Immune genes (19/60) with zero Basal -- the purest
immune factor. F4 and F1 share the Elyada iCAF genes (6/25 and 5/25). The iCAF
program is fragmented across F1, F3, and F4.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.056** | **<0.0001** | **0.779** | **0.003** |
| F2 | 1.098 | 0.292 | **1.153** | **0.027** |
| F3 | 0.910 | 0.271 | 0.933 | 0.196 |
| **F4** | **0.758** | **0.001** | 0.998 | 0.982 |
| **F5** | **1.384** | **0.001** | **1.204** | **0.006** |
| F6 | 1.023 | 0.795 | 1.033 | 0.652 |
| F7 | 1.173 | 0.081 | 1.013 | 0.867 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.044** | **<0.0001** | **0.840** | **0.046** |
| F2 | 1.020 | 0.832 | 1.034 | 0.630 |
| F3 | 1.009 | 0.921 | 0.985 | 0.781 |
| **F4** | **0.701** | **0.0001** | 0.917 | 0.325 |
| F5 | 1.146 | 0.204 | 1.031 | 0.679 |
| F6 | 1.026 | 0.772 | 1.003 | 0.969 |
| F7 | 1.174 | 0.106 | 0.927 | 0.348 |

F1 (iCAF) achieves the most extreme training HR in the entire grid (0.056 unadj,
0.044 adj) and **validates adjusted (HR=0.840, p=0.046)**. Three training factors
are prognostic (F1, F4, F5), but only F1 validates. F4 (HR=0.758 in training)
completely fails validation (p=0.982), confirming severe overfitting of non-iCAF
factors.

**LRT adding F1 (iCAF) to PurIST + DeCAF:** Training p=2.1e-79, **Validation p=0.047**

**Median-split KM for F1 (iCAF):**
- Training (n=273): High 44.4 vs Low 11.1 months, **p<1e-16**
- Validation (n=570): High 30.1 vs Low 17.0 months, **p=2.4e-7**

The validation KM is excellent (13.1-month gap), consistent with the Cox validation.
The training median gap (33.3 months) is the largest in the document, but the
training median OS for the High group (44.4 months) is identical across a055, a075,
a085, and a095, suggesting a ceiling effect.

### Beta structure

F1=**-3.91e-03** (dominant), F2=1.77e-05, F3=-3.60e-06, F4=-2.95e-04, F5=2.32e-04, F6=3.70e-06, F7=3.14e-04

F1 dominates overwhelmingly. F4 and F5 have secondary betas.

---

## Key Questions for K=7

1. **Does K=7 recover the iCAF program?** Yes, surprisingly well. H-cor plateaus
   at ~0.87 for alpha >= 0.55, comparable to K=5's best (0.737 at alpha=0.55).
   K=7 actually recovers the iCAF H-score pattern better than K=5 at matched alpha.

2. **Does validation hold?** Mixed. Unadjusted Cox validates for alpha >= 0.35.
   Adjusted Cox validates at alpha=0.35 (p=0.033), alpha=0.75 (p=0.020), and
   alpha=0.95 (p=0.046). The KM validation is universally strong for alpha >= 0.55
   (all p < 2e-6). The pattern of "strong KM, borderline Cox" suggests the iCAF
   factor captures a real dichotomous signal but its continuous score overlaps
   with PurIST/DeCAF.

3. **How fragmented is the program?** The immune component consistently splits
   across 2-3 factors (e.g., at a055: F3 gets 18/60, F4 gets 23/60, F7 gets
   15/60). The best H-cor factor is not always the factor with the most iCAF gene
   content. This dissociation between H-score pattern and gene content is a K=7
   artifact.

4. **Overfitting evidence?** Severe. Training HRs of 0.056-0.108 (10-18x stronger
   than K=3's 0.344) do not replicate in validation (0.77-0.89). Training KM shows
   44.4 vs 11.1-11.6 months (33-month gap), while validation shows 29-31 vs 17
   months (12-14 month gap). The gap ratio (training/validation) is ~2.5x at K=7
   vs ~1.5x at K=3.

5. **What is the best K=7 fit?** K7_a075 has the strongest validation: adjusted
   Cox p=0.020, LRT p=0.020, KM p=8.0e-9 (14.2-month gap). But its iCAF factor
   (F2) has only 15% DECODER Immune overlap and 12.2% Basal contamination -- it
   is a different kind of factor than K=3's iCAF. K7_a035 is notable for having
   the most modest training HR (0.397) with validated adjusted Cox (p=0.033),
   suggesting less overfitting.

6. **Does n=273 support K=7?** The data suggest it can identify the dominant
   survival axis (iCAF) but overtrain the signal. The pattern of multiple zero
   betas (e.g., at a085: 4 of 7 factors zeroed out) shows the survival penalty
   is effectively reducing K=7 back toward a lower effective rank. The model has
   7 factors but concentrates all survival weight on 1-2.

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
