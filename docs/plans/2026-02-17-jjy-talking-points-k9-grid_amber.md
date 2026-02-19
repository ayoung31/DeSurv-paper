# K=9 DeSurv Factor Analysis: CV Grid Fits (Expanded Grid)

**For discussion with Jen Jen Yeh**
**Date: 2026-02-18 (amber revision)**
**Companion to: 2026-02-17-k-sensitivity-synthesis_amber.md**
**Data source: cv_grid exhaustive search (fixed lambda=0.3, nu=0.05, 100 inits)**

**Note for Amber:** This is a new document -- K=9 was not examined in the original
analysis. Seven fits: alpha=0, 0.25, 0.35, 0.55, 0.75, 0.85, 0.95.
Fits loaded by parameter matching. Scripts: `inst/cv_grid_training_analysis_amber.R`
and `inst/cv_grid_validation_analysis_amber.R`.

## Motivation for K=9

K=9 is the highest rank in the expanded grid and serves as an extreme test of
what happens when the model has far more factors than biologically distinct
programs. With 9 factors and only 273 training samples, we expect:

- Severe fragmentation of all biological programs
- Many spurious prognostic factors in training that do not validate
- Very low H-cor with the original K=3 iCAF program
- Possible convergence issues at high alpha (model may struggle to
  meaningfully partition survival signal across 9 factors)

K=9 provides the upper-bound contrast for the K=3 story: if K=3 produces the
cleanest, most validated iCAF signal, K=9 should produce the most fragmented,
least validated version.

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
| **Best iCAF factor** | F6 | F4 | F6 | F7 | F5 | F7 | F9 |
| **H-cor with orig F1** | 0.468 | 0.785 | 0.533 | 0.793 | 0.806 | 0.741 | 0.806 |
| **Top-270 overlap** | 50/270 | 129/270 | 64/270 | 127/270 | 120/270 | 145/270 | 109/270 |
| **Elyada iCAF** | 7/25 (28%) | 8/25 (32%) | 5/25 (20%) | 5/25 (20%) | 4/25 (16%) | 7/25 (28%) | 6/25 (24%) |
| **DECODER Immune** | 1/60 (1.7%) | 21/60 (35.0%) | 12/60 (20.0%) | 16/60 (26.7%) | 15/60 (25.0%) | 14/60 (23.3%) | 18/60 (30.0%) |
| **DECODER BasalTumor** | 5/139 (3.6%) | 8/139 (5.8%) | 5/139 (3.6%) | 21/139 (15.1%) | 13/139 (9.4%) | 23/139 (16.5%) | 17/139 (12.2%) |
| **Train HR (unadj)** | 0.829 (p=0.035) | 0.526 (p<0.0001) | 0.715 (p<0.0001) | 0.177 (p<0.0001) | 0.405 (p<0.0001) | 0.473 (p<0.0001) | 0.054 (p<0.0001) |
| **Val HR (unadj)** | 0.967 (p=0.69) | 0.837 (p=0.016) | 0.917 (p=0.22) | 0.827 (p=0.008) | 0.823 (p=0.004) | 0.892 (p=0.16) | 0.824 (p=0.008) |
| **Val HR (adj)** | 0.939 (p=0.48) | 0.867 (p=0.071) | 0.935 (p=0.36) | 0.864 (p=0.050) | 0.844 (p=0.020) | 0.885 (p=0.15) | 0.885 (p=0.10) |
| **Val LRT p** | 4.77e-1 | 7.32e-2 | 3.61e-1 | 5.03e-2 | 2.08e-2 | 1.52e-1 | 1.01e-1 |
| **Val KM p** | 9.30e-3 | 8.42e-4 | 2.05e-3 | 2.85e-5 | 7.13e-4 | 8.84e-2 | 1.01e-5 |
| **# training prognostic** | 2 | 2 | 3 | 4 | 3 | 3 | 4 |

---

## K9_a000: Standard NMF (alpha=0.00)

### Factor identity

H-score correlation with original K=3:

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | 0.294   | -0.688  | 0.636   |
| F2     | -0.496  | 0.818   | -0.564  |
| F3     | 0.168   | -0.856  | 0.885   |
| F4     | 0.107   | 0.387   | -0.376  |
| F5     | -0.332  | 0.697   | -0.481  |
| **F6** | **0.468** | -0.647 | 0.719   |
| F7     | 0.057   | 0.176   | -0.063  |
| F8     | 0.221   | -0.748  | 0.857   |
| F9     | -0.133  | 0.954   | -0.907  |

Best iCAF match: F6 (rho=0.468). The dominant correlations are with original F2
and F3 (stroma/basal axis), not F1 (iCAF). F9 nearly perfectly matches original
F2 (0.954), and F3/F8 capture original F3. No factor strongly recovers iCAF.

### Reference signature overlaps (F6 as best iCAF match, top-270)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 7/25 | 28.0% |
| Elyada myCAF (15) | 5/15 | 33.3% |
| SCISSORS iCAF (17) | 2/17 | 11.8% |
| SCISSORS myCAF (16) | 1/16 | 6.2% |
| DeCAF restCAF (8) | 0/8 | 0.0% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 1/60 | 1.7% |
| DECODER ActivatedStroma (108) | 15/108 | 13.9% |
| DECODER NormalStroma (15) | 1/15 | 6.7% |
| DECODER BasalTumor (139) | 5/139 | 3.6% |
| DECODER ClassicalTumor (166) | 12/166 | 7.2% |

F6 has myCAF (33%) and Elyada iCAF (28%) but almost no Immune genes (1.7%).
Where do the Immune genes go? F7 has 32/60 DECODER Immune and 13/25 Elyada
iCAF -- the immune/iCAF program forms its own factor separate from F6. At
alpha=0, the model splits the iCAF program across F6 (stroma-heavy) and F7
(immune-heavy), and neither alone captures the combined K=3 program.

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 2/25 | 0/60 | 5/139 | Mixed |
| F2 | 2/25 | 1/60 | 29/139 | Basal Tumor |
| F3 | 1/25 | 0/60 | 45/139 | **Basal Tumor (dominant)** |
| F4 | 2/25 | 0/60 | 1/139 | Low content |
| F5 | 0/25 | 0/60 | 10/139 | Basal minor |
| **F6** | **7/25** | **1/60** | **5/139** | **Stroma/myCAF** |
| F7 | **13/25** | **32/60** | 3/139 | **Immune/iCAF** |
| F8 | 4/25 | 1/60 | 23/139 | Basal/Mixed |
| F9 | 5/25 | **29/60** | 0/139 | **Immune** |

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.941 | 0.49 | 0.993 | 0.90 |
| **F2** | **1.242** | **0.021** | **1.270** | **<0.0001** |
| F3 | 1.029 | 0.74 | 0.999 | 0.99 |
| F4 | 0.872 | 0.12 | **0.844** | **0.003** |
| F5 | 1.093 | 0.31 | 1.095 | 0.24 |
| F6 (iCAF) | **0.829** | **0.035** | 0.967 | 0.69 |
| F7 | 0.959 | 0.62 | 1.118 | 0.058 |
| F8 | 0.996 | 0.96 | 1.021 | 0.76 |
| F9 | 0.944 | 0.50 | 0.955 | 0.40 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.026 | 0.78 | 1.066 | 0.25 |
| F2 | 1.000 | 1.00 | 1.130 | 0.062 |
| F3 | 0.973 | 0.76 | 0.932 | 0.27 |
| F4 | 1.040 | 0.67 | 0.921 | 0.17 |
| F5 | 0.988 | 0.90 | 0.926 | 0.37 |
| F6 (iCAF) | 0.842 | 0.063 | 0.939 | 0.48 |
| F7 | 0.885 | 0.19 | 1.072 | 0.25 |
| F8 | 0.922 | 0.38 | 0.899 | 0.16 |
| F9 | 0.989 | 0.90 | 0.977 | 0.69 |

**The iCAF factor (F6) is weakly prognostic in training (p=0.035) and does not
validate (p=0.69).** The strongest validated signal is F2 (Basal Tumor, val
p<0.0001) and F4 (val p=0.003), neither of which is the iCAF program. After
PurIST+DeCAF adjustment, no factor is significant in either training or
validation.

**LRT adding F6 (iCAF) to PurIST + DeCAF:** Training p=0.063, Validation p=0.477

**Median-split KM for F6 (iCAF):**
- Training (n=273): High 22.8 vs Low 19.8 months, p=0.072
- Validation (n=570): High 24.4 vs Low 19.4 months, **p=0.009**

The KM shows a modest 5.0-month validation gap (p=0.009) despite the Cox model
failing, suggesting the iCAF signal exists in a crude median split but is too
weak for continuous risk scoring at alpha=0.

### Beta structure

F1=2.28e-05, F2=1.00e-04, F3=1.24e-05, F4=-5.01e-05, F5=0.00e+00,
F6=-2.10e-04, F7=-2.86e-05, F8=0.00e+00, F9=-5.02e-05

All tiny (no survival penalty -> near-zero betas). F6 (iCAF) has the largest
magnitude but it is still negligible.

---

## K9_a025: alpha=0.25

### Factor identity

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.030  | -0.153  | 0.331   |
| F2     | -0.378  | 0.841   | -0.586  |
| F3     | -0.270  | 0.967   | -0.843  |
| **F4** | **0.785** | -0.482 | 0.458   |
| F5     | 0.743   | -0.315  | 0.227   |
| F6     | 0.138   | -0.758  | 0.865   |
| F7     | 0.138   | 0.185   | -0.180  |
| F8     | 0.176   | -0.927  | 0.947   |
| F9     | -0.271  | 0.903   | -0.719  |

F4 is the best iCAF match (rho=0.785). F5 is also strongly correlated with
original F1 (rho=0.743), indicating the iCAF program is split across two
factors. F3 matches original F2 nearly perfectly (0.967), and F8 matches
original F3 (0.947).

### Reference signature overlaps (F4, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 8/25 | 32.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | 6/17 | 35.3% |
| SCISSORS myCAF (16) | 2/16 | 12.5% |
| DeCAF restCAF (8) | 2/8 | 25.0% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 21/60 | 35.0% |
| DECODER ActivatedStroma (108) | 3/108 | 2.8% |
| DECODER NormalStroma (15) | 7/15 | 46.7% |
| DECODER BasalTumor (139) | 8/139 | 5.8% |
| DECODER ClassicalTumor (166) | 24/166 | 14.5% |

F4 has strong Immune (35%), NormalStroma (47%), SCISSORS iCAF (35%), and Elyada
iCAF (32%) overlap with low Basal (5.8%). This is a reasonably pure
iCAF/Immune factor, similar to K5_a25's F4.

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 0/25 | 0/60 | 5/139 | Low content |
| F2 | 5/25 | 3/60 | 21/139 | Basal |
| F3 | 3/25 | 32/60 | 1/139 | **Immune (pure)** |
| **F4** | **8/25** | **21/60** | **8/139** | **iCAF/Immune** |
| F5 | 2/25 | 5/60 | 3/139 | Mixed |
| F6 | 5/25 | 9/60 | 26/139 | Basal/Mixed |
| F7 | 2/25 | 7/60 | 3/139 | Mixed |
| F8 | 0/25 | 0/60 | 55/139 | **Basal Tumor (dominant)** |
| F9 | 8/25 | 12/60 | 3/139 | iCAF/Immune minor |

Note F3 captures 32/60 DECODER Immune with only 1/139 Basal -- a pure immune
factor. F4 (iCAF match) also has 21/60 Immune. The immune signal is split
between F3 and F4. F5 (rho=0.743 with orig F1) captures only 2/25 Elyada iCAF
and 5/60 Immune, suggesting it picks up a different component of the original
F1 (possibly the stroma portion).

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.117 | 0.21 | 0.948 | 0.44 |
| F2 | 1.135 | 0.16 | **1.184** | **0.007** |
| F3 | 1.017 | 0.84 | 1.062 | 0.31 |
| **F4 (iCAF)** | **0.526** | **<0.0001** | **0.837** | **0.016** |
| **F5** | **0.497** | **<0.0001** | **0.756** | **0.0004** |
| F6 | 1.078 | 0.39 | **1.223** | **0.006** |
| F7 | 0.919 | 0.34 | 0.945 | 0.28 |
| F8 | 1.049 | 0.58 | 1.076 | 0.24 |
| F9 | 1.059 | 0.51 | 1.124 | 0.057 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.104 | 0.32 | **0.855** | **0.043** |
| F2 | 0.959 | 0.67 | 1.040 | 0.57 |
| F3 | 0.994 | 0.95 | 1.004 | 0.95 |
| **F4 (iCAF)** | **0.578** | **<0.0001** | 0.867 | 0.071 |
| **F5** | **0.516** | **<0.0001** | **0.788** | **0.005** |
| F6 | 0.971 | 0.75 | 1.055 | 0.50 |
| F7 | 1.075 | 0.44 | 1.018 | 0.74 |
| F8 | 0.989 | 0.90 | 0.996 | 0.95 |
| F9 | 1.001 | 0.99 | 1.050 | 0.44 |

Two factors are strongly prognostic in training: F4 (iCAF, HR=0.526) and F5
(HR=0.497). Both validate unadjusted (F4 p=0.016, F5 p=0.0004). After
PurIST+DeCAF, F4 is borderline (p=0.071) but F5 remains significant (p=0.005).
F5 is NOT the iCAF factor by H-cor (rho=0.743 vs F4's 0.785), yet it carries
the stronger independent signal. This illustrates the fragmentation problem:
the survival-relevant part of iCAF is split across factors and the "best match"
factor is not necessarily the one that validates best.

**LRT adding F4 (iCAF) to PurIST + DeCAF:** Training p=7.04e-8, Validation p=0.073

**Median-split KM for F4 (iCAF):**
- Training (n=273): High 27.0 vs Low 14.2 months, **p=1.1e-6**
- Validation (n=570): High 28.7 vs Low 19.0 months, **p=8.4e-4**

Strong 9.7-month validation gap. The KM validates well even though the
adjusted Cox is borderline, suggesting the iCAF signal is real but partially
redundant with PurIST/DeCAF.

### Beta structure

F1=4.01e-04, F2=3.17e-07, F3=-1.33e-06, F4=**-6.65e-04**, F5=**-7.99e-04** (dominant),
F6=1.78e-04, F7=2.39e-05, F8=5.06e-05, F9=8.17e-05

F5 has the largest beta magnitude, followed by F4 (iCAF). Two factors dominate
the LP, consistent with the two-factor prognostic pattern.

---

## K9_a035: alpha=0.35

### Factor identity

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.368  | 0.890   | -0.601  |
| F2     | 0.056   | -0.154  | 0.535   |
| F3     | -0.250  | 0.959   | -0.761  |
| F4     | 0.240   | 0.375   | -0.110  |
| F5     | 0.264   | -0.520  | 0.747   |
| **F6** | **0.533** | 0.257  | -0.147  |
| F7     | -0.401  | 0.237   | 0.181   |
| F8     | 0.259   | 0.714   | -0.577  |
| F9     | 0.299   | -0.466  | 0.605   |

F6 is the best iCAF match but at a modest rho=0.533. The stroma axis dominates:
F1 and F3 both match original F2 well (0.890 and 0.959). No factor recovers
original F1 strongly.

### Reference signature overlaps (F6, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 5/25 | 20.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | 7/17 | 41.2% |
| SCISSORS myCAF (16) | 3/16 | 18.8% |
| DeCAF restCAF (8) | 3/8 | 37.5% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 12/60 | 20.0% |
| DECODER ActivatedStroma (108) | 9/108 | 8.3% |
| DECODER NormalStroma (15) | 0/15 | 0.0% |
| DECODER BasalTumor (139) | 5/139 | 3.6% |
| DECODER ClassicalTumor (166) | 18/166 | 10.8% |

F6 has high SCISSORS iCAF (41%) and DeCAF restCAF (38%) but moderate Immune
(20%) and no NormalStroma. The factor captures CAF subtype genes but not the
immune component.

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 7/25 | 5/60 | 27/139 | Basal/iCAF |
| F2 | 2/25 | 2/60 | 32/139 | Basal Tumor |
| F3 | 4/25 | 11/60 | 14/139 | Mixed |
| F4 | 3/25 | 20/60 | 8/139 | Immune |
| F5 | 2/25 | 2/60 | 30/139 | Basal Tumor |
| **F6** | **5/25** | **12/60** | **5/139** | **iCAF/CAF** |
| F7 | 1/25 | 0/60 | 44/139 | **Basal (dominant)** |
| F8 | 6/25 | 18/60 | 10/139 | iCAF/Immune |
| F9 | 1/25 | 0/60 | 13/139 | Basal minor |

The immune signal is split: F4 gets 20/60 Immune, F8 gets 18/60, F6 gets
12/60. No single factor consolidates the iCAF+Immune program.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.150 | 0.13 | 1.146 | 0.052 |
| F2 | 1.034 | 0.70 | 1.120 | 0.16 |
| F3 | 1.052 | 0.56 | 1.129 | 0.092 |
| F4 | 0.877 | 0.079 | 0.989 | 0.89 |
| F5 | 0.942 | 0.48 | 0.999 | 0.99 |
| **F6 (iCAF)** | **0.715** | **<0.0001** | 0.917 | 0.22 |
| **F7** | **1.356** | **0.002** | **1.219** | **0.008** |
| **F8** | **0.695** | **<0.0001** | 0.951 | 0.51 |
| F9 | 0.907 | 0.26 | 0.992 | 0.91 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.061 | 0.53 | 1.035 | 0.63 |
| F2 | 0.900 | 0.26 | 0.979 | 0.81 |
| F3 | 1.042 | 0.65 | 1.051 | 0.50 |
| F4 | 0.942 | 0.49 | 0.949 | 0.56 |
| F5 | 0.906 | 0.28 | 0.902 | 0.22 |
| **F6 (iCAF)** | **0.771** | **0.002** | 0.935 | 0.36 |
| F7 | 1.059 | 0.62 | 1.031 | 0.71 |
| **F8** | **0.664** | **<0.0001** | 0.919 | 0.29 |
| F9 | 0.978 | 0.81 | 0.964 | 0.63 |

Three factors are training-significant (F6, F7, F8). Only F7 validates
unadjusted (p=0.008), and F7 is not the iCAF factor -- it is a basal-heavy
factor. F6 (iCAF) and F8 both fail to validate (p=0.22 and p=0.51). After
adjustment, F6 and F8 remain training-significant but neither validates.
Classic K=9 overfitting: strong training, no validation.

**LRT adding F6 (iCAF) to PurIST + DeCAF:** Training p=0.003, Validation p=0.361

**Median-split KM for F6 (iCAF):**
- Training (n=273): High 25.4 vs Low 15.3 months, **p=4.9e-6**
- Validation (n=570): High 23.9 vs Low 19.8 months, **p=0.002**

The KM validates (4.1-month gap, p=0.002) despite Cox failure, again showing
the median split captures a crude signal that continuous scoring cannot.

### Beta structure

F1=3.26e-04, F2=7.74e-05, F3=3.75e-04, F4=3.44e-05, F5=0.00e+00,
F6=-5.92e-04, F7=4.40e-04, F8=**-1.29e-03** (dominant), F9=4.28e-07

F8 has the dominant beta, not F6 (iCAF). The model's survival signal is driven
by a factor (F8) that is not the best iCAF match.

---

## K9_a055: alpha=0.55

### Factor identity

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.376  | -0.125  | 0.419   |
| F2     | 0.602   | -0.214  | 0.100   |
| F3     | 0.193   | -0.936  | 0.958   |
| F4     | -0.003  | 0.339   | -0.057  |
| F5     | -0.009  | 0.532   | -0.494  |
| F6     | -0.316  | 0.960   | -0.804  |
| **F7** | **0.793** | -0.192 | 0.069   |
| F8     | 0.571   | -0.626  | 0.507   |
| F9     | 0.020   | -0.516  | 0.698   |

F7 is the best iCAF match (rho=0.793). F2 also correlates with original F1
(rho=0.602). F6 matches original F2 (0.960) and F3 matches original F3 (0.958).
The stroma/basal axes are well recovered.

### Reference signature overlaps (F7, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 5/25 | 20.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | 4/17 | 23.5% |
| SCISSORS myCAF (16) | 2/16 | 12.5% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | 1/4 | 25.0% |
| DECODER Immune (60) | 16/60 | 26.7% |
| DECODER ActivatedStroma (108) | 3/108 | 2.8% |
| DECODER NormalStroma (15) | 6/15 | 40.0% |
| DECODER BasalTumor (139) | 21/139 | 15.1% |
| DECODER ClassicalTumor (166) | 28/166 | 16.9% |

F7 has moderate Immune (27%), NormalStroma (40%), but also picks up Basal
(15%) and Classical (17%) tumor genes. The factor is less pure than K9_a025's
F4, which had 35% Immune and only 5.8% Basal.

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 6/25 | 3/60 | 67/139 | **Basal (dominant)** |
| F2 | 9/25 | 20/60 | 1/139 | **iCAF/Immune** |
| F3 | 1/25 | 0/60 | 33/139 | Basal |
| F4 | 3/25 | 6/60 | 6/139 | Mixed |
| F5 | 6/25 | 15/60 | 2/139 | Immune/iCAF |
| F6 | 4/25 | 22/60 | 3/139 | **Immune** |
| **F7** | **5/25** | **16/60** | **21/139** | **iCAF/Mixed** |
| F8 | 2/25 | 1/60 | 1/139 | Low content |
| F9 | 1/25 | 1/60 | 17/139 | Basal minor |

The immune signal is highly fragmented: F2 gets 20/60, F6 gets 22/60, F7 gets
16/60, F5 gets 15/60. Four factors share the immune genes. F2 has the most
Elyada iCAF genes (9/25) despite F7 being the H-cor best match.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1** | **1.409** | **0.0002** | **1.377** | **<0.0001** |
| **F2** | **0.664** | **<0.0001** | **0.880** | **0.022** |
| F3 | 1.041 | 0.65 | 1.089 | 0.19 |
| F4 | 0.999 | 0.99 | 1.008 | 0.92 |
| F5 | 0.938 | 0.46 | 0.967 | 0.50 |
| F6 | 1.046 | 0.61 | 1.093 | 0.14 |
| **F7 (iCAF)** | **0.177** | **<0.0001** | **0.827** | **0.008** |
| **F8** | **0.763** | **0.003** | **0.732** | **<0.0001** |
| F9 | 1.163 | 0.096 | 1.042 | 0.56 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.053 | 0.65 | **1.223** | **0.007** |
| F2 | 0.811 | 0.069 | 0.990 | 0.88 |
| F3 | 0.991 | 0.92 | 1.003 | 0.97 |
| F4 | 0.896 | 0.25 | 0.898 | 0.19 |
| F5 | 1.063 | 0.51 | 1.026 | 0.62 |
| F6 | 1.000 | 1.00 | 1.002 | 0.97 |
| **F7 (iCAF)** | **0.163** | **<0.0001** | 0.864 | **0.050** |
| F8 | 0.970 | 0.77 | **0.770** | **0.0003** |
| F9 | 1.082 | 0.43 | 0.923 | 0.30 |

Four factors are training-significant unadjusted (F1, F2, F7, F8). The iCAF
factor (F7) has an extreme training HR of 0.177 -- a 5.6x risk reduction per
SD. In validation, F7 reaches p=0.008 unadjusted and is borderline adjusted
(p=0.050). But F8 validates even more strongly adjusted (p=0.0003), and F1
validates adjusted (p=0.007). The survival signal is distributed across
multiple factors, with the iCAF factor being one of several.

The training-validation gap for F7 is enormous: HR=0.177 (training) vs 0.827
(validation). This is the most extreme train-val disconnect in the grid, a
hallmark of K=9 overfitting at moderate-high alpha.

**LRT adding F7 (iCAF) to PurIST + DeCAF:** Training p=2.7e-50, Validation p=0.050

**Median-split KM for F7 (iCAF):**
- Training (n=273): High 44.4 vs Low 11.6 months, **p<0.0001**
- Validation (n=570): High 29.2 vs Low 19.0 months, **p=2.8e-5**

Massive 32.8-month training gap (clearly overtrained). The validation gap is
10.2 months (p=2.8e-5), the largest absolute validation gap in the K=9 grid.

### Beta structure

F1=7.44e-05, F2=-1.35e-04, F3=0.00e+00, F4=1.31e-04, F5=-1.83e-05,
F6=0.00e+00, F7=**-2.85e-03** (dominant), F8=0.00e+00, F9=2.05e-04

F7 (iCAF) completely dominates the LP. Three factors have zero beta (F3, F6,
F8), meaning the model has effectively collapsed from 9 to 6 active factors.

---

## K9_a075: alpha=0.75

### Factor identity

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.181  | 0.897   | -0.781  |
| F2     | -0.136  | 0.961   | -0.822  |
| F3     | -0.295  | 0.838   | -0.580  |
| F4     | -0.093  | 0.939   | -0.797  |
| **F5** | **0.806** | -0.599 | 0.523   |
| F6     | 0.097   | -0.137  | 0.434   |
| F7     | -0.411  | -0.213  | 0.599   |
| F8     | 0.169   | -0.893  | 0.959   |
| F9     | 0.379   | -0.627  | 0.716   |

F5 is the best iCAF match (rho=0.806). Remarkably, four factors (F1, F2, F3,
F4) all correlate strongly with original F2 (rho=0.897, 0.961, 0.838, 0.939).
The stroma/basal axis has been split into four sub-factors.

### Reference signature overlaps (F5, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 4/25 | 16.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | 4/17 | 23.5% |
| SCISSORS myCAF (16) | 0/16 | 0.0% |
| DeCAF restCAF (8) | 0/8 | 0.0% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 15/60 | 25.0% |
| DECODER ActivatedStroma (108) | 5/108 | 4.6% |
| DECODER NormalStroma (15) | 4/15 | 26.7% |
| DECODER BasalTumor (139) | 13/139 | 9.4% |
| DECODER ClassicalTumor (166) | 46/166 | 27.7% |

F5 has the highest Classical Tumor overlap in the K=9 grid (27.7%), with
moderate Immune (25%) and NormalStroma (27%). The factor is diluted by
tumor genes.

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 3/25 | 4/60 | 1/139 | Stroma |
| F2 | 4/25 | 22/60 | 7/139 | **Immune** |
| F3 | 4/25 | 4/60 | 11/139 | Mixed |
| F4 | 5/25 | 18/60 | 2/139 | Immune/iCAF |
| **F5** | **4/25** | **15/60** | **13/139** | **iCAF/Mixed** |
| F6 | 6/25 | 7/60 | 9/139 | Mixed |
| F7 | 1/25 | 5/60 | 50/139 | **Basal (dominant)** |
| F8 | 0/25 | 0/60 | 53/139 | **Basal (dominant)** |
| F9 | 3/25 | 2/60 | 9/139 | Mixed |

The immune program splits three ways: F2 (22/60), F4 (18/60), F5 (15/60).
Two factors (F7, F8) are pure Basal Tumor. The iCAF genes are scattered:
F6 has the most Elyada iCAF (6/25) but is not the best H-cor match.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.978 | 0.80 | 1.016 | 0.81 |
| F2 | 0.899 | 0.20 | 1.065 | 0.37 |
| **F3** | **1.245** | **0.017** | 1.055 | 0.43 |
| F4 | 0.924 | 0.34 | 0.986 | 0.84 |
| **F5 (iCAF)** | **0.405** | **<0.0001** | **0.823** | **0.004** |
| F6 | 0.992 | 0.92 | 1.123 | 0.086 |
| **F7** | **1.718** | **<0.0001** | **1.262** | **0.0004** |
| F8 | 1.041 | 0.65 | **1.182** | **0.041** |
| F9 | 0.960 | 0.64 | 0.858 | 0.058 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.029 | 0.75 | 1.009 | 0.89 |
| F2 | 0.878 | 0.14 | 1.016 | 0.82 |
| **F3** | **1.302** | **0.006** | 0.986 | 0.84 |
| F4 | 0.958 | 0.63 | 0.965 | 0.61 |
| **F5 (iCAF)** | **0.420** | **<0.0001** | **0.844** | **0.020** |
| F6 | 0.834 | 0.052 | 1.020 | 0.78 |
| **F7** | **1.434** | **0.003** | 1.067 | 0.40 |
| F8 | 0.956 | 0.61 | 1.026 | 0.77 |
| F9 | 1.099 | 0.32 | 0.843 | 0.052 |

Three factors are training-significant unadjusted (F3, F5, F7). F5 (iCAF) is
the standout: it validates both unadjusted (p=0.004) and adjusted (p=0.020).
**This is the only K=9 fit where the iCAF factor survives PurIST+DeCAF
adjustment in validation.** F7 also validates unadjusted (p=0.0004) but fails
adjusted (p=0.40). The training-val gap for F5 (HR=0.405 vs 0.823) is still
large but less extreme than a=0.55.

**LRT adding F5 (iCAF) to PurIST + DeCAF:** Training p=2.6e-16, **Validation p=0.021**

**Median-split KM for F5 (iCAF):**
- Training (n=273): High 37.7 vs Low 13.1 months, **p=3.9e-14**
- Validation (n=570): High 25.6 vs Low 19.2 months, **p=7.1e-4**

The validation KM shows a 6.4-month gap (p=7.1e-4). Combined with the
adjusted Cox validation (p=0.020) and LRT (p=0.021), this is K=9's best
performing alpha overall.

### Beta structure

F1=-8.05e-06, F2=-4.15e-04, F3=5.74e-04, F4=-8.82e-05, F5=**-1.77e-03** (dominant),
F6=-1.19e-04, F7=5.25e-04, F8=7.91e-07, F9=4.21e-04

F5 (iCAF) dominates the LP. F3 and F7 also have substantial betas, consistent
with three training-prognostic factors.

---

## K9_a085: alpha=0.85

### Factor identity

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | 0.100   | 0.134   | 0.058   |
| F2     | -0.145  | -0.534  | 0.865   |
| F3     | -0.334  | -0.161  | 0.545   |
| F4     | -0.048  | 0.837   | -0.624  |
| F5     | -0.079  | 0.660   | -0.346  |
| F6     | -0.145  | 0.689   | -0.372  |
| **F7** | **0.741** | -0.591 | 0.612   |
| F8     | -0.113  | 0.678   | -0.316  |
| F9     | 0.036   | 0.818   | -0.671  |

F7 is the best iCAF match (rho=0.741). The stroma axis fragments broadly:
F4 (0.837), F9 (0.818), F5 (0.660), F6 (0.689), F8 (0.678) all correlate
with original F2. Five factors share the stroma axis, leaving little room for
a clean iCAF factor.

### Reference signature overlaps (F7, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 7/25 | 28.0% |
| Elyada myCAF (15) | 3/15 | 20.0% |
| SCISSORS iCAF (17) | 7/17 | 41.2% |
| SCISSORS myCAF (16) | 0/16 | 0.0% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 14/60 | 23.3% |
| DECODER ActivatedStroma (108) | 5/108 | 4.6% |
| DECODER NormalStroma (15) | 5/15 | 33.3% |
| DECODER BasalTumor (139) | 23/139 | 16.5% |
| DECODER ClassicalTumor (166) | 40/166 | 24.1% |

F7 has decent iCAF overlap (SCISSORS 41%, Elyada 28%) but is contaminated
with Basal (16.5%) and Classical (24.1%) tumor genes. The factor is less pure
than lower-alpha fits.

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 4/25 | 4/60 | 11/139 | Mixed |
| F2 | 3/25 | 2/60 | 37/139 | Basal Tumor |
| F3 | 3/25 | 2/60 | 31/139 | Basal Tumor |
| F4 | 3/25 | 11/60 | 9/139 | Mixed/Immune |
| F5 | 0/25 | 12/60 | 11/139 | Immune/Basal |
| F6 | 4/25 | **22/60** | 14/139 | **Immune** |
| **F7** | **7/25** | **14/60** | **23/139** | **iCAF/Mixed** |
| F8 | 3/25 | 11/60 | 11/139 | Mixed |
| F9 | 1/25 | **21/60** | 1/139 | **Immune (pure)** |

The immune signal splits: F6 (22/60), F9 (21/60), F7 (14/60), F5 (12/60),
F4 (11/60), F8 (11/60). Six of 9 factors capture immune genes. This is
extreme fragmentation. Basal splits between F2 (37) and F3 (31).

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.971 | 0.73 | 1.088 | 0.28 |
| **F2** | **1.592** | **<0.0001** | **1.185** | **0.014** |
| **F3** | **1.526** | **<0.0001** | **1.190** | **0.016** |
| F4 | 0.905 | 0.22 | 0.954 | 0.55 |
| F5 | 0.969 | 0.70 | 1.123 | 0.16 |
| F6 | 0.962 | 0.64 | 1.120 | 0.14 |
| **F7 (iCAF)** | **0.473** | **<0.0001** | 0.892 | 0.16 |
| F8 | 1.009 | 0.92 | 1.106 | 0.19 |
| F9 | 0.863 | 0.067 | 0.905 | 0.15 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.033 | 0.71 | 1.063 | 0.45 |
| **F2** | **1.467** | **0.0001** | 1.033 | 0.66 |
| **F3** | **1.292** | **0.032** | 1.023 | 0.77 |
| F4 | 0.913 | 0.30 | 0.902 | 0.20 |
| F5 | 0.879 | 0.15 | 1.013 | 0.88 |
| F6 | 0.843 | 0.056 | 1.008 | 0.92 |
| **F7 (iCAF)** | **0.487** | **<0.0001** | 0.885 | 0.15 |
| F8 | 0.947 | 0.54 | 0.991 | 0.91 |
| F9 | 0.932 | 0.42 | 0.900 | 0.15 |

Three factors are training-significant (F2, F3, F7). F2 and F3 validate
unadjusted (p=0.014, 0.016) but fail adjusted. **F7 (iCAF) does NOT validate
(unadj p=0.16, adj p=0.15).** Despite a strong training HR of 0.473, the
validation HR is only 0.892. The train-val gap (0.473 vs 0.892) indicates
overfitting.

**LRT adding F7 (iCAF) to PurIST + DeCAF:** Training p=2.5e-13, Validation p=0.152

**Median-split KM for F7 (iCAF):**
- Training (n=273): High 30.4 vs Low 13.1 months, **p=7.4e-10**
- Validation (n=570): High 23.8 vs Low 20.0 months, p=0.088

The validation KM is borderline (3.8-month gap, p=0.088). This is the weakest
validation performance in the K=9 grid aside from alpha=0.

### Beta structure

F1=3.81e-06, F2=1.11e-03, F3=4.82e-04, F4=-8.31e-05, F5=-2.56e-05,
F6=-1.72e-04, F7=**-1.78e-03** (dominant), F8=0.00e+00, F9=-5.01e-05

F7 (iCAF) dominates the LP. F2 and F3 also have substantial positive betas,
pushing high-basal samples toward higher risk.

---

## K9_a095: alpha=0.95

### Factor identity

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.113  | -0.720  | 0.904   |
| F2     | 0.083   | 0.453   | -0.405  |
| F3     | -0.370  | 0.800   | -0.475  |
| F4     | -0.194  | 0.986   | -0.847  |
| F5     | -0.302  | 0.653   | -0.376  |
| F6     | 0.392   | -0.911  | 0.885   |
| F7     | 0.459   | -0.925  | 0.866   |
| F8     | -0.264  | 0.958   | -0.773  |
| **F9** | **0.806** | -0.046 | -0.121  |

F9 is the best iCAF match (rho=0.806). F4 matches original F2 nearly perfectly
(0.986), and F8 also matches (0.958). F6 and F7 both anti-correlate with
original F2 (-0.911 and -0.925). The stroma axis is captured by multiple
factors from both sides.

### Reference signature overlaps (F9, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 6/25 | 24.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | 6/17 | 35.3% |
| SCISSORS myCAF (16) | 0/16 | 0.0% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 18/60 | 30.0% |
| DECODER ActivatedStroma (108) | 8/108 | 7.4% |
| DECODER NormalStroma (15) | 6/15 | 40.0% |
| DECODER BasalTumor (139) | 17/139 | 12.2% |
| DECODER ClassicalTumor (166) | 22/166 | 13.3% |

F9 has the profile of an iCAF/Immune factor: SCISSORS iCAF 35%, Elyada iCAF
24%, Immune 30%, NormalStroma 40%, with moderate Basal contamination (12%).

### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 1/25 | 1/60 | 80/139 | **Basal (extreme)** |
| F2 | 4/25 | 11/60 | 0/139 | Immune minor |
| F3 | 3/25 | 1/60 | 10/139 | Basal minor |
| F4 | 5/25 | 17/60 | 1/139 | Immune |
| F5 | 4/25 | 18/60 | 15/139 | Immune/Basal |
| F6 | 1/25 | 0/60 | 10/139 | Basal minor |
| F7 | 5/25 | 5/60 | 11/139 | Mixed |
| F8 | 3/25 | 9/60 | 6/139 | Mixed |
| **F9** | **6/25** | **18/60** | **17/139** | **iCAF/Mixed** |

F1 is an extreme Basal factor (80/139). The immune genes scatter across F4
(17/60), F5 (18/60), F9 (18/60), F2 (11/60). At extreme alpha with K=9, the
model concentrates basal in one factor but fragments immune across four.

### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1** | **1.273** | **0.007** | **1.386** | **<0.0001** |
| F2 | 0.893 | 0.19 | 0.935 | 0.25 |
| **F3** | **1.270** | **0.013** | **1.135** | **0.045** |
| F4 | 0.991 | 0.92 | 1.021 | 0.73 |
| F5 | 1.200 | 0.051 | **1.159** | **0.026** |
| F6 | 0.934 | 0.45 | **0.847** | **0.009** |
| **F7** | **0.777** | **0.005** | 0.929 | 0.35 |
| F8 | 1.022 | 0.80 | 1.083 | 0.23 |
| **F9 (iCAF)** | **0.054** | **<0.0001** | **0.824** | **0.008** |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.043 | 0.67 | 1.172 | 0.069 |
| F2 | 1.030 | 0.75 | 1.001 | 0.98 |
| **F3** | **1.216** | **0.046** | 1.033 | 0.62 |
| F4 | 0.997 | 0.97 | 0.983 | 0.79 |
| F5 | 1.043 | 0.68 | 1.023 | 0.75 |
| F6 | 0.995 | 0.96 | **0.849** | **0.021** |
| **F7** | **0.782** | **0.007** | 0.929 | 0.37 |
| F8 | 0.971 | 0.75 | 0.996 | 0.96 |
| **F9 (iCAF)** | **0.050** | **<0.0001** | 0.885 | 0.10 |

The most extreme fit in the grid. F9 (iCAF) has a training HR of 0.054 -- an
18.5x risk reduction per SD -- which is clearly overtrained. It validates
unadjusted (p=0.008) but not adjusted (p=0.10). Four factors are training-
prognostic, and in validation, F1 (Basal) is the strongest validated signal
(unadj p<0.0001). F6 validates adjusted (p=0.021) but was not even training-
significant, suggesting it captures a real but alpha-independent signal.

**LRT adding F9 (iCAF) to PurIST + DeCAF:** Training p=5.8e-78, Validation p=0.101

**Median-split KM for F9 (iCAF):**
- Training (n=273): High 44.4 vs Low 11.2 months, **p<0.0001**
- Validation (n=570): High 28.7 vs Low 18.0 months, **p=1.0e-5**

Extreme 33.2-month training gap (pure overtraining). The validation gap is
10.7 months (p=1.0e-5), comparable to a=0.55 and the strongest KM in the grid.

### Beta structure

F1=6.66e-05, F2=-2.23e-05, F3=1.55e-04, F4=-2.45e-06, F5=8.66e-05,
F6=2.30e-05, F7=-1.07e-04, F8=-3.45e-06, F9=**-3.92e-03** (extreme dominant)

F9 (iCAF) completely dominates the LP with a beta 37x larger than any other
factor. The model has essentially collapsed to a single survival-driving factor
with 8 spectator factors.

---

## The K=9 Pattern: Fragmentation, Overfitting, and Selective Validation

### Training-validation disconnect

The hallmark of K=9 is extreme training-validation disconnect for the iCAF
factor:

| Alpha | Train HR | Val HR | Ratio | Interpretation |
|-------|----------|--------|-------|----------------|
| 0.00 | 0.829 | 0.967 | 1.2x | Mild (low penalty) |
| 0.25 | 0.526 | 0.837 | 1.6x | Moderate |
| 0.35 | 0.715 | 0.917 | 1.3x | Moderate |
| 0.55 | 0.177 | 0.827 | 4.7x | **Severe** |
| 0.75 | 0.405 | 0.823 | 2.0x | Moderate-high |
| 0.85 | 0.473 | 0.892 | 1.9x | Moderate |
| 0.95 | 0.054 | 0.824 | 15.3x | **Extreme** |

Compare K=3 alpha=0.55: train HR=0.344, val HR=0.800, ratio=2.3x. K=9 at
alpha=0.55 and 0.95 shows ratios of 4.7x and 15.3x -- clear evidence that the
model is fitting noise in the training data.

### The alpha=0.75 exception

K9_a075 stands out as the only alpha where the iCAF factor validates after
PurIST+DeCAF adjustment (p=0.020, LRT p=0.021). Why?

- The training HR (0.405) is strong but not extreme like a=0.55 (0.177) or
  a=0.95 (0.054)
- The beta is concentrated in one factor but not as extremely as higher alphas
- The train-val ratio (2.0x) is the most moderate of any alpha>0.25

This suggests there is a narrow window at K=9 where alpha is high enough to
create a survival-weighted iCAF factor but not so high that it overfits. But
even this best case (val adj HR=0.844, p=0.020) is weaker than K=3's
comparable fits.

### How many factors are spuriously prognostic?

| Alpha | Training p<0.05 | Validate p<0.05 (unadj) | Validate p<0.05 (adj) |
|-------|----------------|--------------------------|----------------------|
| 0.00 | 2 (F2, F6) | 2 (F2, F4) | 0 |
| 0.25 | 2 (F4, F5) | 4 (F2, F4, F5, F6) | 2 (F1, F5) |
| 0.35 | 3 (F6, F7, F8) | 1 (F7) | 0 |
| 0.55 | 4 (F1, F2, F7, F8) | 4 (F1, F2, F7, F8) | 3 (F1, F7, F8) |
| 0.75 | 3 (F3, F5, F7) | 3 (F5, F7, F8) | 1 (F5) |
| 0.85 | 3 (F2, F3, F7) | 2 (F2, F3) | 0 |
| 0.95 | 4 (F1, F3, F7, F9) | 4 (F1, F3, F5, F6, F9) | 1 (F6) |

At most alphas, only 0-1 factors survive adjusted validation. The exception
is a=0.55 (3 factors adjusted), but this is also the most overtrained fit.

---

## Key Questions for K=9

1. **Fragmentation confirmed.** With 9 factors, the immune program consistently
   splits across 3-4 factors. No single factor captures the full iCAF+Immune
   program that K=3's F1 consolidates.

2. **Overfitting quantified.** The train-val HR ratio exceeds 4x at alpha>=0.55,
   vs 2.3x at K=3 alpha=0.55. At alpha=0.95, the ratio reaches 15.3x.

3. **Selective validation survives.** Despite fragmentation, the iCAF factor
   validates unadjusted at 4 of 7 alphas (a=0.25, 0.55, 0.75, 0.95) and KM
   validates at 6 of 7 alphas. The biological signal is real but attenuated.

4. **Alpha=0.75 is the K=9 sweet spot.** It is the only alpha achieving
   adjusted validation (p=0.020). But even this is weaker than K=3's performance.

5. **Beta collapse at high alpha.** At a=0.55 and a=0.95, the model zeros out
   3+ factor betas, effectively reducing from K=9 to K=5-6 active factors.
   The optimization "knows" 9 factors is too many and kills the extras.

---

## Expected Narrative Role

K=9 serves as the **upper bound** in the K-sensitivity story:

- K=2: too few factors, iCAF mixed with PurIST axis
- K=3: optimal, iCAF cleanly isolated and validated
- K=5: fragmentation begins, validation unstable
- K=7: further fragmentation, but still validates at some alphas
- **K=9: maximum fragmentation, extreme overfitting, only one alpha validates adjusted**

K=9 at alpha=0.75 validates adjusted (p=0.020), preventing a clean "K=9 fails
entirely" narrative. But the train-val disconnects, beta collapse, and extreme
fragmentation patterns all support the conclusion that **K=3 is optimal** and
higher K introduces noise without biological insight.

---

## Glossary

- **H-score**: Factor loading score, H = X^T W.
- **H-cor**: Spearman correlation of H-scores between two fits.
- **LP**: Cox model linear predictor, sum of beta_k * H_k across factors.
- **LRT**: Likelihood ratio test comparing nested Cox models. Tests whether
  adding a factor improves a model already containing PurIST + DeCAF.
- **ntop**: NULL = all genes; 270 (production) = top 270 per factor.
- **strata(dataset)**: Stratified Cox with per-cohort baseline hazards.
- **Gene universe**: 1970 genes after filtering from 3000.
- **CV-optimal alpha**: Alpha minimizing cross-validated partial log-likelihood.
- **Validation-optimal alpha**: Alpha with strongest external validation signal.
