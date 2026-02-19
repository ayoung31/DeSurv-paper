# K=5 DeSurv Factor Analysis: CV Grid Fits (Expanded Grid)

**For discussion with Jen Jen Yeh**
**Date: 2026-02-18 (amber revision)**
**Companion to: 2026-02-17-k-sensitivity-synthesis_amber.md**
**Data source: cv_grid exhaustive search (fixed lambda=0.3, nu=0.05, 100 inits)**

**Note for Amber:** Expanded from 3 alpha values (0, 0.25, 0.55) to 7
(0, 0.25, 0.35, 0.55, 0.75, 0.85, 0.95). Fits loaded by parameter matching.
Scripts: `inst/cv_grid_training_analysis_amber.R` and
`inst/cv_grid_validation_analysis_amber.R`.

## Key Context from Previous K=5 Analysis

The previous analysis found:
- K5_a0: iCAF absent (H-cor=0.313, no validation signal)
- K5_a25: "purest immune factor" (45% DECODER Immune, best KM p=4.4e-6)
  but only H-cor=0.449 with iCAF program
- K5_a55: iCAF partially recovered (H-cor=0.737) but validation failed
  (unadj p=0.13) -- the fragmentation/overfitting problem

The expanded grid tests whether alpha values between 0.25-0.55 or above 0.55
produce more stable validation, and whether the fragmentation pattern persists.

**Training cohort:** TCGA + CPTAC (n=273)
**Validation cohort:** 4 independent datasets (n=570 pooled):
Dijk (n=90), Moffitt GEO array (n=123), Puleo array (n=288), PACA-AU (n=69).
All validation Cox models use `strata(dataset)`.

**iCAF reference**: The "iCAF program" refers to the factor best matching the
original K=3 Factor 1, characterized by overlap with the **Elyada et al. (2019)
iCAF 35-gene signature** (25 of 35 present in the 1970-gene universe).

---

## Summary: K=5 Recovers iCAF But Fragments It Among Competing Factors

At K=5, the iCAF program steadily emerges with increasing alpha (H-cor rises from
0.313 at alpha=0 to 0.893 at alpha=0.95), but is always accompanied by additional
prognostic factors. Unlike K=3 where the iCAF factor is the sole survival signal,
K=5 spreads prognostic weight across 2-3 factors at every alpha level.

**The validation reveals a persistent pattern.** No K=5 iCAF factor survives
PurIST+DeCAF adjustment in validation at any alpha. The best adjusted validation
p-values are borderline (p=0.056 at a=0.85, p=0.069 at a=0.25 and a=0.35).
Meanwhile, unadjusted validation validates at most alpha values (p<0.05 for
a=0.25 through a=0.95), confirming the iCAF signal exists but overlaps with
existing classifiers when K=5 fragments it.

| | a=0.00 | a=0.25 | a=0.35 | a=0.55 | a=0.75 | a=0.85 | a=0.95 |
|--|--------|--------|--------|--------|--------|--------|--------|
| **Best iCAF factor** | F1 | F4 | F1 | F4 | F3 | F3 | F1 |
| **H-cor with orig F1** | 0.313 | 0.449 | 0.667 | 0.737 | 0.816 | 0.823 | 0.893 |
| **Top-270 overlap** | 47/270 | 121/270 | 105/270 | 146/270 | 159/270 | 144/270 | 144/270 |
| **Elyada iCAF** | 3/25 (12%) | 6/25 (24%) | 2/25 (8%) | 6/25 (24%) | 3/25 (12%) | 4/25 (16%) | 4/25 (16%) |
| **DECODER Immune** | 1/60 (2%) | 27/60 (45%) | 11/60 (18%) | 15/60 (25%) | 15/60 (25%) | 14/60 (23%) | 18/60 (30%) |
| **Train HR (unadj)** | 0.969 | 0.652 | 0.641 | 0.529 | 0.438 | 0.418 | 0.288 |
| **Val HR (unadj)** | 0.971 | 0.797 | 0.789 | 0.884 | 0.826 | 0.818 | 0.783 |
| **Val p (unadj)** | 0.651 | **5.0e-4** | **0.002** | 0.127 | **0.032** | **0.014** | **5.0e-4** |
| **Val p (adj)** | 0.568 | 0.069 | 0.070 | 0.167 | 0.090 | 0.056 | 0.082 |
| **Val KM p** | 0.591 | **4.4e-6** | **3.8e-4** | **0.005** | **3.4e-5** | **1.0e-4** | **9.4e-6** |
| **# prognostic factors** | 1 | 2 | 3 | 3 | 3 | 3 | 2 |

---

## Per-Alpha Detailed Analysis

### K5_a0: Standard NMF (alpha=0.00)

#### Factor identity

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

Best iCAF match: **F1** (rho=0.313, top-270 overlap=47/270)

#### Reference signature overlaps (F1 as best iCAF match, top-270)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 3/25 | 12.0% |
| Elyada myCAF (15) | 0/15 | 0.0% |
| SCISSORS iCAF (17) | 0/17 | 0.0% |
| SCISSORS myCAF (16) | 1/16 | 6.2% |
| DeCAF restCAF (8) | 0/8 | 0.0% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 1/60 | 1.7% |
| DECODER ActStroma (108) | 0/108 | 0.0% |
| DECODER NormStroma (15) | 0/15 | 0.0% |
| DECODER BasalTumor (139) | 19/139 | 13.7% |
| DECODER ClassTumor (166) | 35/166 | 21.1% |

F1 is dominated by Classical Tumor genes (21%), not iCAF. Where do the Immune
genes go? F4 has 33/60 DECODER Immune genes (55%), but F4 anti-correlates
with original F1 (rho=-0.133). At K=5 with alpha=0, the immune program forms
its own factor but is not linked to the iCAF program.

#### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 3/25 | 1/60 | 19/139 | Classical/Basal |
| F2 | 3/25 | 1/60 | 39/139 | Basal Tumor |
| F3 | 2/25 | 5/60 | 2/139 | Mixed |
| **F4** | **7/25** | **33/60** | **0/139** | **Immune** |
| F5 | 3/25 | 7/60 | 6/139 | Mixed |

#### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 (best iCAF) | 0.969 | 0.718 | 0.971 | 0.651 |
| **F2** | **1.320** | **0.004** | **1.301** | **<0.0001** |
| F3 | 0.912 | 0.293 | 0.961 | 0.446 |
| F4 | 0.954 | 0.582 | 0.940 | 0.321 |
| F5 | 0.898 | 0.208 | 0.901 | 0.130 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 (best iCAF) | 0.972 | 0.755 | 0.962 | 0.568 |
| F2 | 0.952 | 0.678 | 1.138 | 0.068 |
| F3 | 1.048 | 0.613 | 1.022 | 0.681 |
| F4 | 0.986 | 0.869 | 0.941 | 0.341 |
| F5 | 0.886 | 0.201 | 0.839 | 0.022 |

**The iCAF factor is not prognostic at alpha=0.** Only F2 is significant, and
F2 captures the basal-like tumor axis, not iCAF. In validation, F2 replicates
strongly (p<0.0001), confirming that standard NMF captures tumor subtype but
not the iCAF program.

**LRT adding F1 (iCAF) to PurIST + DeCAF:** Training p=0.755, Validation p=0.568

**Median-split KM for F1 (iCAF):**
- Training (n=273): High 20.5 vs Low 23.4 months, p=0.795 (reversed!)
- Validation (n=570): High 22.5 vs Low 21.4 months, p=0.591

#### Beta structure

F1=0, F2=1.34e-04, F3=-1.37e-05, F4=-4.28e-05, F5=-1.13e-04

All tiny (no survival penalty -> near-zero betas).

---

### K5_a25: Validation-Optimal Alpha (alpha=0.25)

#### Factor identity

H-score correlation with original K=3:

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

Best iCAF match: **F4** (rho=0.449, top-270 overlap=121/270)

#### Reference signature overlaps (F4, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 6/25 | 24.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | **7/17** | **41.2%** |
| SCISSORS myCAF (16) | 1/16 | 6.2% |
| DeCAF restCAF (8) | 2/8 | 25.0% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | **27/60** | **45.0%** |
| DECODER ActStroma (108) | 5/108 | 4.6% |
| DECODER NormStroma (15) | 5/15 | 33.3% |
| DECODER BasalTumor (139) | 0/139 | 0.0% |
| DECODER ClassTumor (166) | 16/166 | 9.6% |

F4 is remarkably pure: **45% DECODER Immune, 0% Basal Tumor**. The SCISSORS
iCAF overlap (41%) and DeCAF restCAF overlap (25%) are the highest across all
fits. This is arguably the purest immune/iCAF factor in the entire grid.

#### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 6/25 | 23/60 | 2/139 | Stroma/Immune |
| F2 | 3/25 | 0/60 | 23/139 | Basal Tumor |
| F3 | 0/25 | 0/60 | 43/139 | Basal Tumor (dominant) |
| **F4** | **6/25** | **27/60** | **0/139** | **iCAF/Immune (purest)** |
| F5 | 7/25 | 4/60 | 13/139 | Mixed |

F1 also has 23/60 DECODER Immune genes and 6/25 Elyada iCAF genes. The
immune/iCAF signal is split between F1 and F4. Both capture parts of it.

#### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.011 | 0.895 | 1.055 | 0.377 |
| **F2** | **1.336** | **0.002** | **1.224** | **0.003** |
| F3 | 1.009 | 0.921 | 1.022 | 0.732 |
| **F4 (iCAF)** | **0.652** | **<0.0001** | **0.797** | **0.0005** |
| F5 | 1.009 | 0.915 | 1.173 | 0.052 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.003 | 0.975 | 1.008 | 0.893 |
| **F2** | **1.256** | **0.017** | 1.093 | 0.208 |
| F3 | 0.955 | 0.613 | 0.941 | 0.378 |
| **F4 (iCAF)** | **0.767** | **0.007** | 0.879 | 0.069 |
| F5 | 0.860 | 0.104 | 1.037 | 0.671 |

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

#### Beta structure

F1=6.02e-05, F2=**5.18e-04**, F3=-6.72e-05, F4=**-6.77e-04** (dominant), F5=-9.62e-05

F4 (iCAF) has the dominant beta, but F2 (Basal) also has substantial weight.

---

### K5_a35: Intermediate Alpha (alpha=0.35)

#### Factor identity

H-score correlation with original K=3:

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| **F1** | **0.667** | -0.709 | 0.643  |
| F2     | 0.453   | 0.043   | 0.190  |
| F3     | -0.247  | -0.491  | **0.742** |
| F4     | -0.219  | **0.985** | -0.851 |
| F5     | 0.218   | -0.574  | 0.696  |

F1 is the best iCAF match (rho=0.667). F4 perfectly matches original F2
(rho=0.985) and F3 matches original F3 (rho=0.742). The stroma/basal axis is
again well preserved; the iCAF program is forming more strongly than at a=0.25
but with a different factor assignment. F2 also correlates with original F1
(rho=0.453), suggesting a partial split.

Best iCAF match: **F1** (rho=0.667, top-270 overlap=105/270)

#### Reference signature overlaps (F1, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 2/25 | 8.0% |
| Elyada myCAF (15) | 1/15 | 6.7% |
| SCISSORS iCAF (17) | 3/17 | 17.6% |
| SCISSORS myCAF (16) | 1/16 | 6.2% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 11/60 | 18.3% |
| DECODER ActStroma (108) | 1/108 | 0.9% |
| DECODER NormStroma (15) | **8/15** | **53.3%** |
| DECODER BasalTumor (139) | 4/139 | 2.9% |
| DECODER ClassTumor (166) | **50/166** | **30.1%** |

F1 has the highest NormalStroma overlap (53.3%) of any fit, but is diluted by
ClassicalTumor genes (30.1%). The Immune overlap (18.3%) is modest compared
to K5_a25's F4 (45.0%). The iCAF program is forming but mixed with tumor genes.

#### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 2/25 | 11/60 | 4/139 | iCAF/Classical (mixed) |
| F2 | 3/25 | 7/60 | 9/139 | Mixed |
| F3 | 1/25 | 0/60 | **83/139** | Basal Tumor (dominant) |
| F4 | 5/25 | 20/60 | 1/139 | Immune |
| F5 | 4/25 | 1/60 | 13/139 | Mixed/Basal |

Note: F4 has 20/60 DECODER Immune genes but only rho=-0.219 with original F1.
The immune signal splits: F1 gets 11/60 Immune + 8/15 NormStroma (the iCAF-like
mix), while F4 gets 20/60 Immune alone (a purer immune factor). The fragmentation
pattern continues.

#### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.641** | **<0.0001** | **0.789** | **0.002** |
| **F2** | **0.716** | **<0.0001** | 0.962 | 0.612 |
| **F3** | **1.348** | **0.001** | **1.336** | **<0.0001** |
| F4 | 1.001 | 0.992 | 1.052 | 0.392 |
| F5 | 1.078 | 0.398 | 0.965 | 0.619 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.754** | **0.006** | 0.860 | 0.070 |
| **F2** | **0.666** | **<0.0001** | 0.896 | 0.179 |
| F3 | 1.042 | 0.704 | 1.121 | 0.177 |
| F4 | 1.005 | 0.954 | 1.013 | 0.826 |
| F5 | 1.198 | 0.049 | 0.961 | 0.587 |

Three training factors are prognostic unadjusted (F1, F2, F3). The iCAF factor
validates unadjusted (HR=0.789, p=0.002) and is borderline adjusted (HR=0.860,
p=0.070). F2 is highly prognostic in training (p<0.0001) but does NOT replicate
(val p=0.612), a clear overfitting signal for that factor. F3 (Basal Tumor)
replicates strongly (val p<0.0001).

**LRT adding F1 (iCAF) to PurIST + DeCAF:** Training p=0.005, Validation p=0.072

**Median-split KM for F1 (iCAF):**
- Training (n=273): High 23.8 vs Low 15.6 months, **p=3.1e-4**
- Validation (n=570): High 24.3 vs Low 19.0 months, **p=3.8e-4**

The KM validates with a 5.3-month gap (p=3.8e-4), comparable to K5_a55 but
weaker than K5_a25.

#### Beta structure

F1=**-7.56e-04**, F2=**-7.58e-04**, F3=3.32e-04, F4=6.45e-05, F5=7.45e-04

Two factors share the dominant negative beta (F1 and F2), confirming the
fragmentation of the survival signal across multiple factors.

---

### K5_a55: Matched to K=3 CV-Optimal Alpha (alpha=0.55)

#### Factor identity

H-score correlation with original K=3:

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

Best iCAF match: **F4** (rho=0.737, top-270 overlap=146/270)

#### Reference signature overlaps (F4, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 6/25 | 24.0% |
| Elyada myCAF (15) | 0/15 | 0.0% |
| SCISSORS iCAF (17) | **7/17** | **41.2%** |
| SCISSORS myCAF (16) | 2/16 | 12.5% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | **2/4** | **50.0%** |
| DECODER Immune (60) | 15/60 | 25.0% |
| DECODER ActStroma (108) | 4/108 | 3.7% |
| DECODER NormStroma (15) | **6/15** | **40.0%** |
| DECODER BasalTumor (139) | 15/139 | 10.8% |
| DECODER ClassTumor (166) | 42/166 | 25.3% |

F4 maintains SCISSORS iCAF overlap (41%) and gains NormalStroma (40%). But
unlike K5_a25 where F4 had 45% DECODER Immune overlap with zero Basal Tumor,
here F4 has only 25% Immune and picks up 11% Basal Tumor and 25% Classical
Tumor. The factor is becoming less pure; it is absorbing more tumor genes as
alpha increases.

#### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 0/25 | 0/60 | 39/139 | Basal/Stroma |
| F2 | 3/25 | 5/60 | 12/139 | Classical Tumor |
| F3 | 5/25 | 5/60 | 31/139 | Mixed/Basal |
| **F4** | **6/25** | **15/60** | 15/139 | **iCAF/Mixed** |
| F5 | 5/25 | **24/60** | 1/139 | Immune/Classical |

At K=5 alpha=0.55, the Immune genes are split: F4 gets 15/60, F5 gets 24/60.
Compare to K=3 alpha=0.55 where F1 (iCAF) gets 16/60 and F3 gets 16/60, a
similar split, but at K=3 the iCAF factor combines Immune + NormalStroma into
a single coherent program, while K=5's F4 dilutes this with Classical Tumor
genes (42/166 = 25%).

#### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1** | **1.318** | **0.003** | **1.166** | **0.039** |
| F2 | 1.008 | 0.922 | 1.119 | 0.128 |
| **F3** | **1.218** | **0.029** | 1.166 | 0.054 |
| **F4 (iCAF)** | **0.529** | **<0.0001** | 0.884 | 0.127 |
| F5 | 0.920 | 0.322 | 0.955 | 0.502 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.178 | 0.091 | 1.001 | 0.991 |
| F2 | 0.941 | 0.508 | 1.014 | 0.853 |
| F3 | 1.047 | 0.638 | 1.014 | 0.866 |
| **F4 (iCAF)** | **0.562** | **<0.0001** | 0.888 | 0.167 |
| F5 | 0.986 | 0.877 | 0.962 | 0.579 |

**Three training factors are prognostic but only F4 survives PurIST+DeCAF
adjustment in training.** In validation, the iCAF factor DOES NOT reach
significance, neither unadjusted (p=0.127) nor adjusted (p=0.167). This is a
critical failure: K=5 at this alpha overtrains the survival signal. The
extremely strong training HR (0.529) does not replicate (val HR=0.884).

**LRT adding F4 (iCAF) to PurIST + DeCAF:** Training p=4.3e-9, Validation p=0.170

**Median-split KM for F4 (iCAF):**
- Training (n=273): High 24.6 vs Low 14.2 months, **p=3.0e-7**
- Validation (n=570): High 25.3 vs Low 19.2 months, **p=0.005**

The KM still validates (6.1-month gap, p=0.005) despite the Cox model failure,
suggesting the iCAF signal exists but is spread across multiple factors that
individually don't reach significance.

#### Beta structure

F1=7.25e-04, F2=-1.12e-04, F3=4.15e-04, F4=**-1.53e-03** (dominant), F5=-4.42e-05

F4 (iCAF) dominates the LP, but F1 and F3 also have substantial betas.

---

### K5_a75: High Alpha (alpha=0.75)

#### Factor identity

H-score correlation with original K=3:

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.218  | **0.773** | -0.547 |
| F2     | 0.018   | -0.733  | **0.945** |
| **F3** | **0.816** | -0.586 | 0.544  |
| F4     | -0.381  | **0.872** | -0.632 |
| F5     | 0.411   | 0.183   | 0.003  |

F3 is the best iCAF match (rho=0.816). F2 matches original F3 (rho=0.945).
F1 and F4 both correlate with original F2 (rho=0.773 and 0.872), continuing
the pattern of the stroma axis splitting across multiple factors.

Best iCAF match: **F3** (rho=0.816, top-270 overlap=159/270)

#### Reference signature overlaps (F3, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 3/25 | 12.0% |
| Elyada myCAF (15) | 1/15 | 6.7% |
| SCISSORS iCAF (17) | 2/17 | 11.8% |
| SCISSORS myCAF (16) | 0/16 | 0.0% |
| DeCAF restCAF (8) | 0/8 | 0.0% |
| DeCAF proCAF (4) | 0/4 | 0.0% |
| DECODER Immune (60) | 15/60 | 25.0% |
| DECODER ActStroma (108) | 3/108 | 2.8% |
| DECODER NormStroma (15) | **8/15** | **53.3%** |
| DECODER BasalTumor (139) | 16/139 | 11.5% |
| DECODER ClassTumor (166) | 38/166 | 22.9% |

F3 has strong NormalStroma (53.3%) and moderate Immune (25.0%), but is diluted
by ClassicalTumor (22.9%) and BasalTumor (11.5%). The Elyada iCAF overlap
drops to 12% and SCISSORS iCAF to 12%, lower than at a=0.25 or a=0.55. The
factor is broadly mixed rather than a coherent iCAF program.

#### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 2/25 | 2/60 | 14/139 | Mixed |
| F2 | 1/25 | 0/60 | **40/139** | Basal Tumor |
| **F3** | **3/25** | **15/60** | 16/139 | **iCAF/Mixed** |
| F4 | 4/25 | 17/60 | 19/139 | Immune/Basal |
| F5 | 5/25 | 5/60 | 4/139 | Mixed |

The Immune genes are split between F3 (15/60) and F4 (17/60). F4 has more
DECODER Immune genes than the iCAF factor itself, plus 19/139 Basal Tumor.
The immune signal is fragmenting across factors with mixed identities.

#### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.041 | 0.639 | 1.065 | 0.398 |
| **F2** | **1.270** | **0.009** | 1.154 | 0.073 |
| **F3 (iCAF)** | **0.438** | **<0.0001** | **0.826** | **0.032** |
| F4 | 1.139 | 0.155 | **1.197** | **0.004** |
| **F5** | **0.812** | **0.007** | 0.896 | 0.133 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 1.058 | 0.527 | 0.996 | 0.959 |
| F2 | 1.165 | 0.105 | 0.996 | 0.959 |
| **F3 (iCAF)** | **0.465** | **<0.0001** | 0.853 | 0.090 |
| F4 | 0.987 | 0.894 | 1.072 | 0.300 |
| F5 | 0.927 | 0.398 | 0.895 | 0.149 |

Three training factors are prognostic (F2, F3, F5). The iCAF factor F3 is
extremely strong in training (HR=0.438, p<0.0001) and validates unadjusted
(HR=0.826, p=0.032) but does not survive adjustment (HR=0.853, p=0.090).
Interestingly, F4 is not prognostic in training but reaches significance in
validation (p=0.004) -- a cross-validation instability signal.

**LRT adding F3 (iCAF) to PurIST + DeCAF:** Training p=1.2e-13, Validation p=0.092

**Median-split KM for F3 (iCAF):**
- Training (n=273): High 30.4 vs Low 13.3 months, **p=1.2e-9**
- Validation (n=570): High 29.0 vs Low 19.0 months, **p=3.4e-5**

The KM validates strongly (10.0-month gap, p=3.4e-5), the largest validation
gap in the entire K=5 grid. The training gap (17.1 months) is very large,
suggesting substantial overfitting of the beta weights.

#### Beta structure

F1=6.55e-05, F2=8.65e-04, F3=**-1.71e-03** (dominant), F4=0, F5=0

F3 (iCAF) dominates the LP. Two factors have exactly zero beta (F4, F5),
meaning the model selected only 3 of 5 factors for the survival component.

---

### K5_a85: Strong Survival Penalty (alpha=0.85)

#### Factor identity

H-score correlation with original K=3:

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| F1     | -0.018  | **0.848** | -0.624 |
| F2     | 0.070   | -0.377  | **0.686** |
| **F3** | **0.823** | -0.600 | 0.549  |
| F4     | -0.337  | -0.310  | 0.683  |
| F5     | -0.131  | **0.942** | -0.788 |

F3 is the best iCAF match (rho=0.823), nearly identical to K5_a75. F5 matches
original F2 extremely well (rho=0.942) and F1 also matches original F2
(rho=0.848). The stroma axis continues splitting across two factors.

Best iCAF match: **F3** (rho=0.823, top-270 overlap=144/270)

#### Reference signature overlaps (F3, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 4/25 | 16.0% |
| Elyada myCAF (15) | 2/15 | 13.3% |
| SCISSORS iCAF (17) | 3/17 | 17.6% |
| SCISSORS myCAF (16) | 1/16 | 6.2% |
| DeCAF restCAF (8) | 0/8 | 0.0% |
| DeCAF proCAF (4) | 1/4 | 25.0% |
| DECODER Immune (60) | 14/60 | 23.3% |
| DECODER ActStroma (108) | 4/108 | 3.7% |
| DECODER NormStroma (15) | **6/15** | **40.0%** |
| DECODER BasalTumor (139) | 17/139 | 12.2% |
| DECODER ClassTumor (166) | 36/166 | 21.7% |

Similar to K5_a75: NormalStroma (40%) and Immune (23%) are present, but
ClassicalTumor (22%) and BasalTumor (12%) dilute the factor. The signature
overlaps are plateauing; the iCAF factor is not becoming purer at higher alpha.

#### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| F1 | 4/25 | 5/60 | 7/139 | Mixed |
| F2 | 4/25 | 2/60 | 21/139 | Basal |
| **F3** | **4/25** | **14/60** | 17/139 | **iCAF/Mixed** |
| F4 | 3/25 | 0/60 | **61/139** | Basal Tumor (dominant) |
| F5 | 4/25 | 10/60 | 2/139 | Immune |

F4 absorbs 61/139 Basal Tumor genes (the largest concentration of any single
factor across the grid). F5 has 10/60 Immune genes but only rho=-0.131 with
original F1. The iCAF program (F3) has only 14/60 Immune, while the immune
signal bleeds into F1 (5/60) and F5 (10/60).

#### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.927 | 0.349 | 1.002 | 0.982 |
| **F2** | **1.266** | **0.011** | 1.079 | 0.334 |
| **F3 (iCAF)** | **0.418** | **<0.0001** | **0.818** | **0.014** |
| **F4** | **1.473** | **0.0001** | **1.339** | **0.0003** |
| F5 | 0.925 | 0.352 | 1.015 | 0.823 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| F1 | 0.951 | 0.568 | 0.959 | 0.568 |
| **F2** | **1.317** | **0.005** | 0.976 | 0.772 |
| **F3 (iCAF)** | **0.447** | **<0.0001** | 0.847 | 0.056 |
| F4 | 1.157 | 0.200 | 1.104 | 0.283 |
| F5 | 0.945 | 0.520 | 0.984 | 0.817 |

Three training factors are prognostic unadjusted. The iCAF factor F3 has the
strongest training HR in the entire K=5 grid below a=0.95 (HR=0.418, p<0.0001).
It validates unadjusted (HR=0.818, p=0.014) and is the closest to surviving
adjustment (HR=0.847, **p=0.056**) of any K=5 fit. F4 replicates unadjusted
(val p=0.0003) but does not survive adjustment.

**LRT adding F3 (iCAF) to PurIST + DeCAF:** Training p=5.7e-15, Validation p=0.058

**Median-split KM for F3 (iCAF):**
- Training (n=273): High 35.3 vs Low 13.1 months, **p=6.4e-13**
- Validation (n=570): High 26.8 vs Low 19.0 months, **p=1.0e-4**

The training KM gap (22.2 months) is enormous, far exceeding the validation gap
(7.8 months). This 2.8x ratio between training and validation gaps is a clear
overfitting indicator; the model is learning training-specific survival patterns
that attenuate in external data.

#### Beta structure

F1=0, F2=9.25e-04, F3=**-1.91e-03** (dominant), F4=3.05e-04, F5=-2.69e-04

F3 (iCAF) has the dominant beta. F1 has exactly zero beta, so only 4 of 5
factors contribute to the survival component.

---

### K5_a95: Near-Maximum Survival Penalty (alpha=0.95)

#### Factor identity

H-score correlation with original K=3:

|        | Orig_F1 | Orig_F2 | Orig_F3 |
|--------|---------|---------|---------|
| **F1** | **0.893** | -0.506 | 0.330  |
| F2     | 0.152   | -0.746  | **0.889** |
| F3     | 0.244   | -0.847  | **0.859** |
| F4     | -0.479  | **0.832** | -0.557 |
| F5     | -0.110  | **0.950** | -0.857 |

F1 matches original F1 the best of any K=5 fit (rho=0.893). F5 matches original
F2 (rho=0.950). F2 and F3 both correlate with original F3 (rho=0.889 and 0.859).
The original F3 axis splits into two K=5 factors.

Best iCAF match: **F1** (rho=0.893, top-270 overlap=144/270)

#### Reference signature overlaps (F1, top-270 genes)

| Signature | Overlap | % |
|-----------|---------|---|
| Elyada iCAF (25) | 4/25 | 16.0% |
| Elyada myCAF (15) | 1/15 | 6.7% |
| SCISSORS iCAF (17) | **6/17** | **35.3%** |
| SCISSORS myCAF (16) | 1/16 | 6.2% |
| DeCAF restCAF (8) | 1/8 | 12.5% |
| DeCAF proCAF (4) | 1/4 | 25.0% |
| DECODER Immune (60) | 18/60 | 30.0% |
| DECODER ActStroma (108) | 3/108 | 2.8% |
| DECODER NormStroma (15) | **8/15** | **53.3%** |
| DECODER BasalTumor (139) | 15/139 | 10.8% |
| DECODER ClassTumor (166) | 30/166 | 18.1% |

F1 has the best SCISSORS iCAF overlap outside of a=0.25/a=0.55 (35.3%), strong
NormalStroma (53.3%), and 30% Immune. But it also carries 18% ClassicalTumor
and 11% BasalTumor. This is the highest H-cor with original F1 (0.893) and
the most iCAF-like factor at K=5, but still more mixed than K=3's F1.

#### Ref overlap ALL factors (Elyada iCAF / DECODER Immune / DECODER BasalTumor)

| Factor | Elyada iCAF | DECODER Immune | DECODER Basal | Identity |
|--------|------------|----------------|--------------|----------|
| **F1** | **4/25** | **18/60** | 15/139 | **iCAF/Mixed** |
| F2 | 3/25 | 0/60 | 32/139 | Basal Tumor |
| F3 | 2/25 | 0/60 | 23/139 | Basal |
| F4 | 5/25 | 6/60 | 22/139 | Mixed/Basal |
| F5 | 2/25 | **20/60** | 0/139 | Immune |

F5 captures 20/60 DECODER Immune genes with zero Basal Tumor, a pure immune
factor (like K5_a25's F4 with 27/60). The immune signal continues to split:
F1 gets 18/60, F5 gets 20/60. At K=3, these would be combined into a single
coherent iCAF+Immune program.

#### Per-factor Cox: Training AND Validation

**Unadjusted**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.288** | **<0.0001** | **0.783** | **0.0005** |
| F2 | 1.010 | 0.912 | 1.146 | 0.077 |
| F3 | 1.094 | 0.312 | 0.958 | 0.480 |
| **F4** | **1.288** | **0.008** | **1.261** | **0.001** |
| F5 | 0.930 | 0.390 | 0.955 | 0.415 |

**Adjusted for PurIST + DeCAF**

| Factor | Train HR | Train p | Val HR | Val p |
|--------|----------|---------|--------|-------|
| **F1 (iCAF)** | **0.286** | **<0.0001** | 0.876 | 0.082 |
| F2 | 0.877 | 0.156 | 0.981 | 0.821 |
| F3 | 1.167 | 0.091 | 0.972 | 0.651 |
| F4 | 1.111 | 0.303 | 1.081 | 0.324 |
| F5 | 0.978 | 0.798 | 0.962 | 0.504 |

F1 (iCAF) has the most extreme training HR in the entire grid (0.288, p<0.0001),
but the validation HR is 0.783 (p=0.0005), a 2.7x attenuation. The adjusted
validation is borderline (p=0.082). Only two factors are prognostic in training
(F1 + F4), fewer than the 3 seen at a=0.55-0.85, suggesting the survival
penalty has concentrated signal more effectively. F4 replicates unadjusted
(val p=0.001).

**LRT adding F1 (iCAF) to PurIST + DeCAF:** Training p=2.0e-24, Validation p=0.083

**Median-split KM for F1 (iCAF):**
- Training (n=273): High 37.7 vs Low 12.5 months, **p<0.0001**
- Validation (n=570): High 27.0 vs Low 18.0 months, **p=9.4e-6**

The training KM gap (25.2 months) is the largest in the entire grid across all
K values, while the validation gap (9.0 months, p=9.4e-6) is the second best
at K=5 after K5_a25's 7.6 months at a stronger p-value. The 2.8x ratio between
training and validation gaps confirms the overfitting pattern.

#### Beta structure

F1=**-1.83e-03** (dominant), F2=-7.90e-05, F3=5.29e-04, F4=1.04e-04, F5=-1.63e-05

F1 (iCAF) has the dominant beta. The second-largest beta is F3 (5.29e-04),
3.5x smaller. At this extreme alpha, the model concentrates on fewer factors.

---

## The K=5 Limitation: Fragmentation and Validation Failure

### The fragmentation problem, confirmed across all 7 alpha values

At K=5, the model has enough capacity to separate biological programs into
finer components. But this creates two problems visible across the full grid:

1. **Immune/iCAF splitting is universal.** At every alpha value, the DECODER
   Immune genes split across at least two factors. The iCAF factor never captures
   more than 27/60 Immune genes (at a=0.25), and at most alpha values gets only
   14-18/60. A separate factor always absorbs the remaining Immune genes. At K=3,
   these are combined into a single coherent program.

2. **No K=5 fit survives PurIST+DeCAF adjustment in validation.** The best
   adjusted validation p-values are: a=0.85 (p=0.056), a=0.25 (p=0.069),
   a=0.35 (p=0.070), a=0.95 (p=0.082). All borderline, none significant.
   This is the defining limitation of K=5: the fragmented iCAF signal overlaps
   with existing classifiers rather than providing independent prognostic value.

3. **Massive training-validation attenuation.** Training HRs range from 0.288
   to 0.641 at alpha >= 0.25, but validation HRs are always 0.783-0.884. The
   attenuation ratio worsens with increasing alpha (1.2x at a=0.25, 2.7x at
   a=0.95), confirming that higher survival penalty at K=5 leads to overfitting.

4. **KM validates robustly despite Cox failure.** Median-split KM p-values are
   significant (p < 0.005) at every alpha >= 0.25, with validation gaps of
   5.3-9.0 months. The iCAF signal exists in the data; it just cannot be
   captured as an independent, PurIST/DeCAF-adjusted signal at K=5.

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
