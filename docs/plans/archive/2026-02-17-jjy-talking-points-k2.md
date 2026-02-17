# K=2 DeSurv Analysis: How Do Conclusions Change?

**For discussion with Jen Jen Yeh**
**Date: 2026-02-17**
**Companion to: 2026-02-16-jjy-talking-points.md (K=3 analysis)**

## Motivation

The K=3 model's central claim is that Factor 1 captures a novel iCAF/B cell
program with independent prognostic value beyond PurIST and DeCAF. A natural
question: **what happens at K=2?** If the third axis is genuinely novel, it
should be absent at K=2. If K=2 recovers the same signal, Factor 1's novelty
is less convincing.

## K=2 Fit Summary

K=2 was fit using the same hyperparameters as K=3 (alpha=0.33, lambda=0.35,
nu=0.056) on the same 273-sample filtered training set (1970 genes). Only
`k` was changed from 3 to 2.

| Metric | K=2 | K=3 |
|--------|-----|-----|
| c-index (training) | **0.697** | **0.786** |
| Non-zero betas | 2/2 | 2/3 (beta_2 = 0) |
| Beta values | [0.000609, -0.000925] | [-0.001003, 0, 9.36e-05] |
| Inter-factor rho | **-0.928** | F2-F3: **-0.971** |

The c-index drops substantially (0.697 vs 0.786), indicating K=2 is a
meaningfully worse model for survival prediction. Both K=2 betas are
non-zero, meaning both factors contribute to the linear predictor --
unlike K=3 where beta_2 = 0 and the LP is dominated by Factor 1.

## The Answer: K=2 Collapses to PurIST/DeCAF

At K=2, DeSurv recovers a **single axis** (PurIST/DeCAF) split into two
anti-correlated factors. The novel iCAF/B cell third axis from K=3 is lost.

| Factor | Biology | DECODER Compartments | PurIST/DeCAF Alignment |
|--------|---------|---------------------|----------------------|
| K2-F1 | Basal-like + activated stroma | Basal Tumor (31), Activated Stroma (24), Classical Tumor (12) | Higher in Basal-like (p=2.6e-18), higher in proCAF (p=5.7e-16) |
| K2-F2 | Classical + mixed stroma | Classical Tumor (32), Activated Stroma (16), Basal Tumor (6) | Higher in Classical (p=1.5e-22), higher in restCAF (p=1.3e-12) |

Both factors track PurIST and DeCAF simultaneously, from opposite ends.
They are near-perfectly anti-correlated (rho = -0.928), meaning they
contain almost the same information with reversed sign. This is analogous
to K=3's Factors 2 and 3 (rho = -0.971), which also collapsed the
Classical/Basal-like axis into two anti-correlated factors.

## Gene-Level Evidence

### K2-F1: Basal-like / Activated Stroma

Top 30 genes: IQGAP1, CD74, MAL2, TUBA1C, P4HB, C1S, B2M, LTBP2, DSP,
CD2AP, COL12A1, RGS5, LAPTM5, CD68, EHF, COL8A1, CDCP1, DCN, JUP, FLNA,
LGALS3BP, SOCS3, C4A, GSN, THBS2, SFRP4, MT2A, SERPINF1, ALOX5, SLC20A1

| Reference Signature | Overlap | Shared Genes |
|---------------------|---------|--------------|
| DeCAF restCAF (9) | 0 | -- |
| DeCAF proCAF (9) | 0 | -- |
| PurIST Classical (8) | 1 | REG4 |
| PurIST Basal-like (8) | **5** | ITGA3, KRT6A, BCAR3, KRT5, S100A2 |
| Elyada iCAF (35) | 2 | SOD2, CCDC80 |
| Elyada myCAF (16) | 2 | POSTN, CTHRC1 |
| SCISSORS iCAF (25) | 0 | -- |
| SCISSORS myCAF (25) | 3 | COL8A1, ADAMTS12, EDIL3 |

K2-F1 is dominated by **Basal Tumor genes (31)** and **Activated Stroma (24)**.
It shares 5/8 PurIST Basal-like genes and 3 SCISSORS myCAF genes. This factor
captures the "bad prognosis" end: Basal-like tumor + myofibroblast stroma.

### K2-F2: Classical / Mixed Stroma

Top 30 genes: ACTN4, HLA-DRB1, FN1, FAM129B, NEAT1, THY1, ASPH, SFRP2,
MXRA7, APOD, VWF, GFPT1, TPM1, PDIA4, CD151, C1R, HSPG2, KRT17, LAMB2,
TGFBI, BAIAP2L1, ACTN1, AXL, IGFBP7, ZFP36, TSPAN8, ATP11A, CD59, MYO6,
PIK3IP1

| Reference Signature | Overlap | Shared Genes |
|---------------------|---------|--------------|
| DeCAF restCAF (9) | **2** | FBLN5, KIAA1217 |
| DeCAF proCAF (9) | 0 | -- |
| PurIST Classical (8) | **5** | ANXA10, SLC40A1, LGALS4, CLDN18, GATA6 |
| PurIST Basal-like (8) | 0 | -- |
| Elyada iCAF (35) | 3 | LMNA, SGK1, FBLN1 |
| Elyada myCAF (16) | **5** | TPM1, IGFBP7, ACTA2, TAGLN, HOPX |
| SCISSORS iCAF (25) | 2 | FBLN5, LIF |
| SCISSORS myCAF (25) | 1 | KIAA1217 |

K2-F2 is dominated by **Classical Tumor genes (32)** and **Activated Stroma (16)**.
It shares 5/8 PurIST Classical genes and 2 DeCAF restCAF genes. This factor
captures the "good prognosis" end: Classical tumor + mixed stroma. Note the
5 Elyada myCAF gene overlap -- both K=2 factors include stroma genes (unlike
K=3 where stroma was concentrated in Factor 2 and Factor 1 had minimal
stroma overlap).

### K2 vs K3: Cross-Factor Gene Overlap

| | K3-F1 (iCAF/B cell) | K3-F2 (Classical/Stroma) | K3-F3 (Basal-like) |
|--|-----|-----|-----|
| **K2-F1** | 22 | 53 | 60 |
| **K2-F2** | 56 | 28 | 34 |

K2-F1 absorbs genes primarily from K3-F2 (53) and K3-F3 (60) -- the tumor
subtype axis. K2-F2 absorbs the most from K3-F1 (56) -- but 71% of K3-F1's
genes go to **neither** K=2 factor:

| K3-F1 iCAF genes | Count | % of 270 |
|-------------------|-------|----------|
| In K2-F1 | 22 | 8.1% |
| In K2-F2 | 56 | 20.7% |
| **In neither** | **192** | **71.1%** |

The iCAF/B cell program from K3 is largely **lost** at K=2. Only 29% of its
genes survive in either K=2 factor, and those are diluted by the dominant
PurIST/DeCAF signal. The 192 lost genes include the core iCAF markers
(CXCL14, CCL19, CCL21 pathway genes) that drove Factor 1's independent
prognostic value at K=3.

### K2 vs K3: H Score Correlations (Model H, Training)

| | K3-F1 | K3-F2 | K3-F3 |
|--|-------|-------|-------|
| **K2-F1** | **-0.858** | 0.494 | -0.333 |
| **K2-F2** | **0.838** | -0.409 | 0.262 |

K2-F1 is strongly anti-correlated with K3-F1 (rho = -0.858): when K3-F1
(iCAF) is high, K2-F1 (Basal/stroma) is low. K2-F2 is strongly positively
correlated with K3-F1 (rho = 0.838): K2-F2 (Classical) partially recovers
the iCAF direction but conflates it with the Classical tumor signal.

The K2 LP correlates moderately with K3 LP (rho = 0.780) -- substantial
overlap but not identical, reflecting the loss of the third axis.

## Per-Factor Prognostic Value (K=2, Training Only)

*Note: K=2 was fit on training data only. No external validation was performed.*

### Unadjusted Cox

| Factor | HR (per unit) | p-value |
|--------|--------------|---------|
| K2-F1 (Basal/stroma) | 1.000268 | **3.6e-06** |
| K2-F2 (Classical) | 0.999711 | **3.3e-08** |

Both factors are prognostic unadjusted. K2-F1 (higher = worse, consistent
with Basal-like prognosis) and K2-F2 (higher = better, consistent with
Classical prognosis).

### Adjusted for PurIST + DeCAF

| Factor | p-value |
|--------|---------|
| K2-F1 | **0.099** (NOT significant) |
| K2-F2 | **0.028** (barely significant) |

**Critical difference from K=3**: At K=3, Factor 1 was highly significant
after PurIST+DeCAF adjustment (training p = 3.1e-13, validation p = 0.004).
At K=2, **neither factor reaches strong significance** after adjustment.
K2-F2's p = 0.028 is marginal and has not been validated externally.

### Joint Model (PurIST + DeCAF + F1 + F2)

| Factor | p-value |
|--------|---------|
| K2-F1 | 0.47 |
| K2-F2 | 0.11 |

Neither factor is significant in the joint model. The near-perfect
anti-correlation (rho = -0.928) causes multicollinearity that absorbs
both signals.

### Median-Split KM (Training)

| Factor Group | Median OS | Log-rank p |
|-------------|-----------|------------|
| K2-F1-High | 15.3 months | 1.3e-05 |
| K2-F1-Low | 24.6 months | |
| K2-F2-High | 25.4 months | 1.2e-07 |
| K2-F2-Low | 14.5 months | |

Both factors show survival separation, but the groups are essentially
mirror images (because F1 and F2 are anti-correlated). F1-High = F2-Low
and vice versa. This is one axis, not two.

### LRT: K2-LP vs PurIST + DeCAF

Adding the K=2 LP to PurIST + DeCAF is marginally significant (LRT p = 0.042),
compared to K=3 where adding the LP was highly significant (LRT p < 2.2e-16
in training, p = 0.001 in validation).

## Comparison: K=2 vs K=3 Conclusions

| Question | K=3 Answer | K=2 Answer |
|----------|-----------|-----------|
| **c-index** | 0.786 | 0.697 (worse) |
| **Number of independent axes** | 2 (tumor subtype + iCAF) | 1 (tumor/stroma combined) |
| **Novel prognostic signal?** | Yes -- Factor 1 (iCAF/B cell) is independently prognostic (adj p = 0.004 in validation) | No -- neither factor significant after PurIST+DeCAF (p = 0.099 and 0.028, training only) |
| **Adds to PurIST + DeCAF?** | Yes (LRT p = 0.001 in validation) | Marginal (LRT p = 0.042, training only) |
| **Inter-factor structure** | F1 is distinct from F2/F3; F2/F3 are anti-correlated | F1 and F2 are anti-correlated (single axis) |
| **iCAF/B cell program** | Captured in Factor 1 (270 genes, dominated by chemokines + B cell markers) | **Lost** -- 71% of K3-F1 genes absent from both K2 factors |
| **PurIST alignment** | F2/F3 track PurIST; F1 does not | Both factors track PurIST |
| **DeCAF alignment** | F2 tracks DeCAF; F1 minimally overlaps | Both factors track DeCAF |
| **External validation** | Yes (n=570, 4 cohorts) | No (training only) |

## What This Means for the K=3 Story

The K=2 analysis **strengthens** the K=3 argument in two ways:

### 1. K=2 demonstrates what you lose without the third axis

At K=2, both factors collapse to the existing PurIST/DeCAF axis. The model
has no room for a separate microenvironmental program. The c-index drops,
neither factor has independent prognostic value after adjustment, and the
LP adds only marginal value to existing classifiers. **K=2 is what you get
if DeSurv merely rediscovers PurIST + DeCAF in continuous form.**

### 2. K=3's extra factor is not just splitting an existing axis

If going from K=2 to K=3 simply subdivided one of the two existing factors,
we'd expect the K=3 factors to have high gene overlap with K=2 factors.
Instead, K=3-Factor 1 has **low overlap** with both K=2 factors (22 and 56
genes out of 270) and 71% of its genes are absent from both K=2 factors
entirely. The third factor brings in genuinely new gene content -- the
iCAF-secreted chemokines and B cell markers that K=2 cannot accommodate.

### The pitch (K=2 edition)

"At K=2, DeSurv recovers the same Basal-like vs Classical axis that PurIST
and DeCAF already capture. Both factors are anti-correlated (rho = -0.93),
neither is independently prognostic after adjustment, and the model's
survival prediction is substantially worse (c-index 0.70 vs 0.79). Going to
K=3 is not merely splitting this axis into three -- it introduces a
genuinely new iCAF-associated program (71% of its genes absent from both
K=2 factors) that is the sole source of independent prognostic information
beyond PurIST + DeCAF. The BO selected K=3 over K=2 because this third
axis exists and improves survival prediction."

## Technical Notes

- **K=2 fit**: `desurv_fit()` with k=2, alpha=0.33, lambda=0.35, nu=0.056,
  ntop=270, ninit=50 (same hyperparameters as K=3 except k). Fit time: 2.5 sec.
- **Saved**: `/tmp/desurv_k2_fit.rds`
- **Training data**: Same 273-sample, 1970-gene filtered set used for K=3.
- **BO context**: BO searched K from 2-12 and selected K=3 as optimal.
  Alpha0 (no survival weighting) selected K=7. K=2 was considered and rejected.
- **No external validation**: K=2 results are training-only. The marginal
  adjusted p-values (0.099, 0.028) have not been tested on held-out data.
  Given that K=3-F2's training-only significance (p=0.006) did not replicate
  in validation (p=0.10), the K=2 results should be interpreted cautiously.
- **Gene overlaps**: All computed using top-270 genes per factor (matching K=3).
  Zero genes shared between K2-F1 and K2-F2 (non-overlapping by construction).
- **DECODER compartments**: From Peng et al. deconvolution -- BasalTumor (394
  genes), ClassicalTumor (372), ActivatedStroma (206), NormalStroma (97),
  Immune (408).
- **Published signatures**: DeCAF (Rashid lab, 18 genes), PurIST (16 genes),
  Elyada CAF (iCAF 35, myCAF 16), SCISSORS CAF (iCAF 25, myCAF 25).
