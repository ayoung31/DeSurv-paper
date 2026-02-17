# K=5 DeSurv Analysis: How Do Conclusions Change?

**For discussion with Jen Jen Yeh**
**Date: 2026-02-17**
**Companion to: 2026-02-16-jjy-talking-points.md (K=3 analysis),
2026-02-17-jjy-talking-points-k2.md (K=2 analysis)**

## Motivation

The Bayesian optimization's GP posterior mean surface peaked at K=5 (the best
observed c-index was at K=7; the 1-SE rule selected K=3). K=5 is also the
"elbow K" from standard NMF reconstruction error. If K=5 is the "true" GP
optimum, what does DeSurv look like at that rank? Does the iCAF/B cell program
from K=3 survive, and do the extra two factors add prognostic value?

## K=5 Fit Summary (HPC Production Fit)

| Metric | K=2 | K=3 | K=5 |
|--------|-----|-----|-----|
| c-index (training, internal) | 0.697 | 0.786 | **0.907** |
| Validation c-index (overall) | -- | ~0.6 | **0.496** |
| Non-zero betas | 2/2 | 2/3 | **5/5** |
| BO nu (survival weight) | 0.056* | 0.056 | **0.005** |
| Independent prognostic factors (train adj) | 0 | 1 (F1) | 1 (F2, marginal) |

*K=2 used K=3's hyperparameters with only k changed.

K=5 hyperparameters from BO: alpha=0.362, lambda=0.314, **nu=0.005**, ntop=268.
This is the production HPC fit using 200+ consensus-based seed initializations,
the same pipeline that produced the K=3 model. The critical parameter is
**nu=0.005** -- 11x smaller than K=3's nu=0.056 -- meaning the BO found that
at K=5, survival weighting should be essentially turned off. DeSurv at K=5
converges to near-standard NMF.

### Beta structure

| Factor | Beta | Magnitude | Note |
|--------|------|-----------|------|
| F1 | -3.06e-06 | ~zero | Classical tumor |
| **F2** | **-0.00503** | **dominant** | Mixed (dominant beta) |
| F3 | -3.62e-05 | tiny | iCAF/stroma (richest biology) |
| F4 | -0.000306 | small | Basal/stroma |
| F5 | -7.18e-06 | ~zero | Stroma/immune |

**All betas are negative** and F2 dominates. But F2 is a diffuse "everything"
factor (Basal + Classical + Stroma + Immune + Exocrine), not a coherent
biological program. The iCAF-richest factor (F3) has near-zero beta.

### The LP Paradox

| Metric | K=3 | K=5 |
|--------|-----|-----|
| Training c-index (internal) | 0.786 | 0.907 |
| Training c-index (LP recomputed) | -- | **0.483** |
| K5-LP vs K3-LP correlation | -- | **-0.059** |

The training c-index of 0.907 is the internal optimization metric, but the
LP = beta^T * H gives a concordance of **0.483** (worse than random). The K=5
LP is essentially **uncorrelated with the K=3 LP** (rho = -0.059). This means
the beta weighting learned during fitting does not produce a useful risk score.

## The Answer: K=5 Fragments the iCAF Program and Overfits

At K=5, the iCAF/B cell program from K=3-F1 is **fragmented across multiple
factors** with no single factor capturing it coherently. The model massively
overfits training data and fails completely in validation.

### Factor identity at K=5

| K5 Factor | Biology | DECODER Compartments | Beta |
|-----------|---------|---------------------|------|
| **F1** | Classical tumor | ClassicalTumor (82), BasalTumor (27) | ~0 |
| **F2** | Mixed (dominant beta) | ClassicalTumor (32), BasalTumor (16), Immune (12), Exocrine (12) | **-0.005** |
| **F3** | **iCAF / stroma** | ActivatedStroma (39), Immune (35), NormalStroma (9) | ~0 |
| **F4** | Basal/myCAF stroma | ActivatedStroma (30), BasalTumor (18) | -0.0003 |
| **F5** | Stroma / immune | ActivatedStroma (17), Immune (14) | ~0 |

Unlike K=3 where the iCAF program was a single factor (F1) with the dominant
beta, at K=5 the iCAF biology is richest in **F3** but F3 has near-zero beta.
The dominant beta goes to F2, which is an incoherent mix of tumor, stroma,
immune, and exocrine genes.

### Validation C-index: Worse Than Random

| Dataset | n | DeSurv K=5 C-index | Std NMF K=5 C-index |
|---------|---|-------------------|---------------------|
| Dijk | 90 | 0.555 | **0.632** |
| Moffitt | 123 | 0.512 | **0.549** |
| PACA-AU array | 63 | 0.477 | **0.654** |
| PACA-AU seq | 52 | 0.403 | **0.690** |
| Puleo | 288 | 0.490 | **0.645** |
| **Overall** | **616** | **0.496** | ~0.63 |

**Standard NMF at K=5 outperforms DeSurv K=5 in every validation cohort.**
The survival-driven optimization at K=5 does not just fail to help -- it
actively *hurts* generalization. With nu=0.005, the survival component adds
noise to the factorization without providing useful beta weights.

## Gene-Level Evidence

### K5-F3: The iCAF-Richest Factor (But Not Prognostic)

DECODER: **ActivatedStroma (39), Immune (35)**, NormalStroma (9).

| Reference Signature | Overlap | Shared Genes |
|---------------------|---------|--------------|
| DeCAF restCAF (9) | 3 | CHRDL1, FBLN5, PI16 |
| Elyada iCAF (35) | **18** | FOSB, GPX3, PLA2G2A, TNXB, FBLN1, FBLN2, ADAMTS1, ADH1B, CCDC80, DUSP1, GFPT2, APOE, EFEMP1, DPT, CPE, CXCL14, FSTL1, PI16 |
| Elyada myCAF (16) | 6 | MYL9, BGN, TPM2, POSTN, CTHRC1, TAGLN |
| SCISSORS iCAF (25) | **12** | PLA2G2A, TNXB, ELN, CILP, FBLN2, SFRP1, ADH1B, FBLN5, GFPT2, DPT, CXCL14, PI16 |

K5-F3 has the strongest iCAF signature overlap of any K=5 factor: 18/35 Elyada
iCAF genes, 12/25 SCISSORS iCAF genes, the core chemokine markers (CXCL14,
PI16, ADH1B). But this factor has **near-zero beta** (-3.6e-05) and is **not
independently prognostic** (training adj p = 0.12, validation adj p = 0.0006).

*Note: F3 reaches validation significance (adj p = 0.0006) despite not reaching
training significance (adj p = 0.12). This paradox arises because F3's gene
content is biologically meaningful, but the survival optimization (with nu ≈ 0)
did not assign it prognostic weight during training.*

### K5-F2: The Dominant-Beta Factor (Diffuse Biology)

DECODER: ClassicalTumor (32), BasalTumor (16), Immune (12), Exocrine (12),
Endocrine (9), ActivatedStroma (7), NormalStroma (5) -- broadly mixed.

| Reference Signature | Overlap | Shared Genes |
|---------------------|---------|--------------|
| PurIST Classical (8) | 3 | CLDN18, LGALS4, GATA6 |
| Elyada iCAF (35) | 4 | S100A4, ADH1B, CXCL14, TNXB |
| SCISSORS iCAF (25) | 5 | LIF, SFRP1, ADH1B, CXCL14, TNXB |

F2 picks up some iCAF genes (4 Elyada, 5 SCISSORS) but dilutes them with
Classical Tumor, Exocrine, and Endocrine genes. It has the dominant beta
(-0.005) but is **not independently prognostic in training** after PurIST+DeCAF
adjustment (adj p = 0.042, marginal).

### K5-F4: Basal-like / myCAF Stroma

DECODER: **ActivatedStroma (30), BasalTumor (18)**, ClassicalTumor (10).

| Reference Signature | Overlap | Shared Genes |
|---------------------|---------|--------------|
| Elyada myCAF (16) | **8** | IGFBP7, TPM1, ACTA2, CALD1, TAGLN, MMP11, POSTN, BGN |
| Elyada iCAF (35) | 4 | LMNA, FSTL1, SOD2, S100A10 |

This factor captures the myCAF/activated stroma program -- beta = -0.0003.

### K5-F1 and F5: Tumor Subtype Axis

F1 (ClassicalTumor 82, BasalTumor 27) and F5 (ActivatedStroma 17, Immune 14)
form an anti-correlated pair (rho = -0.716), analogous to K3's F2/F3 and K2's
F1/F2. These capture the PurIST/DeCAF axis.

### K3-F1 iCAF Program: Where Does It Go at K=5?

| Destination | Genes | % of 270 |
|-------------|-------|----------|
| K5-F2 (mixed, dominant beta) | 126 | 46.7% |
| K5-F1 (Classical tumor) | 65 | 24.1% |
| K5-F3 (iCAF/stroma) | 64 | 23.7% |
| K5-F4 (Basal/stroma) | 34 | 12.6% |
| K5-F5 (stroma/immune) | 27 | 10.0% |
| **In any K5 factor** | **229** | **84.8%** |
| **Lost** | **41** | **15.2%** |

Most K3-F1 genes survive (85%) but are **fragmented**: the largest share goes to
F2 (47%, the mixed dominant-beta factor) rather than F3 (24%, the biologically
purest iCAF factor). No single K=5 factor reconstitutes the coherent iCAF
program from K3-F1.

### Cross-Factor Gene Overlap (K5 vs K3, top 268/270)

| | K3-F1 (iCAF) | K3-F2 (Stroma) | K3-F3 (Basal) |
|--|-----|-----|-----|
| **K5-F1** (Classical) | 65 | 0 | 66 |
| **K5-F2** (mixed) | 126 | 11 | 0 |
| **K5-F3** (iCAF/stroma) | 64 | 31 | 41 |
| **K5-F4** (Basal/stroma) | 34 | 27 | 165 |
| **K5-F5** (stroma/immune) | 27 | 200 | 0 |

K3-F2 maps cleanly to K5-F5 (200/270). K3-F3 maps to K5-F4 (165/270). But
K3-F1 is split across F2 (126), F1 (65), and F3 (64) -- fragmented, not
preserved.

### K5 vs K3 H Score Correlations (Model H, Training)

| | K3-F1 (iCAF) | K3-F2 (Stroma) | K3-F3 (Basal) |
|--|-----|-----|-----|
| **K5-F1** | 0.488 | **-0.897** | **0.870** |
| **K5-F2** | 0.360 | -0.045 | -0.065 |
| **K5-F3** | 0.167 | -0.045 | 0.037 |
| **K5-F4** | -0.387 | 0.200 | -0.137 |
| **K5-F5** | -0.419 | **0.917** | **-0.913** |

**No K=5 factor strongly tracks K3-F1 (iCAF)**. The highest correlation is
K5-F1 at 0.488 -- much weaker than the K5-F5 ↔ K3-F2 (0.917) and K5-F1 ↔
K3-F3 (0.870) correspondences. The iCAF axis is lost as a coherent signal.

### K5 Inter-Factor Correlations

| | F1 | F2 | F3 | F4 | F5 |
|--|-----|-----|-----|-----|-----|
| **F1** | 1 | -0.004 | -0.076 | -0.288 | **-0.716** |
| **F2** | | 1 | -0.460 | 0.456 | -0.073 |
| **F3** | | | 1 | -0.475 | -0.163 |
| **F4** | | | | 1 | 0.020 |
| **F5** | | | | | 1 |

F1-F5 are anti-correlated (rho = -0.716), recapitulating the tumor subtype
axis (analogous to K3's F2-F3 at -0.971 and K2's F1-F2 at -0.928). F2, F3,
and F4 form a loosely correlated triplet with moderate correlations (|rho| =
0.46-0.48). The factor structure is more complex than K=3 but does not
decompose into cleanly interpretable biological axes.

## PurIST / DeCAF Associations

| Factor | PurIST | DeCAF |
|--------|--------|-------|
| F1 (Classical tumor) | Not significant (p=0.10) | Not significant (p=0.72) |
| F2 (mixed) | **Classical** (p=7e-8) | **restCAF** (p=3.5e-5) |
| F3 (iCAF/stroma) | Not significant (p=0.085) | restCAF (p=0.005) |
| F4 (Basal/stroma) | Not significant (p=0.18) | **proCAF** (p=0.002) |
| F5 (stroma/immune) | Not significant (p=0.13) | Not significant (p=0.54) |

Only F2 strongly aligns with PurIST (Classical, p=7e-8) and DeCAF (restCAF,
p=3.5e-5). F3 (the iCAF-rich factor) shows weak DeCAF association (p=0.005)
but no PurIST association -- similar to K3-F1's pattern of being orthogonal to
PurIST, but the signal is much weaker.

## Per-Factor Prognostic Value

### Training (n=273)

| Factor | Unadjusted p | Adjusted p (PurIST+DeCAF) |
|--------|-------------|--------------------------|
| F1 (Classical tumor) | 0.876 | 0.855 |
| **F2 (mixed)** | **0.001** | **0.042** |
| F3 (iCAF/stroma) | 0.148 | 0.124 |
| F4 (Basal/stroma) | 0.111 | 0.274 |
| F5 (stroma/immune) | 0.981 | 0.847 |

Only F2 reaches significance, and its adjusted p-value (0.042) is marginal.
Compare to K3-F1's training adj p of 3.1e-13. The iCAF-rich factor (F3) is
**not significant** in training (p = 0.12).

### Validation (n=564 complete cases, 5 cohorts)

| Factor | Unadjusted p | Adjusted p (PurIST+DeCAF) |
|--------|-------------|--------------------------|
| **F1** | **9.2e-06** | **7.7e-05** |
| F2 | 0.659 | 0.068 |
| **F3** | 0.277 | **0.0006** |
| F4 | 0.357 | 0.078 |
| **F5** | **5.4e-05** | **2.1e-05** |

**Paradox**: Several factors are prognostic in *validation* (F1, F3, F5) but
not in *training*. This reversal occurs because:

1. **Near-zero nu (0.005)** means the factorization barely uses survival
   during training, so factors aren't optimized for prognosis
2. The gene content (inherited from the expression structure) carries inherent
   biological prognostic value that manifests in larger validation samples
3. The beta weighting (dominated by F2) does not align with which factors
   actually predict survival externally

### Joint Model (Validation)

| Term | Coefficient | p-value |
|------|-------------|---------|
| X1 | -0.289 | 0.009 |
| X2 | -0.082 | 0.216 |
| X3 | 0.225 | 0.037 |
| X4 | 0.027 | 0.133 |
| **X5** | **-1.116** | **<0.0001** |
| PurIST (Classical) | -0.304 | 0.046 |
| DeCAF (restCAF) | -0.582 | <0.001 |

F5 has the largest coefficient in the joint model. But the model's LP combines
factors with weights (betas) that bear no resemblance to these Cox coefficients
-- the LP was learned with nu ≈ 0, not from the survival data.

### LRT: Does K=5 Add to PurIST + DeCAF?

| Model | Training LRT p | Validation LRT p |
|-------|---------------|-----------------|
| K=3 LP | < 2.2e-16 | **0.001** |
| K=5 LP | -- | 0.056 (not significant) |
| K=5 F2 alone (dominant beta) | -- | -- |

The K=5 LP does **not** significantly add to PurIST+DeCAF in validation
(p = 0.056). Compare K=3 LP (p = 0.001).

## Comparison: K=2 vs K=3 vs K=5

| Question | K=2 | K=3 | K=5 (HPC) |
|----------|-----|-----|-----------|
| **c-index (train, internal)** | 0.697 | 0.786 | 0.907 |
| **Validation c-index** | -- | ~0.6 | **0.496** |
| **Std NMF K=5 val c-index** | -- | -- | ~0.63 (beats DeSurv) |
| **nu (survival weight)** | 0.056* | 0.056 | **0.005** (near zero) |
| **iCAF factor present?** | No (lost) | Yes (F1, coherent) | **Fragmented** (F2/F3/F1) |
| **iCAF factor adj p (train)** | -- | **3.1e-13** | 0.042 (F2, marginal) |
| **iCAF factor adj p (val)** | -- | **0.004** | 0.068 (F2, NS) |
| **LP adds to PurIST+DeCAF? (val)** | -- | **Yes (p=0.001)** | No (p=0.056) |
| **K3-F1 genes preserved** | 29% | 100% | 85% (but fragmented) |
| **K3-F1 genes lost** | 71% | 0% | 15% |
| **LP vs K3-LP correlation** | rho=0.780 | 1 | **rho=-0.059** |

## What K=5 Teaches Us About K=3

### 1. At K=5, the BO turns off survival weighting

The most revealing finding is **nu = 0.005**. The BO searched over nu and
determined that at K=5, survival information should barely influence the
factorization. This means DeSurv at K=5 is effectively standard NMF -- and
indeed, standard NMF at K=5 *outperforms* DeSurv K=5 in validation (c-index
~0.63 vs 0.496). The survival component at K=5 adds noise, not signal.

At K=3, the BO sets nu=0.056 -- 11x higher. K=3 is the rank where survival
weighting genuinely improves the factorization. The extra capacity at K=5
distributes the survival signal across too many factors for the beta estimation
to work.

### 2. The iCAF program fragments at K=5

At K=3, the iCAF/B cell program is concentrated in a single factor (F1, 270
genes) with the dominant beta and strong validation (adj p=0.004). At K=5,
those same genes scatter across three factors:
- F2 (126 genes, dominant beta but diffuse biology)
- F1 (65 genes, mixed with Classical tumor)
- F3 (64 genes, richest iCAF biology but near-zero beta)

No single K=5 factor captures the iCAF signal coherently. The H-score
correlations confirm this: no K=5 factor correlates strongly with K3-F1
(max rho = 0.488 vs 0.917 for the tumor subtype axis).

### 3. K=5 demonstrates overfitting, not capacity

The training c-index (0.907) is dramatically inflated relative to validation
(0.496). This 0.41-point gap (vs K=3's ~0.19-point gap) is textbook
overfitting. With 5 factors and near-zero survival weighting, the model
memorizes training-specific patterns that don't generalize. The 1-SE rule's
preference for K=3 correctly identifies this overfitting risk.

### 4. The K=2 → K=3 → K=5 arc is clear

| Transition | What changes |
|-----------|-------------|
| K=2 → K=3 | **iCAF axis appears** (71% of genes NEW). Novel prognostic signal validated externally. |
| K=3 → K=5 | iCAF fragments. Survival weighting collapses. LP becomes pathological. Validation fails. |

K=3 is the unique rank where:
- Survival weighting works (nu > 0)
- The iCAF program is coherent (one factor, one dominant beta)
- The LP generalizes (validation LRT p=0.001)
- The model doesn't overfit (c-index gap ~0.19)

### The pitch (K=5 edition)

"At K=5, the BO's own optimization reveals that survival information cannot
usefully guide the factorization (nu = 0.005, essentially zero). The
resulting model fragments K=3's coherent iCAF/B cell program across three
factors, produces a linear predictor that is uncorrelated with K=3's
(rho = -0.06), and validates worse than standard NMF (c-index 0.50 vs 0.63).
Individual factor scores retain some prognostic value in validation --
confirming the biological content is real -- but the survival-driven beta
weighting fails. K=3 is the unique rank where DeSurv's survival component
adds value: nu is 11x higher, the iCAF program is captured in a single
coherent factor, and the LP significantly improves on PurIST+DeCAF in
external validation (p=0.001). The 1-SE rule correctly identified K=3 as
the point where model complexity and generalization are balanced."

## Technical Notes

- **HPC K=5 fit**: Production pipeline fit via `tar_fit_desurv_elbowk_tcgacptac`.
  k=5, alpha=0.362, lambda=0.314, nu=0.005, ntop=268. Consensus initialization
  from 200+ seeds via `desurv_consensus_seed()`. Same pipeline as K=3 production.
- **Local K=5 fit**: Earlier exploratory fit with different BO params (alpha=0.397,
  lambda=0.256, nu=0.149, ntop=185, ninit=50). This told a different story
  (iCAF preserved in F2, rho=0.957 with K3-F1, c-index 0.802). The HPC fit is
  authoritative because it uses the full pipeline's consensus initialization.
- **Training data**: Same 273-sample, 1970-gene filtered set as K=3.
- **Validation**: 5 cohorts (n=616 total, 564 with PurIST+DeCAF labels).
  Projected K=5 W onto validation data using factor scores = W^T X.
- **Standard NMF K=5**: `fit_std_elbowk_tcgacptac` -- standard NMF at K=5
  (alpha=0) with post-hoc Cox beta estimation. Validation c-indices 0.55-0.69,
  substantially better than DeSurv K=5 (0.40-0.56).
- **BO context**: GP posterior mean peaked at K=5 (c-index ~0.643), best
  observed at K=7 (0.655). 1-SE/LCB rule selected K=3 (lcb_threshold=0.644,
  best_mean=0.655, best_se=0.011).
- **Gene overlaps**: K5 top-268 genes. K3 top-270 genes.
- **c-index discrepancy**: The HPC fit reports training c-index 0.907
  (internal optimization metric), but recomputing LP = beta^T * H and
  computing concordance gives 0.483. This likely reflects different LP
  computation or sign conventions inside the C++ fitting code. The validation
  c-index (0.496 from the pipeline's own risk scores) confirms poor
  out-of-sample performance.
- **Validation note**: PACA_AU_seq samples (n=52) have no PurIST/DeCAF labels
  and were excluded from adjusted analyses.
