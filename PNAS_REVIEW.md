# PNAS Publication Review: DeSurv Manuscript

**Review Date:** 2026-01-24 (Updated with editorial strategy analysis)
**Manuscript:** "Survival driven deconvolution (DeSurv) reveals prognostic and interpretable cancer subtypes"
**Target Journal:** Proceedings of the National Academy of Sciences (PNAS)
**Review Type:** Technical compliance + Editorial strategy + Likely reviewer concerns

---

## Executive Summary

The DeSurv manuscript presents a conceptually meaningful framework that integrates nonnegative matrix factorization with Cox proportional hazards modeling for survival-supervised deconvolution of tumor transcriptomes. The key innovation—**supervising W (gene programs) rather than H (sample loadings)**—represents a genuine conceptual reframing that enables interpretable, transferable gene programs.

**Recommendation: Minor-to-Moderate Revision** (not Major)

The core framework is sound, the conceptual contribution is clear, and the results support the claims. Required revisions:

1. **Critical code bugs** must be fixed and analyses re-run
2. **PNAS formatting** (reduce figures 6→4, references 71→50)
3. **Expand underdeveloped sections** (methods, discussion)
4. **Add one harder simulation scenario** (beyond "easy")
5. **Add one comparator method** (LASSO-Cox) with explanation
6. **Complete required sections** (keywords, author contributions)

**What to resist in revision:** Exhaustive benchmarks, deep biological validation, removing bladder transfer result

---

## Overall Editorial Assessment (PNAS Fit)

**Recommendation:** *Minor-to-Moderate Revision—suitable for PNAS with targeted changes.*

### Why It Fits PNAS (Recalibrated)

PNAS judges novelty by **conceptual reframing**, not benchmark performance tables. This paper offers a genuine conceptual advance:

1. **W-level vs H-level supervision distinction is conceptually meaningful.** Existing survival-NMF methods (Le Goff, Huang) supervise through H (sample loadings), limiting utility for subtyping/interpretation. DeSurv supervises through W (gene programs), enabling interpretable programs that can be projected to new data.

2. **The core contribution is structural/methodological:** "Survival supervision at the gene-program level reorganizes latent structure in ways that improve interpretability, stability, and generalization." This is distinct from a prediction benchmark claim or pathway discovery claim.

3. **Cross-cohort generalization** (PDAC → bladder) demonstrates that gene-level supervision enables projection-based transfer—a methodological point, not a universal pan-cancer biology claim.

### Reviewer Risks (Manageable)

| Risk | Concern | Mitigation |
|------|---------|------------|
| **Novelty perception** | Could be misconstrued as incremental | Sharpen W-vs-H distinction earlier in Introduction |
| **"Double-dipping"** | Using survival to select k | Add explicit out-of-sample statement |
| **Modest C-index gains** | ~0.02-0.05 may seem incremental | Lead with structural reorganization, not C-index |
| **Cross-cancer overstated** | Bladder result may seem bold | Minor tempering ("suggests" vs "demonstrates") |
| **Scope creep requests** | Reviewers may ask for exhaustive benchmarks | Use rebuttal language emphasizing conceptual contribution |

---

## Table of Contents

1. [PNAS Submission Compliance](#1-pnas-submission-compliance)
2. [Critical Code Bugs](#2-critical-code-bugs-affecting-paper-claims)
3. [Paper-Code Inconsistencies](#3-paper-code-inconsistencies)
4. [Mathematical Notation Issues](#4-mathematical-notation-issues)
5. [Scientific Content Assessment](#5-scientific-content-assessment)
6. [Narrative Arc Analysis](#6-narrative-arc-analysis)
7. [Editorial Strategy & Reviewer Concerns](#7-editorial-strategy--likely-reviewer-concerns)
8. [Supplement Completeness](#8-supplement-completeness)
9. [Action Checklist](#9-action-checklist)

**Appendices:**
- [Appendix A: Files Reviewed](#appendix-files-reviewed)
- [Appendix B: Editorial Review Analysis](#appendix-b-editorial-review-analysis)
- [Appendix C: Figure Review Analysis](#appendix-c-figure-review-analysis)
- [Appendix D: Recalibrated PNAS Assessment](#appendix-d-recalibrated-pnas-assessment)
- [Appendix E: Panel-by-Panel Figure Analysis](#appendix-e-panel-by-panel-figure-analysis) *(detailed in [FIGURE_ANALYSIS.md](FIGURE_ANALYSIS.md))*

---

## 1. PNAS Submission Compliance

### 1.1 Quantitative Requirements

| Requirement | PNAS Limit | Current Status | Action |
|-------------|------------|----------------|--------|
| Word count | ~4000 words (6 pages) | ~3000-3500 words | ✅ Within limit |
| Maximum pages | 12 pages | Appears reasonable | ✅ OK |
| Figure count | 4 figures (standard) | **6 figures** | ❌ Move 2 to SI |
| Reference count | 50 references | **71 references** | ❌ Reduce by 21+ |
| Abstract | 250 words max | ~120 words | ✅ OK |
| Significance statement | 120 words max | ~110 words | ✅ OK |

### 1.2 Current Figures (6 total - must reduce to 4)

1. `fig-schema` - Model overview schematic
2. `fig-bo` - Bayesian optimization and rank selection
3. `fig-sim-ntop-scenarios` - Simulation results
4. `fig-bio` - Biological interpretation (PDAC)
5. `fig-extval` - External validation
6. `fig-bladder` - Bladder cancer analysis

**Recommendation:** Move `fig-bo` (rank selection details) and `fig-bladder` (cross-cancer analysis) to Supplementary Information.

### 1.3 Required Sections Status

| Section | Status | Current Content | Action Required |
|---------|--------|-----------------|-----------------|
| Title | ✅ Complete | "Survival driven deconvolution (DeSurv)..." | None |
| Abstract | ✅ Complete | ~120 words describing the method | None |
| Significance | ✅ Complete | ~110 words on broader impact | None |
| Keywords | ❌ **Placeholder** | "one, two, optional, optional, optional" | Replace with scientific terms |
| Author contributions | ❌ **Placeholder** | "Please provide details..." | Must complete |
| Conflict of interest | ❌ **Placeholder** | "Please declare any conflict..." | Must complete |
| Acknowledgements | ❌ **Placeholder** | Generic placeholder text | Must complete |

### 1.4 Author Information

| Element | Status | Details |
|---------|--------|---------|
| Authors listed | ✅ | Amber M. Young, Alisa Yurovsky, Didong Li, Naim U. Rashid |
| Affiliations | ❌ **Incomplete** | Contains "Street, City, State, Zip" placeholders |
| Corresponding author | ✅ | ayoung31@live.unc.edu |
| ORCID IDs | ⚠️ Missing | Not specified in YAML header (recommended) |

### 1.5 Manuscript Type

- **Specified:** `pnas_type: pnasresearcharticle` (standard two-column layout) ✅

---

## 2. Critical Code Bugs Affecting Paper Claims

These bugs may invalidate the reported results and **must be fixed before submission**. After fixing, all analyses should be re-run to verify claims still hold.

### 2.1 Bug #1: Variable Assignment Error (CRITICAL)

**File:** `R/compute_metrics.R:9`

**The Bug:**
```r
# Lines 8-9
meanZ = fit$meanZ[,ns]
sdZ = fit$meanZ[,ns]    # BUG: Should be fit$sdZ[,ns]
```

**Impact:**
- Line 12: `XtW = sweep(XtW, 2, sdZ, FUN="/")` divides by the **mean** instead of standard deviation
- Line 13: `lp <- XtW %*% (beta * sdZ)` multiplies by **mean** instead of standard deviation
- **All C-index calculations and survival metrics may be incorrect**
- Paper claims about "higher C-index" and "improved precision" could be based on flawed metrics

**Fix:**
```r
sdZ = fit$sdZ[,ns]
```

### 2.2 Bug #2: Unused Filter Variable (CRITICAL)

**File:** `R/select_best_init.R:5`

**The Bug:**
```r
# Line 5 creates filtered dataframe
keep <- df[!df$flag_nan & is.finite(df$loss), ]

# BUT lines 11-15 still use original `df` instead of `keep`
if(method_select=="surv"){
  best = df %>% dplyr::slice_max(order_by = sloss, n = 1, with_ties = FALSE)
}
```

**Impact:**
- Invalid initializations (NaN or infinite loss) are NOT filtered out
- An initialization with NaN/Inf loss could be selected as "best"
- Paper's claims about stable model selection may be compromised

**Fix:**
Replace all uses of `df` after line 5 with `keep`:
```r
select_best_init <- function(df, method_select="surv") {
  keep <- df[!df$flag_nan & is.finite(df$loss), ]

  if (nrow(keep) == 0) {
    stop("No valid initializations found (all had NaN or infinite loss)")
  }

  keep$method_select = method_select

  if(method_select=="surv"){
    best = keep %>% dplyr::slice_max(order_by = sloss, n = 1, with_ties = FALSE)
  } else if(method_select=="nmf"){
    best = keep %>% dplyr::slice_min(order_by = nloss, n = 1, with_ties = FALSE)
  } else if(method_select=="loss"){
    best = keep %>% dplyr::slice_min(order_by = loss, n = 1, with_ties = FALSE)
  } else {
    stop("This method for selecting the best initialization is not supported")
  }

  return(best)
}
```

### 2.3 Bug #3: Beta Backtracking Logic (CRITICAL - NEEDS VERIFICATION)

**File:** `../DeSurv/src/functions.cpp:472-473`

**The Code:**
```cpp
if (!std::isfinite(surv_loss_candidate) ||
    surv_loss_candidate < surv_loss_prev) {
  // Enter backtracking...
```

**The Issue:**
- `surv_loss` is defined as `cs.loglik * 2.0 / n_event - lambda * penalty_beta`
- Since `loglik` is a Cox partial log-likelihood, **higher values are better** (maximize)
- The condition `surv_loss_candidate < surv_loss_prev` triggers backtracking when the candidate is **lower**
- This means backtracking is triggered when the fit **improves** (if we're maximizing)

**Impact:**
- Unnecessary backtracking when optimization makes progress
- Potentially accepting worse solutions
- Slower convergence
- Convergence proof (Theorem 1) may not hold in practice

**Analysis Required:**
Verify whether `surv_loss` is being maximized or minimized. If maximized, change `<` to `>`.

### 2.4 Bug #4: Debug Statements in Production Code

**File:** `targets_common_pipeline.R`
**Lines:** 1691, 1803, 2075, 2097, 2119, 2141

Six `browser()` statements will halt the pipeline in non-interactive (Slurm) environments.

**Fix:** Remove all `browser()` calls.

---

## 3. Paper-Code Inconsistencies

### 3.1 Convergence Criterion Mismatch

| Aspect | Paper (supp_methods.Rmd) | Code (functions.cpp) |
|--------|--------------------------|----------------------|
| **Criterion** | Relative loss change: `eps = \|lossNew - loss\| / \|loss\|` | Cosine similarity: `eps = max(1-cos(W_prev,W), 1-cos(H_prev,H))` |
| **Location** | Algorithm S1, line 99 | Lines 649-660 |

**Code Reality:**
```cpp
double cosW = (wnorm > 0.0) ? arma::dot(arma::vectorise(Wprev),
                                       arma::vectorise(W)) / wnorm : 1.0;
double cosH = (hnorm > 0.0) ? arma::dot(arma::vectorise(Hprev),
                                       arma::vectorise(H)) / hnorm : 1.0;
epsW = 1.0 - cosW;
epsH = 1.0 - cosH;
eps = std::max(epsW, epsH);
```

**Action:** Either update paper to describe cosine-similarity criterion, or change code to match paper.

### 3.2 W Update Ratio Clamping (Undocumented)

| Aspect | Paper | Code |
|--------|-------|------|
| **W update** | `W^{(t+1)} = max(W^{(t)} ⊙ R^{(t)}, ε_W)` | R clamped to [0.2, 5.0] before application |
| **Location** | supp_methods.Rmd Eq. \ref{eqn:W_update} | functions.cpp lines 209-211 |

**Code Reality:**
```cpp
arma::mat R = num / (denom + eps);
R = arma::clamp(R, 0.2, 5.0);  // UNDOCUMENTED
```

**Action:** Document the clamping in the supplement or remove it from code.

### 3.3 Cox Gradient Scaling Formula Differs

| Aspect | Paper (supp_methods.Rmd Eq. S6) | Code (functions.cpp) |
|--------|--------------------------------|----------------------|
| **Formula** | `δ^{(t)}` as ratio of **squared** Frobenius norms | Uses **unsquared** norms |
| **Clipping** | `δ^{(t)}` clipped to `δ_max = 10^6` | `alpha * (gn/gc)` clipped to `alpha * 1e6` |

**Code Reality:**
```cpp
double gn = arma::norm(grad_nmf, "fro") + 1e-12;
double gc = arma::norm(dW_cox, "fro") + 1e-12;
double cox_scale = alpha * std::min(gn / gc, 1e6);
```

**Action:** Align paper formula with implementation.

### 3.4 H Penalty Scaling (Undocumented)

| Aspect | Paper | Code |
|--------|-------|------|
| **H update denominator** | `W^T W H + λ_H H + ε_H` | `(W^T W H) + ((lambdaH*Xnorm / (n*k)) * H)` |

The implementation scales the H penalty by `Xnorm / (n*k)`, which is not shown in the paper.

**Action:** Document the scaling factor in the supplement.

### 3.5 Missing Function Definitions

**File:** `R/bo_helpers.R:62,66`

Functions called but do not exist:
- `collect_bo_diagnostics()`
- `compute_bo_eval_se()`

The `tryCatch` silently masks these failures. The "one standard error rule" for k-selection may not be properly implemented.

### 3.6 Linear Predictor Computation (Undocumented)

**Paper Claims:**
- Methods: "linear predictor Z^⊤ β"
- Results: "factor scores Z = W^⊤ X... linear predictor Z^⊤ β"

**Code Reality (compute_metrics.R lines 11-13):**
```r
XtW = sweep(XtW, 2, meanZ, FUN="-")    # Centering
XtW = sweep(XtW, 2, sdZ, FUN="/")       # Scaling
lp <- XtW %*% (beta * sdZ)              # Rescaling beta
```

The centering and scaling procedure for validation metrics is not documented in the paper.

---

## 4. Mathematical Notation Issues

### 4.1 High Priority Issues

#### Symbol Overloading: δ (delta)

| Usage | Location | Definition |
|-------|----------|------------|
| Event indicators | Main methods line 2, supp_methods line 34 | `δ ∈ {0,1}^n` |
| Gradient-balancing factor | supp_methods lines 165-176 | `δ^{(t)}` for W update |

**Problem:** Same symbol for fundamentally different quantities creates confusion.

**Fix:** Rename gradient-balancing factor to `γ^{(t)}` or `ρ^{(t)}`.

#### Symbol Overloading: θ (theta)

| Usage | Location | Definition |
|-------|----------|------------|
| Model parameters | supp_methods line 354 | `θ = (W, H, β)` |
| Hyperparameters | supp_methods line 735 | `θ` in Bayesian optimization |

**Fix:** Rename hyperparameter tuple to `Θ` or `ψ`.

#### Incomplete Equation

**Location:** supp_methods.Rmd line 306

**Current:**
```latex
\lambda(\xi||\beta||_1 + \frac{(1-\xi)}{2}
```

**Should be:**
```latex
\lambda(\xi||\beta||_1 + \frac{(1-\xi)}{2}||\beta||_2^2)
```

#### Typo in Gradient Formula

**Location:** supp_methods.Rmd line 646

**Current:**
```latex
\nabla \mathcal{L} = (1-\alpha) \nabla_W \mathcal{L}_{\mathrm{Cox}}(W,H) - \alpha\nabla_W \mathcal{L}_{\mathrm{Cox}}(W,\beta)
```

**Should be:**
```latex
\nabla \mathcal{L} = (1-\alpha) \nabla_W \mathcal{L}_{\mathrm{NMF}}(W,H) - \alpha\nabla_W \mathcal{L}_{\mathrm{Cox}}(W,\beta)
```

### 4.2 Medium Priority Issues

| Issue | Locations | Fix |
|-------|-----------|-----|
| Inconsistent `n_{event}` formatting | Lines 37, 64, 142, 184, 250, 306 | Standardize to `n_{\text{event}}` |
| "Event indicators" vs "censoring indicators" | Main vs supplement | Choose one term |
| Event indicator type | Main: `δ ∈ ℝ^n`, Supp: `δ ∈ {0,1}^n` | Should be binary `{0,1}^n` |
| Nonnegative notation order | Main: `ℝ^{p×n}_{≥0}`, Supp: `ℝ_{≥0}^{p×n}` | Choose one style |
| Inconsistent norm notation | `\|\|` vs `\lVert` | Use `\lVert`/`\rVert` |

### 4.3 Undefined Symbols

The following symbols are used before formal definition:
- `ε_H` and `ε_W` (floor values) - described contextually but not formally defined
- `η̃_i` vs `η_i` - different definitions in different sections

**Recommendation:** Add a notation table at the beginning of the supplement.

---

## 5. Scientific Content Assessment

### 5.1 Strengths

1. **Clear innovation:** Supervision through W (gene signatures) rather than H (sample loadings) is a meaningful distinction that is stated but needs earlier/stronger emphasis
2. **Rigorous theory:** Complete convergence proof with formal theorem and supporting lemmas - methodological depth beyond typical applied genomics
3. **Logical structure:** Results progress naturally from rank selection → simulation → PDAC → validation → bladder
4. **Compelling visualization:** Variance explained vs. survival contribution (Fig. 4C) effectively communicates the key insight - this is the real story and should be emphasized more
5. **Cross-cancer transferability:** PDAC-trained factors retaining signal in bladder cancer is potentially significant and addresses PNAS's preference for generalizable insight
6. **Nested CV design:** Already uses proper out-of-sample validation, just needs explicit statement
7. **Projection advantage:** Closed-form Z = W^⊤X enables stable cross-cohort scoring - needs to be in main text

### 5.2 Weaknesses Requiring Attention

#### Main Methods Section (03_methods.Rmd)

**Problem:** Only ~30 lines - too terse for PNAS which expects self-contained methods.

**Missing details:**
- Specific form of multiplicative updates
- Convergence criteria
- Cox partial likelihood formula
- Sample sizes for simulations
- Bayesian optimization parameters (initial points, acquisition function, total evaluations)

**Action:** Expand to ~1.5 pages with key algorithmic details.

#### Discussion Section (05_discussion.Rmd)

**Problem:** Only 6 lines - completely inadequate for PNAS.

**Missing elements:**
- Limitations discussion (Cox PH assumption, computational cost, small sample sizes)
- Comparison with related approaches (WGCNA, integrative clustering)
- Clinical translation pathway
- Biological interpretation of cross-cancer transfer
- Computational requirements and scalability

**Action:** Expand to 2-3 paragraphs.

#### Simulation Results

**Problem:** Only "R0_easy" scenario shown.

**Missing:**
- Null scenario results (referenced as "R00_null" in code)
- Multiple signal-to-noise ratios
- Different sample sizes
- Various censoring rates

**Action:** Add comprehensive simulation table across parameter ranges.

#### Quantitative Metrics

**Problem:** Underreported throughout.

**Missing:**
- Actual C-index improvement values (DeSurv vs NMF+Cox)
- Hazard ratios with 95% confidence intervals
- Log-rank p-values for Kaplan-Meier comparisons

**Action:** Add numerical values to complement visual presentations.

#### Bladder Cancer Analysis

**Problem:** Underdeveloped compared to PDAC.

**Missing:**
- Comparison with standard NMF on bladder data
- Biological interpretation of transferred factors

**Action:** Either develop to same depth as PDAC or move to supplement.

#### Method Comparisons

**Problem:** Only compared to standard NMF.

**Missing comparisons:**
- Sparse NMF
- PLSR-Cox
- Penalized Cox with all genes
- Prior survival-NMF methods (le2025survnmf, huang2020low)

**Action:** Add at least one additional competing method.

### 5.3 PNAS Suitability Assessment

**Current positioning:** Strong methodological paper with genuine conceptual advance. With targeted editorial revisions, suitable for PNAS.

**What makes it PNAS-worthy:**
1. The W-vs-H supervision distinction represents a conceptual break, not incremental improvement
2. Cross-cohort generalization (PDAC → bladder) addresses PNAS's preference for generalizable insight
3. Variance vs survival realignment (Fig. 4) provides mechanistic interpretability
4. Methodological depth (convergence analysis, BO-based selection) exceeds typical applied genomics

**Key reframing needed:**
1. Lead with **structural reorganization insight**, not predictive gains
2. Emphasize the **interpretability advantage** of W-level supervision
3. Make the **projection argument** (closed-form scoring) explicit in main text
4. Pre-empt "double-dipping" concerns with explicit out-of-sample statements

**Optional enhancements (not required but would strengthen):**
1. Single-cell validation showing discovered programs correspond to real cellular states
2. Testing whether DeSurv signatures predict treatment response (not just prognosis)
3. Deeper biological interpretation of why cross-cancer transfer works

**Most PNAS-suitable aspect:** Cross-cancer transferability of prognostic structure combined with the mechanistic insight that variance-dominant programs ≠ survival-relevant programs.

---

## 6. Narrative Arc Analysis

### 6.1 The Five-Act Story Structure

The DeSurv manuscript employs a **problem-solution-validation** narrative structure with escalating generalization claims:

```
┌─────────────────────────────────────────────────────────────────────┐
│  ACT 1: THE PROBLEM (Introduction)                                  │
│  "Discovery/validation separation misses clinically relevant        │
│   programs; existing survival-NMF supervises H, not W"              │
└──────────────────────────────┬──────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│  ACT 2: THE SOLUTION (Methods + Fig 1)                              │
│  "DeSurv supervises W directly, ensuring programs are               │
│   intrinsically prognostic"                                         │
└──────────────────────────────┬──────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│  ACT 3: PROOF OF CONCEPT (Figs 2-3)                                 │
│  "Standard rank selection fails → DeSurv works"                     │
│  "Simulation validates gene recovery"                               │
└──────────────────────────────┬──────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│  ACT 4: THE INSIGHT (Figs 4-5) ★ TURNING POINT ★                    │
│  "Variance ≠ Prognosis: DeSurv reorganizes programs toward          │
│   survival, deprioritizing variance-dominant noise"                 │
└──────────────────────────────┬──────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│  ACT 5: THE GENERALIZATION (Fig 6) ★ PNAS CLAIM ★                   │
│  "PDAC programs transfer to bladder cancer →                        │
│   This is a general principle, not just PDAC"                       │
└─────────────────────────────────────────────────────────────────────┘
```

### 6.2 Introduction Narrative Beats

The introduction follows a classical five-paragraph problem-solution structure:

| Paragraph | Narrative Function | Key Tension |
|-----------|-------------------|-------------|
| 1 | Establish importance | Molecular subtyping matters for precision oncology |
| 2 | Identify the gap | Discovery/validation separation risks overfitting |
| 3 | Explain why it's hard | Technical barriers (cohort size, reference limitations) |
| 4 | Critique prior art | Existing survival-NMF supervises H, not W |
| 5-6 | Present solution | DeSurv supervises W directly |

**The Core Innovation Claim (Paragraph 4-5):**
> "Existing survival-aware NMF methods integrate the survival objective through the sample-specific loadings (H) rather than the gene-level programs (W). DeSurv integrates survival information directly into the gene signature matrix."

This is the **conceptual hook** the entire paper builds around.

### 6.3 Results Section: The Evidence Cascade

The results follow a strict logical progression:

| Section | Narrative Beat | Figure | Key Claim |
|---------|----------------|--------|-----------|
| Model Overview | Method introduction | Fig 1 | W is shared between NMF and Cox; BO for selection |
| Rank Selection | Problem + solution | Fig 2 | Standard heuristics fail; DeSurv resolves ambiguity |
| Simulation | Controlled validation | Fig 3 | Supervision improves C-index and gene recovery |
| Biological Structure | **TURNING POINT** | Fig 4 | Variance ≠ prognosis; DeSurv reorganizes programs |
| External Validation | Reproducibility | Fig 5 | Structure generalizes across PDAC cohorts |
| Cross-cancer Transfer | **PNAS CLAIM** | Fig 6 | Programs transfer to bladder cancer |

### 6.4 Figure-Narrative Role Mapping

| Figure | Story Function | What Reader Should Take Away |
|--------|----------------|------------------------------|
| **Fig 1** | "Here's what I'm proposing" | DeSurv couples NMF with Cox by sharing W; uses BO for selection |
| **Fig 2** | "Standard methods fail; DeSurv succeeds" | Unsupervised heuristics disagree; survival-guided selection works |
| **Fig 3** | "Under known ground truth, it works" | Supervision improves both prediction AND program recovery |
| **Fig 4** | "Here's the insight that makes this matter" | DeSurv **restructures programs** to isolate survival-relevant biology |
| **Fig 5** | "This isn't overfitting; it generalizes" | Survival-aligned programs generalize to independent datasets |
| **Fig 6** | "This is a general principle, not just PDAC" | Survival structure is shared across cancer types |

### 6.5 The Central Insight: "Variance ≠ Prognosis"

The entire paper builds toward one core insight, visualized in **Figure 4C** and echoed in **Figure 6A**:

> **High-variance transcriptional axes are not necessarily prognostically relevant.** Standard NMF allocates model capacity to variance-dominant signals (like exocrine expression in PDAC) that contribute little to survival. DeSurv redirects this capacity toward survival-aligned programs.

**This insight is what makes the paper PNAS-worthy** - not just "slightly better C-index," but a mechanistic understanding of how transcriptional structure relates to clinical outcomes.

### 6.6 Rhetorical Strategies

1. **Problem-Solution Framing:** Every results section begins by identifying a limitation of standard NMF, then shows DeSurv addresses it

2. **Escalating Generalization Claims:**
   ```
   Simulation (controlled) → PDAC training → PDAC validation → Bladder (different cancer)
   ```

3. **Consistent Baseline Comparison:** DeSurv vs α=0 (standard NMF) runs throughout, providing a consistent ablation

4. **Biological Grounding:** Every claim is tied to known biology (classical/basal-like, CAF subtypes, immune infiltration)

### 6.7 Information Disclosure Hierarchy

| Document | Role | Audience | Detail Level |
|----------|------|----------|--------------|
| Main Methods | Framework + key choices | Everyone | Conceptual |
| Results | Empirical validation | Everyone | Interpretive |
| Supplement | Complete derivations | Methods experts | Technical |
| Discussion | Synthesis + implications | Everyone | Strategic |

### 6.8 Narrative Gaps Identified

| Gap | Current State | Impact |
|-----|---------------|--------|
| W vs H distinction | Buried in Intro para 4 | Novelty may be missed |
| Projection advantage | Mentioned but not highlighted | Generalization mechanism unclear |
| Why cross-cancer works | Biological mechanism implicit | Claim may seem arbitrary |
| Factor biological labels | Fig 4 describes but doesn't label | Interpretation burden on reader |
| Effect sizes on KM | Forest plot has HRs but KM lacks p-values | Statistical rigor questioned |

---

## 7. Editorial Strategy & Likely Reviewer Concerns

This section addresses substantive editorial issues that are likely to trigger reviewer pushback and provides actionable revision guidance.

### 7.0 Reviewer Archetype Framework

When evaluating revision strategies, consider which reviewer types each change helps or might annoy:

| Archetype | Focus | What They Want | What Annoys Them |
|-----------|-------|----------------|------------------|
| **R1: Statistical purist** | Calibration, leakage, rigor | Effect sizes, explicit out-of-sample statements, uncertainty quantification | Sloppy p-values, ambiguous CV procedures, overclaiming |
| **R2: Methods maximalist** | Comparisons, benchmarks | More baselines, ablations, parameter sensitivity | "Why didn't you compare to X?" left unanswered |
| **R3: Biology-first reviewer** | Pathway interpretation, clinical relevance | Factor labels, biological plausibility, therapeutic hooks | Uninterpretable factors, methods-heavy prose |
| **R4: PNAS generalist / AE** | Narrative clarity, broad impact | Clear conceptual hook, visual flow, "why it matters" | Technical jargon, buried insights, weak framing |
| **R5: Skeptical minimalist** | "Less is more," parsimony | Clean figures, focused claims, no scope creep | Over-engineered solutions, too many analyses |

**Strategic principle:** Prioritize changes that satisfy R1 (statistics) and R4 (narrative) without triggering R5 (minimalist). R2 (methods maximalist) requests can often be deflected in rebuttal by emphasizing the conceptual contribution.

### 7.1 Novelty Claim Needs Sharpening and Earlier Placement

**Issue:** The Introduction explains supervised NMF but delays the *key conceptual break*: **supervising W instead of H**.

**Current State:**
- The distinction IS made in the Introduction (paragraph 4): "both integrate the survival objective through the sample-specific loadings rather than the gene-level programs"
- However, this critical point is buried and could be easily missed

**Why Reviewers May Push Back:**
- Supervised factorization + Cox is not new *per se*
- Without early emphasis, reviewers may lump DeSurv with SurvNMF-like approaches

**Actionable Revision:**
In the **first page of the Introduction**, explicitly state:

> "Existing survival-aware NMF methods supervise **sample loadings (H)**, improving prediction but leaving **gene programs (W)** largely unchanged. DeSurv instead embeds survival supervision directly into **W**, forcing the gene programs themselves to align with outcome."

Consider adding a **contrast table** or schematic reference to Fig. 1 listing:
- Supervision target (W vs H)
- Consequence for interpretability
- Portability to new cohorts (closed-form projection vs NNLS)

### 7.2 Outcome-Guided Rank Selection: Anticipate "Double-Dipping" Concerns

**Issue:** Using survival to select k and then reporting survival associations invites skepticism.

**Why This Matters:**
- Statisticians may argue the rank is tuned to maximize survival signal, inflating downstream hazard ratios

**What the Paper Already Does Well:**
- Nested CV with held-out folds
- External validation across cohorts
- Simulation with known true rank

**What's Missing:**
The paper does not explicitly state:
- Rank selection is done *entirely within training folds*
- External cohorts are never used in BO or tuning
- Reported HRs in validation are **out-of-sample**

**Actionable Revision:**
Add **one explicit paragraph** in Methods or Results clarifying:

> "Model selection, including rank (k) and supervision strength (α), was performed entirely within training data using nested cross-validation. Validation cohorts were held out from all tuning procedures. Reported hazard ratios and survival associations in external datasets reflect purely out-of-sample generalization."

Consider adding (if feasible):
- A **sensitivity analysis** showing that nearby k values give similar qualitative biological structure

### 7.3 Reframe Performance Gains Around Structural Reorganization, Not C-Index

**Issue:** Some C-index improvements are modest (~0.02–0.05).

**Why This is Risky:**
- Reviewers may dismiss predictive gains as incremental

**Current Strength:**
- Fig. 4 is excellent: variance vs survival disentanglement is the real story
- Results section (lines 117-123) does discuss structural reorganization

**What's Needed:**
The C-index framing is still prominent. The structural insight should be the lead message.

**Actionable Revision:**
In Results and Discussion, explicitly state:

> "The key gain is not raw discrimination but the **reorganization of transcriptional axes** such that survival-relevant biology is isolated from variance-dominant but prognostically neutral signals."

Add to Discussion:

> "DeSurv is **not optimized primarily for maximal prediction**, but for re-allocation of latent structure toward outcome relevance. This reorganization—not incremental C-index improvement—is the primary contribution."

PNAS reviewers respond well to *mechanistic clarity over marginal AUC/C-index gains*.

### 7.4 Clarify Why Supervising W Improves Generalization

**Issue:** The portability argument is strong but buried in SI.

**Current State:**
- Results mention "applying the learned gene-factor weights to the sample's gene expression profile ($W^TX$)" (line 154)
- But the **closed-form advantage** is not stated explicitly

**Actionable Revision:**
Move or summarize the **projection argument** into the main text:

> "Because DeSurv learns gene-level programs (W), new samples can be scored via closed-form projection Z = W^⊤X without per-sample optimization. This eliminates fitting noise that accumulates in methods requiring sample-wise NNLS or iterative updates, contributing to the stability observed across validation cohorts."

Explicitly connect this to the **PDAC → bladder transfer** result.

### 7.5 Figure-by-Figure Assessment

> **📄 Detailed Analysis:** See [FIGURE_ANALYSIS.md](FIGURE_ANALYSIS.md) for comprehensive panel-by-panel evaluation including visualization effectiveness ratings, proposed improvements, and alternative plotting approaches.

The figures tell the narrative arc correctly (method → why needed → proof it works → biological meaning → generalization → transfer), but several have reviewer vulnerabilities related to **selection transparency** and **effect size reporting**.

#### Figure 1: Model Overview
| Aspect | Status | Issue |
|--------|--------|-------|
| W-sharing visual | ✅ Strong | Communicates novelty hook |
| BO depiction | ✅ Clear | Panel B shows hyperparameter selection |
| Notation consistency | ⚠️ Risk | Z = W^⊤X in caption vs W^TX elsewhere - verify consistency |

**Recommended additions:**
- Micro-callout: "Supervision acts on **W (gene programs)** → improves interpretability + cohort transfer"
- Label "CV C-index objective" visually near evaluation step in Panel B

#### Figure 2: Rank Selection + BO
| Aspect | Status | Issue |
|--------|--------|-------|
| Panels A-C (unsupervised) | ✅ Effective | Shows "no clear elbow / ambiguous criterion" |
| Panel D (GP surface) | ⚠️ Vulnerability | Shows mean only - no uncertainty layer |
| Panel E (simulation) | ⚠️ Vulnerability | Doesn't signal "CV-only" selection |

**Recommended additions:**
- Add thin uncertainty layer (SE contours or evaluated configuration markers) to Panel D
- Add note "selected via CV objective" in Panel E label
- Annotate selected k explicitly on heatmap

#### Figure 3: Simulation Results
| Aspect | Status | Issue |
|--------|--------|-------|
| C-index distribution | ✅ Clear | Distribution shift is convincing |
| Precision comparison | ✅ Stark | Near-zero for α=0 is visually compelling |
| Self-containedness | ⚠️ Missing | Simulation regime (k=3) not in panel |

**Recommended additions:**
- In-panel microlabels showing key simulation parameters (k=3, lethal factors)
- Note "signature size (n_top) tuned via BO" to prevent "unfair baseline" critique

#### Figure 4: PDAC Biological Reorganization (THE KEY FIGURE)
| Aspect | Status | Issue |
|--------|--------|-------|
| Panels A-B (correlations) | ✅ Convincing | Shows organization differs |
| Threshold ">0.2 only" | ⚠️ Risk | Could look like selective reporting |
| Panel C (variance vs survival) | ✅ Strong | Directly supports "variance ≠ prognosis" |
| Survival metric | ⚠️ Unclear | "change in partial log-likelihood" may confuse |
| Panel D (correspondence) | ✅ Good | Shows reorganization not replacement |
| Factor labels | ❌ Missing | Readers must infer factor identity |

**Recommended additions:**
- Provide full correlation matrix in SI (or sensitivity to threshold)
- Add legend sentence in C clarifying survival contribution metric
- **Add factor labels** (e.g., "Exocrine/composition", "Classical tumor", "Activated TME")

#### Figure 5: External Validation
| Aspect | Status | Issue |
|--------|--------|-------|
| Panel A (forest plot) | ✅ Supports claim | Shows more stable HRs |
| Factor selection rule | ⚠️ Not visible | "largest increase in partial log-likelihood in training" not in figure |
| Panels B-C (KM curves) | ⚠️ Missing stats | No log-rank p, HR, or CI on plots |
| Pooling/stratification | ⚠️ Unclear | Not stated whether stratified by cohort |

**Recommended additions:**
- Put factor selection rule in panel title: "Factor selected in training by Δ partial log-lik"
- Add log-rank p and HR annotation directly on KM plots (B and C)
- State "pooled with cohort-stratified Cox" if applicable

#### Figure 6: Bladder + Cross-Cancer Transfer
| Aspect | Status | Issue |
|--------|--------|-------|
| Panel A (variance vs survival) | ✅ Consistent | Matches PDAC message (Fig 4C) |
| Panel B (transfer KM) | ⚠️ Missing stats | No HR/p-value shown |
| Biological interpretation | ⚠️ Missing | No label for what transferred factor represents |
| "Cross-cancer biology" claim | ⚠️ Risk | Reviewers may ask if it's just proliferation/inflammation |

**Recommended additions:**
- Report HR and p-value on Panel B to match Fig 5 standards
- Add caption note about what the transferred factor biologically represents
- Pre-empt "generic signatures" concern in text

### 7.6 Minor Editorial Improvements

| Item | Current State | Recommendation |
|------|---------------|----------------|
| **Significance Statement** | Strong but abstract | Make more concrete: explicitly name *gene programs*, *survival alignment*, and *cross-cancer transfer* |
| **Related Work** | Survival-NMF methods mentioned briefly | Cite survival-supervised topic models (sLDA analogs) and explain *why they fail in genomics* |
| **Optimization Section** | Excellent rigor | Add a "Reader's Guide" sentence explaining why the convergence proof matters for reproducibility |

### 7.7 Likely Reviewer Questions and Pre-emptive Responses

| Reviewer Concern | How to Address |
|------------------|----------------|
| "Isn't this just supervised NMF?" | Emphasize **W-level supervision** and program restructuring in Introduction |
| "Are you overfitting survival?" | Clarify nested CV + external validation; add explicit out-of-sample statement |
| "Why modest C-index gains?" | Reframe around interpretability and structural reorganization |
| "Why PDAC?" | Argue TME heterogeneity + known variance-dominant confounders make it an ideal test case |
| "Will this generalize?" | Highlight bladder transfer + closed-form projection argument |
| "How does this differ from sLDA?" | Explain that topic models assume document-word generating process inappropriate for continuous expression |
| "The figures look cherry-picked" | Add selection protocol microtext, effect sizes, and full matrices in SI |
| "What does the transferred factor represent biologically?" | Add biological interpretation label to Fig 6 |
| "Is bladder transfer just generic proliferation/inflammation?" | Pre-empt in text; discuss what programs are shared |

---

## 8. Supplement Completeness

### 8.1 Algorithms

| Algorithm | Status | Notes |
|-----------|--------|-------|
| Algorithm S1 (BCD for DeSurv) | ✅ Complete | Lines 76-105, properly labeled |
| Algorithm S2 (W update with backtracking) | ✅ Complete | Lines 209-237 |
| Algorithm S3 (Cross-validation) | ✅ Complete | Lines 690-728 |

**Note:** Verify "S" prefix renders correctly in compiled PDF.

### 8.2 Theorem and Proofs

| Element | Status | Location |
|---------|--------|----------|
| Assumptions 1-2 | ✅ | Regularity conditions |
| Lemma 1 | ✅ | Continuity and bounded level sets |
| Lemma 2 | ✅ | Monotone descent |
| Lemma 3 | ✅ | Blockwise optimality |
| Theorem 1 | ✅ | Convergence to stationary point |
| Normalization note | ✅ | Column normalization preserves stationarity |
| Coxnet remark | ✅ | Gap between proof and practice acknowledged |

### 8.3 Dataset Descriptions

| Dataset | Status | Missing Info |
|---------|--------|--------------|
| TCGA-PAAD | ⚠️ Partial | Sample size after filtering |
| CPTAC-PDAC | ⚠️ Partial | Sample size after filtering |
| Dijk | ⚠️ Partial | Sample size after filtering |
| Moffitt | ⚠️ Partial | Sample size after filtering |
| PACA-AU | ⚠️ Partial | Sample size after filtering |
| Puleo | ⚠️ Partial | Sample size after filtering |
| IMVigor210 (bladder) | ⚠️ Partial | Sample size after filtering |

**Missing preprocessing details:**
- Number of genes retained (`ngene` parameter value)
- Low-expression gene filtering threshold
- Log-transformation base (log2 or natural log)
- TPM conversion method

### 8.4 Data Split Inconsistency

| Source | Bladder Data Split |
|--------|-------------------|
| Main methods | 70/30 |
| Supplement (line 806) | 80/20 |

**Action:** Resolve inconsistency.

### 8.5 Hyperparameter Documentation

| Parameter | Documented | Missing |
|-----------|------------|---------|
| α ∈ [0,1) | ✅ | Actual search bounds |
| λ_H ≥ 0 | ✅ | Log-scale bounds used |
| λ ≥ 0 | ✅ | Log-scale bounds used |
| ξ ∈ [0,1] | ✅ | - |
| k (rank) | ❌ | Range (e.g., k ∈ {2,...,15}) |
| n_top | ❌ | Range for top genes per factor |
| ε_H, ε_W | ❌ | Default floor values |
| δ_max | ✅ | 10^6 |
| max_bt | ❌ | Max backtracking iterations |
| ρ | ❌ | Backtracking parameter |
| tol, maxit | ❌ | Convergence tolerance, max iterations |
| BO initial points | ❌ | Number of random initial points |
| BO total evaluations | ❌ | Evaluation budget |
| EI κ value | ❌ | Acquisition function parameter |

### 8.6 Missing Elements

| Element | Status |
|---------|--------|
| Supplementary figure numbering (Figure S1, etc.) | ❌ Not used |
| Supplementary table numbering (Table S1, etc.) | ❌ Not used |
| GitHub repository reference | ❌ Not in supplement |
| Consensus initialization procedure | ❌ Referenced but not described |
| Signature truncation details | ❌ "Details in SI" but not found |

---

## 9. Action Checklist

> **📄 See also:** [FIGURE_ANALYSIS.md](FIGURE_ANALYSIS.md) for detailed figure-specific execution checklist with effort estimates.

### 9.1 Critical (Must Fix Before Submission)

**Code Bugs (blocking):**
| Task | Effort | Risk if Skipped |
|------|--------|-----------------|
| ☐ Fix `sdZ = fit$meanZ[,ns]` → `fit$sdZ[,ns]` in `R/compute_metrics.R:9` | 5 min | All metrics incorrect |
| ☐ Use `keep` variable in `R/select_best_init.R` | 10 min | Invalid models selected |
| ☐ Verify and fix beta backtracking logic in `functions.cpp:472-473` | 1–2 hours | Slow/wrong convergence |
| ☐ Remove all `browser()` statements from `targets_common_pipeline.R` | 10 min | Pipeline hangs on Slurm |
| ☐ **Re-run all analyses** after bug fixes | 4–8 hours | Results may change |

**PNAS Compliance (blocking):**
| Task | Effort | Risk if Skipped |
|------|--------|-----------------|
| ☐ Reduce references from 71 to ≤50 | 1–2 hours | Desk rejection |
| ☐ Move 2 figures to Supplementary Information | 1 hour | Format violation |
| ☐ Replace placeholder keywords | 5 min | Format violation |
| ☐ Complete Author Contributions statement | 15 min | Incomplete submission |
| ☐ Complete Conflict of Interest declaration | 5 min | Incomplete submission |
| ☐ Complete Acknowledgements section | 15 min | Incomplete submission |
| ☐ Add complete institutional addresses | 15 min | Format violation |

### 9.2 High Priority (Editorial Framing)

These address likely reviewer concerns and should be prioritized for PNAS success:

- [ ] **Sharpen W vs H distinction** in Introduction first page (currently buried in para 4)
- [ ] **Add explicit out-of-sample statement** clarifying rank selection is within-fold only
- [ ] **Reframe around structural reorganization**, not C-index gains
- [ ] **Move projection argument** (closed-form Z = W^⊤X) to main text
- [ ] **Add sensitivity analysis** showing nearby k values give similar biological structure
- [ ] Make Significance Statement more concrete (name gene programs, survival alignment, cross-cancer transfer)
- [ ] Add Reader's Guide sentence for convergence proof relevance
- [ ] Cite survival-supervised topic models (sLDA) and explain why they fail in genomics

### 9.2.1 Figure-Specific Improvements (High Leverage)

These are minimal-effort, maximum-reviewer-proofing edits:

**Selection/Evaluation Protocol Transparency:**
- [ ] **Fig 1:** Add micro-callout "Supervision acts on W (gene programs)" and label "CV C-index objective" in Panel B
- [ ] **Fig 2D:** Add uncertainty layer (SE contours or evaluated configuration markers)
- [ ] **Fig 2E:** Add note "selected via CV objective" in panel
- [ ] **Fig 2D:** Annotate selected k explicitly on heatmap
- [ ] **Fig 3:** Add in-panel microlabels (k=3, lethal factors) and note "n_top tuned via BO"

**Effect Size Annotations:**
- [ ] **Fig 5B-C:** Add log-rank p-values and HR annotations directly on KM plots
- [ ] **Fig 5:** State factor selection rule in panel title ("Factor selected by Δ partial log-lik in training")
- [ ] **Fig 6B:** Add HR and p-value to match Fig 5 standards

**Factor Interpretation:**
- [ ] **Fig 4:** Add factor labels ("Exocrine/composition", "Classical tumor", "Activated TME")
- [ ] **Fig 4C:** Add legend sentence clarifying survival contribution metric (univariate Cox Δ log-lik)
- [ ] **Fig 4:** Provide full correlation matrix in SI (address ">0.2 threshold" selectivity concern)
- [ ] **Fig 6:** Add note about what transferred factor biologically represents

**Notation Consistency:**
- [ ] Verify Z = W^⊤X notation is consistent throughout (caption uses W^⊤X, text uses W^TX and X^⊤W)

### 9.3 High Priority (Technical)

- [ ] Fix symbol overloading: rename gradient-balancing `δ` to `γ`
- [ ] Fix symbol overloading: rename hyperparameter `θ` to `Θ`
- [ ] Fix incomplete equation at supp_methods.Rmd line 306
- [ ] Fix NMF→Cox typo at supp_methods.Rmd line 646
- [ ] Resolve 70/30 vs 80/20 data split inconsistency
- [ ] Expand main methods section to ~1.5 pages
- [ ] Expand discussion to 2-3 paragraphs (including limitations, clinical translation)
- [ ] Add quantitative metrics (C-index values, HRs with CIs, p-values)
- [ ] Add simulation robustness across multiple scenarios
- [ ] Document convergence criterion (cosine similarity vs loss change)
- [ ] Document W update clamping [0.2, 5.0]
- [ ] Document H penalty scaling factor

### 9.4 Medium Priority

- [ ] Add ORCID IDs for all authors
- [ ] Add notation table to supplement
- [ ] Standardize `n_{event}` formatting
- [ ] Unify "event indicators" vs "censoring indicators" terminology
- [ ] Add sample sizes table for all datasets
- [ ] Document all hyperparameter search bounds
- [ ] Add GitHub repository reference to supplement
- [ ] Describe consensus initialization procedure
- [ ] Add supplementary figure/table numbering scheme
- [ ] Cite survival-supervised topic models (sLDA) and explain why they fail in genomics

### 9.5 For Enhanced PNAS Suitability

- [ ] Consider single-cell validation of discovered programs
- [ ] Consider treatment response prediction analysis
- [ ] Develop bladder analysis to same depth as PDAC (or move to supplement)
- [ ] Add therapeutic/clinical relevance discussion
- [ ] Add comparison with additional methods (sparse NMF, PLSR-Cox)
- [ ] Explicitly connect PDAC → bladder transfer to closed-form projection advantage

---

## Appendix: Files Reviewed

### Paper Files
- `paper/paper.Rmd` - Main manuscript structure and YAML
- `paper/02_introduction_30102025.Rmd` - Introduction
- `paper/03_methods.Rmd` - Methods section
- `paper/04_results.Rmd` - Results section
- `paper/05_discussion.Rmd` - Discussion section
- `paper/supp_methods.Rmd` - Supplementary methods
- `paper/references_30102025.bib` - Bibliography (71 references)

### Code Files
- `R/compute_metrics.R` - Metric calculations (Bug #1)
- `R/select_best_init.R` - Initialization selection (Bug #2)
- `R/bo_helpers.R` - Bayesian optimization helpers
- `R/run_seed_fits.R` - Seed fitting functions
- `targets_common_pipeline.R` - Main pipeline (browser() statements)
- `../DeSurv/src/functions.cpp` - C++ optimization code

### Related Documents
- `CODE_REVIEW.md` - Comprehensive code review findings
- `BUGS_FOR_STUDENT.md` - Educational bug explanations
- `CLAUDE.md` - Project documentation

---

## Appendix B: Editorial Review Analysis

The following editorial concerns were analyzed for legitimacy against the actual manuscript content:

| Concern | Assessment | Evidence |
|---------|------------|----------|
| W vs H distinction buried | **PARTIALLY CORRECT** | Distinction IS made (Intro para 4) but could be sharper/earlier |
| "Double-dipping" risk | **CORRECT** | Nested CV used but not explicitly stated as within-fold only |
| C-index framing too prominent | **CORRECT** | Structural reorganization story (Fig 4C) should lead |
| Projection argument buried | **CORRECT** | Z = W^⊤X mentioned but closed-form advantage not stated |
| Significance statement abstract | **CORRECT** | Could be more concrete with specific terms |
| Related work incomplete | **CORRECT** | sLDA-type methods not discussed |

## Appendix C: Figure Review Analysis

All six main-text figures were analyzed for (1) whether visual evidence supports caption claims, and (2) whether they advance the narrative arc. Results:

| Figure | Caption Support | Narrative Role | Key Vulnerabilities |
|--------|-----------------|----------------|---------------------|
| **Fig 1** | ✅ Mostly yes | Front door to story | Notation consistency; missing projection callout |
| **Fig 2** | ✅ Yes | "Why BO matters" | No uncertainty on GP surface; no "CV-only" signal |
| **Fig 3** | ✅ Yes | Controlled proof | Simulation details not in-panel |
| **Fig 4** | ✅ Mostly yes | THE turning point | Threshold selectivity; unclear survival metric; no factor labels |
| **Fig 5** | ⚠️ Partially | Reproducibility step | Factor selection rule not visible; no p-values on KM |
| **Fig 6** | ✅ Conceptually yes | Capstone claim | No effect size; no biological label for transferred factor |

**Overall Assessment:** Figures tell the story in the right order, but **selection transparency** is the single biggest gap. The solution is not new experiments—it's **in-figure annotations** that make the evaluation protocol unmistakable.

**Highest-leverage edits (minimal work, maximum reviewer-proofing):**
1. Add selection/evaluation protocol microtext to Fig 2D/E, Fig 5A/B/C, Fig 6B
2. Add effect size annotations (HR, CI, p) directly on KM plots (Figs 5–6)
3. Add factor labels to Fig 4 (reduces interpretability burden)

All legitimate figure concerns have been incorporated into Section 7.5 (Figure-by-Figure Assessment) and Section 9.2.1 (Figure-Specific Improvements).

---

## Appendix D: Recalibrated PNAS Assessment

This section provides a recalibrated assessment that better aligns with PNAS editorial criteria (vs methods-journal standards).

### What Remains Valid (Address These)

1. **Simulation breadth** – Showing only "easy" regime undersells the method. Add one harder scenario (signal aligned with variance, or weak signal) plus the null.

2. **Methods section too terse** – Key algorithmic details need to move from SI to main text.

3. **Discussion underdeveloped** – Six lines is insufficient. Expand to 2-3 paragraphs covering limitations, when to use DeSurv, and computational considerations.

4. **Minor notation/figure issues** – DeSurv vs deSurv, truncated axis labels, cohort characteristics table.

### Where Initial Review Overreached

| Initial Request | Why It's Excessive | What's Actually Needed |
|-----------------|-------------------|------------------------|
| **4-5 method comparisons** | Bioinformatics/JRSS thinking, not PNAS | One comparator (LASSO-Cox) + explanation of why others aren't appropriate (they supervise H, not W) |
| **Pathway analysis, GSEA, single-cell validation** | Treats this as biology discovery paper | This is a **structural/methodological claim**; lightweight enrichment in SI is fine |
| **"Cross-cancer overstated"** | Misses the methodological point | Slight language tempering ("suggests" vs "demonstrates") is sufficient |
| **"Incremental rather than transformative"** | Wrong framing | PNAS judges by **conceptual reframing**, not benchmark tables; W-level supervision is conceptually meaningful |

### Revised Recommendation: Minor-to-Moderate Revision

**Required revisions:**
1. Add 1-2 simulation scenarios beyond "easy" (can be SI with main-text reference)
2. Expand Methods section with key algorithmic details
3. Expand Discussion to 2-3 paragraphs
4. Add one comparator method (e.g., LASSO-Cox) with explanation
5. Minor tempering of cross-cancer language
6. Fix notation inconsistencies and figure labels

**What to resist in revision:**
- Exhaustive benchmark comparisons
- Deep biological validation (pathway analyses, etc.)
- Removing or minimizing the bladder transfer result

### Key Rebuttal Language

If reviewers request scope creep:

> "Our goal is not to outperform all supervised predictors, but to demonstrate that survival supervision applied at the gene-program level reorganizes latent structure in a way that improves interpretability, stability, and generalization—properties not addressed by existing approaches that supervise through sample loadings."

This reorients from "did you beat everything?" to "did you demonstrate the conceptual claim?"—which is what PNAS actually cares about.

### The Core Contribution (Framing Guide)

The paper's contribution is:

> **Survival supervision applied at the gene-program level reorganizes latent structure in ways that improve interpretability, stability, and generalization.**

This is a **structural/methodological claim**, not a pathway discovery claim or a prediction benchmark claim. All revisions should reinforce this framing.

---

## Appendix E: Panel-by-Panel Figure Analysis

A comprehensive panel-by-panel analysis of all six manuscript figures is available in **[FIGURE_ANALYSIS.md](FIGURE_ANALYSIS.md)**.

### Summary

| Figure | Overall Rating | Key Issue | Primary Recommendation |
|--------|----------------|-----------|------------------------|
| **Fig 1** | 7/10 | Panel B dominant, overshadows Cox pathway | Add visual callout for W-level supervision |
| **Fig 2** | 6/10 | Layout inverts problem→solution flow | Reorder: B-C-D above A |
| **Fig 3** | 7/10 | Missing simulation parameters | Add "true k=3, n=100 replicates" annotation |
| **Fig 4** | 6/10 | **Factor labels missing** | Add biological labels (Exocrine, TME, etc.) |
| **Fig 5** | 6/10 | No effect sizes on KM plots | Add HR, 95% CI, p-value to panels B-C |
| **Fig 6** | 6/10 | Same issues as Fig 5; narrative weak | Add effect sizes; explain transfer mechanism |

### High-Impact, Low-Effort Fixes (Prioritized)

1. **Figure 4:** Add factor labels to heatmaps (dramatically improves interpretability)
2. **Figures 5-6:** Add HR, CI, p-value annotations to all KM plots
3. **Figure 2:** Reorder panels to show problem before solution
4. **Figure 3:** Add simulation ground truth annotation
5. **All figures:** Standardize legend terminology

For detailed panel-level analysis, visualization alternatives, and aesthetic recommendations, see the full document.

---

*Review generated: 2026-01-24*
*Updated: 2026-01-24 with editorial strategy analysis*
*Updated: 2026-01-24 with figure-by-figure stress test analysis*
*Updated: 2026-01-24 with narrative arc analysis (Section 6)*
*Updated: 2026-01-24 with recalibrated PNAS assessment (Appendix D)*
*Updated: 2026-01-24 with panel-by-panel figure analysis (Appendix E)*
