# DeSurv Manuscript: Narrative Arc Analysis (Revised)

**Document Purpose:** Detailed analysis of the paper's storytelling structure and how figures/tables support the narrative arc. Revised to reflect the prediction-validation framing proposed for PNAS submission.

**Companion documents:** SUGGESTED_TEXT.md (drop-in language), PNAS_SENIOR_EDITOR_REVIEW.md (review findings), FIGURE_ANALYSIS.md (panel-by-panel review).

---

## Executive Summary

The DeSurv manuscript employs a **prediction-validation** narrative structure. The Introduction establishes a biological problem (variance != prognosis), provides statistical and biological reasons why survival supervision should help, and makes four specific testable predictions. The Results section validates each prediction in order, with escalating generalization: simulation --> PDAC training --> PDAC external validation. Four main-text figures are positioned as evidence for specific introduction claims.

**Core narrative shift (vs. previous arc):** The paper is no longer framed as "here is a method and what it finds" but as "here is why supervision should help, and here is the evidence that it does." This is the hypothesis-driven framing that PNAS expects.

---

## The Prediction-Validation Structure

### Predictions Made in the Introduction

| # | Prediction | Statistical/Biological Basis | Validation Figure |
|---|-----------|------------------------------|-------------------|
| P1 | Variance-dominant factors need not be prognostically relevant | Aran et al. 2015: ~47% of subtyping genes correlate with purity; Bair & Tibshirani 2006: "no guarantee principal components will be associated with survival" | Fig 3C |
| P2 | Supervision should concentrate survival signal into fewer factors | Cook & Forzani 2008: response-guided subspaces target outcome-relevant directions; information bottleneck (Tishby et al. 1999) | Fig 3A-B |
| P3 | Supervision should improve cross-cohort generalization | Sufficient dimension reduction theory: outcome-aligned subspaces are more portable; survClust (Arora et al.): 97% vs 68% accuracy | Fig 4 |
| P4 | Supervision should resolve rank selection ambiguity | Frigyesi & Hoglund 2008: cophenetic correlation can overfit; outcome-driven criterion replaces heuristic disagreement | Fig 2D |

### Results Validate Each Prediction

Each results subsection opens by referencing the corresponding introduction prediction and closes by stating whether the finding confirms, refines, or qualifies it.

---

## The Three-Act Story

### Act I: The Case for Supervision (Introduction)

**Five paragraphs, each with a specific narrative function:**

| Paragraph | Narrative Function | Key Content |
|-----------|-------------------|-------------|
| 1 | The biological problem | Tumor transcriptomes are mixtures; variance != prognosis (cite Aran, Bair, Moffitt). Lead with biology, not methods. |
| 2 | Why the two-step approach creates misalignment | Standard paradigm (unsupervised NMF + post-hoc evaluation) has produced foundational insights (acknowledge Collisson, Moffitt, Bailey, DECODER) but optimizes a different objective (reconstruction) than the evaluation criterion (survival). Frame as tradeoff, not failure. |
| 3 | The case for supervision | Statistical grounding (Cook & Forzani sufficient dimension reduction; Bair & Tibshirani supervised PCA; Arora et al. survClust). Biological argument: supervision concentrates capacity on outcome-relevant programs. Second benefit: simplifies rank selection by providing a clinically meaningful criterion. **This paragraph generates all four predictions.** |
| 4 | What DeSurv does | W-level supervision (Z = W^T X); semi-supervised design (alpha controls balance); Bayesian optimization for model selection. Technical innovation in one paragraph. |
| 5 | What we will show | Explicit statement of the four predictions and the two evaluation settings (simulation, PDAC). Readers know exactly what to expect from Results. |

**Critical tone notes:**
- Paragraph 1 leads with biology (variance != prognosis), not methodology (NMF limitations)
- Paragraph 2 acknowledges the successes of the two-step approach before identifying its structural limitation
- Paragraph 3 grounds the case for supervision in established theory (Cook, Bair, Tishby), not hand-waving
- Paragraph 4 describes DeSurv as a semi-supervised method: both the reconstruction and survival terms are needed (see SUGGESTED_TEXT.md Section H for the argument)
- Paragraph 5 makes predictions explicit, so Results read as confirmation rather than exploration. The focused scope (simulation + PDAC) keeps the paper tightly argued.

### Act II: The Evidence Cascade (Results)

```
+-----------------------------------------------------------------------+
|  MODEL OVERVIEW (Fig 1)                                                |
|  "Here is the DeSurv framework"                                        |
|  Narrative function: Technical foundation for all subsequent claims     |
+-------------------------------+---------------------------------------+
                                |
                                v
+-----------------------------------------------------------------------+
|  SIMULATION + MODEL SELECTION (Fig 2) -- validates P4, P1, P2          |
|  "Under known ground truth, DeSurv recovers programs more reliably.    |
|   Survival-driven BO resolves rank selection ambiguity."                |
+-------------------------------+---------------------------------------+
                                |
                                v
+-----------------------------------------------------------------------+
|  PDAC FACTOR STRUCTURE (Fig 3) -- validates P1, P2 in real data        |
|  *** THE TURNING POINT ***                                             |
|  "Variance != prognosis: exocrine factor dominates variance but not    |
|   survival. DeSurv reorganizes programs toward outcome-relevant axes." |
+-------------------------------+---------------------------------------+
                                |
                                v
+-----------------------------------------------------------------------+
|  EXTERNAL VALIDATION (Fig 4) -- validates P3                           |
|  "Survival-aligned programs generalize across independent PDAC         |
|   cohorts with consistent HRs; NMF factors do not."                    |
+-----------------------------------------------------------------------+
```

**Evidence escalation:** Simulation (controlled) --> PDAC training (real, factor structure) --> PDAC validation (real, independent cohorts).

### Act III: Synthesis and Implications (Discussion)

**Five paragraphs:**

| Paragraph | Narrative Function | Key Content |
|-----------|-------------------|-------------|
| 1 | Summarize contributions | Survival supervision reorganizes the latent landscape --- this is structural reallocation, not just better C-index |
| 2 | PDAC biological context | DeSurv recapitulates known programs (basal/classical, CAF subtypes) --- the contribution is how they are recovered (de novo, without pre-specified signatures), not their identity |
| 3 | Semi-supervised tradeoff | Why DeSurv needs BOTH terms: reconstruction preserves biological interpretability; survival steers toward outcome-relevant decomposition. Alpha controls the balance. When supervision helps vs. when it doesn't (null scenario guidance). |
| 4 | Limitations | Cox PH assumption; computational cost of BO; theory-practice gap in convergence (Theorem 1 vs. implementation); dependence on outcome data quality; validation limited to PDAC |
| 5 | Broader implications | The principle (outcome-guided dimension reduction targets different subspaces than variance-driven reduction) extends beyond cancer. Cook & Forzani, information bottleneck. Natural extensions: multi-omics, spatial transcriptomics, EHR data, cross-cancer transfer. |

---

## Figure-by-Figure Narrative Role

### Figure 1: The Method Schematic
**Narrative Position:** Opening
**Story Function:** "Here is what DeSurv is"
**Prediction Connection:** Sets up the technical basis for all four predictions

| Panel | Content | Narrative Role |
|-------|---------|----------------|
| A | NMF + Cox architecture | Shows W is shared between reconstruction and survival; Z = W^T X |
| B | Bayesian optimization | Shows principled hyperparameter selection (connects to P4) |
| C-E | Model outputs | Shows interpretable outputs (gene programs, factors, coefficients) |

**What reader should take away:** DeSurv couples NMF with Cox by sharing W. Survival gradients act on W (gene programs), not H (sample loadings).

**Accessibility annotation to add:** Callout on panel A: "Survival gradients act on W (gene programs), not H (sample loadings)."

---

### Figure 2: Simulations and Model Selection
**Narrative Position:** First validation (P1, P2, P4)
**Story Function:** "Under controlled conditions DeSurv recovers programs more reliably; survival-driven BO resolves rank selection"

| Panel | Content | Narrative Role |
|-------|---------|----------------|
| A | C-index distributions across scenarios | DeSurv predicts survival better (validates P1, P2) |
| B | Precision distributions across scenarios | DeSurv identifies true genes (validates P1, P2) |
| C | k-selection histograms | DeSurv recovers the true rank more reliably (validates P4) |
| D | BO C-index heatmap (PDAC) | Clear surface --> unambiguous selection in real data (validates P4) |

**What reader should take away:** Supervision improves both prediction and program recovery under controlled conditions; outcome-driven BO resolves rank selection ambiguity in real data.

**Prediction-validation link:** P1 and P2 predicted that supervision concentrates survival signal; panels A-B confirm this under known ground truth. P4 predicted that incorporating outcomes would simplify rank selection; panels C-D confirm this in simulation and real data. Old unsupervised heuristic panels (residuals, cophenetic, silhouette) moved to supplement S4.

**Critical framing:** The null and mixed scenarios must be reported:
- **Null (no survival signal):** DeSurv does not outperform NMF --- confirms the method doesn't overfit when there's nothing to find
- **Mixed (partial overlap between variance and prognosis):** Advantage attenuated --- supervision helps most when variance and prognosis diverge
- **Easy (low-variance lethal programs):** DeSurv advantage is greatest --- the regime where the problem is hardest for unsupervised methods

---

### Figure 3: PDAC Factor Structure (THE KEY FIGURE)
**Narrative Position:** The turning point (P1, P2 in real data)
**Story Function:** "Here is the insight --- variance != prognosis, and DeSurv resolves this"

| Panel | Content | Narrative Role |
|-------|---------|----------------|
| A | DeSurv factor-program heatmap | Shows survival-aligned structure (immune-stromal, basal-like) |
| B | NMF factor-program heatmap | Shows variance-dominant structure (exocrine axis) |
| C | Variance vs survival contribution | **THE PUNCHLINE:** programs explaining the most variance explain the least survival (promoted from old S7) |
| D | W-matrix correlation (DeSurv vs NMF) | Shows reorganization, not replacement (promoted from old S8) |

**What reader should take away:** DeSurv doesn't just predict better --- it **restructures latent programs** to isolate survival-relevant biology from variance-dominant noise. The exocrine factor dominates variance but contributes nothing to survival.

**Prediction-validation link:** P1 predicted this disconnect; P2 predicted that supervision would concentrate survival signal. Panel C confirms P1; panels A-B confirm P2; panel D shows the structural relationship between the two decompositions.

**Accessibility annotation to add:** Label axes on Panel C as "Expression variance explained (%)" vs. "Survival association (concordance index)." Add factor annotations directly on the plot ("Exocrine," "Basal-like," "Activated TME").

---

### Figure 4: External PDAC Validation
**Narrative Position:** Reproducibility checkpoint (P3)
**Story Function:** "Survival-aligned programs generalize; variance-driven programs do not"

| Panel | Content | Narrative Role |
|-------|---------|----------------|
| A | Forest plot of HRs | DeSurv estimates more stable across cohorts (promoted from old S9) |
| B | DeSurv pooled KM curves | Clear survival separation |
| C | NMF pooled KM curves | Weaker separation |

**What reader should take away:** Outcome-aligned programs are more portable across datasets.

**Prediction-validation link:** P3 predicted this, grounded in sufficient dimension reduction theory (Cook & Forzani 2008). The forest plot (now in main text, promoted from supplement) provides the quantitative evidence.

**Accessibility annotation to add:** Add HR, 95% CI, and p-values directly on each KM panel.

---

## The Central Insight: "Variance != Prognosis"

The entire paper builds toward one core insight, visualized in **Figure 3D**:

> **High-variance transcriptional axes are not necessarily prognostically relevant.** Standard NMF allocates model capacity to variance-dominant signals (like exocrine expression in PDAC) that contribute little to survival. DeSurv redirects this capacity toward survival-aligned programs.

In the prediction-validation structure, this insight is **stated as a prediction in the Introduction** (grounded in Aran et al. 2015, Bair & Tibshirani 2006), not revealed as a surprise in the Results. This makes the paper read as hypothesis-driven science rather than exploratory data analysis.

---

## The Semi-Supervised Framing

DeSurv is a **semi-supervised** method. The alpha parameter controls the balance between unsupervised reconstruction (alpha=0) and supervised Cox association (alpha=1). This framing is essential and should be made explicit in the paper.

### Why both terms are needed

| Term | What it ensures | What breaks without it |
|------|----------------|----------------------|
| Reconstruction (1-alpha) | Programs correspond to real expression patterns | Model degenerates into penalized Cox regression on arbitrary projections; no biological interpretability |
| Survival (alpha) | Programs are outcome-aligned | Programs follow variance, not prognosis; the rank selection problem is unsolvable by heuristics |
| Joint optimization | Selects from among comparably reconstructive solutions the one most aligned with outcomes | Neither biological completeness nor clinical relevance is guaranteed |

### When supervision helps vs. doesn't

| Setting | DeSurv advantage | Reason |
|---------|-----------------|--------|
| Variance and prognosis diverge (Fig 2 "easy" scenario) | Large | Unsupervised methods allocate capacity to wrong directions |
| Partial overlap (Fig 2 "mixed" scenario) | Moderate | Some prognostic signal captured by variance-dominant factors |
| No survival signal (Fig 2 "null" scenario) | None | No outcome-relevant subspace to target; DeSurv reduces to NMF |

This simulation-derived guidance should be discussed in both Results (as evidence) and Discussion (as guidance for practitioners).

---

## Rhetorical Strategies

### 1. Prediction-Validation Framing (NEW)
The Introduction makes four specific predictions; each Results subsection opens by referencing the prediction and closes by stating the finding. This creates a hypothesis-driven narrative.

### 2. Tradeoff Framing for Prior Art (NEW)
The two-step approach is acknowledged for its successes (Moffitt, Bailey, DECODER) before identifying the structural misalignment. The tone is "this approach has been foundational, but the objective mismatch creates a specific limitation that can be addressed."

### 3. Escalating Generalization Claims (RETAINED)
Evidence builds in scope:
```
Simulation (controlled) --> PDAC training (factor structure) --> PDAC validation (independent cohorts)
```

### 4. Consistent Baseline Comparison (RETAINED)
DeSurv vs alpha=0 (standard NMF) runs throughout, providing a consistent ablation.

### 5. Biological Grounding (RETAINED + ENHANCED)
Every claim tied to known biology. New: the PDAC subtyping timeline (Collisson 2011 --> consensus 2025) motivates why supervision is needed --- the field took ~14 years of retrospective evaluation to converge.

### 6. Semi-Supervised Argument (NEW)
Explicit framing: DeSurv needs BOTH terms. This prevents the objection "why not just do supervised regression?" and grounds the method in established semi-supervised learning theory (Chapelle et al. 2006).

### 7. Theoretical Grounding (NEW)
Claims grounded in sufficient dimension reduction (Cook & Forzani 2008), supervised PCA (Bair & Tibshirani 2006), and information bottleneck (Tishby et al. 1999). These references give the case for supervision analytical teeth beyond empirical demonstration.

---

## Information Disclosure Strategy

| Document | Role | Audience | New additions |
|----------|------|----------|---------------|
| Main Introduction | Predictions + accessible biological argument | Everyone | Theoretical grounding (Cook, Bair, Tishby); semi-supervised framing |
| Main Methods | Framework + key choices | Everyone | W-vs-H supervision distinction made explicit |
| Main Results | Prediction validation | Everyone | Transitional language linking each finding to its prediction |
| Main Discussion | Synthesis + implications + limitations | Everyone | Semi-supervised tradeoff paragraph; 5 explicit limitations; broader implications |
| Supplement | Complete derivations | Methods experts | Unchanged |

---

## Accessibility for PNAS Readership

PNAS readership spans all scientific disciplines. Key strategies:

1. **Lead every section with intuition, not formalism.** "Tumor transcriptomes are mixtures" before any equation.
2. **Minimize jargon in key paragraphs.** See SUGGESTED_TEXT.md Section K for jargon translation table.
3. **Frame the contribution as a general principle:** "In high-dimensional data where outcome-relevant features explain modest variance, incorporating outcome information during dimensionality reduction yields more portable, more interpretable representations."
4. **Make figures self-explanatory.** Add annotations, label axes in plain language, include HR/CI/p-values directly on KM panels.
5. **Use one-sentence summaries** at the start of each results subsection (bolded).
6. **The Significance Statement is the primary hook.** It must be understandable by an educated scientist outside computational genomics.

---

## Narrative Gaps Still to Address

### Currently Underemphasized
1. **The semi-supervised argument** --- needs explicit treatment (proposed language in SUGGESTED_TEXT.md Section H)
2. **The projection advantage** --- closed-form Z = W^T X enables scoring new samples without their survival data; not shared by H-supervised methods

### Currently Missing
1. **Explicit out-of-sample statement** --- nested CV is used but should be stated clearly and early
2. **Factor biological labels on Figure 3** --- annotations would help general readers
3. **Effect sizes on KM plots** --- HR, 95% CI, and p-values directly on panels
4. **Null scenario discussion** --- the null result is a strength, not a weakness; it validates that DeSurv doesn't overfit when there's no signal

---

## Summary: The Story DeSurv Tells

**Opening Hook (Biology-first):** The programs that explain the most transcriptomic variance are not the programs that drive clinical outcomes. Standard deconvolution methods optimize variance, not prognosis, creating a structural misalignment.

**The Case for Supervision (Theory-grounded):** Sufficient dimension reduction theory and empirical precedent (supervised PCA, survival-weighted clustering) predict that outcome-guided learning should concentrate prognostic signal, improve portability, and simplify model selection.

**Innovation (Technical):** DeSurv supervises through gene programs (W), not sample loadings (H), and uses Bayesian optimization for model selection. The semi-supervised design (alpha parameter) balances biological interpretability with clinical relevance.

**Four Predictions:** (1) variance != prognosis, (2) supervision concentrates survival signal, (3) supervision improves generalization, (4) supervision resolves rank selection.

**Evidence Chain:**
1. Under controlled conditions, supervision improves gene recovery, especially when prognosis diverges from variance; survival-driven BO resolves rank selection (Fig 2) --- validates P1, P2, P4
2. In PDAC, supervision reorganizes the latent landscape: exocrine variance is suppressed, survival signal is concentrated (Fig 3) --- validates P1, P2 in real data
3. Survival-aligned programs generalize across PDAC cohorts with consistent HRs (Fig 4) --- validates P3

**Closing (Broader principle):** Outcome-guided dimension reduction targets different subspaces than variance-driven reduction. This principle extends beyond cancer genomics.

---

## Appendix: Figure-Narrative Mapping Table

| Figure | Section | Prediction Validated | Key Claim | Evidence Type |
|--------|---------|---------------------|-----------|---------------|
| 1 | Model Overview | (Setup) | W is shared; BO for selection | Schematic |
| 2 | Simulation + Model Selection | P1, P2, P4 | Supervision improves recovery; BO resolves rank | Controlled simulation + real data |
| 3 | PDAC Factor Structure | P1, P2 | **Variance != prognosis**; DeSurv reorganizes programs | Real data comparison |
| 4 | External Validation | P3 | Generalizes within-cancer | Independent cohorts |

---

## Change Log

| Date | Change |
|------|--------|
| 2026-01-24 | Initial narrative arc analysis (problem-solution-validation structure) |
| 2026-02-10 | Revised to prediction-validation structure for PNAS. Added: semi-supervised framing, theoretical grounding, accessibility strategies, prediction-figure mapping, Discussion structure, limitation enumeration. Informed by PNAS_SENIOR_EDITOR_REVIEW.md and SUGGESTED_TEXT.md. |
| 2026-03-08 | Updated for figure reorganization: removed P5/Fig 6 (bladder), restructured to 4-prediction 4-figure layout. Promoted S7 (variance-survival), S8 (W-correlation), S9 (forest plot) to main text. |
