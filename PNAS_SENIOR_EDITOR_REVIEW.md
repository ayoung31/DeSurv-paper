# PNAS Senior Editor Review: DeSurv Manuscript

**Manuscript:** "Survival driven deconvolution (DeSurv) reveals prognostic and interpretable cancer subtypes"
**Authors:** Amber M. Young, Alisa Yurovsky, Didong Li, Naim U. Rashid
**Review Date:** February 10, 2026
**Reviewer Role:** Senior Editor, Proceedings of the National Academy of Sciences

---

## I. Summary Statement

This manuscript presents DeSurv, a framework that integrates nonnegative matrix factorization (NMF) with Cox proportional hazards modeling to discover survival-aligned gene programs from bulk tumor transcriptomes. The key innovation is that survival supervision enters through the gene-program matrix W (via factor scores Z = W^T X), rather than through the sample-loading matrix H as in prior survival-NMF methods. This architectural choice yields gene programs that are intrinsically aligned with clinical outcomes and portable to new cohorts via closed-form projection. The method is validated on pancreatic ductal adenocarcinoma (PDAC), simulation studies, and a cross-cancer transfer experiment to bladder cancer.

**Overall Recommendation: Major Revision**

The conceptual contribution---that W-level survival supervision reorganizes latent transcriptional structure toward prognostically relevant axes---is genuine and represents a meaningful advance for translational cancer genomics. However, the manuscript has significant deficiencies in PNAS formatting compliance, methodological reporting, comparator methods, and several technical concerns that must be addressed before publication.

---

## II. Key Strengths

1. **Genuine conceptual advance.** The distinction between supervising W (gene programs) vs. H (sample loadings) is not merely incremental. It fundamentally changes what the model learns: interpretable gene signatures that can be projected onto new data via Z = W^T X, rather than sample-specific coefficients that require per-sample optimization for external validation. This is a structural insight that PNAS readership would value.

2. **The "variance != prognosis" insight is compelling.** Figure 4C powerfully demonstrates that the dominant axes of transcriptomic variation (e.g., exocrine expression in PDAC) are not necessarily prognostically relevant. DeSurv concentrates survival signal into fewer factors while deprioritizing variance-dominant but outcome-neutral programs. This mechanistic understanding of how latent structure relates to clinical outcomes is the paper's strongest contribution.

3. **Cross-cancer transferability.** The demonstration that PDAC-trained DeSurv factors retain prognostic relevance in bladder cancer (Figure 6B) is a striking result that elevates the paper from a single-cancer methods paper to a general framework. This is the kind of generalizable insight PNAS favors.

4. **Rigorous theoretical foundation.** The supplementary materials provide a complete convergence proof (Theorem 1) with supporting lemmas, including careful treatment of the non-standard hybrid W update with backtracking. The block coordinate descent convergence to a stationary point under mild regularity conditions is established with appropriate mathematical rigor.

5. **Bayesian optimization for model selection.** The use of BO to jointly tune rank k, supervision strength alpha, and regularization parameters addresses the well-known ambiguity of rank selection in NMF---a practical contribution that many users of NMF methods will appreciate.

6. **Reproducible research infrastructure.** Software is available as an R package (github.com/ayoung31/DeSurv), with a reproducible analysis pipeline (github.com/ayoung31/DeSurv-paper).

---

## III. Major Comments

### Major 1: PNAS Formatting Non-Compliance (Must Fix)

The manuscript violates several PNAS submission requirements that would result in desk rejection:

| Requirement | PNAS Standard | Current Status | Action |
|-------------|---------------|----------------|--------|
| Figures | 4 (standard 6-page article) | **6 figures** | Move 2 to SI |
| References | ~50 (standard article) | **71 references** | Reduce by 21+ |
| Keywords | Required | **Placeholder text** ("one, two, optional...") | Replace |
| Author contributions | Required | **Placeholder** ("Please provide details...") | Complete |
| Conflict of interest | Required | **Placeholder** ("Please declare any conflict...") | Complete |
| Acknowledgements | Required | **Placeholder text** | Complete |
| Affiliations | Required | **Incomplete** (contains "Street, City, State, Zip") | Complete |
| Abstract | Consistent naming | Uses "deSurv" (lowercase) while rest of paper uses "DeSurv" | Standardize |

**Recommendation:** Move Figures 2 (rank selection details) and 6 (bladder analysis) to Supplementary Information, retaining references in the main text. This preserves the narrative while meeting the 4-figure constraint. Alternatively, merge Figures 5-6 into a single "Generalization" figure.

### Major 2: Discussion Section is Critically Underdeveloped

The Discussion section consists of only 3 short paragraphs (~15 sentences). For a PNAS Research Article, this is inadequate. A proper Discussion should address:

- **Limitations:** The Cox proportional hazards assumption (non-proportional hazards are common in immunotherapy settings); computational cost relative to standard NMF; performance with small sample sizes; sensitivity to censoring rates; the non-convexity of the objective (local minima across initializations).
- **Comparison with related frameworks:** How DeSurv relates to supervised topic models (sLDA), integrative clustering methods (iCluster), WGCNA, and the recent Spectra method (Kotliar et al., Nature Biotechnology 2024) which also performs supervised gene program discovery.
- **When DeSurv should NOT be used:** Settings where unsupervised factorization may be preferred (e.g., exploratory discovery without clinical endpoints, settings with unreliable survival data).
- **Clinical translation pathway:** How the discovered programs could inform therapeutic strategies or patient stratification in clinical trials.
- **Biological interpretation of cross-cancer transfer:** Why do PDAC-trained programs work in bladder cancer? Is this generic proliferation/inflammation signaling, or something more specific? This question must be addressed.

### Major 3: Comparator Methods Are Insufficient

The paper compares DeSurv only to standard NMF (alpha=0). This is the minimum ablation study, not a comparative evaluation. Reviewers will ask:

1. **Why not compare to the two prior survival-NMF methods cited (Huang et al., 2020; Le Goff et al., 2025)?** If these methods are deemed insufficient in the Introduction, the reader expects to see empirical evidence of this insufficiency, or at minimum a clear technical argument for why comparison is unnecessary.

**These methods do not provide software, where and how should we address this in the paper?**

2. **Why not compare to penalized Cox regression on all genes (LASSO-Cox)?** This is the most natural "prediction-only" baseline that doesn't perform factorization. If DeSurv's advantage is structural reorganization rather than prediction, showing that LASSO-Cox achieves similar C-index but without interpretable programs would strengthen the argument.

3. **Why not compare to sparse NMF or semi-NMF?** These are standard alternatives in the NMF literature.

**Recommendation:** Add at minimum LASSO-Cox as a comparator. For Huang and Le Goff, either include empirical comparison or provide a clear, specific technical argument for why comparison is not meaningful (e.g., "these methods supervise through H, so the learned W is not directly outcome-aligned and cannot be projected onto new data").

### Major 4: Simulation Results Are Narrow

Only one simulation scenario ("R0_easy") is shown, despite the code infrastructure containing multiple scenarios (R00_null, R_mixed, etc.). This creates a vulnerability to the charge of cherry-picking.

**Required additions:**
- **Null scenario (R00_null):** What happens when there is no survival signal? DeSurv should not outperform NMF here; showing this confirms that gains reflect genuine signal, not overfitting. Frame as: "Under null signal, DeSurv does not outperform unsupervised NMF, as expected."
- **Mixed/harder scenario (R_mixed):** What happens when lethal programs are partially aligned with variance-dominant programs? The DeSurv advantage should be smaller but still present.
- **Key simulation parameters:** The true rank (k=3), sample size, number of genes, censoring rate, and number of replicates should all be reported in the main text or figure annotation.

### Major 5: Methods Section Lacks Essential Detail

The Materials and Methods section (~30 lines) is too terse for a methods paper at PNAS, where the methods section should be self-contained. Key missing elements:

- **Cox partial likelihood formula** (currently only in supplement)
- **Sample sizes** for each dataset (after filtering)
- **Bayesian optimization specifics:** Number of initial points, total evaluation budget, acquisition function details, search bounds for each hyperparameter
- **Convergence criterion:** What stopping rule is used? (Note: the code uses cosine similarity of W and H between iterations, while the supplement describes relative loss change---these are different criteria; see Technical Comment 1 below)
- **Consensus initialization procedure:** Referenced as "SI Appendix" but not fully described there
- **Runtime estimates:** How long does DeSurv take relative to standard NMF?

### Major 6: Missing Quantitative Metrics Throughout

The paper relies heavily on visual comparisons without reporting key numerical results:

- **C-index values** are not reported (only shown visually in heatmap and boxplots)
- **Hazard ratios with 95% confidence intervals** are shown in the forest plot (Figure 5A) but NOT on the Kaplan-Meier curves (Figures 5B-C, 6B)
- **Log-rank p-values** are missing from all KM plots
- **Precision values** from simulations are not reported numerically
- **Effect sizes** for the DeSurv vs. NMF comparison (e.g., mean C-index difference with standard error) are not provided

This is a significant gap. All survival analyses should report HR (95% CI) and log-rank p-values. All method comparisons should report effect sizes with uncertainty.

---

## IV. Minor Comments

### Minor 1: Paper-Code Inconsistencies

Several discrepancies exist between the mathematical formulations in the paper/supplement and the actual code implementation:

| Aspect | Paper/Supplement | Code (functions.cpp) |
|--------|-----------------|---------------------|
| Convergence criterion | Relative loss change | Cosine similarity of W, H |
| W update ratio | Unbounded | Clamped to [0.2, 5.0] |
| Gradient-balancing factor | Ratio of squared Frobenius norms | Ratio of unsquared norms |
| H penalty scaling | lambda_H * H | lambda_H * (Xnorm/(n*k)) * H |
| Incomplete equation | Eq. S8 in supplement | Missing ||beta||_2^2 term |

These must be resolved. Either update the paper to match the implementation or vice versa. A mathematical reviewer will check these.

### Minor 2: Symbol Overloading

- **delta** is used for both event indicators (delta in {0,1}^n) and the gradient-balancing factor (delta^(t)). Rename the gradient-balancing factor to gamma^(t).
- **theta** is used for both model parameters (W,H,beta) and BO hyperparameters. Rename one.
- **Event indicators** are described as delta in R^n in the main text but delta in {0,1}^n in the supplement. They should be binary.

### Minor 3: Figure Quality Issues

- **Figure 2:** Layout inverts the narrative flow. Panels B-C-D (the "problem": unsupervised heuristics disagree) should appear above Panel A/D (the "solution": DeSurv provides clear guidance). Currently, the solution is shown first.
- **Figure 3:** Internal simulation regime label "R0_easy" is exposed in the strip header. Replace with descriptive text.
- **Figure 4:** Factor columns (F1, F2, F3) lack biological labels. Adding descriptive labels (e.g., "Classical tumor-associated", "Basal-like", "Activated TME") would dramatically improve interpretability.
- **Figure 5B-C:** Number-at-risk tables appear to have formatting issues with overlapping digits (e.g., "0302413 6 2 0"). Verify these render correctly.
- **Figure 6B:** The transferred factor is not identified (which PDAC factor? what does it represent biologically?).

### Minor 4: Bladder Cancer Analysis is Underdeveloped

The bladder analysis (Section "Survival-informed factorization identifies prognostic structure across cancers") shows only:
- Variance vs. survival plot (Figure 6A)
- Cross-cancer transfer KM curve (Figure 6B)

Missing: (a) DeSurv applied directly to bladder data (not just transferred PDAC model); (b) biological interpretation of what the transferred factor represents; (c) comparison with standard NMF applied to bladder data. Either develop to similar depth as PDAC or clearly frame as a "proof of concept" for cross-cancer transfer.

### Minor 5: Data Split Inconsistency

The main methods state the bladder cohort was split "via a 70/30 split" while the supplement states 80/20. Resolve this discrepancy.

### Minor 6: Prose Issues

- The abstract uses "deSurv" (lowercase d) while the rest of the paper uses "DeSurv". Standardize.
- The phrase "the quintessential challenge of rank determination" is unnecessarily flowery. Consider "the well-known challenge of rank selection."
- Several sentences in the Results section are unnecessarily long (>40 words). Consider breaking them up for clarity.
- The Significance Statement is well-written but could be more concrete. Name specific outcomes: "DeSurv identifies three PDAC programs---classical tumor, basal-like, and activated microenvironment---that stratify survival across independent cohorts."

### Minor 7: Missing Notation Table

The supplement should include a notation table defining all symbols at the outset, given the number of parameters (W, H, beta, alpha, lambda, lambda_H, xi, epsilon_H, epsilon_W, delta, k, n_top, etc.).

### Minor 8: Incomplete Supplement

- The supplement text ends with an apparent placeholder or incomplete section
- The consensus initialization procedure is referenced but not fully described
- Signature truncation details are promised but not found
- GitHub repository references should appear in the supplement

---

## V. Technical Concerns

### Technical 1: Convergence Criterion Mismatch

The supplement (Algorithm S1) describes convergence as: eps = |lossNew - loss| / |loss|. The actual code (functions.cpp:649-660) uses cosine similarity: eps = max(1 - cos(W_prev, W), 1 - cos(H_prev, H)). These are fundamentally different criteria. The cosine criterion measures parameter stability, while the loss criterion measures objective progress. The paper must document which criterion is actually used.

### Technical 2: Beta Backtracking Logic

The beta update includes a backtracking step (functions.cpp:472-473) that triggers when surv_loss_candidate < surv_loss_prev. Since surv_loss is computed as loglik * 2/n_event - lambda * penalty_beta, and loglik (Cox partial log-likelihood) is maximized, higher values of surv_loss are better. The condition checks if the candidate is *lower* than the previous value and triggers backtracking---this appears correct (backtracking when the update worsens the objective). However, the paper's existing CODE_REVIEW.md flagged this as potentially inverted. The authors should verify and document the sign convention explicitly.

### Technical 3: Convergence Proof Gap (Coxnet Approximation)

The paper acknowledges (Remark on "Coxnet implementation") that the beta update uses an approximate Hessian rather than the exact second-order information assumed in the idealized convergence proof. While the paper correctly states this is a gap between theory and practice, a mathematical reviewer may ask whether empirical convergence diagnostics (e.g., loss trajectories across iterations) support the claim that the approximation does not materially affect convergence behavior. Consider adding such a diagnostic to the supplement.

### Technical 4: W Update Clamping Not Documented

The code clamps the multiplicative update ratio R to [0.2, 5.0] before applying it to W (functions.cpp:209-211). This is not documented anywhere in the paper or supplement. While likely a numerical stability measure, it changes the theoretical update and should be acknowledged. A reviewer who checks the code will find this.

### Technical 5: Potential Code Bug in Metrics Computation

CODE_REVIEW.md documents a critical bug in R/compute_metrics.R:9 where sdZ is incorrectly assigned from fit$meanZ instead of fit$sdZ. If this bug was present when the reported results were generated, all C-index calculations and survival metrics in the paper may be incorrect. **This must be verified and, if confirmed, all analyses must be re-run.** This is potentially the most serious issue in the manuscript.

---

## VI. Prior Art and Novelty Assessment

### Competing Methods

1. **Huang et al. (2020), CoxNMF** (arXiv:2008.03776): Combines Cox PH regression with NMF. This preprint remains unpublished after 5+ years, suggesting it may not have passed peer review. The DeSurv paper correctly notes that CoxNMF links survival to sample-level factor scores, while DeSurv acts on gene-level programs. However, the precise technical difference (supervision through H vs. W) deserves clearer exposition.

2. **Le Goff et al. (2025), SurvNMF** (PhD thesis, Institut Pasteur Paris / CEA): A survival-supervised NMF formulation from a doctoral dissertation. As an unpublished thesis, its methodological details and reproducibility are unverified. The DeSurv paper's claim that both prior methods are "unpublished and unreviewed" is factually accurate.

3. **DECODER (Peng et al., 2019, Nature Communications)**: Unsupervised compartment deconvolution. DeSurv correctly positions itself as distinct: DECODER discovers compartments without survival supervision, while DeSurv integrates outcomes during discovery. However, the paper should acknowledge that DECODER's compartment ratios are prognostically informative, demonstrating that unsupervised methods can capture survival-relevant structure.

4. **Spectra (Kotliar et al., 2024, Nature Biotechnology)**: Supervised gene program discovery from single-cell data that incorporates prior biological knowledge. While designed for single-cell data (not bulk transcriptomics with survival), this is a conceptually related approach to supervised factorization and should be discussed.

### Novelty Assessment

The novelty claim is sound but needs sharper articulation. The core contribution is NOT "supervised NMF" (which exists) or "NMF + Cox" (which exists). It is:

> **Survival supervision applied at the gene-program level (W) reorganizes latent transcriptional structure in ways that improve interpretability, stability, and cross-cohort generalization---properties not achieved by supervision through sample loadings (H).**

This should be stated explicitly and prominently in the Introduction.

---

## VII. PNAS Suitability Assessment

### Fit for PNAS

PNAS values conceptual reframing over incremental benchmarking. This paper's core insight---that variance-dominant transcriptional structure != prognostically relevant structure, and that W-level supervision reorganizes programs accordingly---is the kind of structural understanding that PNAS publishes.

The cross-cancer transferability result (PDAC to bladder) provides the generalizable insight that PNAS requires. This is not a single-disease observation.

### Concerns for PNAS

1. **Scope:** The paper reads more like a methods paper (suitable for Bioinformatics, JRSS-C, or Biostatistics) than a PNAS paper. To fit PNAS, the framing must emphasize the biological insight (variance != prognosis; survival supervision reveals shared cross-cancer structure) over the methodological machinery.

2. **Breadth of impact:** The paper would be strengthened by discussing implications beyond cancer genomics. Could DeSurv be applied to other high-dimensional survival settings (e.g., electronic health records, financial default prediction)? One paragraph would suffice.

3. **Narrative clarity:** The most important figure (4C: variance vs. survival) is currently the third panel in the fourth figure. Consider making this insight more prominent in the narrative.

---

## VIII. Figures Assessment Summary

| Figure | Rating | Critical Issue | Required Action |
|--------|--------|---------------|-----------------|
| **Fig 1** | 8/10 | Good schematic | Add "Supervision acts on W" callout |
| **Fig 2** | 6/10 | Layout inverts narrative | Reorder panels; add selected k annotation |
| **Fig 3** | 6/10 | Missing context | Add true k, n_replicates; remove "R0_easy" label |
| **Fig 4** | 7/10 | Missing factor labels | Add biological factor labels (critical) |
| **Fig 5** | 5/10 | Missing effect sizes | Add HR, CI, p-value to KM plots (critical) |
| **Fig 6** | 5/10 | Missing effect sizes | Add HR, p-value; identify transferred factor |

---

## IX. Questions for Authors

1. Has the sdZ bug in compute_metrics.R been fixed, and have all results been regenerated with the corrected code? If not, what is the expected impact on reported metrics?

2. What is the actual convergence criterion used in the reported experiments---relative loss change or cosine similarity of parameter matrices?

3. Why was the W update ratio clamped to [0.2, 5.0] in the implementation? Was a sensitivity analysis performed on these bounds?

4. What are the actual C-index values (mean +/- SE) for DeSurv vs. NMF in the PDAC training and validation cohorts?

5. What biologically does the transferred PDAC factor in Figure 6B represent? Is it the classical/basal-like axis, the activated TME factor, or something else?

6. Were the missing functions collect_bo_diagnostics() and compute_bo_eval_se() (referenced in bo_helpers.R but never defined) used in generating the reported results? Does the one-standard-error rule for k-selection function correctly without them?

7. What are the runtimes for DeSurv vs. standard NMF on the PDAC dataset?

---

## X. Action Checklist (Prioritized)

### Blocking (Must Fix for Resubmission)

- [ ] Fix PNAS formatting: reduce figures to 4, references to ~50
- [ ] Complete all placeholder sections (keywords, author contributions, conflicts, acknowledgements, affiliations)
- [ ] Standardize "DeSurv" naming throughout (abstract uses "deSurv")
- [ ] Verify and fix compute_metrics.R sdZ bug; re-run analyses if affected
- [ ] Add HR, CI, p-values to all Kaplan-Meier plots
- [ ] Expand Discussion to 3-4 paragraphs with limitations, related work, clinical implications
- [ ] Add at least one comparator method (LASSO-Cox recommended)
- [ ] Add null and harder simulation scenarios
- [ ] Report all key quantitative metrics numerically (C-index values, effect sizes)
- [ ] Resolve convergence criterion mismatch between paper and code
- [ ] Resolve data split inconsistency (70/30 vs 80/20 for bladder)

### High Priority (Strengthen for PNAS)

- [ ] Expand Methods section with essential algorithmic details
- [ ] Add factor labels to Figure 4 heatmaps
- [ ] Sharpen W-vs-H novelty claim in Introduction (first page)
- [ ] Add explicit out-of-sample statement (nested CV, held-out validation)
- [ ] Move projection advantage (Z = W^T X) argument to main text
- [ ] Add simulation parameter annotations to Figure 3
- [ ] Document W update clamping and other code-paper discrepancies
- [ ] Fix symbol overloading (delta, theta)
- [ ] Add notation table to supplement

### Medium Priority (Improve Quality)

- [ ] Restructure Figure 2 layout (problem before solution)
- [ ] Add sensitivity analysis showing nearby k values give similar structure
- [ ] Describe consensus initialization procedure completely
- [ ] Add dataset characteristics table (sample sizes, platforms, censoring rates)
- [ ] Discuss relationship to Spectra and supervised topic models
- [ ] Add runtime comparison

---

## XI. Concluding Assessment

This manuscript presents a conceptually meaningful framework with genuine methodological innovation. The core insight---that gene-level survival supervision reorganizes transcriptional latent structure toward outcome-relevant biology---is well-supported by the presented analyses and represents an advance over existing survival-NMF approaches. The cross-cancer transfer result adds generalizability that elevates the work beyond a single-application methods paper.

However, in its current form, the manuscript has critical formatting deficiencies, incomplete reporting of quantitative results, an underdeveloped Discussion, insufficient comparator methods, and several paper-code inconsistencies that undermine confidence in the reported results. The potential code bug in compute_metrics.R is particularly concerning and must be resolved with certainty.

With the revisions outlined above---particularly PNAS compliance, expanded Discussion, additional comparators, broader simulations, and complete quantitative reporting---this manuscript would be suitable for publication in PNAS.

**Recommendation: Major Revision**

---

*Review conducted: February 10, 2026*
*This review synthesizes analysis of the main manuscript (11 pages), supplementary methods (10 pages), source code (R and C++), and existing internal review documents.*
