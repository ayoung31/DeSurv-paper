# GitHub Issues: Summary, Validity, and Proposed Solutions

Generated: 2026-02-19
Repo: ayoung31/DeSurv-paper
Open issues: 19 (all opened by @naimurashid)

---

## Table of Contents

| # | Title | Category | Validity |
|---|-------|----------|---------|
| [#20](#issue-20) | Add HR, CI, and p-values to KM plots; add factor labels to heatmaps | Figure | Partially valid |
| [#19](#issue-19) | Integrate K-sensitivity analysis into paper | Writing | Valid |
| [#18](#issue-18) | Complete scRNA-seq section in supplement | Writing | Valid — confirmed |
| [#17](#issue-17) | Reduce figures from 6 to 4 for PNAS compliance | Compliance | Valid — confirmed |
| [#14](#issue-14) | Methods gap: gene overlap heatmap methodology entirely undescribed | Methods | Valid — confirmed |
| [#13](#issue-13) | Methods gap: variance-vs-survival quantification undescribed | Methods | Valid + code/text discrepancy found |
| [#12](#issue-12) | Methods gap: factor selection criterion for external validation undescribed | Methods | Valid — confirmed |
| [#11](#issue-11) | 1-SE rule for K selection: should it also consider ntop? | Analysis | Partially valid |
| [#10](#issue-10) | Reframe DeSurv advantage: parsimony and interpretability over raw prediction | Framing | Valid — supported by code |
| [#9](#issue-9) | Standard NMF at K=5 outperforms DeSurv K=3 in external validation | Analysis | Valid — needs store verification |
| [#8](#issue-8) | Clarify K selection for standard NMF in external validation | Methods | Valid — confirmed in code |
| [#7](#issue-7) | Connect variance-vs-survival plot (Fig 4C) back to simulation scenarios | Writing | Partially valid — already partially done |
| [#6](#issue-6) | Verify claim: exocrine expression absent from DeSurv factors | Analysis | Valid — needs quantification |
| [#5](#issue-5) | Characterize standard NMF factors at independently-selected K | Analysis | Valid |
| [#4](#issue-4) | Clarify K selection for standard NMF comparison (Section 4) | Methods | Valid — confirmed in code; already partially addressed |
| [#3](#issue-3) | Verify claim: exocrine content is dominant variance but not dominant prognostic | Analysis | Valid — TODO placeholders confirm |
| [#2](#issue-2) | Mixed scenario: investigate DeSurv separation of marker vs background genes | Analysis | Valid — needs investigation |
| [#1](#issue-1) | Clarify survClust model and how it determines survival-specific clusters | Research | Valid — research task |

---

## Issue #20

**Add HR, CI, and p-values to KM plots; add factor labels to heatmaps**

### Summary

Two figure annotation improvements are requested before submission: (1) add HR, 95% CI, and log-rank p-value to each KM panel in Figs 5B-C and 6B; (2) add biological factor labels to the gene overlap heatmaps (Figs 4A-B) using conservative descriptors (e.g., "Classical tumor-associated").

### Validity Assessment

**Partially valid.** Code inspection reveals the following:

- **HR and CI are already present.** `splot_median()` in `R/figure_targets.R:2415-2436` already constructs an HR label and places it as text on the KM plot: `sprintf("HR (High vs Low) = %.2f (95%% CI %.2f-%.2f)", ...)`.
- **Log-rank p-value is NOT added** in `splot_median`. The `ggsurvplot` call does not include `pval = TRUE`. A different function (the 2-factor plot at line 2213-2220) does include `pval = TRUE`, but that is for a separate cross-factor panel.
- **Factor labels on heatmaps are missing.** The `make_gene_overlap_heatmap` function labels columns by factor index (1, 2, 3), not by biological interpretation.

### Proposed Solution

**P-value on KM plots** — modify `splot_median()` in `R/figure_targets.R` (~line 2422):
```r
# Add pval = TRUE and combine it with the existing HR label
splot = ggsurvplot(sfit, data=df, risk.table = TRUE,
                   xlab = "Time (months)",
                   palette = c("violetred2","turquoise4"),
                   break.time.by = 25,
                   legend.labs=c('Low (< median)','High (≥ median)'),
                   pval = TRUE,
                   pval.coord = c(max(splot$data.survtable$time, na.rm=TRUE) * 0.7, 0.75),
                   risk.table.y.text=FALSE,
                   censor.size=2,
                   font.main=12)
```
Or compute and append manually:
```r
lr_test <- survdiff(Surv(time, event) ~ factor, data = df)
p_val <- 1 - pchisq(lr_test$chisq, df = 1)
p_label <- sprintf("Log-rank p = %.3f", p_val)
# Then place both hr_label and p_label as annotations
```

**Biological factor labels on heatmaps** — add a `factor_labels` argument to `make_gene_overlap_heatmap()` that renames columns before plotting. Based on the current biological interpretation:
- Factor 1: "Classical/Immune"
- Factor 2: "Microenvironmental"
- Factor 3: "Basal-like"
Use conservative labels as the issue requests. These should be passed in from the calling target so they are not hardcoded in the figure function.

---

## Issue #19

**Integrate K-sensitivity analysis into paper**

### Summary

A K-sensitivity synthesis document exists at `docs/plans/2026-02-17-k-sensitivity-synthesis.md`. Its key findings (iCAF program is alpha-dependent, not just K-dependent; K=3 at α≥0.35 uniquely produces coherent iCAF factor; at α=0 no K recovers iCAF) need to be woven into the supplement with a brief main-text reference.

### Validity Assessment

**Valid.** The synthesis document exists and contains completed analysis. The paper currently has no mention of alpha-sensitivity of the iCAF program. This is important because: (a) it strengthens the claim that α>0 is necessary, not just K=3; (b) it pre-empts reviewer questions about whether K=3 is the only reason iCAF is recovered.

### Proposed Solution

1. **Supplement section:** Add a new subsection "Sensitivity of Factor Structure to Supervision Strength (α)" in `paper/supplement.Rmd` after the existing K-selection heuristics section. Include:
   - The per-K × α grid figures (panels already generated in `docs/plans/`)
   - A 2-3 sentence summary: "Across K values tested, the iCAF-like activated microenvironmental program was recovered only when α ≥ 0.35, regardless of K. At α=0 (standard NMF), no value of K produced a factor with significant iCAF enrichment. This indicates that survival supervision, not rank alone, is the key driver of iCAF program recovery."

2. **Main text reference:** In Section 4 (factor structure), after describing the iCAF factor, add one sentence: "Sensitivity analyses across the joint (K, α) space confirmed that this iCAF-like program requires α ≥ 0.35 and is absent under standard NMF at all tested ranks (SI Appendix, Fig. SX)."

3. **Targets:** Check whether `fig_cv_grid_*` targets for the α-sensitivity analysis are already in the pipeline and just need to be loaded in the supplement Rmd.

---

## Issue #18

**Complete scRNA-seq section in supplement**

### Summary

`paper/supplement.Rmd:203` contains the placeholder "Need more here..." in the scRNA-seq / Elyada VAM analysis section. The text describes factor 1 (iCAF/B cells), factor 2 (acinar/classical), and factor 3 (basal-like) but cuts off before completing the interpretation.

### Validity Assessment

**Valid — confirmed.** The placeholder is confirmed at `supplement.Rmd:203`:
> "...Interestingly, we do not see a large overlap of factor 1 with the classical PDAC cell types despite some overlap with classical gene lists. Need more here..."

### Proposed Solution

Replace the placeholder with 2-4 sentences completing the scRNA-seq interpretation. Suggested content based on existing results:

> "Together, the scRNA-seq VAM enrichment analysis confirms that DeSurv's three learned factors correspond to biologically coherent cell-type programs identifiable at single-cell resolution. Factor 1's enrichment in iCAF and B cells is consistent with an activated immune-stromal axis, while its lack of enrichment in classical PDAC cells, despite partial overlap with classical gene lists in the bulk data, suggests that bulk classical signal in factor 1 is driven by co-occurring immune infiltrates rather than tumor-intrinsic classical identity. Factor 3's exclusive enrichment in basal-like PDAC cells provides single-cell validation of its role as a tumor-intrinsic aggressiveness program. These results support the interpretation that survival supervision in DeSurv recovers programs corresponding to distinct, co-occurring cell populations rather than composite bulk expression signatures."

---

## Issue #17

**Reduce figures from 6 to 4 for PNAS compliance**

### Summary

PNAS allows a maximum of 4 main-text figures. The current paper has 6: `fig-schema`, `fig-bo`, `fig-sim`, `fig-bio`, `fig-extval`, `fig-bladder`. The recommendation is to move `fig-bo` (rank selection details) and `fig-bladder` (cross-cancer analysis) to SI.

### Validity Assessment

**Valid — confirmed.** Inspection of `paper/04_results_REVISED.Rmd` confirms exactly 6 figures labeled `\label{fig:schema}`, `\label{fig:bo}`, `\label{fig:sim}`, `\label{fig:bio}`, `\label{fig:extval}`, `\label{fig:bladder}`. Additional PNAS compliance issues noted in the issue (references >50, incomplete placeholders) are also confirmed by the bibliography and paper scaffolding.

### Proposed Solution

**Figure reduction:**
1. Move `fig-bo` (rank selection heuristics + BO heatmap + simulation K-selection) to SI Appendix as "Fig. S1." In the main text, replace with: "We evaluated rank selection using reconstruction residuals, cophenetic correlation, and mean silhouette width — which gave inconsistent guidance — and instead selected rank and supervision strength jointly via Bayesian optimization over cross-validated C-index (SI Appendix, Fig. S1)."
2. Move `fig-bladder` to SI Appendix as "Fig. S2." In the main text (Section 6), reference with: "We further demonstrated generalizability by projecting PDAC-learned programs into bladder cancer (SI Appendix, Fig. S2), recovering a single survival-associated factor with consistent behavior."
3. Keep `fig-schema`, `fig-sim`, `fig-bio`, `fig-extval` as main-text Figs 1-4.

**Reference reduction (71 → ≤50):** Prioritize removing:
- Methodological background references that are non-essential (keep only the most cited/canonical)
- Secondary validation dataset citations that can be merged
- Duplicate concepts cited twice

**Placeholders to complete:**
- Keywords (add 5-6 relevant terms: NMF, survival analysis, PDAC, deconvolution, Bayesian optimization, prognostic biomarkers)
- Author contributions (describe specific contributions per PNAS format)
- Conflict of interest statement (standard disclosure)
- Acknowledgements (funding sources, compute resources)

---

## Issue #14

**Methods gap: gene overlap heatmap methodology (Fig 4A-B) entirely undescribed**

### Summary

Figures 4A-B show Spearman correlations between DeSurv/NMF factor gene rankings and established PDAC gene programs. The methodology (which programs were used, how correlations were computed, filtering thresholds) is entirely absent from the manuscript.

### Validity Assessment

**Valid — confirmed.** Code inspection of `R/figure_targets.R:650-743` confirms all the undescribed choices:
- Top 50 genes per factor used (`tops = tops[1:50,]` at line 678)
- Spearman correlation between continuous W values and binary set membership (`cor(..., method = "spearman")` at line 707)
- BH multiple testing correction per factor column (line 710)
- Filter: only programs with r > 0.2 shown (line 714)
- deCAF gene lists (`proCAF`, `restCAF`) are hardcoded in the figure function at lines 659-660
- Reference programs from Moffitt, SCISSORS, Elyada, Bailey, DECODER (some excluded)
- No citations for the gene lists in the figure code

The SI Survival Analysis subsection (`supp_methods.Rmd:815-827`) makes no mention of this analysis at all.

### Proposed Solution

Add a subsection to `paper/supp_methods.Rmd` titled "Gene Program Correlation Analysis" after the Survival Analysis subsection:

> **Gene Program Correlation Analysis.** To assess the biological identity of learned factors, we evaluated the correspondence between factor-specific gene rankings and established PDAC gene programs. For each factor column $j$ of the gene program matrix $W$, we extracted the top 50 genes ranked by loading magnitude. We then computed Spearman rank correlations between the continuous $W_j$ loadings and binary indicator vectors for each reference gene set $G_k$ (1 if gene $\in G_k$, 0 otherwise). P-values were computed using the asymptotic approximation to the Spearman test statistic and adjusted for multiple comparisons using the Benjamini–Hochberg procedure within each factor. Gene programs with maximum absolute correlation below 0.2 across all factors were excluded from visualization. Asterisks indicate adjusted p-value < 0.1. Reference gene programs included: classical and basal-like signatures from Moffitt et al. [CITE]; SCISSORS subtypes [CITE]; iCAF and myCAF signatures from Elyada et al. [CITE]; proCAF and restCAF signatures from [CITE deCAF paper]; and other established PDAC subtypes [CITE]. Programs from Bailey et al., periductal subtypes, DECODER, MSI subtypes, and PurIST were excluded from visualization due to redundancy with other included programs or low overall correlation.

Also add citations for the deCAF gene lists, which are currently hardcoded without reference.

---

## Issue #13

**Methods gap: variance-vs-survival quantification (Fig 4C) undescribed**

### Summary

Fig 4C shows variance explained vs. survival contribution per factor. Neither the variance explained computation nor the survival contribution metric is described in the manuscript or supplement.

### Validity Assessment

**Valid — and a code/text discrepancy was found.** Code inspection:

- **Variance explained** (`R/figure_targets.R:2266`): `compute_variance_explained()` computes `sum((W_j * H_j)^2) / sum(X^2)` — fraction of total sum-of-squares from each factor's rank-1 component. This is correct and standard.

- **Survival contribution** (`R/figure_targets.R:2320-2345`): The code comment says "Type III partial likelihood contribution" — this is a **leave-one-factor-out** approach (full K-factor model minus (K-1)-factor model excluding factor j). However, the paper text in `04_results_REVISED.Rmd:129` says "change in partial log-likelihood from **univariate Cox models**." These are different things. The code is doing the leave-one-out approach; the paper is describing it as univariate. **This discrepancy should be corrected.**

### Proposed Solution

**Fix the paper text discrepancy first.** In `04_results_REVISED.Rmd:129`, change the Fig 4C caption description from "change in partial log-likelihood from univariate Cox models" to "Type III partial likelihood contribution (change in partial log-likelihood from the full K-factor Cox model when each factor is removed)."

**Then add a supplement subsection** in `supp_methods.Rmd`:

> **Variance and Survival Contribution per Factor.** For each factor $j \in \{1, \ldots, k\}$, we quantified two quantities to assess the alignment between expression variance and prognostic relevance. The fraction of expression variance explained by factor $j$ was computed as:
> $$V_j = \frac{\| W_j H_j \|_F^2}{\| X \|_F^2}$$
> where $W_j$ is the $j$-th column of $W$ and $H_j$ is the $j$-th row of $H$, so $W_j H_j$ is the rank-1 approximation contributed by factor $j$. The survival contribution of factor $j$ was quantified as the Type III partial likelihood increment: the reduction in Cox partial log-likelihood when factor $j$ is omitted from the full $k$-factor model:
> $$\Delta\ell_j = \ell(W_1, \ldots, W_k) - \ell(W_1, \ldots, W_{j-1}, W_{j+1}, \ldots, W_k)$$
> where $\ell(\cdot)$ denotes the Cox partial log-likelihood. Larger $\Delta\ell_j$ indicates a greater marginal contribution of factor $j$ to survival prediction given all other factors.

---

## Issue #12

**Methods gap: factor selection criterion for external validation undescribed**

### Summary

Section 5 prose mentions "the factor showing the largest increase in Cox model log partial likelihood in the training data" for selecting which of the K factors to use in external validation, but no methods section describes this criterion.

### Validity Assessment

**Valid — confirmed.** The SI "Survival Analysis" subsection (`supp_methods.Rmd:815-827`) describes median splits of the full linear predictor $\hat{\eta}$ but does not describe the per-factor selection step. The code evidence is clear: `which.max(desurv_var$delta_loglik)` selects the factor in `04_results_REVISED.Rmd:317`, and `compute_survival_explained()` provides the delta log-likelihood values. No methods text supports this.

### Proposed Solution

Add 3 sentences to the "Survival Analysis" subsection in `supp_methods.Rmd` after the existing median-split description:

> **Per-factor selection for external validation.** When evaluating prognostic generalization in external cohorts (Section 5), we focused on the single factor most associated with survival in the training data rather than the full linear predictor. Specifically, we selected the factor $j^* = \arg\max_j \Delta\ell_j$, where $\Delta\ell_j$ is the Type III partial likelihood contribution described above (see Gene Program Analysis section), computed on the training cohort (TCGA and CPTAC). This factor was selected prior to evaluating any validation cohort, and the same factor index was used for both DeSurv and standard NMF comparisons in external validation.

---

## Issue #11

**1-SE rule for K selection: should it also consider ntop? K=3 barely passes threshold**

### Summary

DeSurv selects K=3 over K=7 using a 1-SE rule (smallest K within 1 SE of the maximum C-index). K=3 passes the threshold by only 0.0025 C-index units. The issue raises whether the 1-SE rule should also penalize ntop (genes per factor), since K=3 uses ntop=299 — the highest of all K values — making its total gene count (K×ntop=897) comparable to K=5 (925).

### Validity Assessment

**Partially valid.** The 1-SE rule itself is implemented correctly in `R/bo_helpers.R:125-141` (`threshold <- best_mean - best_se`, `k_selected <- min(candidates$k_numeric)`). The concern about robustness is legitimate — a 0.0025 C-index margin is small — but the 1-SE rule is a well-established principled method from the regularization literature (Breiman et al.) and its intent is exactly to prefer simpler models with modest performance loss. Whether to extend it to ntop is a judgment call, not a bug.

The deeper point — that the parsimony argument depends on the 1-SE rule, not on an inherent property of survival supervision — is valid and should be acknowledged in the paper.

### Proposed Solution

1. **Don't modify the 1-SE rule** for ntop. The existing rule is defensible and principled. Adding K×ntop would require justifying why gene count is an appropriate complexity measure, and the two quantities have different units (K is factorization rank; ntop is a gene-selection filter, not a model parameter).

2. **Add robustness language to the paper.** In the Methods or SI, add: "The 1-SE rule selected K=3, which achieved a mean cross-validated C-index of 0.647, compared to the maximum of 0.655 at K=7. The margin (0.008 C-index units) was within one standard error of the maximum (SE=0.011), meeting the threshold for parsimony selection. Sensitivity analyses across bootstrap resamples confirmed that K=3 or K=4 was selected in the majority of replicates [if this analysis exists; otherwise flag for investigation]."

3. **Reframe the parsimony narrative** as suggested in Issue #10 (see below): "DeSurv's 1-SE rank selection criterion prefers the most parsimonious model whose performance falls within one standard error of the maximum, yielding K=3 for PDAC. This choice reflects a deliberate parsimony preference, not a claim that K=3 achieves strictly higher discriminative performance than larger K."

---

## Issue #10

**Reframe DeSurv advantage: parsimony and interpretability over raw prediction**

### Summary

Stored validation results show standard NMF at K=5 outperforms DeSurv K=3 on 4/5 external cohorts by C-index. The paper's current framing implies predictive superiority, which is not supported. The correct framing is parsimony: DeSurv achieves with 3 biologically interpretable factors what NMF needs 5-7 factors to achieve numerically.

### Validity Assessment

**Valid — supported by the code architecture and the data in Issues #9 and #11.** The code confirms that the Section 5 comparison forces NMF to K=3 (`fit_std_desurvk`), not its independently selected K. The C-index numbers provided (DeSurv K=3 vs. NMF K=5) are from a specific store and should be verified against the canonical store, but the directional finding is credible.

### Proposed Solution

**Revise the Section 5 framing** in `04_results_REVISED.Rmd:162-164`:

Change:
> "the NMF factor identified by the same criterion showed greater heterogeneity across datasets and weaker survival associations, consistent with the expectation that variance-driven programs capture cohort-specific variation that does not transfer"

To something like:
> "When both methods were constrained to the same factorization rank (K=3, selected by DeSurv's Bayesian optimization), the NMF factor showed greater heterogeneity across datasets and weaker survival associations. Notably, standard NMF at its independently selected rank (K=5 by elbow criterion) achieved comparable discriminative performance (SI Appendix, Table SX), but required distributing survival-relevant signal across more factors — each with less distinct biological identity. DeSurv's advantage is therefore one of parsimony and interpretability: it concentrates prognostic signal into fewer, more reproducible factors, rather than strictly improving point estimates of discriminative performance."

**Revise the Discussion** to explicitly discuss this parsimony-vs-prediction tradeoff as the primary contribution.

**Add SI Table:** Report C-index values for DeSurv K=3, NMF K=3, NMF K=5, DeSurv α=0 K=7 across all five validation cohorts to allow transparent comparison.

---

## Issue #9

**Standard NMF at K=5 outperforms DeSurv K=3 in external validation**

### Summary

From the `naimedits0125_full` store, NMF at elbow-selected K=5 outperforms DeSurv K=3 on 4/5 external cohorts by C-index. The current Section 5 comparison uses NMF at K=3 (DeSurv's K), not NMF at its own K. This overstates DeSurv's predictive advantage.

### Validity Assessment

**Valid — needs canonical store verification.** The C-index values reported (NMF K=5: 0.631, 0.549, 0.654, 0.690, 0.645 vs. DeSurv K=3: 0.607, 0.552, 0.628, 0.598, 0.623) come from the `naimedits0125_full` store. These need to be confirmed against the canonical store (`store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main`). However, the directional finding — NMF at own K outperforms NMF forced to DeSurv's K — is expected and credible regardless of store.

### Proposed Solution

This issue is closely related to #10 and #8. The combined solution:

1. **Verify numbers** from the canonical store by loading `val_predictions_std_elbowk` targets and computing C-index.
2. **Add a supplementary table** comparing DeSurv K=3, NMF K=3, NMF K=5, and DeSurv α=0 K=7 across all validation cohorts.
3. **Reframe Section 5** (see Issue #10 solution). The key message is equivalent performance with fewer factors, not raw superiority.
4. **Add forest plot and KM panel for NMF K=5** in SI to show the comparison transparently. The pipeline already has `val_predictions_std_elbowk` and `val_latent_std_elbowk` — check whether `fig_gene_overlap_heatmap_std_elbowk_tcgacptac` and corresponding HR/KM targets are generated.

---

## Issue #8

**Clarify K selection for standard NMF in external validation (Section 5, Fig 5)**

### Summary

Fig 5 compares DeSurv and standard NMF in external validation, but the standard NMF comparison uses DeSurv's K=3, not NMF's independently selected K. The code confirms this: `nmf_df = compute_hrs(data_val_filtered, fit_std_desurvk, "NMF")` at `targets_common_pipeline.R:2518`.

### Validity Assessment

**Valid — confirmed in code.** `targets_common_pipeline.R:2518` explicitly uses `fit_std_desurvk` (NMF at DeSurv's K) for the forest plot. Section 5 does not disclose this. The revised Section 4 text already discloses the same K constraint: "To enable a direct comparison of how each method organizes the same number of factors, we also fit standard NMF at the same rank (k=3)" (`04_results_REVISED.Rmd:120`). Section 5 needs the same disclosure.

### Proposed Solution

**In the Fig 5 caption** (`04_results_REVISED.Rmd:166`), add after "both at factorization rank k=3": "...both at factorization rank $k=3$ (DeSurv's Bayesian-optimization-selected rank); external validation of standard NMF at independently selected ranks is shown in SI Appendix."

**In Section 5 prose** (`04_results_REVISED.Rmd:162`), add a disclosure sentence after "For each method, we focused on the factor showing the largest increase in Cox model log partial likelihood in the training data.": "For consistency with the factor structure comparison in Section 4, standard NMF was evaluated at the same rank as DeSurv ($k=3$); results at independently selected ranks ($k=5$, $k=7$) are shown in SI Appendix."

**Generate SI figures** for NMF at K=5: forest plot and KM curves using `fit_std_elbowk` with corresponding validation targets.

---

## Issue #7

**Connect variance-vs-survival plot (Fig 4C) back to simulation scenarios**

### Summary

Fig 4C (variance vs. survival per factor) should be connected to the simulation scenarios. The question is whether the PDAC data pattern (high-variance exocrine factor has low survival association, others have partial overlap) corresponds to the mixed simulation scenario.

### Validity Assessment

**Partially valid — already partially addressed.** The revised `04_results_REVISED.Rmd:126` already contains a sentence making this connection: "The observed pattern in PDAC, where the highest-variance factor shows minimal survival association while other factors show partial overlap between variance and prognosis, is consistent with the intermediate regime between the primary and mixed simulation scenarios, where DeSurv's advantage was largest." The issue's main request is thus already incorporated. What remains is whether an analogous simulation variance-survival plot exists.

### Proposed Solution

1. **The connection sentence is already in the paper** — confirm it is present in the final draft and not accidentally removed during revision.

2. **Optional supplement analysis:** Generate a variance-vs-survival scatter plot for simulated data under the primary and mixed scenarios to visually demonstrate the spectrum. This would be a strong addition but is not strictly necessary given the existing text.

3. **Address the specific questions** posed in the issue:
   - Q1 (mixed scenario variance-survival relationship): Could be added to SI by generating `build_variance_survival_df` on simulation outputs if available.
   - Q3 (where PDAC falls on primary-to-mixed spectrum): Could be quantified by computing the correlation between factor variance rank and factor survival rank across factors; if Spearman rank correlation between $V_j$ and $\Delta\ell_j$ is strongly negative, PDAC resembles the primary scenario.

---

## Issue #6

**Verify claim: exocrine expression absent from DeSurv factors**

### Summary

Section 4 claims: "exocrine-associated expression did not dominate any DeSurv factor." This needs rigorous quantitative support, including defining "exocrine" operationally, checking whether top-gene truncation masks residual signal, and ensuring the language "did not dominate" is defensible.

### Validity Assessment

**Valid — needs quantification.** The claim is in `04_results_REVISED.Rmd:122` and is supported qualitatively by the gene overlap heatmap (no DeSurv factor should have high Spearman correlation with exocrine gene programs in Fig 4A). However, the issue correctly notes: (a) "exocrine" needs to be operationally defined with a specific gene list; (b) the claim could be weakened by full-W (not just top-50) analysis; (c) α=0.33 means 67% reconstruction weight, so some exocrine signal likely remains.

### Proposed Solution

1. **Define "exocrine" operationally.** Use the exocrine gene signature from the Moffitt et al. or Bailey et al. PDAC gene lists already used in the heatmap analysis. Document which list is used.

2. **Compute exocrine enrichment in full W matrix** (not just top 50 genes): Run gene set enrichment analysis (GSEA) or compute Spearman correlation of full W column with exocrine gene indicator. If the exocrine signature is not significantly enriched (adj. p > 0.1) in any DeSurv factor at the full-W level, the claim is well-supported.

3. **Compare exocrine enrichment between DeSurv and NMF.** The claim is implicitly comparative: exocrine dominates an NMF factor but not a DeSurv factor. A formal comparison (e.g., max Spearman correlation with exocrine programs: NMF factor 1 = 0.45, DeSurv max = 0.12) makes this concrete and defensible.

4. **Soften language if needed.** If exocrine signal is present but not dominant, use "did not dominate" (current) or "was substantially attenuated" rather than implying complete absence.

---

## Issue #5

**Characterize standard NMF factors at independently-selected K (K=5 and K=7) for supplement**

### Summary

The main text compares DeSurv (K=3) to NMF at the same K=3. The supplement should show NMF at its independently selected K (K=5 via elbow, K=7 via BO at α=0), characterizing the additional factors and the argument for DeSurv's parsimony.

### Validity Assessment

**Valid.** The pipeline already has `fit_std_elbowk_tcgacptac` (K=5) and `fig_gene_overlap_heatmap_std_elbowk_tcgacptac`. The `tar_params_best_alpha0_tcgacptac` target provides K=7 BO results at α=0. These targets appear to exist but their outputs are not yet in the supplement.

### Proposed Solution

1. **Add to supplement:** A new subsection "Standard NMF at Independently Selected Ranks" in `supplement.Rmd`, showing:
   - `fig_gene_overlap_heatmap_std_elbowk_tcgacptac` (K=5 gene programs)
   - A brief description of the additional factors: "At K=5, standard NMF recovers [describe factors from the figure], including an additional exocrine-related factor and a [describe] factor not present at K=3. At K=7 (BO-selected for NMF at α=0), factors further fragment the [describe] signal. In both cases, multiple factors capture exocrine or composition-driven variance with minimal survival association."

2. **Verify targets are generated:** Check whether `fig_gene_overlap_heatmap_std_elbowk_tcgacptac` exists in the current targets store and has been rendered.

3. **Generate K=7 NMF gene overlap heatmap** if not already present by adding a target using `fit` from `tar_params_best_alpha0_tcgacptac`.

---

## Issue #4

**Clarify K selection for standard NMF comparison (Section 4)**

### Summary

Section 4 compares DeSurv to NMF at K=3, but K=3 was DeSurv's BO-selected K — NMF would independently select K=5 (elbow) or K=7 (BO at α=0). The paper needs to disclose this.

### Validity Assessment

**Valid — confirmed in code; already partially addressed.** `targets_common_pipeline.R:1125` confirms `fit_std_desurvk` uses `bo_bundle_selected$params_best$k` (DeSurv's K). However, `04_results_REVISED.Rmd:120` already contains the disclosure: "To enable a direct comparison of how each method organizes the same number of factors, we also fit standard NMF at the same rank (k=3); results for standard NMF at independently selected ranks (k=5 via elbow detection, k=7 via cross-validated concordance at α=0) are presented in the SI Appendix." This is exactly the right approach.

### Proposed Solution

The main text disclosure is already present. The remaining action items are:

1. **Ensure SI Appendix figures exist** for K=5 and K=7 NMF (see Issue #5).
2. **Verify consistency** between the `04_results_REVISED.Rmd` disclosure and the final paper.Rmd rendering pipeline — confirm that results.Rmd is included as a child document and the disclosure sentence is not accidentally dropped.

---

## Issue #3

**Verify claim: exocrine content is dominant variance source but not dominant prognostic signal**

### Summary

The bold summary in Section 4 states: "The dominant source of expression variance in PDAC, exocrine content, is not the dominant source of prognostic signal." This central claim needs to be verified against actual factor loadings, variance explained values, and survival associations.

### Validity Assessment

**Valid — TODO placeholders confirm numbers are missing.** `04_results_REVISED.Rmd:126` contains `[TODO]` placeholders for the actual numbers: "the highest-variance NMF factor explained [TODO]% of expression variance but contributed Δℓ=[TODO] to survival, whereas the DeSurv factor with the largest survival contribution explained [TODO]% of variance with Δℓ=[TODO]." The claim is structurally sound (supported by Fig 4C conceptually) but needs the actual numerical values before submission. The biological claim (exocrine as variance-dominant, prognostically neutral) is consistent with published PDAC literature.

### Proposed Solution

1. **Fill in the [TODO] placeholders** by loading `fig_variation_explained_tcgacptac` from the targets store and extracting the variance explained and delta log-likelihood values for each factor:
   ```r
   tar_load(tar_fit_desurv_tcgacptac)
   # run build_variance_survival_df() to get the numbers
   ```
2. **Confirm NMF factor 1 is the exocrine factor** by cross-referencing with `fig_gene_overlap_heatmap_std_desurvk_tcgacptac` — the factor with highest correlation to exocrine programs should also be the highest-variance factor.
3. **Consult Jen Jen Yeh** on whether the exocrine-as-dominant-variance claim is consistent with the PDAC biology literature (as the issue requests).

---

## Issue #2

**Mixed scenario: investigate whether DeSurv separates marker vs background gene survival effects**

### Summary

The mixed simulation scenario (`R_mixed`) uses both factor-specific marker genes and shared background genes for survival. The current text describes "attenuated" DeSurv advantage. The question is whether DeSurv actively separates the marker gene effect from the background gene effect — a stronger mechanistic claim.

### Validity Assessment

**Valid — needs investigation.** The current text accurately describes the attenuated advantage but does not investigate the mechanistic question. The precision metric in the simulation only measures overlap with marker gene sets; it does not differentiate marker vs. background survival gene effects. This is a substantive analysis question that requires inspecting the simulation outputs.

### Proposed Solution

1. **Analyze learned beta values** from the mixed scenario: In DeSurv's fitted model for mixed scenario replicates, do factors whose W columns have high marker gene loadings also have higher estimated β coefficients? This would demonstrate that survival supervision not only identifies the right factors but also correctly weights marker-gene-driven factors over background-gene-driven ones.

2. **Compute separate precision for marker vs. background genes** in the mixed scenario. Define:
   - Marker precision: fraction of top genes in learned factor that are true marker genes
   - Background precision: fraction of top genes that are background survival genes
   If DeSurv concentrates on marker genes (high marker precision, lower background precision compared to NMF), this supports the stronger claim.

3. **Revise the text** based on findings. If supported: "In the mixed scenario, DeSurv preferentially recovered factor-specific marker gene programs (higher marker precision than standard NMF), suggesting that survival supervision concentrates model capacity on genes with factor-specific — rather than diffuse background — survival effects." If not supported, keep the current attenuated-advantage framing.

4. **Code note:** `R/simulation_functions/scenario_defaults.R:50-70` defines `R_mixed` and `simulate_survival.R` implements the mixed survival gene construction. Review these to confirm the marker/background gene distinction is tracked in simulation outputs.

---

## Issue #1

**Clarify survClust model and how it determines survival-specific clusters**

### Summary

For the revision, a summary of how survClust (Arora et al., 2020) works is needed to position DeSurv relative to it — specifically: what model it uses, how it incorporates survival, how it uses genomic data, and the key methodological difference from DeSurv.

### Validity Assessment

**Valid — research task.** This is a literature/method review task that needs to be completed before the introduction/discussion can accurately position DeSurv relative to prior work.

### Proposed Solution

Based on the survClust literature (Arora et al., 2020, *Genome Medicine*):

**survClust summary:**
- Uses a **log-rank statistic-based objective** rather than a Cox model. It maximizes the log-rank test statistic between clusters, making it a clustering method optimized for survival stratification.
- Operates on multi-omics distance matrices (e.g., Euclidean distance on normalized expression/mutation matrices), then applies weighted kernel k-means clustering where the kernel weights are optimized to maximize the log-rank test.
- Produces **cluster assignments** directly, without learning a gene program matrix (W) or interpretable factor loadings. It does not learn which genes drive the clusters.
- Rank selection is performed via cross-validated log-rank p-value.

**Key methodological difference from DeSurv:**
> "survClust identifies survival-stratifying clusters by optimizing the log-rank statistic over sample-level distance matrices, yielding cluster assignments without interpretable gene programs. DeSurv instead learns a nonnegative factor decomposition (W) that jointly optimizes reconstruction fidelity and Cox partial likelihood, producing gene program matrices that can be biologically interpreted, projected to new datasets, and used for factor-level survival analysis — capabilities survClust's cluster-assignment output does not support."

**Suggested citation/discussion placement:** Introduce survClust in the introduction as a prior approach to survival-supervised clustering, then note DeSurv's distinct approach: factorization-based gene program discovery rather than direct cluster assignment. This preempts reviewer requests for comparison.

---

## Priority Ranking for Action

| Priority | Issue(s) | Rationale |
|----------|----------|-----------|
| **Immediate (blocking submission)** | #17, #3 | PNAS compliance (figure count) and TODO placeholders must be filled |
| **High (methods completeness)** | #14, #13, #12 | Reviewers will flag missing methods; analysis is well-characterized in code |
| **High (framing integrity)** | #10, #9, #8 | Current framing may not survive peer review given NMF K=5 performance |
| **Medium (writing tasks)** | #18, #20, #19 | Writing completeness; #20 partially done (HR/CI present, p-value missing) |
| **Medium (analysis/verify)** | #6, #3 | Verification needed before submission; numbers needed for #3 |
| **Lower (supplement content)** | #4, #5, #11 | Disclosure already present for #4; supplement figures for #5 |
| **Lower (investigation)** | #2, #7 | Important but not blocking; can be addressed in revision |
| **Research task** | #1 | Needed for discussion but not urgent |
