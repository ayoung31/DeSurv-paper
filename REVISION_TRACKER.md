# DeSurv Manuscript: Revision Tracker

**Purpose:** Maps each action item from the PNAS Senior Editor Review to its resolution status in the REVISED manuscript files. Created to track progress toward resubmission.

**Source review:** [PNAS_SENIOR_EDITOR_REVIEW.md](PNAS_SENIOR_EDITOR_REVIEW.md)
**Revised files:** `paper/02_introduction_REVISED.Rmd`, `paper/03_methods_REVISED.Rmd`, `paper/04_results_REVISED.Rmd`, `paper/05_discussion_REVISED.Rmd`
**Supporting documents:** [SUPPLEMENTARY_UPDATES.md](SUPPLEMENTARY_UPDATES.md), [NARRATIVE_ARC.md](NARRATIVE_ARC.md), [SUGGESTED_TEXT.md](SUGGESTED_TEXT.md)

---

## Status Key

- DONE = Addressed in REVISED files
- PLANNED = Documented in SUPPLEMENTARY_UPDATES.md but not yet implemented
- TODO = Not yet addressed anywhere
- PARTIAL = Partially addressed

---

## I. Blocking Items (Must Fix for Resubmission)

### 1. PNAS Formatting

| Item | Status | Location | Notes |
|------|--------|----------|-------|
| Reduce figures from 6 to 4 | TODO | — | NARRATIVE_ARC.md identifies Fig 4 and Fig 6 as essential. Senior editor suggested moving Figs 2 and 6 to SI, but Fig 6 is the cross-cancer claim. Consider merging Figs 5+6 into one "Generalization" figure, and moving Fig 2 (rank selection details) to SI. |
| Reduce references from 71 to ~50 | TODO | — | **Worsened**: 8 new references added. Net target: cut ~29 references. Strategy needed. |
| Keywords placeholder | DONE | paper.Rmd YAML | Set to: nonnegative matrix factorization, survival analysis, semi-supervised learning, tumor deconvolution, pancreatic cancer |
| Author contributions placeholder | TODO | — | Not addressed in any revision file |
| Conflict of interest placeholder | TODO | — | Not addressed in any revision file |
| Acknowledgements placeholder | TODO | — | Not addressed in any revision file |
| Affiliations incomplete | TODO | — | Not addressed in any revision file |
| Standardize "DeSurv" naming | DONE | paper.Rmd abstract | Fixed "deSurv" → "DeSurv" in abstract |

### 2. Discussion Expansion

| Item | Status | Location | Notes |
|------|--------|----------|-------|
| Expand from 3 to 5 paragraphs | DONE | 05_discussion_REVISED.Rmd | 5 paragraphs: contributions, PDAC biology, semi-supervised tradeoff, limitations, broader implications |
| Limitations paragraph | DONE | 05_discussion_REVISED.Rmd para 4 | 5 explicit limitations: Cox PH assumption, BO cost, single transfer pair, theory-practice gap, outcome data quality |
| When DeSurv should NOT be used | DONE | 05_discussion_REVISED.Rmd para 3 | "DeSurv is therefore not intended to replace unsupervised NMF in all settings" |
| Clinical translation pathway | TODO | — | Not discussed. Consider adding one sentence to discussion para 5 |
| Biological interpretation of cross-cancer transfer | DONE | 05_discussion_REVISED.Rmd para 2 | "consistent with prior evidence that basal-like transcriptional programs are shared across epithelial cancers" |
| Comparison with related frameworks (iCluster, Spectra, WGCNA) | TODO | — | Not addressed. Senior editor specifically asked for Spectra comparison |

### 3. Comparator Methods

| Item | Status | Location | Notes |
|------|--------|----------|-------|
| Add LASSO-Cox as comparator | TODO | — | Requires new analysis. Most natural "prediction-only" baseline |
| Empirical comparison with CoxNMF/SurvNMF or clear technical argument | PARTIAL | 02_introduction_REVISED.Rmd line 17 | Technical argument present ("a property not shared by methods that supervise through H") but no empirical comparison |

### 4. Simulation Results

| Item | Status | Location | Notes |
|------|--------|----------|-------|
| Add null scenario (R00_null) | DONE | 04_results_REVISED.Rmd (simulation section) | Text added: describes null result (no advantage, C-index ~0.5, BO selects low alpha). References SI Appendix Fig. S\ref{fig:sim-null-mixed} |
| Add mixed scenario (R_mixed) | DONE | 04_results_REVISED.Rmd (simulation section) | Text added: describes attenuated advantage. Establishes dose-response across three scenarios |
| Report simulation parameters (k, n, G, replicates) | PARTIAL | 03_methods_REVISED.Rmd lines 31-34, supplement.Rmd | Three scenarios described in methods; full params (G=3000, N=200, K=3, gamma params) now in supplement text. SUPPLEMENTARY_UPDATES.md 1D plans additional detail |
| Include null/mixed figures in supplement | DONE | supplement.Rmd (fig:fig-sim-null-mixed) | 5-panel figure: null (A-B: C-index, k-hist) and mixed (C-E: C-index, precision, k-hist). Precision omitted for null (undefined when beta=0). Uses tar_load(sim_figs_by_scenario) |

### 5. Missing Quantitative Metrics

| Item | Status | Location | Notes |
|------|--------|----------|-------|
| Report C-index values numerically | TODO | — | Still only shown visually (heatmap, boxplots) |
| Add HR, CI, p-values to KM plots | TODO | — | Requires figure regeneration. Flagged as critical in senior editor review |
| Report precision values from simulations | TODO | — | Only shown in boxplots |
| Report effect sizes for DeSurv vs NMF | TODO | — | No mean C-index difference with SE reported |

### 6. Other Blocking Items

| Item | Status | Location | Notes |
|------|--------|----------|-------|
| Verify sdZ bug in compute_metrics.R | TODO | — | CODE_REVIEW.md:9 flags `sdZ` assigned from `fit$meanZ`. If bug present during reported results, all metrics may be affected. **Highest-priority code fix.** |
| Resolve convergence criterion mismatch | TODO | — | Paper: relative loss change. Code: cosine similarity of W, H. Must document which was actually used |
| Resolve bladder data split (70/30 vs 80/20) | TODO | — | Main text says 70/30; supplement says 80/20 |

---

## II. High Priority Items (Strengthen for PNAS)

| Item | Status | Location | Notes |
|------|--------|----------|-------|
| Expand Methods with algorithmic detail | PARTIAL | 03_methods_REVISED.Rmd | Semi-supervised framing and out-of-sample statement added. Still missing: sample sizes per dataset, BO specifics, convergence criterion, runtime estimates, Cox PL formula |
| Add factor labels to Fig 4 heatmaps | TODO | — | Requires figure regeneration |
| Sharpen W-vs-H novelty claim | DONE | 02_introduction_REVISED.Rmd para 4, 03_methods_REVISED.Rmd line 21 | Explicit: "the survival gradient acts explicitly on the gene program matrix W" and "sample-level loadings H enter only through the reconstruction term" |
| Add explicit out-of-sample statement | DONE | 03_methods_REVISED.Rmd line 26 | "All model selection was performed entirely within training data using nested cross-validation; no validation cohort data were used during any stage of tuning" |
| Add projection advantage to main text | DONE | 02_introduction_REVISED.Rmd line 17, 04_results_REVISED.Rmd line 58 | "new samples can be scored by simple projection (Z_new = W^T X_new) without requiring their survival data" |
| Add simulation parameter annotations | TODO | — | Requires figure update |
| Document W update clamping | PLANNED | SUPPLEMENTARY_UPDATES.md 1C | Theory-to-practice gap paragraph planned for supplement |
| Fix symbol overloading (delta, theta) | TODO | — | Not addressed |
| Add notation table to supplement | TODO | — | Not addressed |

---

## III. Medium Priority Items

| Item | Status | Location | Notes |
|------|--------|----------|-------|
| Restructure Fig 2 layout | TODO | — | Requires figure regeneration |
| Add sensitivity analysis for k | TODO | — | Requires new analysis |
| Describe consensus initialization | TODO | — | Referenced as "SI Appendix" but incomplete |
| Add dataset characteristics table | TODO | — | Sample sizes, platforms, censoring rates |
| Discuss Spectra and supervised topic models | TODO | — | |
| Add runtime comparison | TODO | — | |

---

## IV. Items Added by Revision (Not in Original Review)

These are new elements introduced in the REVISED files that need verification.

| Item | Status | Notes |
|------|--------|-------|
| 8 new bib entries | DONE | All 8 entries added to references_30102025.bib. `consensus2025pdac` placeholder replaced with `collisson2019molecular` (already in bib). Duplicate `brunet2004metagenes` removed. |
| Cook citation verification | DONE | Using `cook2007fisher` (2007 Fisher Lecture on dimension reduction in regression). This is the correct SDR theory reference. |
| Semi-supervised framing paragraph in supplement | PLANNED | SUPPLEMENTARY_UPDATES.md 1A |
| Formal dL_Cox/dH = 0 derivation in supplement | PLANNED | SUPPLEMENTARY_UPDATES.md 1B |
| SDR theory connection in supplement | PLANNED | SUPPLEMENTARY_UPDATES.md 1F |
| Complete scRNA-seq text in supplement | PLANNED | SUPPLEMENTARY_UPDATES.md 2A — currently "Need more here..." at line 41 of supplement.Rmd |
| Significance statement revision | PLANNED | SUPPLEMENTARY_UPDATES.md 4A |
| Abstract update | PLANNED | SUPPLEMENTARY_UPDATES.md 4C — needs semi-supervised language, cross-cancer transfer, simulation validation |
| Word count check | TODO | Revised introduction appears ~1200 words; PNAS body limit is 4,500 words total. Full word count needed |

---

## V. Simulation Scenario Reference

For convenience, the three scenarios from `R/simulation_functions/scenario_defaults.R`:

| Parameter | R00_null | R0_easy | R_mixed |
|-----------|----------|---------|---------|
| G (genes) | 3000 | 3000 | 3000 |
| N (samples) | 200 | 200 | 200 |
| K (factors) | 3 | 3 | 3 |
| Markers/factor | 150 | 150 | 150 |
| Background genes | 500 | 500 | 500 |
| **beta** | **(0, 0, 0)** | **(2, 0, 0)** | **(2, 0, 0)** |
| **Survival genes (Factor 1)** | **N/A (beta=0)** | **150 markers** | **150 markers + 150 background** |
| survival_marker_frac | — | — | 0.5 |
| Marker gamma | shape=3, rate=0.8 | shape=3, rate=0.8 | shape=3, rate=0.8 |
| Background gamma | shape=2, rate=1.0 | shape=2, rate=1.0 | shape=2, rate=1.0 |
| Expected DeSurv advantage | **None** | **Large** | **Moderate** |

**What each scenario tests:**
- **R00_null:** beta = 0 means survival is independent of gene expression. DeSurv should not outperform NMF. Tests specificity — does the method hallucinate structure?
- **R0_easy:** One factor drives survival through its 150 marker genes. Background genes dominate variance but are outcome-neutral. Tests the core claim: supervision helps when variance and prognosis diverge.
- **R_mixed:** Survival depends on 300 genes — half markers (factor-specific), half background (shared across factors). Partial overlap between variance-dominant and survival-relevant structure. Tests whether DeSurv's advantage attenuates when the assumption of divergence is relaxed.

**Figure files already exist** in `figures/sim/` (see SUPPLEMENTARY_UPDATES.md Section 5 for complete listing).

---

## VI. Priority Sequence for Remaining Work

### Phase 1: Unblock rendering (effort: low) — COMPLETE
1. ~~Add 8 bib entries to `paper/references_30102025.bib`~~ DONE
2. ~~Replace YAML placeholders (keywords)~~ DONE; significance statement already present
3. ~~Verify/fix `consensus2025pdac` citation key~~ DONE — replaced with `collisson2019molecular`

### Phase 2: Critical content gaps (effort: medium)
4. Write null/mixed scenario text for Results (replace HTML comments at 04_results_REVISED.Rmd:114-123)
5. Verify sdZ bug in `R/compute_metrics.R` — re-run if needed
6. Resolve convergence criterion mismatch (document which was used)
7. Resolve bladder data split discrepancy

### Phase 3: Quantitative reporting (effort: medium)
8. Add numeric C-index values to results prose
9. Add HR, CI, p-values to KM plot figure code
10. Report simulation precision values

### Phase 4: PNAS compliance (effort: high)
11. Figure consolidation plan (6 to 4)
12. Reference reduction plan (71+ to ~50)
13. Complete all remaining YAML placeholders (author contributions, conflicts, acknowledgements, affiliations)
14. Word count verification

### Phase 5: Supplement completion (effort: medium)
15. Items from SUPPLEMENTARY_UPDATES.md Sections 1A-1F, 2A-2C

### Phase 6: Additional analyses (effort: high, optional for initial resubmission)
16. LASSO-Cox comparator
17. Runtime comparison
18. Spectra/iCluster discussion

---

*Created: 2026-02-10*
*Source: Peer review analysis of REVISED files against PNAS_SENIOR_EDITOR_REVIEW.md*
