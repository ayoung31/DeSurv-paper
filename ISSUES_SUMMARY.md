# GitHub Issues Summary — DeSurv-paper

Generated: 2026-03-02

---

## Status Overview

| # | Title | Status |
|---|-------|--------|
| [#24](#24-citation-errors-and-unresolved-references) | Citation errors and unresolved references | ❌ Open — not addressed |
| [#23](#23-bladder-cross-cancer-transfer-non-significant-result) | Bladder cross-cancer transfer: non-significant result | ❌ Open — not addressed |
| [#17](#17-reduce-figures-from-6-to-4-for-pnas-compliance) | Reduce figures from 6 to 4 for PNAS compliance | ⚠️ Partially addressed |
| [#10](#10-reframe-desurv-advantage-parsimony-and-interpretability) | Reframe DeSurv advantage: parsimony and interpretability | ⚠️ Discussed, not implemented |
| [#9](#9-standard-nmf-at-k5-outperforms-desurv-k3-in-external-validation) | Standard NMF K=5 outperforms DeSurv K=3 in external validation | ❌ Open — not addressed |
| [#8](#8-clarify-k-selection-for-standard-nmf-in-external-validation) | Clarify K selection for standard NMF in external validation | ❌ Open — not addressed |
| [#7](#7-connect-variance-vs-survival-plot-back-to-simulation-scenarios) | Connect variance-vs-survival plot back to simulation scenarios | ❌ Open — not addressed |
| [#5](#5-characterize-standard-nmf-factors-at-independently-selected-k) | Characterize standard NMF factors at K=5 and K=7 for supplement | ⚠️ Discussed, not implemented |
| [#4](#4-clarify-k-selection-for-standard-nmf-comparison-section-4) | Clarify K selection for standard NMF comparison (Section 4) | ✅ Resolved (decision made, disclosure needed in text) |
| [#2](#2-mixed-scenario-investigate-desurv-marker-vs-background-gene-effects) | Mixed scenario: marker vs background gene survival effects | ✅ Resolved — fully implemented |
| [#1](#1-clarify-survclust-model) | Clarify survClust model | ✅ Resolved — questions answered |

---

## Issue Summaries

---

### #24: Citation errors and unresolved references

**Status: ❌ Not addressed** (no comments, no commits referencing this issue)

**Summary:**
Seven citation errors, placeholder references, and broken cross-references identified during manuscript review:

1. Two literal `[REF]` placeholders remain in `paper/04_results_REVISED.Rmd:175` (iCAF biology claims)
2. "33 cancer types" wrongly attributed to DECODER (Peng et al.) — DECODER was PDAC-specific, not pan-cancer
3. Schwarzová 2023 (a bladder cancer paper) incorrectly cited for PDAC subtyping
4. Seung & Lee 2001 NMF paper has reversed author names in the bib entry
5. Three `\ref{fig:...}` cross-references in main text point to supplement figures — they render as "??" in the PDF
6. Same broken `\ref{fig:varsurvival}` in `supp_methods.Rmd:885`
7. Huang 2020 still listed as arXiv preprint (6 years old); Le Goff 2025 listed as PhD thesis — worth checking for published versions

**Proposed solutions:**

1. **`[REF]` placeholders** (`04_results_REVISED.Rmd:175`): Replace with:
   - iCAF-classical co-occurrence → `@ohlund2017distinct` or `@elyada2019cross` (both in bib)
   - iCAF outcomes → `@Maurer2019` (in bib)

2. **33 cancer types / DECODER** (`02_introduction_REVISED.Rmd:15`): Change to "compartment-specific deconvolution in pancreatic cancer (5)" or add the correct pan-cancer reference alongside DECODER.

3. **Schwarzová mis-citation** (`02_introduction_REVISED.Rmd:17`): Replace ref 14 with a PDAC-subtyping ref such as `@puleo2018stratification` or `@Bailey2016`.

4. **Seung/Lee bib fix** (`references_30102025.bib`):
   ```bibtex
   author={Lee, Daniel D and Seung, H Sebastian},
   ```

5. **Broken `\ref{}`**: Replace all three cross-document refs with hardcoded supplement figure numbers, e.g.:
   - `\ref{fig:forest}` → `SI Appendix, Fig. S10`
   - `\ref{fig:varsurvival}` → `SI Appendix, Fig. SX`
   - `\ref{fig:corr}` → `SI Appendix, Fig. SY`
   Add a comment block at the top of `04_results_REVISED.Rmd` mapping SI figure numbers.

6. **Huang/Le Goff**: Search for published peer-reviewed versions. Update bib entries if found.

---

### #23: Bladder cross-cancer transfer: non-significant result, text/figure mismatch

**Status: ❌ Not addressed** (comment from naimurashid adds detail, but no resolution decision or code change)

**Summary:**
The bladder cross-cancer section has three compounding problems:

1. **Non-significant result with overclaiming prose**: P = 0.68, HR = 0.66, CI 0.09–4.81 in inline text. The text says "clear separation," "suggests prognostic signal," "transferred to bladder cancer" — none of these are supported.

2. **Sample size collapse**: Previous version used full cohort (n = 260). Current version uses only the held-out 30% test split from the bladder pipeline (n = 51), destroying statistical power. The issue notes that this is arguably the wrong split to use: the point is cross-cancer transfer from a *PDAC-trained* model; bladder samples were never in the PDAC training set, so all bladder samples are genuinely held-out.

3. **Cutpoint mismatch between inline text and figure**: Inline stats use log-rank-optimized cutpoint (P = 0.68); Figure 5B uses C-index-optimized cutpoint (P = 0.477). The reader sees different statistics for the same analysis. Additionally, the figure caption says "median split" but the code uses an optimized cutpoint.

**Proposed solutions:**

The most defensible option (Option 3 from the issue) is to **project the PDAC W matrix onto the full bladder cohort (n ≈ 260)**. Rationale: the claim is about transfer from a PDAC-trained model; the bladder samples were never used to train PDAC-DeSurv. The bladder pipeline's 70/30 split is irrelevant to evaluating PDAC→bladder transfer.

**Implementation steps:**
1. In `paper/04_results_REVISED.Rmd`, load `raw_data_bladder_nobcg` (or equivalent full cohort) rather than `data_val_filtered_bladder_nobcg`
2. Preprocess and project PDAC W matrix onto all ~260 samples
3. Unify cutpoint method between inline stats and figure — pick one approach (median split for simplicity, or document the cutpoint selection method)
4. Update the figure caption to match the cutpoint method used
5. Update all prose claims to match actual statistics:
   - If P < 0.05 after using full cohort: update language to honestly report the value
   - If still non-significant: change to "non-significant trend (n = 260, P = X)" and reframe as preliminary cross-cancer evidence
6. Remove "clear separation," "transferred prognostic signal" language unless it becomes statistically justified
7. Update `02_introduction_REVISED.Rmd:23` and `05_discussion_REVISED.Rmd:7,9` to match

---

### #17: Reduce figures from 6 to 4 for PNAS compliance

**Status: ⚠️ Partially addressed**

**Comment from @ayoung31:** "Moved parts of figure 5 to supplement and combined key panels from 4 and 5 into one figure (figure 4)"

**Remaining items not confirmed addressed:**
- References still need reduction from 71 to ≤50 (21+ to cut)
- Placeholder sections (keywords, author contributions, conflict of interest, acknowledgements) need completion
- Need to verify actual figure count after reorganization is ≤4

**Proposed remaining steps:**
1. Count current main-text figures after reorganization — confirm ≤4
2. Audit reference list: identify and remove the least-cited or most-redundant references. Good candidates: review articles that can be combined into a single representative citation, preprints with published versions, or topic-adjacent refs not strictly needed for the argument.
3. Complete manuscript metadata:
   - Keywords (5–10 terms relevant to PNAS audience)
   - Author contributions (CRediT taxonomy)
   - Conflict of interest statement
   - Acknowledgements

---

### #10: Reframe DeSurv advantage: parsimony and interpretability over raw prediction

**Status: ⚠️ Discussed but not implemented in paper text**

**Summary:**
Analysis of stored validation results shows that standard NMF at K=5 (independently selected by elbow) outperforms DeSurv K=3 on 4/5 external validation cohorts. The key insight is that DeSurv's advantage is **parsimony** (achieves with K=3 what NMF needs K=7–11 to achieve), not raw predictive superiority.

From BO landscape data (comment by @naimurashid):
- DeSurv K=3 achieves CV C-index 0.647
- Standard NMF doesn't reach 0.647 until K=11
- DeSurv max (K=7, C-index 0.655) beats NMF max (K=11, C-index 0.647) with fewer factors
- Alpha decreases as K increases for DeSurv — at high K, DeSurv converges toward NMF behavior

**The three-layer argument identified in comments:**
1. At same K (K=3): DeSurv dramatically outperforms (0.647 vs 0.554) — supervision allocates factors to prognostic biology
2. At each method's optimal K: DeSurv still wins (0.655 at K=7 vs 0.647 at K=11) with fewer factors
3. The 1-SE rule selects K=3, where DeSurv's advantage is largest

**Proposed text changes (Section 5 and Discussion):**
- Replace: "DeSurv outperforms standard NMF in generalization"
- With: "DeSurv achieves equivalent or superior discriminative performance with a more parsimonious factorization — concentrating prognostic signal into K=3 factors where standard NMF requires K=5–11 to reach comparable performance"
- Add to Discussion: explicit parsimony-vs-prediction tradeoff paragraph
- Note that alpha decreasing with K is consistent with the interpretability story

---

### #9: Standard NMF at K=5 outperforms DeSurv K=3 in external validation

**Status: ❌ Not addressed**

**Summary:**
Same finding as #10 but focused specifically on the external validation section. At K=5 (NMF's independently selected rank), NMF outperforms DeSurv K=3 on 4/5 external cohorts. The current Section 5 framing implies DeSurv generalizes better, which is only true at K=3.

**Proposed solutions:**

1. **Section 5 text revision**: Add explicit disclosure that the comparison uses K=3 for both methods (referencing Section 4 disclosure). Add a statement like:
   > "At the same factorization rank (K=3), DeSurv factors showed more consistent survival associations across cohorts. When standard NMF was permitted to select its own rank (K=5), discriminative performance in external validation was comparable (Table SX), though distributed across more factors."

2. **Add supplementary table**: Show external validation C-index for DeSurv K=3, Std NMF K=3, Std NMF K=5 across all cohorts side by side.

3. **Check HR consistency and KM separation at NMF K=5**: C-index alone may not capture the full picture. If NMF K=5 factors show inconsistent hazard ratios across cohorts (some protective in one cohort, harmful in another), this reinforces the interpretability/reliability argument for DeSurv even when C-indices are similar.

---

### #8: Clarify K selection for standard NMF in external validation (Section 5, Fig 5)

**Status: ❌ Not addressed**

**Summary:**
The forest plot (Fig 5A) and KM curves (Fig 5B-C) compare DeSurv to standard NMF using `fit_std_desurvk` — standard NMF forced to DeSurv's K. Section 5 does not disclose this. The `val_predictions_std_elbowk` and `val_latent_std_elbowk` targets exist but are not used in the main figures.

**Proposed solutions:**

1. Add a disclosure sentence to Section 5 (parallel to any Section 4 disclosure): "Standard NMF factors were evaluated at K=3, matching DeSurv's BO-selected rank, to enable direct comparison of factor structure."

2. Either:
   - Add supplementary forest plot and KM curves for NMF at K=5/7 using `val_predictions_std_elbowk`, or
   - Reference Issue #5 supplement characterization figures here

3. Minimum: ensure Section 5 and Section 4 are internally consistent in how they describe the K=3 constraint.

---

### #7: Connect variance-vs-survival plot (Fig 4C) back to simulation scenarios

**Status: ❌ Not addressed** (no comments)

**Summary:**
Fig 4C shows that the highest-variance NMF factor has minimal survival association, while DeSurv concentrates survival signal. The issue asks whether this can be explicitly tied to the primary and mixed simulation scenarios, noting that real PDAC data likely falls between the two.

**Proposed solutions:**

1. **Add one sentence to Section 4** connecting the real-data observation to the simulation framework:
   > "This pattern — where the highest-variance factor shows minimal survival association while other factors show partial overlap — is consistent with the intermediate regime between the primary and mixed simulation scenarios characterized above."

2. **Optional additional analysis**: Compute the fraction of survival-associated genes that also load on variance-dominant NMF factors in PDAC. A value near 0% would place PDAC in the "primary" scenario; values near 50% in the "mixed" scenario. This could be reported as a single sentence or added as a supplement panel.

3. **In the mixed scenario simulation supplement**: Add a variance-vs-survival analog figure (if simulatable) to show that DeSurv also reorganizes factor allocation in the mixed case, paralleling Fig 4C. This visual parallel would strengthen the simulation↔real-data narrative bridge.

---

### #5: Characterize standard NMF factors at independently-selected K (K=5 and K=7) for supplement

**Status: ⚠️ Discussed but not implemented**

**Comment from @naimurashid** raises additional question: K=5 (elbow, purely unsupervised) vs K=7 (BO at α=0, CV C-index) represent different "standard NMF" baselines. Both should likely be shown.

**Summary:**
`fit_std_elbowk_tcgacptac` (K=5) and the α=0 BO result (K=7) exist in the pipeline but haven't been characterized biologically or added to the supplement.

**Proposed solutions:**

1. **Add supplement section** "Standard NMF factor structure at independently selected K":
   - Load `fig_gene_overlap_heatmap_std_elbowk_tcgacptac` (K=5) and describe which factors map to exocrine, classical, and microenvironmental programs
   - Show that extra factors (beyond K=3) capture variance-dominant programs (likely exocrine subtypes or composition signals) with weak survival association
   - Note that K=7 (BO at α=0) selects higher K than elbow because it's optimizing survival prediction — but still higher than DeSurv K=3

2. **Add supplement table** showing factor biological annotations at K=3, K=5, K=7 for standard NMF

3. **Closing supplement claim**: "When standard NMF selects its own rank — whether by variance elbow (K=5) or cross-validated concordance (K=7) — it requires more factors than DeSurv (K=3) to achieve comparable discriminative performance, and those additional factors capture variance-dominant programs with minimal survival relevance."

4. Note in the supplement which K-selection criterion is more appropriate for each comparison:
   - K=5 (elbow) for "purely unsupervised" comparison
   - K=7 (BO at α=0) for "same selection criterion, only difference is α" comparison

---

### #4: Clarify K selection for standard NMF comparison (Section 4)

**Status: ✅ Resolved**

**Comment from @ayoung31:** "I believe this has been discussed. Comparing the methods at the same k provides a better illustration of factor reorganization. other k to be explored in supplement"

**Remaining action**: The decision is documented in this issue, but the paper text should include an explicit disclosure statement so reviewers aren't surprised. A one-sentence clarification in Section 4 ("Standard NMF was evaluated at K=3, matching DeSurv's BO-selected rank") would close this out.

---

### #2: Mixed scenario — DeSurv separates marker vs background gene survival effects

**Status: ✅ Fully implemented**

**Comment from @ayoung31** documents a complete implementation:
- `_targets_sims.R`: Added `marker_only_metrics` columns (`marker_only_lethal_factor_metrics`, `marker_only_ari`) via a second `purrr::pmap` call evaluating precision against `truth$marker_sets` (150 genes) specifically
- `sim_figs.R`: Added `extract_matched_factor_beta()`, `plot_mixed_precision_breakdown()`, `plot_matched_factor_beta()`, and updated `build_sim_figs()` to append these for R_mixed scenario
- `paper/supplement.Rmd`: Updated `fig-sim-null-mixed` to 3-row layout with new panels; updated mixed scenario prose with four-point mechanistic narrative

The scientific claim is now supported: DeSurv's marker-specific precision (150-gene set) substantially exceeds NMF even when overall survival-gene precision (300-gene set) improvement is modest, and the matched factor receives a substantially larger |β| under DeSurv.

**No further action needed** — issue can be closed.

---

### #1: Clarify survClust model and how it determines survival-specific clusters

**Status: ✅ Resolved**

**Comment from @ayoung31** answers all four questions:
1. **Model**: Cox model
2. **Survival incorporation**: Fits univariate Cox per genomic feature; weights Euclidean distance by a diagonal matrix of log(HR) values
3. **Genomic data**: Computes weighted distance matrix across all data types, averages, and clusters
4. **Key difference from DeSurv**: survClust learns discrete cluster assignments; DeSurv learns continuous latent factor programs (W matrix). survClust cannot decompose survival into independent biological axes — two patients with different cell cycle activity and immune infiltration would be collapsed into one or two clusters; DeSurv would learn two separate factors.

**Remaining action**: If this summary hasn't been incorporated into the introduction or discussion prose, add 1–2 sentences contrasting survClust with DeSurv. Suggested text:
> "survClust [ref] incorporates survival by weighting gene-gene distances via univariate Cox hazard ratios, then clustering on these survival-weighted distances. Unlike DeSurv, which learns continuous latent factor programs (W), survClust produces discrete cluster assignments and cannot decompose survival signal into independent biological axes."

Issue can be closed once prose is added.

---

## Priority Action List

### P0 — Fix before submission (blocking)
1. **#24**: Replace `[REF]` placeholders with actual citations
2. **#24**: Fix broken `\ref{fig:...}` cross-references (renders as "??" in PDF)
3. **#23**: Resolve bladder section: decide on full-cohort vs drop, unify cutpoint, fix overclaiming prose
4. **#24**: Fix Schwarzová mis-citation in PDAC subtyping claim

### P1 — Important for scientific integrity
5. **#9/#10**: Add disclosure in Section 5 that comparison uses K=3; reframe DeSurv advantage as parsimony
6. **#8**: Add K disclosure in Section 5 (parallel to Section 4)
7. **#17**: Verify figure count ≤4; reduce references to ≤50

### P2 — Strengthens narrative and preempts reviewers
8. **#5**: Add supplement section characterizing NMF at K=5 and K=7
9. **#7**: Add one sentence in Section 4 connecting Fig 4C to simulation scenarios
10. **#10**: Add Discussion paragraph on parsimony-vs-prediction tradeoff

### P3 — Minor cleanup (quick)
11. **#24**: Fix Seung/Lee author order in bib
12. **#24**: Check publication status of Huang 2020 and Le Goff 2025
13. **#4**: Add K disclosure sentence to Section 4 text
14. **#1**: Add survClust contrast prose to Introduction or Discussion
15. Close issues #2, #1, and #4 on GitHub once confirmatory text edits are made
