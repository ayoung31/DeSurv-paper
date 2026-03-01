# DeSurv Manuscript Review — 2026-02-28

Comprehensive review of recent changes (commits `044257d`–`a5db7f9`), rendered
paper.pdf, supplement.pdf, supp_methods.pdf, source Rmd files, and bibliography.

---

## 1. Statistical Claims vs Evidence

### 1a. [P0] Bladder cancer transfer result is non-significant — text overclaims

**Rendered PDF (p. 5):** n = 51, 34 events; log-rank P = 0.68; HR 0.66; 95% CI 0.09–4.81

The CI spans nearly two orders of magnitude and easily crosses 1.0, yet the text states:

- Results: "demonstrated **clear separation** of survival curves" (p. 5)
- Results: "**suggests** that survival supervision captures cross-cancer biology" (p. 5)
- Intro: "we show that a PDAC-trained DeSurv factor **retains prognostic signal**" (p. 2)
- Discussion ¶1: "prognostic signal **transferred** to bladder cancer" (line 7)

None of these claims are supported by P = 0.68.

**Root cause:** The old code used the full bladder cohort (n = 260) with a simple median
split on projected PDAC factor scores. The new code uses only the held-out test split
(n = 51) from `_targets_bladder.R` with an optimized cutpoint, destroying statistical
power.

**Options:**
1. Temper language: "showed a non-significant trend", "exploratory", "hypothesis-generating"
2. Use the full bladder cohort if the train/test split is too aggressive for this analysis
3. Add a power analysis showing n = 51 is underpowered for the observed effect size
4. Remove the claim from the abstract/intro if it cannot be supported

**Files:** `paper/04_results_REVISED.Rmd` (lines 398–406), `paper/02_introduction_REVISED.Rmd`
(line 23), `paper/05_discussion_REVISED.Rmd` (lines 7, 9)

### 1b. [P0] Inline text statistics vs Figure 5B use different cutpoints

The `bladder-precompute` chunk computes inline stats using the **log-rank cutpoint**:
```r
# 04_results_REVISED.Rmd:383
sdf$factor <- as.integer(z_val > desurv_lp_stats_bladder_nobcg$optimal_z_cutpoint)
```
Inline result: P = 0.68, HR = 0.66, CI 0.09–4.81

The figure uses the **C-index cutpoint** via a different target:
```r
# 04_results_REVISED.Rmd:413
tar_load(fig_cindex_cutpoint_survival_desurv_bladder_nobcg)
```
Figure annotation: P = 0.477, HR = 0.78, CI 0.39–1.56

The reader sees different statistics in the text and the figure for the same analysis.

**Fix:** Pick one cutpoint method and use it consistently in both the inline computation
and the figure. Then verify the caption matches the actual stratification method.

**Files:** `paper/04_results_REVISED.Rmd` (lines 368–395 vs 413–414)

### 1c. [P2] Validation cohort counts need clarification

- Page 4 text: "Across **5** independent PDAC cohorts (n = **616**)"
- K×α sensitivity section: references `n_val_total` patients across `prod_sum$n_val_cohorts` cohorts
- Supplement sensitivity analysis: "**570** independent validation patients across **4** cohorts"

Are these referring to different analyses (e.g., 4 cohorts for iCAF factor vs 5 for pooled)?
If so, clarify in text. If not, reconcile the counts.

---

## 2. Unresolved References and Cross-References

### 2a. [P1] Two [REF] placeholders in main text

`paper/04_results_REVISED.Rmd:175`:
> "consistent with the known co-occurrence of classical tumor identity and iCAF-enriched
> stroma in PDAC **[REF]** and with evidence that iCAF activation is independently
> associated with patient outcomes **[REF]**"

These render directly into the PDF on page 3.

**Suggested citations:**
- iCAF-classical co-occurrence: `@ohlund2017distinct` (Öhlund et al. 2017, already in bib)
  or `@elyada2019cross`
- iCAF and outcomes: `@Maurer2019` (Maurer et al. 2019, already in bib)

### 2b. [P1] Undefined `\ref{fig:forest}` in main text

`paper/04_results_REVISED.Rmd:307` references `Supplementary Fig. \ref{fig:forest}`.
The forest plot is defined in `supplement.Rmd:439`. LaTeX `\ref{}` does not resolve
across separate documents.

`paper.log:1170`: "LaTeX Warning: Reference `fig:forest' on page 4 undefined"

**Fix:** Replace with hardcoded "SI Appendix, Fig. 10" (or whatever the final number is).

### 2c. [P1] Undefined `\ref{fig:varsurvival}` and `\ref{fig:corr}` in main text

Both referenced in `04_results_REVISED.Rmd` (lines 179, 181) but defined in
`supplement.Rmd` (lines 409, 425). These render as "??" in the compiled PDF (page 3).

Also: `supp_methods.Rmd:885` references `\ref{fig:varsurvival}` — this will also show
as "??" since supp_methods is compiled separately.

**Fix:** Replace all cross-document `\ref{fig:...}` with explicit text references
(e.g., "SI Appendix, Fig. S8"). Consider adding a mapping comment in the Rmd listing
supplement figure numbers for easy maintenance.

---

## 3. Caption and Figure Consistency

### 3a. [P2] Figure 5B caption says "median split" but uses optimized cutpoint

Caption: "stratified by **median split** of projected factor scores from a PDAC-trained
DeSurv model"

Code: uses `fig_cindex_cutpoint_survival_desurv_bladder_nobcg` which applies a
C-index-optimized z-score cutpoint, not a median split.

**Fix:** Update caption to describe the actual cutpoint method used.

**File:** `paper/04_results_REVISED.Rmd:408` (the `fig.cap` parameter)

### 3b. [P3] Figure 5B legend says "Strata" instead of "Risk group"

The PDAC KM plots in Figure 4C–D use "Risk group" in the legend. Figure 5B uses
"Strata." Standardize across all KM panels.

### 3c. [P3] Figure 2 panels A–C shared legend may confuse

The shared legend shows 4 measure types (Basis, Best fit, Coefficients, Consensus).
Panel A (residuals) only has one line, but the legend implies 4 series. Panels B–C
(cophenetic, silhouette) show the multiple distance metrics. Consider annotating
or splitting the legend to avoid confusion about what Panel A shows.

### 3d. [P3] Dead code: `splot_bl` object is constructed but never used

`04_results_REVISED.Rmd:416-428` builds `splot_bl` from inline data but the final
`plot_grid` call (line 430) uses `km_bladder` instead. Remove the dead code to
avoid confusion.

---

## 4. Citation Errors

### 4a. [P2] "33 cancer types" attributed to DECODER (wrong reference)

Intro, paragraph 1: "compartment-specific deconvolution across 33 cancer types (5)"

Reference 5 is Peng et al. 2019 (DECODER), which was applied to **pancreatic cancer**
specifically, not 33 cancer types. The 33-cancer-type work is Hoadley et al. 2018
(ref 21) or Peng et al. 2019 mentions TCGA pan-cancer but DECODER itself was PDAC-focused.

**Fix:** Either change to "compartment-specific deconvolution in pancreatic cancer"
or cite the correct pan-cancer reference.

### 4b. [P2] Schwarzová 2023 mis-cited for PDAC subtyping claim

Intro, paragraph 2: "Over a decade of PDAC subtyping efforts proposed between two
and six subtypes (14, 13)"

Reference 14 (Schwarzová et al. 2023) is about **stroma-rich bladder cancers**,
not PDAC subtyping. This appears to be an incorrect citation.

**Fix:** Replace with an appropriate PDAC subtyping review (e.g., Collisson 2019
is already ref 13 and covers this; or cite Puleo 2018, Bailey 2016).

### 4c. [P2] Seung & Lee 2001 bib entry has incorrect author names

`references_30102025.bib:9-17`:
```bibtex
author={Seung, D and Lee, L},
```

The actual paper is **Lee, D.D. and Seung, H.S.** "Algorithms for non-negative matrix
factorization," NeurIPS 2001. Author names are abbreviated incorrectly and the
standard first-author ordering is Lee, not Seung.

### 4d. [P3] Huang et al. 2020 — check publication status

Reference 18 (`huang2020low`) is still listed as "arXiv preprint arXiv:2008.03776"
from 2020. After 6 years, check whether a peer-reviewed version was published.

### 4e. [P3] Le Goff 2025 — check for published paper version

Reference 19 (`le2025survnmf`) is cited as a PhD thesis. If a journal paper version
of SurvNMF exists, it should be cited instead.

---

## 5. Compilation and Build Issues

### 5a. [P1] Supplement PDF appears stale — needs recompilation

The rendered `supplement.pdf` still shows:
- The scRNA-seq section with raw citation keys (`[@elyada2019cross]`)
- Placeholder text: "Need more here..."
- Content that was commented out in `supplement.Rmd` via HTML comments

The supplement.Rmd changes (commenting out the scRNA-seq section) have not been
reflected in the compiled PDF.

**Fix:** Re-render `supplement.Rmd` after the HTML comment changes.

### 5b. [P3] paper.log overfull hbox warnings

Several `Overfull \hbox` warnings in `paper.log`. These are cosmetic but worth
addressing before final submission.

---

## 6. PNAS Formatting

### 6a. [P3] Acknowledgments are still placeholder text

`paper.Rmd` renders: "Please include your acknowledgments here, set in a single paragraph."

### 6b. [P3] Versioning block in rendered PDF

The git commit hashes at the end of the paper are useful for reproducibility but
should be removed (or moved to SI) for journal submission.

### 6c. [INFO] Reference count is within PNAS limits

Current: 26 references. PNAS limit: 50. No action needed.

### 6d. [INFO] Figure count

Currently 5 main-text figures. PNAS limit is typically 6 (was 4 for some article types).
Verify against the target article type.

---

## 7. Code Quality (from diff review)

### 7a. [DONE] `browser()` statements removed

Two debug `browser()` calls were correctly removed from
`targets_common_pipeline.R:2332` and `:2448`. No further action needed.

### 7b. [INFO] R0k6 simulation scenario commented out

`_targets_sims.R:125-133` comments out the "R0 scenario with k = 6". Verify no
references to "R0k6" or "k = 6 scenario" remain in paper or supplement text.

### 7c. [INFO] Silhouette legend logic removed from `make_nmf_metric_plot`

`R/figure_targets.R:549-560` removed the if/else block that handled per-metric
legend display. The shared legend is now extracted externally from `fit_std_tcgacptac`
and placed below panels A–C. This is a cleaner approach but see item 3c about
potential confusion.

### 7d. [INFO] New `fig_cindex_cutpoint_survival_desurv` target added

`targets_common_pipeline.R:2851-2854` adds a new figure target using
`cutpoint_field = "optimal_z_cutpoint_cindex"`. This is the target used in Figure 5B.
The existence of two cutpoint methods (log-rank and C-index) is by design, but the
inconsistency in the bladder section (issue 1b) needs resolution.

---

## Priority Summary

| Priority | Count | Category |
|----------|-------|----------|
| **P0**   | 2     | Statistical claims vs evidence |
| **P1**   | 4     | Unresolved references, stale compilation |
| **P2**   | 5     | Caption consistency, citation errors |
| **P3**   | 7     | Minor polish, formatting, dead code |
| **INFO** | 4     | No action needed, for awareness |
