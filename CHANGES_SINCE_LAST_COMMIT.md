# Changes Since Last Commit

Detailed summary of all modifications made since commit `f024de7` ("finalize formatting for figures and supplementary sections").

---

## 1. Fixed NMF Pooled Hazard Ratio (NA in paper.pdf)

**File:** `paper/04_results_REVISED.Rmd` (lines 105–151)

**Problem:** The NMF pooled HR displayed as "NA" in the rendered paper because it relied on an external RDS file (`results/nmf_pooled_hr.rds`) that was either missing or contained nonsensical values (HR = 137.76). The root cause was a scale mismatch: DeSurv risk scores had SD ≈ 0.65 while NMF risk scores had SD ≈ 0.01, making per-unit HR comparisons meaningless.

**Fix:** Replaced the fragile RDS-based approach with inline computation using SD-standardized risk scores for both methods:
- Added `extract_pool_df()` helper to extract pooled risk score data from validation latent targets
- Added `pooled_cox_stats()` helper to fit stratified Cox on SD-standardized risk scores
- Now loads `val_latent_std_desurvk_tcgacptac` directly for NMF stats (same target infrastructure as DeSurv)
- Both methods' HRs are now reported "per SD" for fair comparison

**Computed values:** DeSurv HR per SD = 1.45 (95% CI 1.29–1.63); NMF HR per SD = 1.05 (95% CI 0.95–1.17)

---

## 2. Updated Abstract HR to Per-SD Scale

**File:** `paper/paper.Rmd` (line 45)

**Change:** Updated abstract from "pooled validation HR per unit risk score = 1.78, 95% CI 1.48–2.13" to "pooled validation HR per SD = 1.45, 95% CI 1.29–1.63, P < 0.001" to match the SD-standardized computation.

---

## 3. Corrected False Claims About NMF C-index Performance

Table S4 clearly shows NMF at k=7 matches or exceeds DeSurv k=3 C-index in most cohorts. Three locations contained false claims that have been corrected:

### 3a. Main text results (`paper/04_results_REVISED.Rmd`, ~line 349)

**Before:** "standard NMF at elbow-selected k=5 and BO-selected k=7 showed no systematic concordance advantage over NMF at k=3"

**After:** Honestly states that NMF k=5 and k=7 substantially improve over NMF k=3, with k=7 matching or exceeding DeSurv in most cohorts, but at the cost of fragmentation and overlap with existing classifiers. Added "HR per SD" labels throughout.

### 3b. Supplement prose below Table S4 (`paper/supplement.Rmd`, ~line 530)

**Before:** "However, neither consistently surpasses DeSurv at k=3: DeSurv achieves comparable or higher pooled concordance with a more parsimonious factorization"

**After:** "At k=7, NMF matches or exceeds DeSurv's per-cohort C-index in most cohorts, demonstrating that unsupervised factorization can approach supervised concordance given enough factors. However, this comes at the cost of a less parsimonious and less interpretable factorization..." — pivots to the fragmentation and classifier-overlap argument.

### 3c. Discussion (`paper/05_discussion_REVISED.Rmd`, line 10)

**Before:** "concordance generally increased with k"

**After:** "concordance increased steadily with k, approaching DeSurv's level by k=7 (SI Appendix, Table S4) but requiring more than twice as many factors and producing a fragmented structure whose prognostic content overlaps with existing molecular classifiers (SI Appendix, Tables S1–S2)."

---

## 4. Added Pooled Mean Row to Table S4

**File:** `paper/supplement.Rmd` (lines 510–518)

Added a "Pooled (mean)" row to the C-index comparison table (Table S4) averaging across all cohorts for each method column. Uses `kableExtra::kable_styling(latex_options = "hold_position")` to prevent the table from splitting across pages.

---

## 5. Expanded W-Matrix Correlation Paragraph

**File:** `paper/04_results_REVISED.Rmd` (~line 253)

**Change:** Added observation that DeSurv's most prognostic factor (D1) shows only weak, roughly uniform correlations with all three NMF factors, with mechanistic explanation: D1 couples classical tumor and restCAF-associated stromal signatures into a single survival-aligned program, so its signal is distributed across multiple NMF factors, none of which individually captures the same tumor–stroma coupling.

---

## 6. Supplement Rendering Fixes

### 6a. Figure S1 Missing / "Fig. ??" Broken Reference

**File:** `paper/supplement.Rmd` (line 12, line 125)

**Problem:** The convergence figure (S1) was not rendering, causing "Fig. ??" in Section 1 text and the figure being absent from the List of Figures.

**Fix:**
- Added `\usepackage{float}` to the YAML `header-includes`
- Changed `fig.pos='H'` to `fig.pos='htbp'` on the convergence figure chunk to avoid float placement issues
- Added `\\label{fig:fig-converge}` to the `fig.cap` string (was missing)

### 6b. Figure S2 Caption Overflow (Page 10)

**File:** `paper/supplement.Rmd` (line 173)

**Fix:** Reduced figure dimensions from `fig.height=8, out.height="8in"` to `fig.height=7, out.height="6.5in"` to prevent the long caption from overflowing into the page number.

### 6c. Figure S6 Caption Overflow (Page 13)

**File:** `paper/supplement.Rmd` (line 402)

**Fix:** Reduced from `fig.height=14` (no out.height) to `fig.height=11, out.height="8in"` so the 3-row KM figure plus caption fits on a single page.

### 6d. Section 3 Consistency Fix

**File:** `paper/supplement.Rmd` (line 258)

**Before:** "DeSurv maintains a comparable or improved C-index relative to NMF"

**After:** "DeSurv at k=3 achieves concordance comparable to NMF at higher ranks (e.g., k=7; Table S4), demonstrating that survival supervision concentrates prognostic information into fewer factors"

---

## 7. CAF Terminology Clarification

**Files:** `paper/04_results_REVISED.Rmd` (~line 249), `paper/05_discussion_REVISED.Rmd` (line 8)

**Change:** Standardized terminology from "iCAF-associated (restCAF)" to "restCAF-associated (iCAF)" and "restCAF-associated stromal signature" throughout, for consistency with the biological literature where restCAF and iCAF are distinct but related concepts.

---

## 8. Minor Text Edits

### 8a. Introduction wording (`paper/02_introduction_REVISED.Rmd`)
- Changed "without requiring their survival data" to "without requiring re-optimization" (more precise description of the projection property)

### 8b. Methods hyperparameter list (`paper/03_methods_REVISED.Rmd`)
- Changed $(k,\alpha,\lambda_H,\lambda,\xi)$ to $(k,\alpha,\lambda,\xi,n_{\text{top}})$ to match the actual BO search space

### 8c. Discussion opening (`paper/05_discussion_REVISED.Rmd`)
- Removed redundant clause "Aligning the discovery objective with the evaluation criterion ---" from the opening sentence

### 8d. Results word choice (`paper/04_results_REVISED.Rmd`)
- "per-factor survival contribution" → "per-factor contribution" (removed redundancy)
- "strongly and significantly correlated" → "strongly correlated" (simplified)

---

## 9. Paper.Rmd Structural Changes

**File:** `paper/paper.Rmd`

- **Co-authors commented out** (lines 7–16, 22–27): All co-authors and their affiliations are now commented out with `#`, leaving only Amber M. Young. This appears to be for single-author draft review.
- **Corresponding author code:** Changed from `code: 2` to `code: 1` (matches the single remaining author)
- **Section header fix:** Changed "Versioning" to "References" (line 157) — the versioning code block is `eval=FALSE` so this section heading was incorrectly labeled

---

## 10. Supplementary Methods Fixes

**File:** `paper/supp_methods.Rmd`

- **"program loadings" → "factor scores"**: Changed terminology for Z = W^T X (line 43)
- **Algorithm float positioning:** Changed `\begin{algorithm}[H]` to `\begin{algorithm}[htbp]` (line 84) to avoid float package dependency issues
- **Citation format fixes in convergence proof:** Changed `[@simon2011regularization]` and `[@seung2001algorithms; @pascualmontano2006nonsmooth; @lin2007convergence]` to inline author-year format within the LaTeX `enumerate` environment where pandoc citation processing doesn't work (lines 392, 401)
- **Citation format fixes in PDAC Datasets section:** Changed `\cite{...}` LaTeX commands to `[@...]` pandoc format for consistency (lines 822–826, 901–907)
- **Table font size:** Changed cohort table from `\small` to `\footnotesize` (line 838) for better fit
- **Table total rows:** Removed `\textit{}` from "Training total" and "Validation total" (lines 856–857)

---

## 11. Bibliography Cleanup

**File:** `paper/references_30102025.bib`

- **Removed duplicate entry:** Deleted `bailey2016genomic` which duplicated `Bailey2016` (both cited Bailey et al. 2016 Nature paper)
- **Updated SurvClust reference:** Changed `arora2020survclust` from bioRxiv preprint to the published Genome Medicine (2020) version with full journal metadata

---

## Files Modified (Summary)

| File | Type of Change |
|------|---------------|
| `paper/04_results_REVISED.Rmd` | Pooled HR computation, W-matrix paragraph, NMF C-index claims, terminology |
| `paper/05_discussion_REVISED.Rmd` | NMF concordance claim, CAF terminology, opening sentence |
| `paper/paper.Rmd` | Abstract HR, co-authors commented out, section header fix |
| `paper/supplement.Rmd` | Table S4 pooled row, NMF C-index prose, figure overflow fixes, Section 3 consistency, `\usepackage{float}`, convergence figure label |
| `paper/supp_methods.Rmd` | Terminology, algorithm float, citation format, table formatting |
| `paper/02_introduction_REVISED.Rmd` | Minor wording fix |
| `paper/03_methods_REVISED.Rmd` | Hyperparameter list correction |
| `paper/references_30102025.bib` | Removed duplicate, updated preprint to published |
