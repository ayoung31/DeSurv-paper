# Co-Author Revisions Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Address all 39 co-author comments from DL, Alisa, and Peng, including combining the two supplement files into one SI Appendix document and updating all cross-references.

**Architecture:** Changes span the main paper (4 REVISED child .Rmd files + paper.Rmd YAML), the two supplement files (merge into one), and `supp_methods.Rmd` (formatting fixes). The main paper is rendered from `paper/paper.Rmd` which includes `_REVISED` child documents. The supplements are standalone rendered documents.

**Tech Stack:** R Markdown, LaTeX, `targets` (read-only; do NOT re-run pipeline)

**Key constraint:** Do NOT re-run any pipeline targets. All figure/table objects are pre-computed in the targets store. Only edit `.Rmd` prose, YAML, LaTeX, and figure captions.

---

## Phase 1: Combine Supplements into Single SI Appendix

This is the highest-impact structural change. Currently there are two separate files:
- `paper/supp_methods.Rmd` (951 lines) — math derivations, algorithms, proofs, dataset details
- `paper/supplement.Rmd` (689 lines) — figures S1-S10, tables S1-S5

The PI directed combining them into one file. The combined document should place Methods first (since the main text references "SI Appendix" for both), then Results.

### Task 1: Create combined `paper/si_appendix.Rmd`

**Files:**
- Create: `paper/si_appendix.Rmd`
- Reference: `paper/supp_methods.Rmd` (source for Part 1)
- Reference: `paper/supplement.Rmd` (source for Part 2)

**Step 1:** Create `paper/si_appendix.Rmd` with:
- YAML header titled "SI Appendix" (not "Supplementary Methods" or "Supplementary Results")
- Combine `header-includes` from both files (union of LaTeX packages)
- Include `bibliography: references_30102025.bib` from supp_methods
- Add `\renewcommand{\thefigure}{S\arabic{figure}}` and `\renewcommand{\thetable}{S\arabic{table}}` from supplement.Rmd
- Part I heading: "Supplementary Methods" containing all content from `supp_methods.Rmd` (sections 1-12)
- Part II heading: "Supplementary Results" containing all content from `supplement.Rmd`
- Single TOC, single List of Figures, single List of Tables

**Step 2:** Verify the combined file renders:
```bash
cd /home/naimrashid/Downloads/DeSurv-paper
Rscript -e 'rmarkdown::render("paper/si_appendix.Rmd", knit_root_dir = getwd())'
```
Expected: PDF output with both methods and results sections, continuous figure/table numbering.

**Step 3:** Commit
```bash
git add paper/si_appendix.Rmd
git commit -m "combine supp_methods + supplement into single SI Appendix"
```

**IMPORTANT:** Do NOT delete the original `supplement.Rmd` or `supp_methods.Rmd` yet — keep them as reference until rendering is verified.

### Task 2: Renumber supplement algorithms as S1, S2, S3

**Files:**
- Modify: `paper/si_appendix.Rmd` (the algorithms section inherited from supp_methods)

DL flagged that Algorithm 1, 2, 3 should be S1, S2, S3 since they are in the supplement.

**Step 1:** Add to the LaTeX preamble of `si_appendix.Rmd`:
```latex
\renewcommand{\thealgorithm}{S\arabic{algorithm}}
```

**Step 2:** Verify algorithm labels render as "Algorithm S1", "Algorithm S2", "Algorithm S3".

**Step 3:** Commit.

---

## Phase 2: Main Paper — Introduction Edits

All edits to `paper/02_introduction_REVISED.Rmd`.

### Task 3: Replace "Sufficient" with "Supervised" (DL comment, p1)

**File:** `paper/02_introduction_REVISED.Rmd:19`

DL: "Our work is technically not a sufficient dimension reduction approach"

**Step 1:** In the paragraph starting "Sufficient dimension reduction theory...", change the opening to avoid claiming DeSurv IS sufficient dimension reduction. The paragraph should still reference the SDR theory as motivation but frame DeSurv as "supervised" or "outcome-guided":

Change:
```
Sufficient dimension reduction theory establishes that response-guided subspace estimation targets directions most relevant to the outcome, whereas variance-maximizing projections can miss outcome-relevant structure entirely
```
To:
```
Supervised dimension reduction theory establishes that response-guided subspace estimation targets directions most relevant to the outcome, whereas variance-maximizing projections can miss outcome-relevant structure entirely
```

Also update `paper/paper.Rmd` line 48 (Significance Statement) — change "nuisance variation" to "noise variation" per Alisa's comment. (Note: "Sufficient" does not appear in the significance statement in the REVISED version, only "nuisance".)

**Step 2:** Commit.

### Task 4: Define X ≈ WH before discussing W vs H (DL + Alisa)

**File:** `paper/02_introduction_REVISED.Rmd:21`

DL: "We haven't defined W or H, how about we briefly define X as data matrix and X ≈ WH?"
Alisa: "Insert a sentence about Mixture = W * H"

**Step 1:** In paragraph 4 (line 21), before "The key architectural choice is that the survival gradient acts on the gene program matrix $W$...", insert:

```
DeSurv approximates the expression matrix as $X \approx WH$, where columns of $W$ define nonnegative gene programs and rows of $H$ encode sample-level loadings.
```

The full sentence flow becomes:
"Here we present DeSurv, a survival-supervised deconvolution framework that integrates NMF with Cox proportional hazards modeling. DeSurv approximates the expression matrix as $X \approx WH$, where columns of $W$ define nonnegative gene programs and rows of $H$ encode sample-level loadings. The key architectural choice is that the survival gradient acts on..."

**Step 2:** Consider removing the α/k technical detail from the intro per DL's suggestion. The sentence "A supervision parameter $\alpha$ and the factorization rank $k$ are jointly selected via cross-validated concordance (Methods)." can be simplified to just "(Methods)" or a brief forward reference. **Decision for PI:** Accept or reject removing this detail.

**Step 3:** Commit.

### Task 5: Add missing citations per Peng (5 comments)

**File:** `paper/02_introduction_REVISED.Rmd`, `paper/references_30102025.bib`

Peng's 5 citation requests (all on Introduction or Results):

1. **Line 15**, after "identification of molecular subtypes [@collisson2011subtypes]" — Add PurIST [@rashid2020purity], DeCAF [@peng2024determination], and Moffitt [@moffitt2015virtual]. Check: `rashid2020purity` and `peng2024determination` already appear later in the intro. Adding them here acknowledges this work earlier. Also add PMID:30718832 (Maurer 2019 — check if already in .bib as `Maurer2019`).

2. **Line 15**, after "key biological insights, including the identification of molecular subtypes [@collisson2011subtypes]" — Also cite Bailey [@Bailey2016] and Moffitt [@moffitt2015virtual]. Check: both are already cited on line 17. Move their first appearance earlier.

3. **Results section** (`04_results_REVISED.Rmd:250`), where "restCAF-associated (iCAF) stromal signatures" and "proCAF-associated" are mentioned — cite DeCAF paper [@peng2024determination] for restCAF/proCAF terminology.

4. **Results section**, where "iCAF" is mentioned — cite Elyada et al. PMID:31197017 (likely `elyada2019cross` in .bib).

5. **Results section**, where "myCAF" or "SCISSORS panCAF and myCAF" appears — cite SCISSORS paper PMID:37498558. Need to check .bib for this entry; if missing, add it.

**Step 1:** Verify each citation key exists in `references_30102025.bib`:
```bash
grep -c "Maurer2019\|elyada2019cross\|Bailey2016\|rashid2020purity\|peng2024determination" paper/references_30102025.bib
```
Add any missing entries.

**Step 2:** Add citations at the locations specified above.

**Step 3:** Check that SCISSORS paper (PMID:37498558, Leary et al. 2023) is in the .bib. If not, add the BibTeX entry.

**Step 4:** Commit.

---

## Phase 3: Main Paper — Results Edits

All edits to `paper/04_results_REVISED.Rmd`.

### Task 6: Define BO acronym at first use (DL, 3 comments)

**File:** `paper/04_results_REVISED.Rmd`

BO is first used implicitly in the Results. DL wants it defined. The acronym "Bayesian optimization (BO)" should appear at its first use.

**Step 1:** Find the first mention of "Bayesian optimization" in the Results. It appears on line 201: "Applying the 1-SE rule to the Bayesian optimization surface...". Change to "Applying the 1-SE rule to the Bayesian optimization (BO) surface..." and use "BO" for all subsequent occurrences in the Results.

Actually, "Bayesian optimization" first appears in the Methods section (line 26 of `03_methods_REVISED.Rmd`): "Bayesian optimization (BO)". It IS already defined there. Check if Methods appears before Results in the rendered document — in PNAS format, Methods comes AFTER Discussion. So DL is right: BO is used in Results before Methods defines it.

**Fix:** Define BO at its first use in Results (line 201), and verify the Methods section doesn't redundantly redefine it — or keep both definitions since PNAS readers may jump to Methods first.

**Step 2:** Replace subsequent uses of "Bayesian optimization" with "BO" in Results where appropriate (e.g., line 248: "Bayesian optimization selected" → "BO selected"; line 602: any remaining full form).

**Step 3:** Commit.

### Task 7: Fix Fig 3 caption — asterisks and correlation threshold (DL, 2 comments)

**File:** `paper/04_results_REVISED.Rmd:257`

DL: "I couldn't find any asterisks" and "Do we really mean absolute value > 0.2?"

**Step 1:** In the Fig 3 caption (line 257, `fig.cap` string), find:
```
Asterisks denote significance after multiple testing correction; only correlations $> 0.2$ shown.
```
Two options:
- **(a)** If asterisks truly don't appear in the rendered heatmaps, REMOVE the asterisk sentence entirely.
- **(b)** If the threshold is meant to be |correlation| > 0.2, change to `only $|r| > 0.2$ shown`.

**Decision for PI:** Verify whether asterisks appear in the actual rendered Fig 3A-B heatmaps. If not, delete the asterisk clause. Also confirm whether the filter is absolute correlation or signed.

**Recommended edit** (assuming no asterisks, absolute value filter):
```
(A--B) Spearman correlations between factor gene rankings and established PDAC gene programs for DeSurv (A) and standard NMF (B) at $k = 3$; only $|r| > 0.2$ shown.
```

**Step 2:** Commit.

### Task 8: Add summary sentences to Fig 3 legend (Alisa)

**File:** `paper/04_results_REVISED.Rmd:257`

Alisa: "just looking at Figure 3, A and B, this is not immediately clear. Maybe add a sentence in Figure3 legend for A and B to summarize the main points."

**Step 1:** After the existing Fig 3 caption text for panels A-B, add a one-sentence summary:
```
DeSurv concentrates survival signal into D1 (classical tumor + restCAF stroma) while suppressing exocrine-compositional variation; standard NMF devotes one factor to exocrine signals, distributing prognostic content across all three factors.
```

**Step 2:** Commit.

### Task 9: Reorder Fig 4A forest plot panels (DL)

**File:** `paper/04_results_REVISED.Rmd:404-408`

DL: "Any reason for this order? Why not D1 D2 D3?"

The current factor ordering in the forest plot is D3, D2, D1 (bottom to top) and N3, N2, N1. DL wants D1 at top.

**Step 1:** Change the factor level ordering in lines 404-408:
```r
desurv_dat$factor_name <- factor(desurv_dat$factor_name,
  levels = c("D1", "D2", "D3"))
...
nmf_dat$factor_name <- factor(nmf_dat$factor_name,
  levels = c("N1", "N2", "N3"))
```

**Decision for PI:** The current bottom-to-top ordering (D3→D1) follows the convention of most prognostic factor at bottom with HR to the left. Top-to-bottom (D1→D3) is more intuitive for reading order. Accept DL's suggestion?

**Step 2:** Commit.

### Task 10: Emphasize ~100-patient reclassification (Alisa)

**File:** `paper/04_results_REVISED.Rmd:352`

Alisa: "not just a better p-value but a difference of 100 patients being classified into a different risk group"

**Step 1:** After the sentence about DeSurv vs NMF stratification, add a sentence quantifying the reclassification difference. The KM plots show different group sizes between DeSurv and NMF cutpoints. Extract the numbers from the KM data (already loaded as `fig_median_survival_desurv_tcgacptac` and `fig_median_survival_std_desurvk_tcgacptac`).

Add something like:
```
The DeSurv cutpoint classified [X] of [N] validation patients as high-risk versus [Y] by NMF, reclassifying approximately [Z] patients between risk groups.
```

This requires extracting numbers from the KM plot objects. Check whether the risk table data is accessible from the stored `survfit` objects.

**Step 2:** Commit.

### Task 11: Add Fig S4 legend clarification (Alisa)

**File:** `paper/supplement.Rmd` (now `si_appendix.Rmd`) — Fig S4 caption

Alisa: "For the figure S4 legend, maybe spell out explicitly why the figures suggest a different k"

**Step 1:** In the Fig S4 caption, after "These criteria yield inconsistent guidance", add:
```
Residuals decrease smoothly without a clear elbow, cophenetic correlation fluctuates without a distinct transition, and mean silhouette width favors very low ranks, providing no consensus on a single optimal $k$.
```

**Step 2:** Commit.

---

## Phase 4: Main Paper — Introduction/Methods/Discussion Cross-Cutting Edits

### Task 12: Update all "SI Appendix" references to be specific (DL)

**Files:**
- `paper/04_results_REVISED.Rmd` — lines 199, 241, 248, 254, 350, 602
- `paper/05_discussion_REVISED.Rmd` — lines 10, 14
- `paper/03_methods_REVISED.Rmd` — lines 23, 26, 28, 37, 40

DL: "Which file, which section?" and "Do we need better names?"

Now that the supplements are combined into one SI Appendix, all "SI Appendix" references are unambiguous as to FILE. But DL also wants section specificity.

**Step 1:** Replace vague "SI Appendix" with section-specific references where possible:
- "convergence to a stationary point under mild conditions (SI Appendix)" → "(SI Appendix, Section 3)"
- "Complete derivations and algorithmic details are provided in the SI Appendix" → "...in SI Appendix, Sections 1–2"
- "consensus-based initialization... (SI Appendix)" → "(SI Appendix, Section 5)"
- "details in SI Appendix" for gene truncation → "(SI Appendix, Section 5)"
- "Further training, validation, and runtime details appear in the SI Appendix" → "...SI Appendix, Section 7"
- "SI Methods" in line 40 → "SI Appendix, Section 7" (since there is no separate SI Methods file anymore)

**Step 2:** Ensure "Fig. S" and "Table S" references are preceded by "SI Appendix" consistently:
- Line 350: "Fig. S6C--F" → "SI Appendix, Fig. S6C--F"

**Step 3:** Commit.

### Task 13: Change \tilde{W} to \widetilde{W} throughout (DL, 3 comments)

**Files:**
- `paper/02_introduction_REVISED.Rmd:21` — `$\tilde{W}$`
- `paper/04_results_REVISED.Rmd:348` — `$\tilde{W}$`
- `paper/03_methods_REVISED.Rmd:26, 28, 34` — `$\tilde{W}$`
- `paper/si_appendix.Rmd` — any `\tilde{W}` instances

**Step 1:** Global find-and-replace `\tilde{W}` → `\widetilde{W}` across all `.Rmd` files in `paper/`.

**Step 2:** Commit.

### Task 14: Fix re-expanded acronyms (DL — NMF, PDAC, BO)

**Files:** `paper/03_methods_REVISED.Rmd`

DL: "NMF has been defined" (p5) and "PDAC... already defined in intro" (p6)

**Step 1:** In `03_methods_REVISED.Rmd`:
- Line 12: "Nonnegative Matrix Factorization (NMF)" → "NMF" (already defined in intro)
- Line 37: "pancreatic ductal adenocarcinoma (PDAC)" → "PDAC" (already defined in intro)

**Step 2:** Commit.

### Task 15: Mention convergence theorem explicitly (DL)

**File:** `paper/03_methods_REVISED.Rmd:23`

DL: "Shall we explicitly say we have a theorem for it?"

**Step 1:** Change:
```
these updates converge to a stationary point under mild conditions (SI Appendix)
```
To:
```
these updates converge to a stationary point under mild regularity conditions (SI Appendix, Theorem 1)
```

Verify that the convergence result in `supp_methods.Rmd` is labeled "Theorem 1". Check line ~300+ of supp_methods.Rmd.

**Step 2:** Commit.

---

## Phase 5: Supplement Formatting Fixes (DL comments on supp_methods)

All changes to `paper/si_appendix.Rmd` (the combined file).

### Task 16: Fix 5 missing periods in derivations

**File:** `paper/si_appendix.Rmd` — inherited from supp_methods page 5

DL flagged 5 locations on page 5 of the Supp Methods PDF (Section 2.6, derivation of W update). These are likely missing terminal periods after displayed equations.

**Step 1:** Review the W-update derivation section (around lines 350-450 of supp_methods content). Add periods after any equations that end sentences.

**Step 2:** Commit.

### Task 17: Fix equation numbering — remove (8)-(15) if unreferenced

**File:** `paper/si_appendix.Rmd`

DL: "If (8)-(15) are never referred later, we may not number them."

**Step 1:** Check whether equations (8)-(15) are referenced anywhere. If not, convert them from `\begin{equation}` to `\begin{equation*}` (unnumbered).

**Step 2:** Commit.

### Task 18: Fix Matérn accent, β spacing, indicator typos

**File:** `paper/si_appendix.Rmd`

1. Line 758 (inherited): "Matern" → "Mat\\'ern" (add accent)
2. Missing space before β (supp methods p8)
3. Indicator symbol typos (p10, 2 instances of "⊮") — likely should be `\mathbb{1}` or `\mathbbm{1}`

**Step 1:** Fix each of the three issues.

**Step 2:** Commit.

### Task 19: Consolidate Table S2/S3 α headers and remove duplicate Fig S3 legend

**File:** `paper/si_appendix.Rmd`

DL: "Do we really need 2 rows for α values?" (Tables S2/S3) and "same legend repeated" (Fig S3)

**Step 1:** In Tables S2 and S3, consolidate the α header to a single row if possible.

**Step 2:** For Fig S3, remove the duplicated legend that appears for both panels.

**Step 3:** Commit.

### Task 20: Standardize k/K notation (DL)

**File:** `paper/si_appendix.Rmd`

DL: "sometimes K=3 in plain text, sometimes $k=3$ in math"

**Step 1:** Global convention: use `$k = 3$` in math mode everywhere. Find plain-text "K=3", "K=5", "K=7" etc. and convert to `$k = 3$`, `$k = 5$`, `$k = 7$`.

Exception: When referring to the sensitivity analysis grid labels (e.g., "K=3" as a row label in tables), keep as-is since those are categorical labels.

**Step 2:** Commit.

---

## Phase 6: Significance Statement and Abstract Edits

### Task 21: Word choice — "dominated" and "nuisance" (Alisa)

**File:** `paper/paper.Rmd:44-48`

Alisa flagged "dominated" and "nuisance" as words to replace with synonyms.

These appear in the abstract (line 44) and significance statement (line 48).

**Step 1:**
- "dominated by variance-dominant" → "driven by high-variance" (abstract line 44)
- "nuisance variation" → "noise variation" (significance statement line 48)

Note: "nuisance" also appears in `05_discussion_REVISED.Rmd:6` — change there too if desired. **Decision for PI:** "nuisance variation" is standard statistical terminology. Accept or reject this change?

**Step 2:** Commit.

---

## Phase 7: Figure 1 Schematic Fixes (DL)

### Task 22: Fix Fig 1 multiplication symbol and β alignment

**File:** `figures/model_schematic_final.pdf` (source needed)

DL: "Used × for X=WH but · elsewhere — be consistent. Also move β vector up for alignment."

**Decision for PI:** This requires editing the original figure source (likely PowerPoint, Illustrator, or tikz). If the schematic was generated by code, identify the source file and fix. If it was created externally, flag for Amber to fix.

**Step 1:** Identify the source of `figures/model_schematic_final.pdf`.

**Step 2:** Fix the multiplication symbol inconsistency and β alignment.

**Step 3:** Commit.

### Task 23: Add (α) after "strength" in Fig 2 caption (DL)

**File:** `paper/04_results_REVISED.Rmd:203`

DL: "Add (α) after 'strength'"

**Step 1:** In the Fig 2 caption, find "supervision strength" and change to "supervision strength ($\\alpha$)".

**Step 2:** Commit.

---

## Phase 8: Additional Alisa Comments

### Task 24: Add citation near "survival modeling" (Alisa)

**File:** `paper/02_introduction_REVISED.Rmd:15`

Alisa: "Maybe add a citation or two here" near "survival modeling"

**Step 1:** Add relevant survival modeling citations. Candidates: Bair & Tibshirani 2004 [@bair2004semi] is already cited in P2. Consider adding a PDAC-specific survival reference here, e.g., Aung/COMPASS [@aung2018compass].

**Decision for PI:** Which citation(s) to add?

**Step 2:** Commit.

### Task 25: Add "this work provides some evidence" sentence (Alisa)

**File:** `paper/02_introduction_REVISED.Rmd:19`

Alisa: "Maybe add sentence to the effect that this work provides some evidence..."

After the paragraph ending "remains an open empirical question", the next paragraph ("Here we present DeSurv...") answers it. The transition is already implicit. **Decision for PI:** Add an explicit bridge sentence, or is the current flow sufficient?

If adding, insert at end of paragraph 3:
```
The present work provides empirical evidence that survival-guided NMF recovers prognostic programs more reliably than unsupervised factorization in this setting.
```

**Step 2:** Commit.

### Task 26: Add methods cross-reference for simulation details (Alisa)

**File:** `paper/04_results_REVISED.Rmd:237`

Alisa: "Make a reference to the methods section and/or briefly re-define X, δ, y"

**Step 1:** The simulation paragraph already says "Methods" at the end of the first sentence. Consider adding "(see Materials and Methods, Simulation Studies)" for clarity. The variables $p$, $n$, $k$ are already defined inline.

**Step 2:** Commit.

---

## Render and Verify

### Task 27: Render all documents and verify

**Step 1:** Render the combined SI Appendix:
```bash
Rscript -e 'rmarkdown::render("paper/si_appendix.Rmd", knit_root_dir = getwd())'
```

**Step 2:** Render the main paper:
```bash
Rscript -e 'rmarkdown::render("paper/paper.Rmd", knit_root_dir = getwd())'
```

**Step 3:** Verify:
- All "SI Appendix" references are consistent
- Algorithm numbers show S1, S2, S3
- Figure/table numbering is S1, S2, etc.
- No LaTeX compilation errors
- All citations resolve

**Step 4:** Final commit.

---

## Comment Disposition Summary

| # | Reviewer | Comment | Disposition | Task |
|---|----------|---------|-------------|------|
| 1 | DL | Capitalize N? | **Reject** — lowercase n is standard for sample size | — |
| 2 | DL | "Sufficient" → "Supervised" | **Accept** | 3 |
| 3 | DL | \widetilde vs \tilde | **Accept** | 13 |
| 4 | DL | Define X ≈ WH before W/H discussion | **Accept** | 4 |
| 5 | DL | Remove α/k from intro | **PI decision** | 4 |
| 6 | DL | SI Appendix naming | **Accept** — merge solves this | 1, 12 |
| 7 | DL | Define BO acronym | **Accept** | 6 |
| 8 | DL | Use BO consistently | **Accept** | 6 |
| 9 | DL | \widetilde (2nd instance) | **Accept** | 13 |
| 10 | DL | NMF already defined | **Accept** | 14 |
| 11 | DL | SI — which file/section | **Accept** | 12 |
| 12 | DL | SI Appendix — which file/section | **Accept** | 12 |
| 13 | DL | BO should be defined earlier | **Accept** | 6 |
| 14 | DL | "converge" → reference theorem | **Accept** | 15 |
| 15 | DL | PDAC already defined | **Accept** | 14 |
| 16 | DL | Fig 1 × vs · and β alignment | **Accept** — needs source file | 22 |
| 17 | DL | Fig 2 add (α) after "strength" | **Accept** | 23 |
| 18 | DL | Fig 3 — no asterisks found | **Accept** — remove clause | 7 |
| 19 | DL | |correlation| > 0.2? | **Accept** — clarify as absolute | 7 |
| 20 | DL | Fig 4 — D1 D2 D3 order | **PI decision** | 9 |
| 21 | DL | Equation into 1-2 lines | **Accept** | 16 |
| 22 | DL | Algorithm 1 → S1 | **Accept** | 2 |
| 23 | DL | Algorithm 2 → S2 | **Accept** | 2 |
| 24 | DL | 5 missing periods | **Accept** | 16 |
| 25 | DL | Don't number eqs (8)-(15) | **Accept** | 17 |
| 26 | DL | Missing space before β | **Accept** | 18 |
| 27 | DL | Indicator typos (×2) | **Accept** | 18 |
| 28 | DL | Algorithm 3 → S3 | **Accept** | 2 |
| 29 | DL | Matérn accent | **Accept** | 18 |
| 30 | DL | Table S2/S3 α header rows | **Accept** | 19 |
| 31 | DL | Same as Table S2 | **Accept** | 19 |
| 32 | DL | k vs K consistency | **Accept** | 20 |
| 33 | DL | Duplicate Fig S3 legend | **Accept** | 19 |
| 34 | Alisa | "dominated" → synonym | **PI decision** | 21 |
| 35 | Alisa | "nuisance" → "noise" | **PI decision** | 21 |
| 36 | Alisa | Add citations near modeling | **Accept** | 24 |
| 37 | Alisa | Add evidence sentence | **PI decision** | 25 |
| 38 | Alisa | Define X, W, H (same as DL #4) | **Accept** | 4 |
| 39 | Alisa | Or methods reference | **Accept** | 26 |
| 40 | Alisa | Fig S4 legend clarity | **Accept** | 11 |
| 41 | Alisa | Fig 3 legend summary | **Accept** | 8 |
| 42 | Alisa | ~100 patient reclassification | **Accept** | 10 |
| 43 | Peng | Cite PurIST, DeCAF, PMID:30718832 | **Accept** | 5 |
| 44 | Peng | Cite Moffitt, Bailey | **Accept** | 5 |
| 45 | Peng | Cite DeCAF for Maurer etc. | **Accept** | 5 |
| 46 | Peng | Cite Elyada (iCAF) | **Accept** | 5 |
| 47 | Peng | Cite SCISSORS (myCAF) | **Accept** | 5 |

**Items requiring PI decision (5):**
1. Remove α/k technical detail from intro? (Task 4, DL #5)
2. Reorder Fig 4 forest panels D1→D3? (Task 9, DL #20)
3. Replace "dominated"? (Task 21, Alisa #34)
4. Replace "nuisance" with "noise"? (Task 21, Alisa #35) — note: "nuisance" is standard stats terminology
5. Add explicit bridge sentence about "this work provides evidence"? (Task 25, Alisa #37)
