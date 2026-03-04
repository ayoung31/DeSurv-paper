# DeSurv Paper — Open Narrative & Structural Decisions

*Created 2026-03-04. Revisit after Amber addresses issues #26–#30.*

## Completed

- Laura Peng added as middle author (Pharmacology + Lineberger)
- Jen Jen Yeh added as 2nd-to-last author (Pharmacology + Surgery + Lineberger)
- Author contributions updated
- GitHub issue #29: iCAF/myCAF vs restCAF/proCAF discussion (for Laura/Jen Jen)
- GitHub issue #30: NMF rank selection procedure (for Amber)

---

## Decision 1: D3 "Basal-like-associated" label

**Problem**: D3's W-matrix is most similar to NMF's *classical* factor (N1), not a basal factor. D3 absorbs residual exocrine variation + basal programs, contributes minimally to survival. Label may mislead reviewers into expecting basal-like correlations that won't appear.

**Options**:
- (a) Soften to descriptive: "D3 (Residual tumor/exocrine)"
- (b) Keep but explain anti-correlation: absence of classical ≈ basal-like
- (c) Wait for Laura/Jen Jen input (issue #29) — coordinates with restCAF/proCAF nomenclature

**Recommendation**: (c) — decide after issue #29 discussion

**Files**: `targets_figure_configs.R`, `paper/04_results_REVISED.Rmd` (lines 191–197), `paper/supplement.Rmd` (line 408+)

---

## Decision 2: Missing NMF HR/p for "heterogeneity" claim

**Problem**: Line 323 says NMF showed "greater heterogeneity across datasets and weaker survival associations" but only reports DeSurv's continuous HR/CI/p. NMF's values are missing.

**Fix**: Add NMF continuous risk score pooled HR, CI, and p-value inline. Values should come from the validation analysis target (ask Amber which object contains these).

**File**: `paper/04_results_REVISED.Rmd` (line 323)

---

## Decision 3: K=3 vs K=7 paragraph contradiction

**Problem**: Paragraph 16 (line 374) says "Without this rule, DeSurv would have selected k=7" — implying k=7 is better. Paragraphs 18–19 (lines 378–380) then show k=7 fails validation at ntop=270. Reads as contradictory.

**Options**:
- (a) Remove paragraph 16 entirely
- (b) Merge into paragraph 17 opening: "BO selected k=3 via the 1-SE rule, favoring it over the BO global max k=7 (margin Y). Whether this parsimony is warranted..."
- (c) Keep both but add bridging sentence

**Recommendation**: (b) — merge eliminates the standalone question-without-answer

**File**: `paper/04_results_REVISED.Rmd` (lines 374–376)

---

## Decision 4: Results section ordering

**Current order**: Model overview → Rank selection + simulations → PDAC factor characterization → External validation → K-sensitivity

**Finding**: PNAS NMF/subtyping precedents (Brunet 2004, etc.) use logical progression. Current order follows this convention.

**Problem**: External validation (strongest result) is only 293 words vs K-sensitivity at 902 words. The payoff is compressed.

**Options**:
- (a) Keep order but rebalance: expand validation, move K-sensitivity detail to SI
- (b) Move validation before factor characterization: sims → validation → factor interpretation
- (c) No change

**Recommendation**: (a) — rebalance without reordering. Move K-sensitivity tables/detail to SI, keep one-paragraph summary in main text.

---

## Decision 5: Figure placement (supplement → main)

**Current main figures (4, at PNAS limit)**: Fig 1 (schema), Fig 2 (rank selection), Fig 3 (simulations), Fig 4 (gene programs + KM)

**Supplement figures that could strengthen main text**:
- **Fig S9 (forest plot)**: Per-cohort HR heterogeneity — direct visual proof of generalization
- **Fig S7 (variance vs survival)**: The paper's core thesis visualized
- **Fig S3 (C-index vs k)**: Concordance plateau evidence

**Options**:
- (a) Swap Fig S9 for Fig 2A–C (move NMF diagnostics to SI, bring forest plot to main)
- (b) Swap Fig S7 for Fig 3 (move simulation plots to SI, bring variance-vs-survival to main)
- (c) Combine panels — merge S9 into Fig 4 as panel E
- (d) No change

**Recommendation**: Discuss with Amber. (a) is strongest — forest plot is more impactful than cophenetic/silhouette for PNAS reviewers.

---

## Decision 6: Cherry-picking concern (NMF at k=3)

**Concern**: Comparing DeSurv k=3 to NMF k=3 when NMF's own BO selects k=7 could look unfair.

**Already mitigated**: Paragraph 16 states NMF also selects k=7; Table S4 has NMF at independent ranks; K-sensitivity shows k=7 fails with gene focusing.

**Additional options**:
- Add sentence to validation section noting NMF at BO-selected k was also evaluated (→ Table S4)
- Emphasize k=3 comparison is *most favorable* to NMF (fewer factors = less overfitting)
- Show NMF k=7 validation results in main text alongside NMF k=3

---

## Decision 7: Narrative compellingness

**Diagnosis**: Paper reads as technical report rather than story. "Variance ≠ prognosis" thesis is powerful but diluted by methodological defense.

**Suggestions**:
- Open each results subsection with the *finding*, not the *method*
- Reduce K-sensitivity in main text to one paragraph
- Consider a take-home sentence at the end of each subsection
- Use Discussion more aggressively (already well-written)

---

## Related GitHub issues

| Issue | Topic | Assigned to |
|-------|-------|-------------|
| #26 | Figure rendering fixes | Amber |
| #27 | `tar_params_best` α bug | Amber |
| #28 | Render + review revised PDFs | Amber |
| #29 | iCAF/myCAF nomenclature | Amber (for Laura/Jen Jen) |
| #30 | NMF rank selection procedure | Amber |
