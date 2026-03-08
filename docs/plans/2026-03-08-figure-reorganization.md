# Figure Reorganization Plan for PNAS Submission

**Date:** 2026-03-08
**Status:** Implemented (2026-03-08) — on `figure-reorg` branch, pending render verification
**PNAS constraint:** Maximum 4 main-text figures

---

## Motivation

Three structural issues weaken the current figure layout:

1. **The paper's central insight is buried in the supplement.** Fig S7 (variance vs. survival contribution per factor) is the visual embodiment of the thesis — "the highest-variance factors are not the most prognostic." The main text devotes an entire paragraph to describing this result and the discussion opens with it, yet readers must navigate to the supplement to see it.

2. **Fig 2 mixes setup with results and crosses data domains.** Panels A–C (unsupervised heuristics) are *motivation*. Panel D (BO heatmap) is *method performance* on PDAC. Panel E (simulation k-histograms) is from *simulated data*. Mixing motivation, real data, and simulations in one figure dilutes impact and forces the reader to track three different contexts.

3. **The forest plot (Fig S9) — the most direct evidence for "reproducible" — is in the supplement.** The title now says "reproducible," but per-cohort consistency is only visible via inline statistics or in the supplement. The pooled KM curves (current Fig 4C–D) are a weaker visual argument than the forest plot showing consistent HRs across 5 independent cohorts.

---

## Current Layout

| Main fig | Label | Content | Panels |
|----------|-------|---------|--------|
| Fig 1 | `fig:schema` | Model schematic | A: DeSurv model, B: Bayesian optimization |
| Fig 2 | `fig:bo` | Rank selection | A: residuals, B: cophenetic, C: silhouette, D: BO heatmap, E: sim k-histogram |
| Fig 3 | `fig:sim` | Simulation (primary) | A: C-index boxes, B: precision boxes |
| Fig 4 | `fig:bio` | Biology + validation | A: DeSurv heatmap, B: NMF heatmap, C: DeSurv KM, D: NMF KM |

### Current supplement figures

| Supp fig | Label | Content |
|----------|-------|---------|
| S1 | `fig:fig-converge` | Convergence trajectories |
| S2 | `fig:sim-null-mixed` | Null + mixed simulation scenarios (6 panels) |
| S3 | `fig:fig-cindex-by-k` | C-index vs K curves |
| S4 | `fig:k3-k7-nesting` | K=3 vs K=7 factor correspondence |
| S5 | `fig:fig-cutpoint-km` | Cutpoint selection + per-cohort KM (6 panels) |
| S6 | `fig:fig-subtype-overlap` | Subtype composition of risk groups |
| S7 | `fig:varsurvival` | Variance vs survival contribution (k=3) |
| S8 | `fig:corr` | NMF vs DeSurv W-matrix correlation |
| S9 | `fig:forest` | Per-cohort HR forest plot |
| S10 | `fig:fig-nmf-k5-heatmap` | NMF gene overlap at k=5 |
| S11 | `fig:fig-nmf-k7-heatmap` | NMF gene overlap at k=7 |
| S12 | `fig:varsurvival-k5` | Variance-survival at k=5 |

---

## Proposed Layout

### New main figures

| Main fig | New label | Content | Panels | Narrative role |
|----------|-----------|---------|--------|----------------|
| **Fig 1** | `fig:schema` | Model schematic | A: DeSurv model, B: Bayesian optimization | Method introduction |
| **Fig 2** | `fig:sim` | Simulations + model selection | A: C-index boxes, B: precision boxes, C: k-selection histograms, D: BO C-index heatmap | Method validation: does it work? |
| **Fig 3** | `fig:pdac` | PDAC factor structure | A: DeSurv heatmap, B: NMF heatmap, C: W-matrix correlation, D: variance vs survival | Biological insight: what's different and why |
| **Fig 4** | `fig:val` | External validation | A: forest plot, B: DeSurv pooled KM, C: NMF pooled KM | Generalization: does it reproduce? |

### Panels promoted from supplement to main text

| Panel | Source | Destination | Rationale |
|-------|--------|-------------|-----------|
| Variance vs survival scatterplot | Fig S7 (`fig:varsurvival`) | Fig 3D | Paper's central thesis visualized; discussed in main text paragraph but figure was in supplement |
| W-matrix correlation heatmap | Fig S8 (`fig:corr`) | Fig 3C | Completes structural comparison; main text devotes full paragraph to describing it |
| Per-cohort HR forest plot | Fig S9 (`fig:forest`) | Fig 4A | Direct visual evidence for "reproducible" claim in title |

### Panels demoted from main text to supplement

| Panel | Source | Destination | Rationale |
|-------|--------|-------------|-----------|
| Unsupervised heuristics (residuals, cophenetic, silhouette) | Fig 2A–C (`fig:bo`) | New Fig S-new | Motivational setup; one text paragraph suffices; these are standard NMF diagnostics familiar to the field |

### New Fig 3 rationale (4-panel story)

| Panel | Content | Question answered |
|-------|---------|-------------------|
| A | DeSurv gene program heatmap | What does DeSurv find? |
| B | NMF gene program heatmap | What does NMF find? |
| C | W-matrix correlation heatmap | How do the two factorizations relate? |
| D | Variance vs survival scatterplot | Why does the structural difference matter? |

A–B establish the factor identities. C shows they aren't just relabelings — N2/exocrine has no DeSurv counterpart. D delivers the punchline: variance and prognosis are misaligned, and DeSurv concentrates prognostic signal where NMF doesn't.

---

## Updated supplement figure numbering

After reorganization, the supplement figures become:

| New # | Old # | Label | Content |
|-------|-------|-------|---------|
| S1 | S1 | `fig:fig-converge` | Convergence trajectories |
| S2 | S2 | `fig:sim-null-mixed` | Null + mixed simulation scenarios |
| S3 | S3 | `fig:fig-cindex-by-k` | C-index vs K curves |
| S4 | **NEW** | `fig:nmf-diagnostics` | Standard NMF rank selection heuristics (former Fig 2A–C) |
| S5 | S4 | `fig:k3-k7-nesting` | K=3 vs K=7 factor correspondence |
| S6 | S5 | `fig:fig-cutpoint-km` | Cutpoint selection + per-cohort KM |
| S7 | S6 | `fig:fig-subtype-overlap` | Subtype composition of risk groups |
| S8 | S10 | `fig:fig-nmf-k5-heatmap` | NMF gene overlap at k=5 |
| S9 | S11 | `fig:fig-nmf-k7-heatmap` | NMF gene overlap at k=7 |
| S10 | S12 | `fig:varsurvival-k5` | Variance-survival at k=5 |

Note: Former S7 (variance-survival), S8 (W-matrix correlation), and S9 (forest plot) are now in main text. Former Fig 2A–C becomes new S4. Net change: supplement shrinks from 12 to 10 figures.

---

## Implementation steps

### Phase 1: Main text figure code (`04_results_REVISED.Rmd`)

**Step 1.1 — Restructure Fig 2 (simulations + model selection)**

- Remove panels A–C (unsupervised heuristics) from the `fig-bo` chunk
- Keep panel D (BO heatmap) and panel E (k-selection histogram)
- Add simulation C-index boxes (from current `fig-sim-ntop-scenarios` chunk) and precision boxes
- Relabel: A = C-index boxes, B = precision boxes, C = k-selection histograms, D = BO heatmap
- Update `\label{fig:sim}` (or keep `fig:bo` and remap — decision needed)
- Update figure caption

**Step 1.2 — Build new Fig 3 (PDAC factor structure, 4 panels)**

- Panel A: DeSurv gene program heatmap (from current `fig-bio` chunk, `fig_gene_overlap_heatmap_desurv_tcgacptac`)
- Panel B: NMF gene program heatmap (from current `fig-bio` chunk, `fig_gene_overlap_heatmap_std_desurvk_tcgacptac`)
- Panel C: W-matrix correlation heatmap (from supplement's `fig-corr` chunk — `tar_load(fig_desurv_std_correlation_tcgacptac)`)
- Panel D: Variance vs survival scatterplot (from supplement's `fig-varsurvival` chunk — `tar_load(fig_variation_explained_tcgacptac)`)
- Assemble with `cowplot::plot_grid()` in 2×2 or 1×4 layout
- New `\label{fig:pdac}` (or similar)
- Write new caption covering all 4 panels

**Step 1.3 — Build new Fig 4 (external validation, 3 panels)**

- Panel A: Per-cohort forest plot (from supplement's `fig-forest` chunk — `tar_load(fig_hr_forest_tcgacptac)`)
- Panel B: DeSurv pooled KM (from current `fig-bio` chunk, `fig_median_survival_desurv_tcgacptac`)
- Panel C: NMF pooled KM (from current `fig-bio` chunk, `fig_median_survival_std_desurvk_tcgacptac`)
- Assemble with `plot_grid()`: forest plot on top or left, KM curves below or right
- New `\label{fig:val}` (or similar)
- Write new caption

**Step 1.4 — Remove old chunks**

- Delete the old `fig-sim-ntop-scenarios` chunk (content merged into Fig 2)
- Delete the old `fig-bio` chunk (content split across Fig 3 and Fig 4)

### Phase 2: Cross-reference updates (`04_results_REVISED.Rmd`, `05_discussion_REVISED.Rmd`, `02_introduction_REVISED.Rmd`)

All `\ref{fig:...}` references need updating. Key changes:

| Old reference | New reference | Locations |
|---------------|---------------|-----------|
| `Fig. \ref{fig:bo}A-C` | `SI Appendix, Fig. S4` (or text-only description) | Results line 157 |
| `Fig. \ref{fig:bo}D` | `Fig. \ref{fig:sim}D` | Results line 159, Discussion line 10 |
| `Fig. \ref{fig:bo}E` | `Fig. \ref{fig:sim}C` | Results line 207 |
| `Fig. \ref{fig:sim}A` | `Fig. \ref{fig:sim}A` | Results line 205 (unchanged if label reused) |
| `Fig. \ref{fig:sim}B` | `Fig. \ref{fig:sim}B` | Results line 205 (unchanged if label reused) |
| `Fig. \ref{fig:bio}A--B` | `Fig. \ref{fig:pdac}A--B` | Results lines 221, 223, 225 |
| `Fig. \ref{fig:bio}C` | `Fig. \ref{fig:val}B` | Results line 383 |
| `Fig. \ref{fig:bio}D` | `Fig. \ref{fig:val}C` | Results line 383 |
| `SI Appendix, Fig. S7` | `Fig. \ref{fig:pdac}D` | Results line 223, 227 |
| `SI Appendix, Fig. S8` | `Fig. \ref{fig:pdac}C` | Results lines 229 |
| `SI Appendix, Fig. S9` | `Fig. \ref{fig:val}A` | Results line 381 |

### Phase 3: Supplement updates (`supplement.Rmd`)

**Step 3.1 — Add demoted unsupervised heuristics figure**

- Create new chunk with the residuals, cophenetic, and silhouette plots
- Load targets: `fig_residuals_tcgacptac`, `fig_cophenetic_tcgacptac`, `fig_silhouette_tcgacptac`, `fit_std_tcgacptac`
- Assemble as 3-panel figure with shared legend
- Place in logical order (after convergence, before simulations — new Fig S4)
- Write caption: "Standard NMF rank selection heuristics yield inconsistent guidance..."

**Step 3.2 — Remove promoted figures**

- Delete `fig-varsurvival` chunk (now main Fig 3D) and its prose section
- Delete `fig-corr` chunk (now main Fig 3C) and its prose section
- Delete `fig-forest` chunk (now main Fig 4A) and its prose section

**Step 3.3 — Update supplement cross-references**

- Update any internal supplement references to S7, S8, S9 numbering
- Add reference to main text figures where supplement previously had its own

**Step 3.4 — Renumber supplement figures**

- Verify LaTeX auto-numbering handles this (it should, via `\renewcommand{\thefigure}{S\arabic{figure}}`)
- Update any hardcoded "Fig. S7" etc. references in supplement prose

### Phase 4: Prose adjustments

**Step 4.1 — Results section rank-selection paragraph (lines 155–159)**

Current text references `Fig. \ref{fig:bo}A-C` for the unsupervised heuristics. Options:
- (a) Change to `SI Appendix, Fig. S4` and keep the prose as-is
- (b) Shorten the prose since the figure is no longer visible in main text — summarize in 1–2 sentences rather than describing each heuristic individually

Recommendation: (b) — trim lines 157 to something like: "Standard NMF diagnostics (reconstruction residuals, cophenetic correlation, mean silhouette width) yielded inconsistent guidance across candidate ranks (SI Appendix, Fig. S4)." Then proceed directly to the BO heatmap description.

**Step 4.2 — Factor structure section (lines 227–229)**

Change "SI Appendix, Fig. S7" to direct figure reference. Change "SI Appendix, Fig. S8" to direct figure reference. These become main-text figure panel references, making the prose more immediate.

**Step 4.3 — Validation section (line 381)**

Change "per-cohort forest plot in SI Appendix, Fig. S9" to main-text reference. This is now a stronger claim since the reader can see the evidence directly.

### Phase 5: Abstract update

The abstract currently says "SI Appendix" nowhere, so no changes needed. But verify that the figure numbers referenced in the significance statement (if any) are correct.

### Phase 6: Verification

- [ ] Render `paper.Rmd` and verify all 4 main figures appear correctly
- [ ] Render `supplement.Rmd` and verify renumbered figures
- [ ] Grep for any remaining references to old figure labels (`fig:bo`, `fig:bio`) or old supplement numbers (`Fig. S7`, `Fig. S8`, `Fig. S9`)
- [ ] Verify discussion reference to concordance plateau now points to correct panel
- [ ] Check that no `tar_load()` calls were lost in the reorganization

---

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| Fig 3 becomes visually crowded with 4 panels | Use 2×2 grid with adequate spacing; heatmaps (A–B) on top, correlation + scatterplot on bottom |
| Fig 4 forest plot may need resizing to fit alongside KM curves | Forest plot as full-width panel A on top; KM curves side-by-side below as B–C |
| Supplement renumbering breaks existing cross-references | Systematic grep for all `Fig. S[0-9]` patterns after implementation |
| `tar_load()` calls duplicated or missing after moving chunks | Track which targets each figure needs; verify no target is loaded twice unnecessarily |

---

## Decision log

- **2026-03-08**: Plan drafted. Awaiting PI review.
- **2026-03-08**: Implemented on `figure-reorg` branch. All Rmd source changes complete. Render verification pending (requires targets store on student's machine).
- Previous analysis in `NARRATIVE_DECISIONS.md` (Decision 5) independently identified S7 and S9 as promotion candidates.
