# Figure Reorganization — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Restructure manuscript figures from 4+12 supplement → 4+10 supplement, promoting three supplement figures to main text and demoting unsupervised heuristics to supplement, per the design in `docs/plans/2026-03-08-figure-reorganization.md`.

**Architecture:** All changes are to Rmd prose and R code chunks. No new analyses or targets are required — every plot object already exists in the targets store. The work is: (1) rearrange `tar_load()` + `plot_grid()` calls in `04_results_REVISED.Rmd`, (2) add a new chunk and delete three chunks in `supplement.Rmd`, (3) update all `\ref{fig:...}` and "SI Appendix, Fig. S*" cross-references across 4 Rmd files, and (4) render both documents and verify.

**Tech Stack:** R, rmarkdown, cowplot, ggplot2, targets (read-only — no pipeline runs)

**Key files:**
- `paper/04_results_REVISED.Rmd` — main results (primary edit target)
- `paper/supplement.Rmd` — supplement (add heuristics, delete 3 promoted chunks)
- `paper/05_discussion_REVISED.Rmd` — discussion (cross-ref updates)
- `paper/02_introduction_REVISED.Rmd` — introduction (verify only, no figure refs expected)
- `paper/paper.Rmd` — parent document (render target)

**Targets store:** `store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main` (read-only via `paper/_targets.yaml`)

---

## Pre-Implementation: Create a Working Branch

### Task 0: Branch and snapshot

**Step 1: Create feature branch**

```bash
cd /home/naimrashid/Downloads/DeSurv-paper
git checkout -b figure-reorg
```

**Step 2: Commit the design plan (already untracked)**

```bash
git add docs/plans/2026-03-08-figure-reorganization.md
git add docs/plans/2026-03-08-figure-reorganization-implementation.md
git commit -m "add figure reorganization design and implementation plans"
```

---

## Phase 1: Restructure Main-Text Figures

### Task 1: Build new Fig 2 — Simulations + Model Selection

This merges the old `fig-sim-ntop-scenarios` chunk (simulation C-index + precision) with the BO heatmap and k-histogram from the old `fig-bo` chunk. The unsupervised heuristics (residuals, cophenetic, silhouette) are removed from this figure.

**Files:**
- Modify: `paper/04_results_REVISED.Rmd` — `fig-bo` chunk (~lines 161–197) and `fig-sim-ntop-scenarios` chunk (~lines 209–215)

**Step 1: Replace the `fig-bo` chunk**

Replace the entire `fig-bo` code chunk (lines 161–198) with:

````r
```{r fig-sim, fig.width = 6.5, fig.height = 6, fig.cap = "Method validation in simulation and PDAC data. (A--B) Primary simulation scenario ($p = 3{,}000$ genes, $n = 200$ samples, true $k = 3$, 100 replicates), where prognostic programs explain low variance relative to outcome-neutral programs. Both methods were tuned via Bayesian optimization over cross-validated concordance index. (A) Test-set concordance index across simulation replicates. (B) Precision (fraction of top-weighted genes in a learned factor that overlap with a true prognostic program's gene set). (C) Distribution of selected $k$ values across simulation replicates for DeSurv versus standard NMF ($\\alpha = 0$). (D) Gaussian process predicted mean cross-validated C-index from Bayesian optimization over the joint space of factorization rank ($k$) and supervision strength ($\\alpha$) in PDAC training data (TCGA + CPTAC). \\label{fig:sim}", fig.env='figure*', fig.pos='t', out.height= "6in", out.width='\\textwidth'}

tar_load(fig_bo_heat_tcgacptac)
fig_bo_heat_tcgacptac <- set_fig_font(fig_bo_heat_tcgacptac, size = 10)

# Top row: simulation C-index + precision (from alt_plots, loaded in setup chunk)
sim_row <- plot_grid(
  alt_plots$cindex_box + labs(title = NULL) + scale_y_continuous(limits = c(.5, 1)),
  alt_plots$precision_box + labs(title = NULL),
  ncol = 2, labels = c("A", "B")
)

# Bottom row: k-selection histogram + BO heatmap
k_hist <- plot_grid(alt_plots$k_hist, NULL, nrow = 2, rel_heights = c(4, 1))
bo_row <- plot_grid(
  k_hist,
  fig_bo_heat_tcgacptac,
  ncol = 2, labels = c("C", "D"), rel_widths = c(2, 3)
)

plot_grid(sim_row, bo_row, nrow = 2, rel_heights = c(1, 1.2))
```
````

Note: the old `tar_load` calls for `fig_residuals_tcgacptac`, `fig_cophenetic_tcgacptac`, `fig_silhouette_tcgacptac`, and `fit_std_tcgacptac` are removed from this chunk (they move to supplement in Task 4). Only `fig_bo_heat_tcgacptac` remains.

**Step 2: Delete the old `fig-sim-ntop-scenarios` chunk**

Delete the entire chunk at ~lines 209–215 (the standalone simulation figure), since its content is now incorporated into the new `fig-sim` chunk above. Keep the simulation prose sections above and below it — they still apply.

**Step 3: Update the figure label reference in prose**

In the rank-selection prose paragraph (~line 161 caption area), the label changes from `\label{fig:bo}` to `\label{fig:sim}`. The old `fig-sim-ntop-scenarios` chunk had `\label{fig:sim}` — since we deleted that chunk, the label is now free for the new combined figure.

**Step 4: Verify no duplicate labels**

Run:
```bash
grep -n '\\\\label{fig:sim}' paper/04_results_REVISED.Rmd
```
Expected: exactly 1 match (the new chunk).

**Step 5: Commit**

```bash
git add paper/04_results_REVISED.Rmd
git commit -m "fig reorg: build new Fig 2 (simulations + model selection)"
```

---

### Task 2: Build new Fig 3 — PDAC Factor Structure (4 panels)

This creates a new 4-panel figure combining:
- A: DeSurv gene program heatmap (from old `fig-bio` chunk)
- B: NMF gene program heatmap (from old `fig-bio` chunk)
- C: W-matrix correlation (from supplement `fig-corr`)
- D: Variance vs survival scatterplot (from supplement `fig-varsurvival`)

**Files:**
- Modify: `paper/04_results_REVISED.Rmd` — replace `fig-bio` chunk (~lines 230–374)

**Step 1: Replace the `fig-bio` chunk**

The old `fig-bio` chunk contains panels A-B (heatmaps) + C-D (KM curves). We keep A-B and replace C-D with the promoted supplement panels. The KM curves move to the new Fig 4 (Task 3).

Replace the entire `fig-bio` code chunk (lines 230–374) with:

````r
```{r fig-pdac, fig.cap="Survival supervision produces different factor structures in PDAC. (A--B) Spearman correlations between factor gene rankings and established PDAC gene programs for DeSurv (A) and standard NMF (B) at $k = 3$. Asterisks denote significance after multiple testing correction; only correlations $> 0.2$ shown. (C) Pairwise correspondence between NMF and DeSurv factors, measured by Spearman rank correlation of W-matrix gene loadings across all genes. (D) Fraction of expression variance explained versus survival contribution ($\\Delta$ Cox partial log-likelihood upon factor removal) for each factor. DeSurv concentrates survival signal into D1 while the highest-variance NMF factor (N2, exocrine-associated) contributes negligibly to survival. \\label{fig:pdac}", fig.env='figure*', fig.pos='t', out.height= "7in", out.width='\\textwidth'}

library(cowplot)
library(ggplot2)

# Panels A-B: gene program heatmaps (already loaded targets in setup or load here)
tar_load(fig_gene_overlap_heatmap_desurv_tcgacptac)
tar_load(fig_gene_overlap_heatmap_std_desurvk_tcgacptac)

# Panel C: W-matrix correlation (promoted from supplement)
tar_load(fig_desurv_std_correlation_tcgacptac)

# Panel D: variance vs survival (promoted from supplement; plot object already loaded
# in computed-stats chunk for statistics — reload here for the plot)
tar_load(fig_variation_explained_tcgacptac)

# ── Panels A-B: heatmaps with shared legend ──
plot_3a <- fig_gene_overlap_heatmap_desurv_tcgacptac$plot
plot_3b <- fig_gene_overlap_heatmap_std_desurvk_tcgacptac$plot
legend_ab <- fig_gene_overlap_heatmap_desurv_tcgacptac$legend
legend_ab <- gtable::gtable_add_padding(
  legend_ab, padding = unit(c(0, 10, 0, 0), "pt")
)

top_row <- plot_grid(
  plot_3a, plot_3b, cowplot::ggdraw(legend_ab),
  ncol = 3, labels = c("A", "B", ""), align = "hv",
  label_size = 12, rel_widths = c(2.5, 2.5, 0.6)
)

# ── Panel C: W-matrix correlation ──
plot_3c <- plot_grid(
  fig_desurv_std_correlation_tcgacptac$plot,
  plot_grid(
    NULL,
    cowplot::ggdraw(fig_desurv_std_correlation_tcgacptac$legend),
    nrow = 2, rel_heights = c(0.08, 0.92)
  ),
  ncol = 2, rel_widths = c(4, 1)
)

# ── Panel D: variance vs survival scatterplot ──
plot_3d <- set_fig_font(fig_variation_explained_tcgacptac, size = 10) +
  theme(
    legend.position    = c(1, 0.5),
    legend.justification = c(1, 0),
    legend.background  = element_rect(color = "black"),
    axis.title         = element_text(size = 8)
  )

bottom_row <- plot_grid(
  plot_3c, plot_3d,
  ncol = 2, labels = c("C", "D"),
  label_size = 12, rel_widths = c(1.2, 1)
)

plot_grid(top_row, bottom_row, nrow = 2, rel_heights = c(1.4, 1))
```
````

**Step 2: Remove the helper functions that were only used for old KM panels**

The `km_extract_label`, `km_parse_stats`, `stack_surv`, and the KM-related code blocks (~lines 241–374 of the old chunk) are no longer in this figure. They move to the new Fig 4 chunk (Task 3). Keep them available by placing them in the new Fig 4 chunk or a shared setup chunk.

**Step 3: Verify the `tar_load` in the `computed-stats` chunk still works**

The `computed-stats` chunk (~line 88) already loads `fig_variation_explained_tcgacptac` for statistics. The new `fig-pdac` chunk also loads it for the plot. Loading twice is safe (`tar_load` is idempotent).

**Step 4: Commit**

```bash
git add paper/04_results_REVISED.Rmd
git commit -m "fig reorg: build new Fig 3 (PDAC factor structure, 4 panels)"
```

---

### Task 3: Build new Fig 4 — External Validation (3 panels)

This creates a 3-panel figure:
- A: Per-cohort forest plot (promoted from supplement `fig-forest`)
- B: DeSurv pooled KM (from old `fig-bio` panels C-D)
- C: NMF pooled KM (from old `fig-bio` panels C-D)

**Files:**
- Modify: `paper/04_results_REVISED.Rmd` — add new chunk after the validation prose section (~after line 383)

**Step 1: Add new `fig-val` chunk**

Insert a new code chunk after the validation prose paragraph (after current line 383, before the "Factorization rank K=3 is robust" section). This chunk uses the KM helper functions previously in `fig-bio` and adds the forest plot.

````r
```{r fig-val, fig.cap="Survival-aligned programs generalize across independent PDAC cohorts. (A) Forest plot of per-cohort hazard ratios (95\\% CIs) for the latent factor most predictive of survival in the training data. DeSurv (blue) shows more consistent associations across cohorts than standard NMF (red). (B--C) Kaplan--Meier curves for pooled external validation samples stratified by the cross-validated optimal cutpoint: DeSurv linear predictor (B) and standard NMF (C). \\label{fig:val}", fig.env='figure*', fig.pos='t', out.height= "6in", out.width='\\textwidth'}

library(cowplot)
library(ggplot2)

# ── Panel A: forest plot (promoted from supplement) ──
tar_load(fig_hr_forest_tcgacptac)

base_size <- 8
theme_pnas <- theme_classic(base_size = base_size) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.margin = margin(6, 10, 6, 15),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(6, "pt"),
    legend.key.width = unit(12, "pt")
  )

plot_forest <- fig_hr_forest_tcgacptac +
  scale_color_discrete(
    labels = c(
      "Puleo_array"       = "Puleo",
      "Dijk"              = "Dijk",
      "PACA_AU_seq"       = "PACA seq",
      "PACA_AU_array"     = "PACA array",
      "Moffitt_GEO_array" = "Moffitt"
    )
  ) +
  theme_pnas +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = base_size),
    legend.title = element_text(size = base_size)
  )

# ── Panels B-C: KM curves (moved from old fig-bio) ──
tar_load(fig_median_survival_desurv_tcgacptac)
tar_load(fig_median_survival_std_desurvk_tcgacptac)

km_text_size    <- base_size - 3
risk_title_size <- 5

stack_surv <- function(surv_obj, title, keep_legend = FALSE) {
  p <- surv_obj$plot +
    theme_pnas +
    theme(
      legend.position = if (keep_legend) "bottom" else "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = km_text_size),
      axis.text.x  = element_text(size = km_text_size),
      axis.text.y  = element_text(size = km_text_size),
      plot.margin  = margin(10, 10, 2, 10),
      plot.title   = element_text(size = risk_title_size + 3)
    ) +
    labs(title = title)

  t <- surv_obj$table +
    theme_pnas +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.x  = element_text(size = km_text_size),
      axis.text.y  = element_text(size = km_text_size),
      text = element_text(size = km_text_size),
      axis.title.x = element_text(size = km_text_size),
      plot.title   = element_text(size = risk_title_size),
      plot.margin  = margin(10, 10, 6, 10)
    ) +
    labs(x = "Time (months)")
  t$layers[[1]]$aes_params$size <- km_text_size / ggplot2::.pt
  plot_grid(p, t, ncol = 1, rel_heights = c(2.5, 1.5), align = "v", axis = "lr")
}

km_b <- stack_surv(fig_median_survival_desurv_tcgacptac, "DeSurv")
km_c <- stack_surv(fig_median_survival_std_desurvk_tcgacptac, "NMF")

# Shared KM legend
legend_plot <- ggplot(
  data.frame(x = 1:2, y = 1:2,
             group = factor(c("Low", "High"), levels = c("Low", "High"))),
  aes(x = x, y = y, colour = group)
) +
  geom_line() +
  scale_colour_manual(values = c("Low" = "violetred2", "High" = "turquoise4"),
                      name = "Risk group") +
  theme_pnas +
  theme(legend.position = "bottom",
        legend.text = element_text(size = km_text_size),
        legend.title = element_text(size = km_text_size)) +
  guides(colour = guide_legend(nrow = 1))

legend_grob <- ggplotGrob(legend_plot)
km_legend <- legend_grob$grobs[
  sapply(legend_grob$grobs, function(x) x$name) == "guide-box"
][[1]]

km_block <- plot_grid(
  km_b, km_c, ggdraw(km_legend),
  nrow = 3, labels = c("B", "C", ""),
  label_size = 12, rel_heights = c(5, 5, 1)
)

# ── Full figure: forest (left/top) + KM curves (right/bottom) ──
plot_grid(plot_forest, km_block, ncol = 2, rel_widths = c(1.2, 1))
```
````

**Step 2: Commit**

```bash
git add paper/04_results_REVISED.Rmd
git commit -m "fig reorg: build new Fig 4 (external validation, 3 panels)"
```

---

## Phase 2: Cross-Reference Updates (Main Text)

### Task 4: Update all figure references in results prose

**Files:**
- Modify: `paper/04_results_REVISED.Rmd`

**Reference mapping (apply all changes):**

| Old reference | New reference | Location in prose |
|---------------|---------------|-------------------|
| `Fig. \ref{fig:bo}A-C` | `SI Appendix, Fig. S4` | Rank-selection paragraph (~line 157) |
| `Fig. \ref{fig:bo}D` | `Fig. \ref{fig:sim}D` | Rank-selection paragraph (~line 159) |
| `Fig. \ref{fig:bo}E` | `Fig. \ref{fig:sim}C` | Simulation prose (~line 207) |
| `Fig. \ref{fig:sim}A` | `Fig. \ref{fig:sim}A` | Simulation prose (~line 205) — label reused, no change needed |
| `Fig. \ref{fig:sim}B` | `Fig. \ref{fig:sim}B` | Simulation prose (~line 205) — no change needed |
| `Fig. \ref{fig:bio}A--B` | `Fig. \ref{fig:pdac}A--B` | Factor structure section (~lines 221, 223, 225) |
| `Fig. \ref{fig:bio}C` | `Fig. \ref{fig:val}B` | Validation section (~line 383) |
| `Fig. \ref{fig:bio}D` | `Fig. \ref{fig:val}C` | Validation section (~line 383) |
| `SI Appendix, Fig. S7` | `Fig. \ref{fig:pdac}D` | Factor structure section (~lines 223, 227) |
| `SI Appendix, Fig. S8` | `Fig. \ref{fig:pdac}C` | Factor structure section (~line 229) |
| `SI Appendix, Fig. S9` | `Fig. \ref{fig:val}A` | Validation section (~line 381) |

**Step 1: Apply each substitution**

Use find-and-replace carefully. Some references appear multiple times. Key replacements:

1. Line 157: Change `(Fig. \ref{fig:bo}A-C)` → `(SI Appendix, Fig. S4)`
2. Line 159: Change `(Fig. \ref{fig:bo}D)` → `(Fig. \ref{fig:sim}D)`
3. Line 207: Change `(Fig. \ref{fig:bo}E)` → `(Fig. \ref{fig:sim}C)`
4. All occurrences of `\ref{fig:bio}` → `\ref{fig:pdac}` (for panels A, B) or `\ref{fig:val}` (for panels C, D)
5. Line 223: Change `(SI Appendix, Fig. S7)` → `(Fig. \ref{fig:pdac}D)`
6. Line 227: Change `(SI Appendix, Fig. S7; SI Methods)` → `(Fig. \ref{fig:pdac}D; SI Methods)`
7. Line 229: Change `(SI Appendix, Fig. S8)` → `(Fig. \ref{fig:pdac}C)`
8. Line 381: Change `SI Appendix, Fig. S9` → `Fig. \ref{fig:val}A`

**Step 2: Shorten the heuristics paragraph (optional but recommended)**

Consider trimming lines 157 to:
> "Standard NMF diagnostics (reconstruction residuals, cophenetic correlation, mean silhouette width) yielded inconsistent guidance across candidate ranks (SI Appendix, Fig. S4)."

This replaces the 3-sentence description of each heuristic. The figure is now in the supplement, so a brief summary suffices.

**Step 3: Verify no stale references remain**

```bash
grep -n 'fig:bo' paper/04_results_REVISED.Rmd
grep -n 'fig:bio' paper/04_results_REVISED.Rmd
grep -n 'Fig\. S7' paper/04_results_REVISED.Rmd
grep -n 'Fig\. S8' paper/04_results_REVISED.Rmd
grep -n 'Fig\. S9' paper/04_results_REVISED.Rmd
```

Expected: all return 0 matches (or only appear inside comments/SI numbering header).

**Step 4: Commit**

```bash
git add paper/04_results_REVISED.Rmd
git commit -m "fig reorg: update all cross-references in results"
```

---

### Task 5: Update discussion cross-references

**Files:**
- Modify: `paper/05_discussion_REVISED.Rmd`

**Step 1: Find and replace `fig:bo` references**

Line 10 references `Fig. \ref{fig:bo}D` (the concordance plateau). Change to `Fig. \ref{fig:sim}D`.

```bash
grep -n 'fig:bo' paper/05_discussion_REVISED.Rmd
```

Apply each match.

**Step 2: Check for any supplement figure references that changed**

```bash
grep -n 'Fig\. S[0-9]' paper/05_discussion_REVISED.Rmd
```

Update any S7/S8/S9 references if present.

**Step 3: Verify**

```bash
grep -n 'fig:bo\|fig:bio' paper/05_discussion_REVISED.Rmd
```

Expected: 0 matches.

**Step 4: Commit**

```bash
git add paper/05_discussion_REVISED.Rmd
git commit -m "fig reorg: update discussion cross-references"
```

---

## Phase 3: Supplement Updates

### Task 6: Add demoted NMF diagnostics figure to supplement

The unsupervised heuristics (residuals, cophenetic, silhouette) that were removed from main-text Fig 2 become new supplement Fig S4.

**Files:**
- Modify: `paper/supplement.Rmd` — insert new chunk after the `fig-cindex-by-k` section (~after line 249)

**Step 1: Insert new section and chunk**

Insert after the "Cross-validated C-index across factorization rank" section (after the `fig-cindex-by-k` chunk, before the "K-sensitivity analysis" section):

````markdown
## Standard NMF rank selection heuristics {.unnumbered}

Standard unsupervised rank selection criteria — reconstruction residuals, cophenetic correlation coefficient, and mean silhouette width — yielded inconsistent guidance for selecting the factorization rank $k$ in the PDAC training cohort (TCGA + CPTAC). Reconstruction residuals decreased smoothly without a clear elbow (Fig. \ref{fig:nmf-diagnostics}A). Cophenetic correlation fluctuated without a distinct transition point (Fig. \ref{fig:nmf-diagnostics}B). Mean silhouette width favored very low ranks that conflicted with other criteria (Fig. \ref{fig:nmf-diagnostics}C). These contradictions motivated the survival-supervised rank selection approach described in the main text.

```{r fig-nmf-diagnostics, fig.width = 6, fig.height = 3.5, fig.cap = "Standard unsupervised NMF rank selection heuristics for the PDAC training cohort (TCGA + CPTAC). (A) Reconstruction residuals as a function of factorization rank $k$. (B) Cophenetic correlation coefficient as a function of $k$. (C) Mean silhouette width across multiple distance metrics as a function of $k$. These criteria yield inconsistent guidance, motivating the survival-supervised model selection approach (main text Fig. 2D). \\label{fig:nmf-diagnostics}", fig.env='figure*', fig.pos='t', out.height= "3.5in", out.width='\\textwidth'}
library(cowplot)
library(ggplot2)

tar_load(fig_residuals_tcgacptac)
tar_load(fig_cophenetic_tcgacptac)
tar_load(fig_silhouette_tcgacptac)
tar_load(fit_std_tcgacptac)

fig_residuals_tcgacptac  <- set_fig_font(fig_residuals_tcgacptac, size = 10)
fig_cophenetic_tcgacptac <- set_fig_font(fig_cophenetic_tcgacptac, size = 10)
fig_silhouette_tcgacptac <- set_fig_font(fig_silhouette_tcgacptac, size = 10)

# Shared legend from fit_std plot
legend_plot <- plot(fit_std_tcgacptac) +
  ggplot2::theme(legend.position = "bottom")
legend_grob <- ggplotGrob(legend_plot)
legend <- legend_grob$grobs[
  sapply(legend_grob$grobs, function(x) x$name) == "guide-box"
][[1]]

panels <- plot_grid(
  fig_residuals_tcgacptac  + theme(legend.position = "none"),
  fig_cophenetic_tcgacptac + theme(legend.position = "none"),
  fig_silhouette_tcgacptac + theme(legend.position = "none"),
  ncol = 3, labels = c("A", "B", "C")
)

plot_grid(panels, legend, nrow = 2, rel_heights = c(1, 0.15))
```
````

**Step 2: Commit**

```bash
git add paper/supplement.Rmd
git commit -m "fig reorg: add demoted NMF diagnostics as supplement Fig S4"
```

---

### Task 7: Delete promoted figures from supplement

Three supplement sections are now in the main text and must be removed.

**Files:**
- Modify: `paper/supplement.Rmd`

**Step 1: Delete the variance-survival section (~lines 419–433)**

Delete from `## Variance explained versus survival contribution per factor` through the end of the `fig-varsurvival` chunk (up to but not including `## DeSurv and NMF factor correspondence`).

**Step 2: Delete the W-matrix correlation section (~lines 435–451)**

Delete from `## DeSurv and NMF factor correspondence` through the end of the `fig-corr` chunk (up to but not including `## Per-cohort hazard ratios`).

**Step 3: Delete the forest plot section (~lines 453–480)**

Delete from `## Per-cohort hazard ratios for the most prognostic factor` through the end of the `fig-forest` chunk (up to but not including `## Standard NMF factor structure at independently selected rank`).

**Step 4: Verify the deleted labels are gone**

```bash
grep -n 'fig:varsurvival}' paper/supplement.Rmd
grep -n 'fig:corr}' paper/supplement.Rmd
grep -n 'fig:forest}' paper/supplement.Rmd
```

Expected: 0 matches for `fig:varsurvival` (note: `fig:varsurvival-k5` should still exist). 0 for `fig:corr`. 0 for `fig:forest`.

**Step 5: Commit**

```bash
git add paper/supplement.Rmd
git commit -m "fig reorg: remove promoted figures (S7, S8, S9) from supplement"
```

---

### Task 8: Update supplement internal cross-references

After removing S7-S9 and inserting new S4, the internal supplement references need updating.

**Files:**
- Modify: `paper/supplement.Rmd`

**Step 1: Inventory all hardcoded supplement figure references**

```bash
grep -n 'Fig\. S[0-9]' paper/supplement.Rmd
grep -n 'Fig\.~S[0-9]' paper/supplement.Rmd
grep -n 'Fig\. \\\\ref' paper/supplement.Rmd
```

LaTeX auto-numbering via `\renewcommand{\thefigure}{S\arabic{figure}}` handles most references *if* they use `\ref{fig:...}`. Only hardcoded "Fig. S7" etc. strings need manual updating.

**Step 2: Update hardcoded references**

Key patterns to fix:
- Any remaining `Fig. S7` that referred to variance-survival → change to `main text Fig. 3D`
- Any remaining `Fig. S8` that referred to W-matrix corr → change to `main text Fig. 3C`
- Any remaining `Fig. S9` that referred to forest plot → change to `main text Fig. 4A`
- References to `Fig. \ref{fig:varsurvival}` within supplement prose (e.g., in the NMF k=5 section, line 547) → change to `main text Fig. \ref{fig:pdac}D` or `main text Fig. 3D`

Specific known references:
1. Supplement line ~437 (in `fig-corr` prose): References `Fig. \ref{fig:varsurvival}` → This prose section is now deleted (Task 7), so no action needed.
2. Supplement line ~547: "The variance--survival decomposition at $k = 5$ (Fig. \ref{fig:varsurvival-k5})" → This label is unchanged (varsurvival-k5 stays in supplement). No action.
3. Supplement line ~545: "consistent with the mixed simulation scenario described earlier in this supplement (Fig. \ref{fig:sim-null-mixed})" → Unchanged. No action.

**Step 3: Verify**

```bash
grep -n 'Fig\. S[789]' paper/supplement.Rmd
grep -n 'varsurvival}' paper/supplement.Rmd   # should show only varsurvival-k5
grep -n 'fig:corr' paper/supplement.Rmd       # should show 0
grep -n 'fig:forest' paper/supplement.Rmd     # should show 0
```

**Step 4: Commit**

```bash
git add paper/supplement.Rmd
git commit -m "fig reorg: update supplement internal cross-references"
```

---

### Task 9: Update the SI figure numbering comment header in results

**Files:**
- Modify: `paper/04_results_REVISED.Rmd` — the comment block at top (~lines 5–15)

**Step 1: Replace the SI numbering comment**

Replace the existing comment block with:

```html
<!-- SI Appendix figure numbering (supplement.Rmd order):
     S1  fig-converge         convergence trajectories
     S2  fig-sim-null-mixed   null and mixed simulation scenarios
     S3  fig-cindex-by-k      C-index vs K sensitivity
     S4  fig-nmf-diagnostics  standard NMF rank selection heuristics (NEW — former main Fig 2A-C)
     S5  fig-k3-k7-nesting    K=3 vs K=7 factor correspondence
     S6  fig-cutpoint-km      cutpoint selection and KM curves
     S7  fig-subtype-overlap  subtype composition of risk groups
     S8  fig-nmf-k5-heatmap   NMF gene overlap at k=5
     S9  fig-nmf-k7-heatmap   NMF gene overlap at k=7
     S10 fig-varsurvival-k5   variance-survival at k=5

     Former S7 (variance-survival), S8 (W-corr), S9 (forest) promoted to main Figs 3-4.
-->
```

**Step 2: Commit**

```bash
git add paper/04_results_REVISED.Rmd
git commit -m "fig reorg: update SI figure numbering comment header"
```

---

## Phase 4: Verification

### Task 10: Grep verification — no stale references

**Step 1: Check for stale figure labels across all Rmd files**

```bash
grep -rn 'fig:bo' paper/*.Rmd
grep -rn 'fig:bio' paper/*.Rmd
```

Expected: 0 matches in active files (matches in `paper/old/` are fine).

**Step 2: Check for stale hardcoded supplement numbers**

```bash
grep -rn 'Fig\. S7' paper/04_results_REVISED.Rmd paper/05_discussion_REVISED.Rmd paper/02_introduction_REVISED.Rmd
grep -rn 'Fig\. S8' paper/04_results_REVISED.Rmd paper/05_discussion_REVISED.Rmd paper/02_introduction_REVISED.Rmd
grep -rn 'Fig\. S9' paper/04_results_REVISED.Rmd paper/05_discussion_REVISED.Rmd paper/02_introduction_REVISED.Rmd
```

Expected: 0 matches in main-text files.

**Step 3: Check for orphaned `tar_load` calls**

Verify every `tar_load()` in the results file loads a target that is actually used in a figure or computation:

```bash
grep -n 'tar_load' paper/04_results_REVISED.Rmd
```

Cross-check: `fig_residuals_tcgacptac`, `fig_cophenetic_tcgacptac`, `fig_silhouette_tcgacptac`, `fit_std_tcgacptac` should NOT appear (moved to supplement). `fig_desurv_std_correlation_tcgacptac` and `fig_hr_forest_tcgacptac` SHOULD appear (promoted to main text).

**Step 4: Check for duplicate figure labels**

```bash
grep -oh '\\\\label{fig:[^}]*}' paper/04_results_REVISED.Rmd | sort | uniq -c | sort -rn
```

Expected: each label appears exactly once. Labels should be: `fig:schema`, `fig:sim`, `fig:pdac`, `fig:val`.

**Step 5: Commit verification notes (optional)**

No commit needed — this is a check step.

---

### Task 11: Render main paper

**Step 1: Render paper.Rmd**

```bash
cd /home/naimrashid/Downloads/DeSurv-paper
Rscript -e 'rmarkdown::render("paper/paper.Rmd", knit_root_dir = getwd())'
```

**Step 2: Check for errors**

If rendering fails, read the error output. Common issues:
- Missing `tar_load` target → check target name spelling
- Duplicate label → grep for the label
- Undefined function → check that helper functions (like `stack_surv`) are defined before use in the chunk where they're called

**Step 3: Visual inspection of output PDF**

Open `paper/paper.pdf` and verify:
- Fig 1: Model schematic (unchanged)
- Fig 2: 4 panels — A (C-index boxes), B (precision boxes), C (k-histogram), D (BO heatmap)
- Fig 3: 4 panels — A (DeSurv heatmap), B (NMF heatmap), C (W-matrix correlation), D (variance-survival scatter)
- Fig 4: 3 panels — A (forest plot), B (DeSurv KM), C (NMF KM)
- No Fig 5 or Fig 6
- All cross-references in prose point to correct figure/panel
- No "SI Appendix, Fig. S7/S8/S9" references remain in main text

**Step 4: Commit rendered PDF (if policy allows)**

```bash
git add paper/paper.pdf
git commit -m "fig reorg: re-render paper with new figure layout"
```

---

### Task 12: Render supplement

**Step 1: Render supplement.Rmd**

```bash
Rscript -e 'rmarkdown::render("paper/supplement.Rmd", knit_root_dir = getwd())'
```

**Step 2: Visual inspection**

Open `paper/supp_methods.pdf` (or `paper/supplement.pdf`) and verify:
- S1: Convergence trajectories
- S2: Null + mixed simulation
- S3: C-index vs K
- **S4: NMF diagnostics (NEW — residuals, cophenetic, silhouette)**
- S5: K=3 vs K=7 nesting (was S4)
- S6: Cutpoint + KM curves (was S5)
- S7: Subtype overlap (was S6)
- S8: NMF k=5 heatmap (was S10)
- S9: NMF k=7 heatmap (was S11)
- S10: Variance-survival k=5 (was S12)
- **S7 (variance-survival), S8 (W-corr), S9 (forest) are GONE**

**Step 3: Verify figure count**

10 supplement figures total (down from 12).

**Step 4: Commit**

```bash
git add paper/supplement.pdf paper/supp_methods.pdf
git commit -m "fig reorg: re-render supplement with updated figure numbering"
```

---

## Phase 5: Documentation Cleanup

### Task 13: Update NARRATIVE_ARC.md

The narrative arc document still describes a 5-prediction, 6-figure structure. It needs to reflect the current paper (4 predictions, 4 figures, no bladder).

**Files:**
- Modify: `NARRATIVE_ARC.md`

**Key changes:**
1. Remove P5 (cross-cancer transfer) from the prediction table (line 27)
2. Update "three evaluation settings (simulation, PDAC, cross-cancer)" → "two evaluation settings (simulation, PDAC)" (line 47)
3. Remove the entire "Figure 6: Cross-Cancer Transfer" section (lines 211–225)
4. Update the evidence cascade (line 103) to end at PDAC validation
5. Remove bladder from the Three-Act Story diagram (lines 96–100)
6. Update the summary predictions (lines 338–345) — remove item 5
7. Update the Figure-Narrative Mapping Table (lines 353–361) — remove row 6
8. Update figure numbers throughout to match new 4-figure layout:
   - Old Fig 2 → new Fig 2 (content changed)
   - Old Fig 3 → merged into new Fig 2
   - Old Fig 4 → new Fig 3 (with promoted S7, S8)
   - Old Fig 5 → new Fig 4 (with promoted S9)
9. Update the "Narrative Gaps" section — remove cross-cancer items
10. Remove "single transfer pair (PDAC→bladder)" from limitations (line 114)

**Step 1: Apply all changes listed above**

**Step 2: Commit**

```bash
git add NARRATIVE_ARC.md
git commit -m "update narrative arc to reflect 4-prediction, 4-figure structure"
```

---

### Task 14: Update the design plan status

**Files:**
- Modify: `docs/plans/2026-03-08-figure-reorganization.md`

**Step 1: Change status from "Proposed" to "Implemented"**

Update line 4:
```
**Status:** Implemented — 2026-03-XX
```

**Step 2: Add decision log entry**

Append to the decision log:
```
- **2026-03-XX**: Plan implemented. Bladder results (Fig 6) previously removed from paper; plan covers Figs 1-4 only. NARRATIVE_ARC.md updated to 4-prediction structure.
```

**Step 3: Commit**

```bash
git add docs/plans/2026-03-08-figure-reorganization.md
git commit -m "mark figure reorganization plan as implemented"
```

---

## Verification Checklist (Final)

Run after all tasks are complete:

- [ ] `grep -rn 'fig:bo' paper/*.Rmd` → 0 matches (excluding `old/` and comments)
- [ ] `grep -rn 'fig:bio' paper/*.Rmd` → 0 matches
- [ ] `grep -rn 'Fig\. S[789][^0-9]' paper/04_results_REVISED.Rmd paper/05_discussion_REVISED.Rmd` → 0 matches
- [ ] `grep -rn 'fig:varsurvival[^-]' paper/supplement.Rmd` → 0 matches (only `varsurvival-k5` survives)
- [ ] `grep -rn 'fig:corr' paper/supplement.Rmd` → 0 matches
- [ ] `grep -rn 'fig:forest' paper/supplement.Rmd` → 0 matches
- [ ] `paper/paper.pdf` renders without error
- [ ] `paper/supplement.pdf` renders without error
- [ ] Main text has exactly 4 figures (schema, sim, pdac, val)
- [ ] Supplement has exactly 10 figures (S1–S10)
- [ ] No `tar_load()` calls reference targets that aren't used
- [ ] Discussion references `fig:sim` not `fig:bo`
- [ ] All KM panels in main text have HR/CI/p-value annotations (from existing `km_extract_label` code)
- [ ] NARRATIVE_ARC.md matches the current 4-prediction, 4-figure structure

---

## Risk Mitigation Notes

| Risk | Mitigation |
|------|------------|
| `tar_load` fails because target not in store | All targets used already exist in the active store — no new targets needed. Test with `Rscript -e 'targets::tar_read("target_name")'` before rendering. |
| ggplot objects from HPC fail to render locally | Known issue (see MEMORY.md). If `fig_variation_explained_tcgacptac` errors, rebuild from raw data using the `computed-stats` chunk approach. |
| Fig 3 visually crowded with 4 panels | Use `rel_heights = c(1.4, 1)` to give heatmaps more space. Iterate on dimensions after first render. |
| Duplicate `\label{}` causes LaTeX error | The grep checks in Task 10 catch this before rendering. |
| Supplement renumbering confuses readers of earlier drafts | The SI numbering comment header (Task 9) documents the mapping for future reference. |
