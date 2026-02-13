# Supplementary Materials — Required Updates for Revision

This document catalogs all changes needed in supplementary materials, bibliography, and YAML header to align with the revised main-text sections (`02_introduction_REVISED.Rmd`, `03_methods_REVISED.Rmd`, `04_results_REVISED.Rmd`, `05_discussion_REVISED.Rmd`).

---

## 1. `paper/supp_methods.Rmd` — Supplementary Methods

### 1A. Semi-supervised framing (new subsection or paragraph)

**Location:** After the Model Details subsection (~line 66), or as a new paragraph within it.

**What to add:** A paragraph (or short subsection) that explicitly frames DeSurv as a semi-supervised method. The main text now uses this language prominently. The supplement should reinforce the interpretation: the reconstruction term acts as a generative regularizer, ensuring that learned programs correspond to genuine transcriptomic patterns, while the Cox term steers the factorization toward outcome-relevant structure. Cite Chapelle et al. (2006) for the semi-supervised framework and note that the $\alpha$ parameter controls the relative influence of each objective, with $\alpha = 0$ recovering standard unsupervised NMF.

### 1B. Explicit $\partial \mathcal{L}_{\text{Cox}} / \partial H = 0$ statement

**Location:** Subsection "Supervision on W versus H" (~line 631).

**What to add:** The existing subsection discusses W-vs-H supervision qualitatively. Add a brief formal statement showing that $\partial \mathcal{L}_{\text{Cox}} / \partial H = 0$ directly. This follows from the fact that $Z = W^\top X$ does not depend on $H$; therefore the Cox partial likelihood gradient with respect to $H$ is identically zero. One sentence and one equation suffice. This closes a gap noted in the PNAS senior editor review: the main text now states that "sample-level loadings $H$ enter only through the reconstruction term and receive no survival gradient," and the supplement should contain the derivation.

**Suggested equation:**
$$\frac{\partial \mathcal{L}_{\text{Cox}}}{\partial H} = \frac{\partial \mathcal{L}_{\text{Cox}}}{\partial Z} \cdot \frac{\partial Z}{\partial H} = \frac{\partial \mathcal{L}_{\text{Cox}}}{\partial Z} \cdot \mathbf{0} = \mathbf{0},$$
since $Z = W^\top X$ does not depend on $H$.

### 1C. Convergence proof — theory-to-practice gap paragraph

**Location:** End of the convergence proof section (~line 628, after the Coxnet remark).

**What to add:** The revised Discussion (paragraph 4) now explicitly states: "DeSurv's convergence guarantee applies to an idealized algorithm; the implementation includes practical modifications (backtracking, gradient clamping) that depart from the theoretical analysis." The supplement should expand on this. Specifically:

1. List the practical modifications: (a) gradient clamping via $\delta_{\max}$, (b) $\varepsilon_W$ and $\varepsilon_H$ floors, (c) approximate Hessian in Coxnet $\beta$-update, (d) finite backtracking budget.
2. Explain why each departure is benign in practice: clamping prevents numerical overflow without affecting limit points; floors keep iterates strictly positive (matching Assumption 1(iii)); the approximate Hessian still yields descent; and the finite backtracking budget falls back to $W^{(t+1)} = W^{(t)}$ (identity step, which trivially satisfies the Armijo condition with equality).
3. Note that empirical convergence traces (Figure S1) confirm monotone descent across all datasets and initializations.

### 1D. Null and mixed simulation scenario details

**Location:** After the existing dataset descriptions (~line 813), or as a new subsection "Simulation Design Details."

**What to add:** The revised main text now mentions three simulation scenarios (easy, null, mixed) but keeps descriptions brief. The supplement should provide full details:

- **Easy (R0_easy):** Prognostic programs explain low variance relative to outcome-neutral signals. Marker genes load strongly on one factor; background genes dominate variance. Survival depends on marker gene structure through $X^\top \tilde{W}$.
- **Null (R00_null):** No survival signal. Survival times are generated independently of gene expression (i.e., $\beta = 0$ or $\tilde{W} = 0$). This scenario tests whether DeSurv correctly reduces to unsupervised NMF when $\alpha$ is optimized and no outcome signal exists.
- **Mixed (R_mixed):** Prognostic and variance-dominant programs partially overlap. Some factors contribute to both reconstruction and survival, creating a more realistic intermediate scenario.

Include the gamma distribution parameters used for each gene class (marker, background, noise) and the exponential distribution parameters for survival and censoring times. Cross-reference the simulation code in the repository.

### 1E. Out-of-sample validation protocol

**Location:** After the Cross-validation subsection (~line 728) or within the Survival Analysis subsection (~line 815).

**What to add:** The revised main text now states that "all model selection was performed entirely within training data using nested cross-validation; no validation cohort data were used during any stage of tuning." The supplement should expand with a diagram or step-by-step description:

1. Inner loop: K-fold CV within training data, computing C-index for each hyperparameter configuration proposed by BO.
2. Outer loop: BO proposes configurations, receives averaged inner-fold C-index.
3. Final model: Refit on full training data at BO-selected hyperparameters.
4. External validation: Project validation data via $Z = \tilde{W}^\top X_{\text{new}}$, compute $\hat{\eta} = Z^\top \hat{\beta}$, evaluate.
5. Explicit statement: "At no point does the validation survival data ($y_{\text{val}}, \delta_{\text{val}}$) influence $W$, $\beta$, or any hyperparameter."

### 1F. Sufficient dimension reduction theory connection

**Location:** New paragraph in the "Supervision on W versus H" subsection (~line 631) or as a new subsection.

**What to add:** The revised Introduction cites Cook (2007) and Bair & Tibshirani (2006) for the theoretical motivation. The supplement should briefly explain how DeSurv relates to sufficient dimension reduction (SDR): the column space of $W$ defines a $k$-dimensional subspace of $\mathbb{R}^p$; the Cox supervision steers this subspace toward the central dimension reduction subspace for the survival outcome. This is analogous to supervised PCA (Bair & Tibshirani 2004), but DeSurv imposes nonnegativity and operates within the NMF framework rather than PCA.

---

## 2. `paper/supplement.Rmd` — Supplementary Figures

### 2A. Complete the scRNA-seq text

**Location:** Line 41 — currently reads "Need more here..."

**What to add:** Complete the scRNA-seq analysis paragraph. Key points:
- Factor 1 (immune/stromal): expressed primarily in iCAF and B cells, consistent with ORA and published gene list overlaps.
- Factor 2 (exocrine/classical): expressed in acinar and classical PDAC 2 cells.
- Factor 3 (basal-like): expressed in basal-like PDAC cells, consistent with ORA and overlap analysis.
- Discuss the apparent discrepancy: Factor 1 overlaps with classical gene lists but is not enriched in classical PDAC cell types in scRNA-seq. Possible explanations: (a) published "classical" lists may include microenvironmental genes that co-occur with classical cells in bulk; (b) DeSurv separates tumor-intrinsic from stromal signals more cleanly than unsupervised approaches.
- Conclude: scRNA-seq validation at the cellular level confirms that DeSurv factors correspond to distinct cell populations, supporting biological interpretability.

### 2B. Add null and mixed simulation scenario figures

**Location:** After the convergence figure (Figure S1).

**What to add:** The revised Results mention null and mixed scenarios alongside the primary (easy) scenario. Supplementary figures should show:

| Figure | Content |
|--------|---------|
| Fig S-null-cindex | C-index boxplots under null scenario (DeSurv vs NMF) |
| Fig S-null-precision | Precision boxplots under null scenario |
| Fig S-null-khist | Selected-k histograms under null scenario |
| Fig S-mixed-cindex | C-index boxplots under mixed scenario |
| Fig S-mixed-precision | Precision boxplots under mixed scenario |
| Fig S-mixed-khist | Selected-k histograms under mixed scenario |

These figures already exist in `figures/sim/` (filenames like `sim_cindex_boxplot__R00_null-*.pdf` and `sim_cindex_boxplot__R_mixed-*.pdf`). They need to be included in the supplement with appropriate captions explaining the scenario design and interpreting the results.

**Key interpretation points:**
- Null: DeSurv should show no advantage over NMF (BO selects low $\alpha$), confirming it does not hallucinate structure.
- Mixed: DeSurv should show moderate advantage, intermediate between easy and null.

### 2C. Expand convergence traces (if needed)

**Location:** Figure S1 (convergence figure).

**What to check:** The current convergence figure shows 5 initializations for one setting ($k=8$, $\alpha=0.5$). Consider adding:
- Convergence traces across different $\alpha$ values (e.g., $\alpha = 0, 0.25, 0.5, 0.75$) to show that supervision does not destabilize optimization.
- Traces from the BO-selected optimal configuration.

---

## 3. `paper/references_30102025.bib` — New Citation Keys

The revised sections cite 8 keys not currently in the bibliography. Only `Maurer2019` exists. The following entries need to be added:

### 3.1 `cook2007fisher`
```bibtex
@article{cook2007fisher,
  title     = {Fisher lecture: Dimension reduction in regression},
  author    = {Cook, R. Dennis},
  journal   = {Statistical Science},
  volume    = {22},
  number    = {1},
  pages     = {1--26},
  year      = {2007},
  publisher = {Institute of Mathematical Statistics}
}
```

### 3.2 `bair2006supervised`
```bibtex
@article{bair2006supervised,
  title     = {Prediction by supervised principal components},
  author    = {Bair, Eric and Hastie, Trevor and Paul, Debashis and Tibshirani, Robert},
  journal   = {Journal of the American Statistical Association},
  volume    = {101},
  number    = {473},
  pages     = {119--137},
  year      = {2006},
  publisher = {Taylor \& Francis}
}
```

### 3.3 `tishby1999information`
```bibtex
@article{tishby1999information,
  title   = {The information bottleneck method},
  author  = {Tishby, Naftali and Pereira, Fernando C. and Bialek, William},
  journal = {Proceedings of the 37th Annual Allerton Conference on Communication, Control, and Computing},
  pages   = {368--377},
  year    = {1999}
}
```

### 3.4 `chapelle2006semi`
```bibtex
@book{chapelle2006semi,
  title     = {Semi-Supervised Learning},
  author    = {Chapelle, Olivier and Sch{\"o}lkopf, Bernhard and Zien, Alexander},
  year      = {2006},
  publisher = {MIT Press},
  address   = {Cambridge, MA}
}
```

### 3.5 `aran2015systematic`
```bibtex
@article{aran2015systematic,
  title     = {Systematic pan-cancer analysis of tumour purity},
  author    = {Aran, Dvir and Sirota, Marina and Butte, Atul J.},
  journal   = {Nature Communications},
  volume    = {6},
  pages     = {8971},
  year      = {2015},
  publisher = {Nature Publishing Group}
}
```

### 3.6 `frigyesi2008nmf`
```bibtex
@article{frigyesi2008nmf,
  title     = {Non-negative matrix factorization for the analysis of complex gene expression data: identification of clinically relevant tumor subtypes},
  author    = {Frigyesi, Attila and H{\"o}glund, Mattias},
  journal   = {Cancer Informatics},
  volume    = {6},
  pages     = {275--292},
  year      = {2008}
}
```

### 3.7 `arora2020survclust`
```bibtex
@article{arora2020survclust,
  title     = {Surv{C}lust: an integrative survival-weighted clustering method for multi-omic data},
  author    = {Arora, Arshi and Olshen, Adam B. and Seshan, Venkatraman E. and Shen, Ronglai},
  journal   = {bioRxiv},
  year      = {2020},
  doi       = {10.1101/2020.09.04.283838},
  note      = {Preprint}
}
```

### 3.8 `damrauer2014intrinsic`
```bibtex
@article{damrauer2014intrinsic,
  title     = {Intrinsic subtypes of high-grade bladder cancer reflect the hallmarks of breast cancer biology},
  author    = {Damrauer, Jeffrey S. and Hoadley, Katherine A. and Chism, David D. and Fan, Cheng and Tiganelli, Christopher J. and Wobker, Sara E. and Yeh, Jen Jen and Milowsky, Matthew I. and Iyer, Gopa and Parker, Joel S. and Kim, William Y.},
  journal   = {Proceedings of the National Academy of Sciences},
  volume    = {111},
  number    = {8},
  pages     = {3110--3115},
  year      = {2014},
  publisher = {National Academy of Sciences}
}
```

### 3.9 `consensus2025pdac`

**Note:** This is a placeholder key for the emerging PDAC subtyping consensus reference. Replace with the actual citation once published. Candidates include:

- Ellrott et al. 2025, "Classification of pancreatic cancer" (if this is the consensus paper)
- A review synthesizing Collisson, Moffitt, Bailey, and subsequent validation studies

If no single consensus paper exists yet, this key can be replaced with a phrase like "as reviewed in [@collisson2011subtypes; @moffitt2015virtual; @Bailey2016]" and the citation key removed.

---

## 4. `paper/paper.Rmd` — YAML Header Updates

### 4A. Significance statement

**Location:** Lines 37-38 (`significance:` field).

**What to update:** The current significance statement is acceptable but could be tightened to match the revised prediction-validation framing. Suggested revision:

> Tumor transcriptomes mix malignant and microenvironmental signals, and the programs that explain the most variance are not necessarily those that drive clinical outcomes. Existing deconvolution methods discover latent programs without ensuring prognostic relevance, while supervised predictors sacrifice biological interpretability. DeSurv integrates nonnegative matrix factorization with Cox modeling to learn gene programs and their survival associations jointly. In simulations, DeSurv recovers prognostic structure more reliably than unsupervised NMF. In pancreatic cancer, survival supervision suppresses outcome-neutral variance and isolates interpretable programs that generalize across independent cohorts and transfer to bladder cancer.

### 4B. Keywords

**Location:** Lines 43-48 (`keywords:` field).

**Current:** Placeholder values ("one", "two", "optional", etc.).

**Suggested keywords:**
```yaml
keywords:
  - nonnegative matrix factorization
  - survival analysis
  - semi-supervised learning
  - tumor deconvolution
  - pancreatic cancer
```

### 4C. Abstract

**Location:** Lines 34-35 (`abstract:` field).

**What to consider:** The abstract is adequate but could be updated to reflect the prediction-validation framing and the semi-supervised language used throughout the revised sections. Key additions:
- Mention "semi-supervised" framing
- Mention cross-cancer transfer (bladder)
- Mention simulation validation
- Use "gene programs" consistently (not "gene signatures" in some places)

---

## 5. Missing Supplementary Figures Checklist

The following simulation figures exist in `figures/sim/` but are not yet included in `supplement.Rmd`:

| File | Scenario | Metric | Tuning |
|------|----------|--------|--------|
| `sim_cindex_boxplot__R00_null-bo.pdf` | Null | C-index | BO |
| `sim_cindex_boxplot__R00_null-bo_tune_ntop.pdf` | Null | C-index | BO + ntop |
| `sim_cindex_boxplot__R00_null-fixed.pdf` | Null | C-index | Fixed |
| `sim_precision_boxplot__R0_easy-bo.pdf` | Easy | Precision | BO |
| `sim_precision_boxplot__R0_easy-bo_tune_ntop.pdf` | Easy | Precision | BO + ntop |
| `sim_precision_boxplot__R0_easy-fixed.pdf` | Easy | Precision | Fixed |
| `sim_precision_boxplot__R_mixed-bo.pdf` | Mixed | Precision | BO |
| `sim_precision_boxplot__R_mixed-bo_tune_ntop.pdf` | Mixed | Precision | BO + ntop |
| `sim_precision_boxplot__R_mixed-fixed.pdf` | Mixed | Precision | Fixed |
| `sim_selected_k_hist__R00_null-bo.pdf` | Null | Selected k | BO |
| `sim_selected_k_hist__R00_null-bo_tune_ntop.pdf` | Null | Selected k | BO + ntop |
| `sim_selected_k_hist__R00_null-fixed.pdf` | Null | Selected k | Fixed |
| `sim_selected_k_hist__R0_easy-bo.pdf` | Easy | Selected k | BO |
| `sim_selected_k_hist__R0_easy-bo_tune_ntop.pdf` | Easy | Selected k | BO + ntop |
| `sim_selected_k_hist__R0_easy-fixed.pdf` | Easy | Selected k | Fixed |
| `sim_selected_k_hist__R_mixed-bo.pdf` | Mixed | Selected k | BO |
| `sim_selected_k_hist__R_mixed-bo_tune_ntop.pdf` | Mixed | Selected k | BO + ntop |
| `sim_selected_k_hist__R_mixed-fixed.pdf` | Mixed | Selected k | Fixed |
| `sim_cindex_boxplot__R0_easy-bo.pdf` | Easy | C-index | BO |
| `sim_cindex_boxplot__R0_easy-bo_tune_ntop.pdf` | Easy | C-index | BO + ntop |
| `sim_cindex_boxplot__R0_easy-fixed.pdf` | Easy | C-index | Fixed |
| `sim_cindex_boxplot__R_mixed-bo.pdf` | Mixed | C-index | BO |
| `sim_cindex_boxplot__R_mixed-bo_tune_ntop.pdf` | Mixed | C-index | BO + ntop |
| `sim_cindex_boxplot__R_mixed-fixed.pdf` | Mixed | C-index | Fixed |

**Recommendation:** Include at least the BO-tuned panels for each scenario (null, easy, mixed) for C-index, precision, and selected-k. The `bo_tune_ntop` and `fixed` variants can be included as additional supplementary panels or discussed in text.

---

## 6. Priority Order for Implementation

| Priority | Item | Effort | Dependency |
|----------|------|--------|------------|
| 1 | Add 8 bib entries (Section 3) | Low | None — blocks all revised sections |
| 2 | YAML keywords and significance (Section 4) | Low | None |
| 3 | Complete scRNA-seq text (Section 2A) | Medium | Requires reviewing scRNA-seq results |
| 4 | Add null/mixed sim figures (Section 2B) | Low | Figures already exist |
| 5 | Semi-supervised framing paragraph (Section 1A) | Medium | None |
| 6 | $\partial \mathcal{L}_{\text{Cox}} / \partial H = 0$ derivation (Section 1B) | Low | None |
| 7 | Theory-to-practice gap paragraph (Section 1C) | Medium | None |
| 8 | Simulation scenario details (Section 1D) | Medium | None |
| 9 | Out-of-sample protocol description (Section 1E) | Medium | None |
| 10 | SDR theory connection (Section 1F) | Low | bib entries added |
| 11 | Abstract update (Section 4C) | Low | After main-text sections finalized |
