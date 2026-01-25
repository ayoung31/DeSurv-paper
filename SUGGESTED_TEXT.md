# DeSurv Manuscript: Suggested Replacement Text

**Purpose:** Drop-in replacement language implementing the editorial recommendations from PNAS_REVIEW.md.

**Alignment:** This text aligns with the recalibrated PNAS assessment (Appendix D) by:
- Acknowledging prior work in prose (not requiring new benchmarks)
- Tempering cross-cancer claims appropriately
- Reframing novelty as methodological, not biological
- Providing the W-vs-H distinction that reviewers will expect

---

## A. Replacement Introduction Language

*(Replace current paragraphs discussing prior survival-NMF work, discovery/validation framing, and novelty)*

### Revised Introduction (Drop-in Replacement)

> Nonnegative matrix factorization (NMF) is widely used in cancer genomics because its nonnegativity constraints yield additive, interpretable molecular programs that often correspond to cell types or transcriptional states. However, standard NMF is variance-driven: it prioritizes dominant sources of expression variation, which may reflect tissue composition or differentiation rather than clinically relevant biology. As a result, factors learned by unsupervised NMF often require extensive downstream filtering and retrospective evaluation to assess their prognostic or therapeutic relevance.
>
> Several prior efforts have sought to incorporate survival outcomes into matrix factorization frameworks. For example, proportional-hazards–based NMF formulations have been proposed in which the Cox partial likelihood is combined with a low-rank reconstruction objective, enabling joint optimization of latent factors and survival associations (e.g., Huang et al., 2020; Le Goff et al., 2025). In these approaches, survival is typically linked to sample-level factor scores derived from the decomposition, while gene-level programs are updated indirectly through their coupling to reconstruction error. Although these methods demonstrate that outcome-aware factorization is feasible, their practical behavior, reproducibility, and biological interpretability remain incompletely characterized, and formal strategies for model selection and validation are often limited.
>
> Importantly, outcome-guided discovery introduces a fundamental tradeoff. Unsupervised discovery followed by independent validation is designed to reduce outcome-specific overfitting, but it can preferentially capture cohort-specific variance that is weakly related to clinical endpoints. Conversely, jointly incorporating outcomes during discovery can better align learned structure with prognosis, at the cost of increased risk of outcome-driven overfitting. Addressing this tradeoff requires careful regularization, cross-validation, and external validation across independent cohorts.
>
> Here we introduce **DeSurv**, a survival-supervised deconvolution framework that integrates nonnegative matrix factorization with Cox proportional hazards modeling in a way that directly couples survival gradients to gene-level programs. In DeSurv, survival supervision enters through factor scores defined as linear projections of expression onto the gene program matrix, (Z = W^⊤X), so that the Cox partial likelihood depends explicitly on the gene programs (W) and their associated coefficients. Sample-level loadings are updated solely through the reconstruction objective, while gene programs are shaped by both reconstruction fidelity and survival association. This design yields gene-level programs that are intrinsically aligned with clinical outcomes and are portable across datasets via closed-form projection.
>
> To mitigate overfitting and ambiguity in rank selection, DeSurv employs Bayesian optimization to tune model complexity and supervision strength using cross-validated concordance, followed by validation in independent cohorts. We apply this framework to pancreatic ductal adenocarcinoma and additional cancer types, demonstrating that survival-informed factorization concentrates prognostic signal into a small number of interpretable programs while suppressing variance-dominant but outcome-neutral structure.

---

## B. New "Related Work" Paragraph

*(Place at end of Introduction or beginning of Discussion)*

### Related Work and Context

> DeSurv builds on a substantial body of prior work on transcriptomic deconvolution and molecular subtyping in pancreatic cancer. Virtual microdissection approaches first demonstrated that bulk PDAC expression profiles can be decomposed into tumor- and stroma-specific compartments with distinct prognostic associations, including basal-like and classical tumor programs and activated versus normal stroma (Moffitt et al., 2015). Subsequent NMF-based frameworks, including DECODER, formalized unsupervised deconvolution of tumor and microenvironmental compartments and showed that derived compartment ratios, such as basal-to-classical expression, are associated with survival across cohorts (Peng et al., 2019).
>
> In contrast to these approaches, which discover latent structure without outcome supervision and assess prognostic relevance post hoc, DeSurv explicitly incorporates survival information during factor discovery. The goal is not to redefine established PDAC biology, but to evaluate whether outcome-guided learning can improve the stability, portability, and prognostic alignment of gene programs relative to variance-driven factorization. Accordingly, the biological programs recovered by DeSurv overlap with known tumor and microenvironmental states, while the methodological contribution lies in how these programs are learned, selected, and validated.

**Note:** This paragraph defuses claims of rediscovery, acknowledges UNC/DECODER work explicitly, and correctly frames novelty as methodological.

---

## C. Targeted Fixes to Specific Claims

### Fix 1: Huang et al. (CoxNMF) claim

**Current:**
> "both integrate the survival objective through the sample-specific loadings rather than the gene-level programs"

**Replace with:**
> "prior survival-aware NMF formulations typically link survival to sample-level factor scores, with gene programs updated indirectly through their coupling to reconstruction error. In contrast, DeSurv defines the survival objective directly in terms of gene-program projections, so that Cox gradients act explicitly on the gene program matrix."

**Why:** Mathematically accurate and defensible.

---

### Fix 2: Discovery-validation overfitting claim

**Current:**
> "Separating discovery from validation can risk overfitting and limits biological and clinical generalizability."

**Replace with:**
> "Unsupervised discovery followed by independent validation reduces outcome-specific overfitting but may prioritize cohort-specific variance that is weakly related to clinical endpoints. Outcome-guided discovery offers an alternative tradeoff, potentially improving prognostic alignment while requiring careful regularization and external validation to control overfitting."

**Why:** More nuanced; doesn't dismiss standard practice.

---

### Fix 3: Cross-cancer transfer novelty claim

**Current:**
> "demonstrating cross-cancer transfer of prognostic programs"

**Replace with:**
> "demonstrating that survival-aligned programs learned in PDAC recapitulate known basal-like transcriptional structure that generalizes across epithelial cancers, consistent with prior reports."

**Why:** Tempers novelty claim; acknowledges prior knowledge.

---

### Fix 4: W vs H supervision (add to Methods or SI)

**Add:**
> "Because the Cox partial likelihood depends on W^⊤X rather than on H, the gradient of the survival objective with respect to H is identically zero in DeSurv, whereas the gradient with respect to W is nonzero, yielding direct outcome supervision of gene programs."

**Why:** Preempts mathematical reviewers who will check this.

---

### Fix 5: Out-of-sample statement (add to Methods or Results)

**Add:**
> "Model selection, including rank (k) and supervision strength (α), was performed entirely within training data using nested cross-validation. Validation cohorts were held out from all tuning procedures. Reported hazard ratios and survival associations in external datasets reflect purely out-of-sample generalization."

**Why:** Pre-empts "double-dipping" concerns.

---

### Fix 6: Structural reorganization framing (add to Results or Discussion)

**Add:**
> "The key gain is not raw discrimination but the reorganization of transcriptional axes such that survival-relevant biology is isolated from variance-dominant but prognostically neutral signals. DeSurv is not optimized primarily for maximal prediction, but for re-allocation of latent structure toward outcome relevance."

**Why:** Reframes away from C-index comparison toward conceptual contribution.

---

## D. Significance Statement (Revised)

If the current significance statement is too abstract, consider:

> Molecular subtyping of tumors typically discovers gene programs through unsupervised methods and evaluates their clinical relevance retrospectively. We introduce DeSurv, a framework that integrates survival outcomes directly into gene program discovery, yielding programs that are intrinsically aligned with prognosis. Applied to pancreatic cancer, DeSurv concentrates survival signal into interpretable programs while suppressing variance-dominant but outcome-neutral structure. The resulting gene signatures are portable across independent cohorts and cancer types, demonstrating that outcome-guided learning can improve the stability and clinical utility of transcriptomic deconvolution.

---

## Usage Notes

1. **Introduction replacement:** Can be used as-is or adapted to match author voice
2. **Related Work paragraph:** Critical for addressing prior art concerns; place strategically
3. **Claim fixes:** Apply individually to specific sentences in the manuscript
4. **W vs H gradient statement:** Add to Methods supplement if space-constrained

---

*Generated: 2026-01-24*
*Complements: PNAS_REVIEW.md, FIGURE_ANALYSIS.md*
