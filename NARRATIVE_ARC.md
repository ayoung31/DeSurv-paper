# DeSurv Manuscript: Narrative Arc Analysis

**Document Purpose:** Detailed analysis of the paper's storytelling structure and how figures/tables support the narrative arc.

---

## Executive Summary

The DeSurv manuscript employs a **problem-solution-validation** narrative structure with escalating generalization claims. The story progresses from controlled validation (simulation) through single-cancer analysis (PDAC) to cross-cancer transfer (bladder), building credibility at each step. Six main-text figures are strategically positioned to provide evidence at each narrative beat.

---

## The Five-Act Story Structure

### Act 1: The Problem (Introduction)
**Narrative Question:** Why do current methods fail to identify clinically relevant gene programs?

### Act 2: The Solution (Methods + Fig 1)
**Narrative Question:** How does DeSurv address this problem?

### Act 3: The Proof of Concept (Simulation + Figs 2-3)
**Narrative Question:** Does the method work under controlled conditions?

### Act 4: The Biological Insight (PDAC + Figs 4-5)
**Narrative Question:** What does DeSurv reveal about real tumor biology?

### Act 5: The Generalization (Cross-cancer + Fig 6)
**Narrative Question:** Is this a general principle or a PDAC-specific finding?

---

## Detailed Narrative Beat Analysis

### Introduction: Setting the Stakes

The introduction follows a classical five-paragraph problem-solution structure:

| Paragraph | Narrative Function | Key Tension |
|-----------|-------------------|-------------|
| 1 | Establish importance | Molecular subtyping matters for precision oncology |
| 2 | Identify the gap | Discovery/validation separation risks overfitting |
| 3 | Explain why it's hard | Technical barriers (cohort size, reference limitations) |
| 4 | Critique prior art | Existing survival-NMF supervises H, not W |
| 5-6 | Present solution | DeSurv supervises W directly |

**The Core Innovation Claim (Paragraph 4-5):**
> "Existing survival-aware NMF methods integrate the survival objective through the sample-specific loadings (H) rather than the gene-level programs (W). DeSurv integrates survival information directly into the gene signature matrix."

This is the **conceptual hook** that the entire paper builds around.

---

### Results Section: The Evidence Cascade

The results section follows a strict logical progression:

```
┌─────────────────────────────────────────────────────────────────────┐
│  MODEL OVERVIEW (Fig 1)                                             │
│  "Here's what DeSurv is"                                            │
└──────────────────────────────┬──────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│  RANK SELECTION PROBLEM (Fig 2)                                     │
│  "Standard methods can't pick k; DeSurv can"                        │
└──────────────────────────────┬──────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│  SIGNATURE QUALITY (Fig 3)                                          │
│  "Supervision identifies correct genes"                             │
└──────────────────────────────┬──────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│  BIOLOGICAL STRUCTURE (Fig 4)                                       │
│  "DeSurv reorganizes programs toward prognosis"                     │
│  *** THE TURNING POINT - This is why you should care ***            │
└──────────────────────────────┬──────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│  EXTERNAL VALIDATION (Fig 5)                                        │
│  "Structure generalizes across PDAC cohorts"                        │
└──────────────────────────────┬──────────────────────────────────────┘
                               │
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│  CROSS-CANCER TRANSFER (Fig 6)                                      │
│  "Programs transfer to bladder cancer"                              │
│  *** THE PNAS-LEVEL CLAIM - General principle, not just PDAC ***    │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Figure-by-Figure Narrative Role

### Figure 1: The Method Schematic
**Narrative Position:** Opening
**Story Function:** "Here's what I'm proposing"

| Panel | Content | Narrative Role |
|-------|---------|----------------|
| A | NMF + Cox architecture | Shows W is shared between reconstruction and survival |
| B | Bayesian optimization | Shows principled hyperparameter selection |
| C-E | Model outputs | Shows interpretable outputs (gene programs, factors, coefficients) |

**What reader should take away:** DeSurv is a principled framework that couples NMF with Cox by sharing the gene program matrix W.

**Connection to story arc:** Sets up the technical foundation for all subsequent claims.

---

### Figure 2: The Rank Selection Problem
**Narrative Position:** Problem establishment + first solution
**Story Function:** "Standard methods fail; DeSurv succeeds"

| Panel | Content | Narrative Role |
|-------|---------|----------------|
| A | Reconstruction residuals | No clear elbow → ambiguous |
| B | Cophenetic correlation | Fluctuates → ambiguous |
| C | Silhouette width | Conflicts with A-B → ambiguous |
| D | DeSurv C-index heatmap | Clear surface → unambiguous selection |
| E | Simulation k recovery | Ground-truth validation |

**What reader should take away:** Unsupervised heuristics disagree; survival-guided selection resolves ambiguity.

**Connection to story arc:** Establishes that DeSurv can be trusted to select appropriate models before showing what those models reveal.

---

### Figure 3: Simulation Validation
**Narrative Position:** Controlled proof
**Story Function:** "Under known ground truth, DeSurv works better"

| Panel | Content | Narrative Role |
|-------|---------|----------------|
| A | C-index distributions | DeSurv predicts survival better |
| B | Precision distributions | DeSurv identifies true genes |

**What reader should take away:** Survival supervision improves both prediction AND program recovery.

**Connection to story arc:** Validates the method before applying to real (uncontrolled) data.

---

### Figure 4: PDAC Biological Reorganization (THE KEY FIGURE)
**Narrative Position:** Turning point
**Story Function:** "Here's the insight that makes this matter"

| Panel | Content | Narrative Role |
|-------|---------|----------------|
| A | NMF factor-program correlations | Shows variance-dominant structure (exocrine axis) |
| B | DeSurv factor-program correlations | Shows survival-aligned structure (immune-stromal) |
| C | Variance vs survival contribution | THE PUNCHLINE: variance ≠ prognosis |
| D | NMF-DeSurv correspondence | Shows reorganization, not replacement |

**What reader should take away:** DeSurv doesn't just predict better; it **restructures latent programs** to isolate survival-relevant biology from variance-dominant noise.

**Connection to story arc:** This is the "aha moment" that transforms the paper from a methods contribution to a biological insight.

---

### Figure 5: External PDAC Validation
**Narrative Position:** Reproducibility checkpoint
**Story Function:** "This isn't overfitting; it generalizes"

| Panel | Content | Narrative Role |
|-------|---------|----------------|
| A | Forest plot of HRs | DeSurv estimates more stable across cohorts |
| B | DeSurv KM curves | Clear survival separation |
| C | NMF KM curves | Weaker separation |

**What reader should take away:** The survival-aligned programs identified in training generalize to independent datasets.

**Connection to story arc:** Addresses the "overfitting" concern raised in the introduction by showing real external validation.

---

### Figure 6: Cross-Cancer Transfer (THE PNAS CLAIM)
**Narrative Position:** Culmination
**Story Function:** "This is a general principle, not just PDAC"

| Panel | Content | Narrative Role |
|-------|---------|----------------|
| A | Bladder variance vs survival | Same pattern as PDAC (variance ≠ prognosis) |
| B | PDAC→Bladder transfer KM | Programs learned in PDAC work in bladder |

**What reader should take away:** Survival-relevant transcriptional structure is shared across cancer types; DeSurv captures it.

**Connection to story arc:** Elevates the paper from "a better method for PDAC" to "a general framework for cancer transcriptomics."

---

## The Information Disclosure Strategy

The paper employs **hierarchical information disclosure**:

| Document | Role | Audience |
|----------|------|----------|
| Main Methods | Framework + key choices | Everyone |
| Results | Empirical validation | Everyone |
| Supplement | Complete derivations | Methods experts |
| Discussion | Synthesis + implications | Everyone |

This allows the main text to remain accessible while providing full rigor in the supplement.

---

## Rhetorical Strategies

### 1. Problem-Solution Framing
Every results section begins by identifying a limitation of standard NMF, then shows DeSurv addresses it:
- "Standard heuristics yield inconsistent guidance... To resolve this ambiguity, we applied DeSurv"
- "In the standard NMF solution, factors largely reflected dominant sources of variance... By contrast, DeSurv produced..."

### 2. Escalating Generalization Claims
The evidence builds in scope:
```
Simulation (controlled) → PDAC training → PDAC validation → Bladder (different cancer)
```

### 3. Consistent Baseline Comparison
DeSurv vs α=0 (standard NMF) runs throughout, providing a consistent ablation.

### 4. Biological Grounding
Every claim is tied to known biology (classical/basal-like, CAF subtypes, immune infiltration) to establish plausibility.

---

## The Central Insight: "Variance ≠ Prognosis"

The entire paper builds toward one core insight, visualized in **Figure 4C** and echoed in **Figure 6A**:

> **High-variance transcriptional axes are not necessarily prognostically relevant.** Standard NMF allocates model capacity to variance-dominant signals (like exocrine expression in PDAC) that contribute little to survival. DeSurv redirects this capacity toward survival-aligned programs.

This insight is what makes the paper PNAS-worthy rather than just "slightly better C-index."

---

## Narrative Gaps and Opportunities

### Currently Underemphasized
1. **The W vs H distinction** - buried in Introduction paragraph 4
2. **The projection advantage** - closed-form Z = W^⊤X is mentioned but not highlighted
3. **Why cross-cancer transfer works** - biological mechanism is implicit

### Currently Missing
1. **Explicit out-of-sample statement** - nested CV is used but not clearly stated
2. **Factor biological labels** - Figure 4 describes but doesn't label
3. **Effect sizes on KM plots** - Forest plot has HRs but KM plots lack p-values

---

## Summary: The Story DeSurv Tells

**Opening Hook:** Current molecular subtyping separates discovery from validation, missing clinically relevant programs.

**Innovation:** DeSurv supervises through gene programs (W), not sample loadings (H), ensuring discovered programs are intrinsically prognostic.

**Evidence Chain:**
1. Standard rank selection fails; DeSurv's BO-guided selection works (Fig 2)
2. Under controlled conditions, supervision improves gene recovery (Fig 3)
3. In PDAC, supervision reorganizes programs toward prognosis (Fig 4)
4. This structure generalizes across PDAC cohorts (Fig 5)
5. This structure even transfers to bladder cancer (Fig 6)

**Closing:** DeSurv provides a general framework for deriving actionable, survival-aligned insights from heterogeneous tumor transcriptomes.

---

## Appendix: Figure-Narrative Mapping Table

| Figure | Section | Narrative Beat | Key Claim | Evidence Type |
|--------|---------|----------------|-----------|---------------|
| 1 | Model Overview | Method introduction | W is shared, BO for selection | Schematic |
| 2 | Rank Selection | Problem + solution | Standard fails, DeSurv works | Real data + simulation |
| 3 | Simulation | Controlled validation | Supervision improves recovery | Simulation |
| 4 | Biological Structure | **Turning point** | Variance ≠ prognosis | Real data comparison |
| 5 | External Validation | Reproducibility | Generalizes within-cancer | Independent cohorts |
| 6 | Cross-cancer | **PNAS claim** | Generalizes across-cancer | Transfer experiment |

---

*Analysis generated: 2026-01-24*
