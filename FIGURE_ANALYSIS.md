# DeSurv Figure Analysis: Panel-by-Panel Evaluation

**Analysis Date:** 2026-01-24
**Purpose:** Evaluate whether each panel optimally supports its narrative point, propose improvements, and assess overall figure structure.

---

## Executive Summary

The figures effectively tell the DeSurv story but have several areas for improvement:

1. **Figure 1** (Schematic): Strong conceptual communication, minor notation improvements needed
2. **Figure 2** (Rank Selection): Layout issues - key message (Panel A) is excellent but relegated to bottom; NMF diagnostics dominate
3. **Figure 3** (Simulation): Good clarity but missing key context (true k, simulation parameters)
4. **Figure 4** (Biological Structure): THE KEY FIGURE - strong content but missing factor labels and variance vs survival scatter plot not shown
5. **Figure 5** (External Validation): Forest plot excellent; KM curves need effect sizes
6. **Figure 6** (Bladder): Variance vs survival plot strong; transfer KM needs annotation

**Overall Narrative Arc Impact:** Improvements would strengthen the "variance ≠ prognosis" central insight and make the evidence cascade more compelling.

---

## Reviewer Archetype Reference

When evaluating recommendations, consider which reviewer types each change helps or might annoy:

| Archetype | Focus | Priority |
|-----------|-------|----------|
| **R1: Statistical purist** | Calibration, leakage, rigor | Effect sizes, double-dipping prevention |
| **R2: Methods maximalist** | Comparisons, benchmarks | Simulation breadth, baselines |
| **R3: Biology-first reviewer** | Pathway interpretation, clinical relevance | Factor labels, biological plausibility |
| **R4: PNAS generalist / AE** | Narrative clarity, broad impact | Visual flow, conceptual hooks |
| **R5: Skeptical minimalist** | "Less is more," parsimony | Clean figures, no scope creep |

---

## Backfire Risk Assessment

Not all recommendations are risk-free. This section stress-tests each major proposed change.

### Summary: Backfire Risk Table

| Change | Backfire Risk | Who It Helps | Who Might Push Back | Verdict |
|--------|---------------|--------------|---------------------|---------|
| HR / p-values on KM | **None** | R1, R3, R4 | Nobody | **Must do** |
| Factor labels | **Low** (if labels overclaim) | R3, R4 | R1 (if speculative) | **Must do (label conservatively)** |
| Null/harder simulation | **Low** (if framed poorly) | R1, R2 | R5 (paper length) | **Strongly recommended** |
| Fig 2 restructuring | **None** | R4, R5 | Nobody | **Must do** |
| "variance ≠ prognosis" annotation | **Very low** | R4, R3 | R1 (if slogan-like) | **Do carefully** |
| Color consistency | **None** (time cost only) | Everyone | Nobody | Optional |
| I²/meta-analysis | **Moderate** | R1 | R4, R5 | **Skip** |

### Detailed Risk Analysis

**A. Adding HRs / p-values to KM plots (Figs 5–6)**
- **Risk:** Essentially none. This is standard survival reporting.
- **Edge case:** If HRs are modest (e.g., ~1.2 with wide CI), frame as "evidence of generalization," not clinical readiness.
- **Verdict:** Do it. No downside.

**B. Adding factor labels to Figure 4 heatmaps**
- **Risk:** Over-interpretive labels (e.g., "Immunosuppressive CAF-like state") invite "prove it" demands.
- **Mitigation:** Use descriptive, conservative labels:
  - ✅ "Classical tumor–associated"
  - ✅ "Immune/stromal–enriched"
  - ❌ "Immunosuppressive CAF-mediated invasion program"
- **Verdict:** Do it, but keep labels cautious.

**C. Adding "variance ≠ prognosis" annotation (Fig 4C)**
- **Risk:** A pedantic statistician might bristle at oversimplification.
- **Mitigation:** Phrase precisely: "Variance explained is not necessarily aligned with survival relevance"
- **Verdict:** Do it, but phrase carefully.

**D. I² or pooled meta-analytic diamonds (Fig 5)**
- **Risk:** Raises expectations of full meta-analysis rigor. May trigger random vs fixed effects questions.
- **Verdict:** Skip unless trivial. Not needed.

---

## Figure 1: Model Schematic (`model_schematic_final.pdf`)

### Current State
- **Panel A:** NMF + Cox architecture showing W is shared
- **Panel B:** Bayesian optimization cycle with hyperparameter sliders

### Panel-by-Panel Evaluation

#### Panel A: The DeSurv Model
| Aspect | Current | Assessment |
|--------|---------|------------|
| W-sharing emphasis | Red "W is shared" callout | ✅ Strong - this IS the innovation |
| Data flow | X → W,H decomposition → Z = X^T W → Cox | ✅ Clear left-to-right |
| α interpretation | Endpoints labeled (α=0 NMF, α→1 Cox) | ✅ Intuitive |
| Matrix dimensions | All labeled (p×n, p×k, k×n, k×1) | ✅ Complete |

**What's Working:** The W-sharing callout is prominent and well-positioned. The visual flow clearly shows how NMF and Cox share the gene program matrix.

**Proposed Improvements:**
1. **Add micro-callout:** Near W, add text: "Supervision acts on **gene programs (W)**, not sample loadings (H) → enables projection-based transfer"
2. **Notation consistency:** Ensure Z = W^⊤X matches paper text (currently shows Z = (X^T W) which is transposed)
3. **Add interpretability hook:** Small annotation showing "W columns = interpretable gene signatures"

#### Panel B: Bayesian Optimization
| Aspect | Current | Assessment |
|--------|---------|------------|
| Cycle diagram | Propose → Evaluate → Update → repeat | ✅ Clear |
| Hyperparameters | k, α, λ with sliders | ✅ Intuitive |
| Objective | "cross-validated C-index" | ⚠️ Should emphasize "survival-guided" |

**What's Working:** The BO cycle is intuitive; slider metaphor for hyperparameters is accessible.

**Proposed Improvements:**
1. **Label the objective explicitly:** Add "CV C-index objective" in bold near evaluation box
2. **Add "out-of-sample" annotation:** Near evaluation box, add "(held-out folds)" to pre-empt double-dipping concerns
3. **Consider:** Adding 1-standard-error rule annotation to output

### Overall Figure 1 Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Clarity** | 9/10 | Excellent conceptual communication |
| **Aesthetic balance** | 8/10 | A is denser than B; could equalize |
| **Narrative support** | 8/10 | Strong but W-level innovation could be more prominent |
| **PNAS fit** | 9/10 | Clean, professional, appropriate complexity |

**Structural Recommendation:** Consider making Panel A slightly larger (it's the conceptual core) and adding the micro-callouts above. The current 50/50 split slightly underweights the model architecture.

---

## Figure 2: Rank Selection (`fig_bo_tcgacptac.pdf`)

### Current Layout
```
Row 1: [Panel A - CV C-index curves (DeSurv vs α=0)]
Row 2: [Panel B - Cophenetic] [Panel C - Residuals] [Panel D - Silhouette]
```

### Panel-by-Panel Evaluation

#### Panel A: CV C-index by Rank (DeSurv vs NMF)
| Aspect | Current | Assessment |
|--------|---------|------------|
| Comparison | Two lines with error bars | ✅ Clear method comparison |
| Message | DeSurv (green) higher and flatter | ✅ Supports claim |
| Legend | "alpha = 0" vs "DeSurv" | ⚠️ Inconsistent terminology |
| Y-axis range | 0.52-0.65 | ⚠️ Starts at 0.52, could exaggerate differences |

**What's Working:** The comparison is stark - DeSurv achieves higher C-index at all k values and is more stable.

**Issues:**
1. **Legend inconsistency:** "alpha = 0" is technical; should be "Standard NMF" to match text
2. **Y-axis starting point:** Starting at 0.52 makes the gap look larger than it is (0.05-0.10 difference)
3. **Missing:** Selected k annotation - which k was actually chosen?

**Proposed Improvements:**
1. Change legend to "Standard NMF (α=0)" vs "DeSurv (α>0)"
2. Annotate the selected k* with vertical line or marker
3. Add standard error bands (ribbons) instead of error bars for better visual
4. Consider starting y-axis at 0.5 for honesty (or add axis break indicator)

#### Panels B-D: Unsupervised Heuristics
| Panel | Metric | Message | Assessment |
|-------|--------|---------|------------|
| B | Cophenetic correlation | Declines, fluctuates at high k | ✅ Shows ambiguity |
| C | Residuals | Smooth decrease, no elbow | ✅ Shows ambiguity |
| D | Silhouette (multiple metrics) | Monotonic decrease | ✅ Shows ambiguity |

**What's Working:** Each panel shows a different heuristic giving different guidance, effectively demonstrating the "ambiguous rank selection" problem.

**Issues:**
1. **No selected k indicated:** Reader can't see what each method would choose
2. **X-axis labels:** "rank" (lowercase) vs "Rank k" in Panel A - inconsistent
3. **Y-axis labels:** Scientific notation in Panel C (1.00e+10) is hard to parse
4. **Panel D colors:** Three lines (green, red, purple) but no legend visible in PDF

**Proposed Improvements:**
1. Add vertical dashed lines showing "would select k=X" for each heuristic
2. Standardize axis labels across all panels
3. Use relative residuals (% of max) instead of absolute values
4. Add legend to Panel D or remove multi-metric display

### Layout Issue: Narrative Hierarchy

**Problem:** The most important message (Panel A: DeSurv outperforms NMF) is given equal weight with the "problem demonstration" panels (B-D). But the narrative is:
1. FIRST: Show the problem (NMF heuristics disagree)
2. THEN: Show the solution (DeSurv provides clear guidance)

**Current layout inverts this** - solution on top, problem below.

**Proposed Restructure:**
```
Option 1: Reorder
Row 1: [B - Cophenetic] [C - Residuals] [D - Silhouette]
Row 2: [A - CV C-index comparison, spanning full width]

Option 2: Size hierarchy
Row 1: [A - larger, ~60% height]
Row 2: [B] [C] [D] - smaller, ~40% height
```

Option 1 better matches the narrative flow (problem → solution).

### Overall Figure 2 Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Clarity** | 7/10 | Individual panels clear, but flow confusing |
| **Aesthetic balance** | 6/10 | A is cramped into same space as B-D |
| **Narrative support** | 6/10 | Layout contradicts problem→solution narrative |
| **PNAS fit** | 7/10 | Acceptable but could be tighter |

**Key Recommendation:** Restructure to show the problem (B-D) BEFORE the solution (A), either through reordering or size hierarchy.

---

## Figure 3: Simulation Results

### Current Layout
```
[Panel A - C-index boxplots] [Panel B - Precision boxplots]
```

### Panel-by-Panel Evaluation

#### Panel A: Test C-index Distribution
| Aspect | Current | Assessment |
|--------|---------|------------|
| Comparison | Boxplots with jittered points | ✅ Shows distribution shape |
| Message | DeSurv (blue) higher than α=0 (red) | ✅ Clear separation |
| Strip header | "R0_easy" | ❌ Technical jargon exposed |
| Y-axis | 0.6-0.9, labeled "C-index" | ✅ Appropriate |

**What's Working:** The distribution separation is stark and convincing. Jittered points show n and spread.

**Issues:**
1. **"R0_easy" label exposed:** This is internal simulation regime naming - meaningless to reader
2. **No true k indicated:** Reader doesn't know ground truth
3. **Sample size unclear:** How many replicates?
4. **Title redundant:** "Test C-index" in title AND y-axis label

**Proposed Improvements:**
1. Remove or replace "R0_easy" with descriptive text: "Simulation: k=3 true factors, clear signal"
2. Add in-panel annotation: "n=100 replicates, true k=3"
3. Add statistical comparison: p-value or effect size annotation
4. Remove title (redundant with y-axis)

#### Panel B: Precision Distribution
| Aspect | Current | Assessment |
|--------|---------|------------|
| Comparison | Boxplots with jittered points | ✅ Shows distribution shape |
| Message | DeSurv ~0.45, α=0 ~0.07 | ✅ DRAMATIC difference |
| Strip header | "R0_easy" | ❌ Same issue |
| Y-axis | "Precision" | ⚠️ Needs definition |

**What's Working:** This is the STRONGEST evidence in the simulation - near-zero precision for unsupervised NMF vs ~45% for DeSurv. This directly supports the "supervision improves gene program recovery" claim.

**Issues:**
1. **Precision undefined:** What's the denominator? (Answer from code: precision among top n_top genes per factor)
2. **Near-zero baseline not emphasized:** The fact that α=0 gets ~7% precision is devastating to unsupervised methods but not called out
3. **Title too long:** "Best precision (mean across lethal factors)" - truncate

**Proposed Improvements:**
1. Add caption text explaining precision metric
2. Add annotation: "NMF: ~7% precision (near-random)"
3. Consider adding recall or F1 as secondary metric
4. Simplify title to "Gene Recovery Precision"

### Missing Panel: k Selection Accuracy

The k-selection histogram (`sim_selected_k_hist__R0_easy-bo_tune_ntop.pdf`) is shown in Figure 2 as Panel E, but it could be more effectively integrated here with the other simulation results.

**The histogram shows:**
- DeSurv: Mode at k=3 (true value), reasonable spread
- α=0: Strong mode at k=2 (under-selection)

This is important evidence that belongs with other simulation validation.

### Overall Figure 3 Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Clarity** | 7/10 | Boxplots clear, but missing context |
| **Aesthetic balance** | 8/10 | Two equal panels works well |
| **Narrative support** | 6/10 | Missing "ground truth k=3" context |
| **PNAS fit** | 7/10 | Needs polish on labels |

**Key Recommendation:** Add in-panel annotations for simulation parameters (k=3, n replicates) and define precision metric in caption.

---

## Figure 4: Biological Structure Reorganization

**Note:** The PDF I viewed (`fig_bio_tcgacptac.pdf`) appears to be an older version with GO enrichment dotplots, not the heatmap comparison described in the paper. The paper describes panels A-D as:
- A: NMF gene-program correlations
- B: DeSurv gene-program correlations
- C: Variance vs survival scatter
- D: NMF-DeSurv factor correspondence

I'll analyze based on the paper description and the full panels I viewed.

### Panel A-B: Gene-Program Correlation Heatmaps

Based on panel D which shows the full heatmap:

| Aspect | Current | Assessment |
|--------|---------|------------|
| Rows | Reference gene signatures (Puleo, Moffitt, SCISSORS, etc.) | ✅ Comprehensive |
| Columns | Factors (F1, F2, F3...) | ⚠️ No interpretive labels |
| Color scale | Blue-white-red (-0.4 to 0.4) | ✅ Appropriate |
| Clustering | Rows clustered, columns not | ✅ Groups related signatures |

**What's Working:** The heatmap clearly shows which factors align with which biological programs. The blue-red diverging palette is appropriate for correlations.

**Issues:**
1. **Factor columns unlabeled:** Reader must infer biological meaning from row correlations
2. **No significance indicators:** Asterisks mentioned in caption but hard to see
3. **Threshold filtering (>0.2):** May look like selective reporting

**Proposed Improvements:**
1. **ADD FACTOR LABELS:** Below factor numbers, add interpretive labels:
   - F1: "Classical tumor" (if positively correlated with Classical signatures)
   - F2: "Basal-like/aggressive"
   - F3: "Activated TME" (if correlated with immune/stromal)
2. Add statistical significance markers (*, **, ***)
3. Provide full matrix in SI to address threshold concern
4. Consider side-by-side layout for NMF vs DeSurv instead of stacked

### Panel C: Variance vs Survival Contribution (THE KEY PANEL)

This is the "variance ≠ prognosis" visualization - the paper's core insight.

**Based on code analysis:** This should be a scatter plot with:
- X-axis: Fraction of variance explained by factor
- Y-axis: Survival contribution (Δ partial log-likelihood from univariate Cox)
- Points: Individual factors, colored by method (NMF vs DeSurv)

| Aspect | Expected | Assessment |
|--------|----------|------------|
| Message | High-variance factors ≠ high survival | ✅ Core insight |
| Comparison | NMF factors (high var, low surv) vs DeSurv (reorganized) | ✅ Critical |
| Metric clarity | "Change in partial log-likelihood" | ⚠️ Technical |

**Proposed Improvements:**
1. **Add annotations:** Label specific factors (e.g., "Exocrine-associated" for NMF high-variance factor)
2. **Add regression lines:** Show slope difference between methods
3. **Clarify metric:** In legend, explain "Δ log-lik from univariate Cox = survival contribution"
4. **Emphasize punchline:** Add text box: "NMF: variance ≠ prognosis; DeSurv: realigned"

### Panel D: NMF-DeSurv Correspondence

| Aspect | Current | Assessment |
|--------|---------|------------|
| Purpose | Show reorganization, not replacement | ✅ Important |
| Format | Correlation heatmap (NMF cols vs DeSurv rows) | ✅ Clear |
| Message | Tumor structure preserved, microenvironment refined | ✅ Supports claim |

**What's Working:** Shows DeSurv doesn't invent new structure - it reorganizes existing structure toward survival relevance.

**Proposed Improvement:**
- Add annotations for the key correspondences (e.g., "Classical tumor axis preserved")

### Overall Figure 4 Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Clarity** | 6/10 | Missing factor labels is a major gap |
| **Aesthetic balance** | 7/10 | Four-panel layout reasonable |
| **Narrative support** | 7/10 | Core insight present but underemphasized |
| **PNAS fit** | 6/10 | Needs interpretive annotations |

**Key Recommendation:** ADD FACTOR LABELS. This single change would dramatically improve interpretability and reduce reader burden. Also add "variance ≠ prognosis" annotation to Panel C.

---

## Figure 5: External Validation

### Current Layout (from paper code)
```
[Panel A - Forest plot] [Panel B - DeSurv KM] [Panel C - NMF KM]
```

### Panel A: Forest Plot of Hazard Ratios

| Aspect | Current | Assessment |
|--------|---------|------------|
| Cohorts | Puleo, Dijk, PACA seq/array, Moffitt | ✅ 5 independent cohorts |
| Comparison | DeSurv (blue) vs NMF (red) | ✅ Clear color coding |
| Effect display | Point estimates + 95% CI | ✅ Standard |
| Y-axis labels | Abbreviated cohort names | ✅ Space-efficient |

**What's Working:** Forest plot is the gold standard for displaying effect heterogeneity across studies. Color-coded comparison is intuitive.

**Issues:**
1. **Factor selection criterion not visible:** How was the "best" factor chosen? (Answer: largest Δ log-lik in training)
2. **No pooled estimate:** Missing summary diamond for meta-analytic effect
3. **X-axis reference line:** Is HR=1 (null) clearly marked?

**Proposed Improvements:**
1. Add subtitle: "Factor selected by max Δ log-lik in training data"
2. Add pooled HR with 95% CI at bottom (if defensible)
3. Ensure HR=1 reference line is prominent
4. Add I² heterogeneity statistic (optional but PNAS-friendly)

### Panels B-C: Kaplan-Meier Curves

| Aspect | Panel B (DeSurv) | Panel C (NMF) | Assessment |
|--------|------------------|---------------|------------|
| Stratification | Median split | Median split | ✅ Consistent |
| Risk table | Present | Present | ✅ Good |
| Legend | High/Low | High/Low | ✅ Clear |
| Colors | violetred2/turquoise4 | Same | ⚠️ Hard to distinguish |

**What's Working:** Risk tables are present (important for survival plots). Median split is simple and reproducible.

**Issues:**
1. **NO P-VALUES OR HRs:** This is the biggest gap - no effect size annotation on the plots
2. **Colors similar:** violetred2 and turquoise4 may be hard to distinguish for colorblind readers
3. **Pooled vs stratified unclear:** Are these pooled across cohorts or stratified?
4. **Time axis units:** "Time (months)" - verify units are correct

**Proposed Improvements:**
1. **ADD EFFECT SIZES:** On each plot, add:
   - Log-rank p-value
   - HR (95% CI)
   - Median survival per group
2. Use more divergent colors (e.g., blue vs orange)
3. Add note: "Pooled validation cohorts" or "Stratified by cohort" in caption
4. Consider adding median survival lines (horizontal dashed)

### Overall Figure 5 Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Clarity** | 7/10 | Missing effect sizes on KM plots |
| **Aesthetic balance** | 8/10 | Three-panel layout works |
| **Narrative support** | 7/10 | Good but incomplete statistical reporting |
| **PNAS fit** | 6/10 | PNAS expects complete statistical annotation |

**Key Recommendation:** ADD p-values and HRs to KM plots. This is standard practice and their absence may invite reviewer criticism.

---

## Figure 6: Cross-Cancer Transfer (Bladder)

### Current Layout
```
[Panel A - Variance vs survival (bladder)] [Panel B - Transfer KM]
```

### Panel A: Variance vs Survival in Bladder Cancer

This mirrors Figure 4C but for bladder cancer - showing the same "variance ≠ prognosis" pattern holds in a different cancer type.

| Aspect | Current | Assessment |
|--------|---------|------------|
| Message | Same pattern as PDAC | ✅ Supports generalization |
| Factors | DeSurv F1, F2 + NMF factors | ✅ Method comparison |
| Axes | Variance explained vs Δ log-lik | ⚠️ Same clarity issues as Fig 4 |

**What's Working:** Demonstrates the central insight is not PDAC-specific.

**Proposed Improvements:**
1. Add factor labels/annotations
2. Explicitly connect to Fig 4C: "Same pattern observed in bladder cancer"
3. Consider overlaying PDAC and bladder points on same plot (in SI)

### Panel B: Cross-Cancer Transfer KM

This is the "PNAS claim" - PDAC-trained factors work in bladder cancer.

| Aspect | Current | Assessment |
|--------|---------|------------|
| Transfer direction | PDAC → Bladder | ✅ Clear in caption |
| Stratification | Median split on PDAC factor scores | ✅ |
| Risk table | Present | ✅ |

**Issues:**
1. **NO HR or p-value:** Same issue as Figure 5
2. **Which factor transferred?** Not specified on plot
3. **Biological interpretation missing:** What does this factor represent?

**Proposed Improvements:**
1. **ADD HR and p-value** directly on plot
2. **Label the transferred factor:** "PDAC Factor X (Classical/Activated TME/etc.)"
3. **Add biological interpretation:** In caption, explain what the transferred factor represents and why transfer is biologically plausible

### Overall Figure 6 Assessment

| Criterion | Rating | Notes |
|-----------|--------|-------|
| **Clarity** | 6/10 | Missing statistical annotations |
| **Aesthetic balance** | 8/10 | Two-panel layout appropriate |
| **Narrative support** | 7/10 | Capstone claim needs more emphasis |
| **PNAS fit** | 6/10 | Biological interpretation needed |

**Key Recommendation:** This is the PNAS-level claim (cross-cancer generalization). It needs to be bulletproof with complete statistical annotation and biological interpretation of what's being transferred.

---

## Consolidated Recommendations

### High-Impact Changes (Minimal Effort, Maximum Reviewer-Proofing)

| Figure | Change | Impact |
|--------|--------|--------|
| **1** | Add "Supervision acts on W (gene programs)" micro-callout | Emphasizes innovation |
| **2** | Reorder panels: B-C-D above A | Matches problem→solution narrative |
| **3** | Add "true k=3, n=100 replicates" annotation | Provides missing context |
| **4** | ADD FACTOR LABELS to heatmaps | Dramatically improves interpretability |
| **5** | Add HR, CI, p-value to KM plots | Standard practice, prevents criticism |
| **6** | Add HR, p-value, factor label to transfer KM | Completes the PNAS claim |

### Aesthetic/Structural Improvements

| Issue | Current | Recommendation |
|-------|---------|----------------|
| Color consistency | Variable across figures | Standardize: DeSurv=green, NMF=blue/red |
| Legend terminology | "alpha=0" vs "DeSurv" | Use "Standard NMF (α=0)" vs "DeSurv (α>0)" |
| Axis label consistency | "rank" vs "Rank k" | Capitalize consistently |
| Panel labeling | Mixed (A. vs A) | Standardize to bold "A" without period |

### Narrative Arc Improvements

The figures tell this story:
```
Fig 1: Here's DeSurv
Fig 2: NMF rank selection fails → DeSurv works
Fig 3: Simulation validates supervision
Fig 4: THE INSIGHT - variance ≠ prognosis
Fig 5: Generalizes within PDAC
Fig 6: Generalizes across cancers
```

**To strengthen:**
1. Fig 2: Restructure to emphasize solution (Panel A), not just show it
2. Fig 4: Make "variance ≠ prognosis" punchline more prominent with text annotation
3. Fig 5-6: Complete statistical reporting elevates from "looks good" to "is proven"

---

## Execution Checklist (Delegatable)

This checklist is formatted for delegation to a student/postdoc with effort estimates.

### 🔴 HIGH-PRIORITY (Reviewer-blocking)

| Task | Effort | Risk Reduced | Owner |
|------|--------|--------------|-------|
| ☐ Add HR, 95% CI, log-rank p to all KM plots (Figs 5–6) | 1–2 hours | "Incomplete survival reporting" | R/plotting |
| ☐ Add factor labels to Fig 4 heatmaps | 1–2 hours | "Uninterpretable latent factors" | R/plotting |
| ☐ Add simulation context to Fig 3 (true k, # replicates, define precision) | 1 hour | "Toy simulation / unclear setup" | R/plotting |

**Guidance for factor labels:** Use conservative descriptors:
- "Classical tumor–associated" (not "Classical tumor phenotype driver")
- "Basal-like–associated" (not "Aggressive basal-like invasion program")
- "Immune/stromal–enriched" (not "Immunosuppressive CAF-mediated state")

### 🟠 STRONGLY RECOMMENDED

| Task | Effort | Risk Reduced | Owner |
|------|--------|--------------|-------|
| ☐ Add null or harder simulation regime (SI + main text reference) | 3–6 hours | "Cherry-picked simulation" | Analysis |
| ☐ Restructure Fig 2 (problem panels above solution panel) | 2–3 hours | "Confusing rank-selection story" | R/plotting |
| ☐ Add "variance ≠ prognosis" annotation to Fig 4C | 30 min | Missed conceptual insight | R/plotting |
| ☐ Annotate selected k* on Fig 2A | 30 min | Selection transparency | R/plotting |

**Key text for null simulation:** "As expected, under null signal DeSurv does not outperform unsupervised NMF, confirming that observed gains reflect genuine signal recovery."

### 🟡 OPTIONAL / TIME-PERMITTING

| Task | Effort | Risk Reduced | Owner |
|------|--------|--------------|-------|
| ☐ Standardize colors across main figures (DeSurv=green, NMF=red) | 2–4 hours | Visual inconsistency | R/plotting |
| ☐ Clarify biological plausibility in Fig 6 caption (1–2 sentences) | 30 min | "Generic immune signal" critique | Writing |
| ☐ Standardize legend terminology ("Standard NMF (α=0)" vs "DeSurv") | 1 hour | Terminology inconsistency | R/plotting |

### 🔵 SKIP UNLESS REQUESTED

- I² heterogeneity statistics
- Meta-analytic pooling visuals / summary diamonds
- New biological validation analyses
- Additional benchmark methods in figures

---

## Appendix: Alternative Visualization Approaches

### Figure 2 Alternative: Integrated Heatmap

Instead of separate panels for cophenetic/residuals/silhouette, consider:
- Single heatmap with rows = metrics, columns = k
- Color = "recommendation strength" (normalized)
- Shows disagreement more directly

### Figure 4 Alternative: Sankey Diagram

For Panel D (NMF-DeSurv correspondence):
- Sankey/alluvial diagram showing how NMF factors map to DeSurv factors
- Width = correspondence strength
- More intuitive than heatmap for showing "reorganization"

### Figure 5-6 Alternative: Integrated Validation

Consider combining into single "Generalization" figure:
- Panel A: Forest plot (all cohorts including bladder)
- Panel B: PDAC KM
- Panel C: Bladder transfer KM
- Shows generalization story in one figure (useful if need to reduce to 4 figures for PNAS)

---

*Analysis generated: 2026-01-24*
