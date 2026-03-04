# DeSurv PNAS Final Draft — Comprehensive Revision Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Revise the DeSurv manuscript to address all issues identified across three independent AI reviews, producing a final draft ready for PNAS submission.

**Architecture:** Edit-in-place on the four REVISED child Rmd files + paper.Rmd wrapper. No rendering locally — student renders after push. All changes are prose/formatting only; no R code chunks are modified (figures are generated from the targets store and cannot be re-rendered locally).

**Tech Stack:** R Markdown (.Rmd files), BibTeX (.bib), PNAS LaTeX class

**Key Constraint:** We cannot render locally. The student will render using her targets store after we push. Therefore: (1) do NOT modify any R code chunks, (2) do NOT change tar_load/tar_read calls, (3) do NOT alter figure chunk options, (4) only edit prose, LaTeX, YAML, and .bib formatting.

---

## Review Consensus Summary

Three independent reviews converged on these priorities:

| Priority | Issue | Agreement |
|----------|-------|-----------|
| 1 | K-sensitivity Results section (~800 words) is redundant with SI — compress to ~2-3 sentences | All 3 reviews |
| 2 | Discussion opens with Results re-narration — restructure: biology first, broader implications early | All 3 reviews |
| 3 | Discussion limitations over-enumerated (7 items) — consolidate to 3 thematic paragraphs | All 3 reviews |
| 4 | Introduction P3 too cautious / self-objecting | All 3 reviews |
| 5 | Model Overview Results subsection — cut or fold into intro | Reviews 2 & 3 |
| 6 | Prose hedging throughout ("suggesting" → "demonstrating") | Reviews 2 & 3 |
| 7 | Abstract opening sentence is weak (technical, not clinical hook) | Review 2 |
| 8 | Significance Statement final sentence too vague | Reviews 1 & 3 |
| 9 | Title could be more conceptual (less method-name-centric) | Review 3 |
| 10 | PNAS formatting: P not p, "Fig." not "Figure", hyphenation, SI S-prefix | All 3 reviews |
| 11 | Discussion needs more PDAC biology (iCAF-classical D1 implications) | Review 2 |
| 12 | Methods: add simulation parameter values, compress model subsection | Review 2 |
| 13 | Reference formatting inconsistencies | Review 2 |

---

## Task 1: Revise paper.Rmd — Title, Abstract, Significance Statement, Keywords

**Files:**
- Modify: `paper/paper.Rmd` (lines 2, 44-48, 47-48, 53-58)

### Step 1: Revise the title

Current:
```
title: Survival driven deconvolution (DeSurv) reveals prognostic and interpretable cancer subtypes
```

Revised:
```
title: Survival-guided matrix factorization reveals prognostic transcriptional programs in pancreatic cancer
```

Rationale: PNAS titles emphasize conceptual advance over method names. "DeSurv" can be introduced in the abstract and intro. Hyphenate "Survival-guided". Anchor in PDAC biology. All 3 reviews flagged that the title should be more conceptual.

### Step 2: Revise the abstract opening and add quantitative result

Current opening: "Bulk tumor transcriptomes reflect a mixture of malignant cells and microenvironmental components..."

Revised opening: "Identifying the gene programs that drive patient outcomes from bulk tumor expression data remains a central challenge because the transcriptional programs explaining the most variance are often not those most associated with survival."

Also add one quantitative validation result per Review 1's recommendation. Near the end, change "consistent hazard ratios and clearer survival separation" to include: "consistent hazard ratios (pooled validation HR = 1.93, 95% CI 1.55–2.41, P < 0.001) and clearer survival separation".

The abstract should remain under 250 words. Current is ~230; with these edits it should stay within limit.

### Step 3: Revise the Significance Statement

Current final sentence: "This advances tumor deconvolution and provides a general tool for identifying prognostically informative drivers of disease progression."

Revised: "This provides a general framework for identifying prognostically informative gene programs without requiring post hoc filtering or pre-specified signatures, with potential applications wherever high-dimensional nonnegative measurements contain both signal and nuisance variation."

### Step 4: Verify keywords are complete

Current keywords look fine. No change needed.

### Step 5: Commit
```
git add paper/paper.Rmd
git commit -m "revise title, abstract, and significance statement for PNAS"
```

---

## Task 2: Revise Introduction (02_introduction_REVISED.Rmd)

**Files:**
- Modify: `paper/02_introduction_REVISED.Rmd`

### Step 1: Tighten Paragraph 1 (lines 15)

Current P1 is solid. Minor edit: ensure it opens with the clinical problem more directly. The current opening "Identifying transcriptional signatures associated with clinical outcomes has become central to molecular subtyping..." is acceptable but could be slightly stronger.

Revised opening sentence: "Identifying the transcriptional programs that predict clinical outcomes is central to molecular subtyping and prognostic stratification in cancer."

Remove "has become" (throat-clearing). Rest of P1 stays.

### Step 2: Tighten Paragraph 2 (line 17)

P2 is the longest paragraph (~200 words). Two changes:
1. The sentence beginning "Over a decade of PDAC subtyping efforts proposed between two and five subtypes..." is important context but runs long. Compress to: "Over a decade of PDAC subtyping ultimately converged on a robust basal/classical dichotomy confirmed across independent cohorts [@rashid2020purity; @puleo2018stratification], single-cell analyses [@chansengyue2020transcription; @werba2023single], and prospective profiling [@aung2018compass]."
2. Change "Identifying the prognostically relevant subset then requires extensive downstream filtering that is ad hoc, cohort-specific, and difficult to reproduce" — keep as is, it's strong.

### Step 3: Compress Paragraph 3 (line 19) — the self-objection paragraph

This is the paragraph all 3 reviews flagged. Current is ~150 words of self-objection.

Compress to ~60 words:
"Whether this advantage extends to NMF-based deconvolution is less clear: nonnegative constraints restrict the factor space geometry, censored survival outcomes provide a noisier supervisory signal than continuous responses, and survival gradients that reshape gene programs too aggressively could compromise the mixture-coefficient interpretation of sample loadings. These tensions have not been empirically resolved."

Cut the opening throat-clearing sentence about SDR ("Whether incorporating survival information during factorization improves generalizability..."). Move the SDR citation to the preceding sentence in P2 where it fits naturally.

### Step 4: Trim Paragraph 4 (line 21) — DeSurv description

The sentence "DeSurv operates as a semi-supervised method: a parameter $\alpha$ balances..." through "preventing the survival term from overwhelming the factorization" is repeated in Methods. Cut it here; the Methods section already covers it. Keep only: the architectural W-vs-H choice, the projection property, and the BO selection sentence.

Revised P4 (~120 words):
"Here we present DeSurv, a survival-supervised deconvolution framework that integrates NMF with Cox proportional hazards modeling. The key architectural choice is that the survival gradient acts on the gene program matrix $W$ but not on the sample loadings $H$, directing gene programs toward outcome-relevant structure while preserving the interpretation of sample loadings as mixture coefficients. Because $W$ defines shared transcriptomic programs, new samples can be scored by projection ($Z_{\text{new}} = \tilde{W}^\top X_{\text{new}}$) without requiring their survival data — a property not shared by methods that route supervision through $H$ [@huang2020low; @le2025survnmf]. A supervision parameter $\alpha$ and the factorization rank $k$ are jointly selected via cross-validated concordance (Methods)."

### Step 5: Tighten Paragraph 5 (line 23) — evaluation roadmap

Minor edits only:
- Change "These results suggest that the prognostically relevant structure..." to "These results demonstrate that..."
- Change "suggesting that outcome supervision during factorization yields" to "indicating that outcome supervision during factorization yields"

### Step 6: Commit
```
git add paper/02_introduction_REVISED.Rmd
git commit -m "tighten introduction: compress P3 self-objections, trim P4 redundancy with Methods"
```

---

## Task 3: Revise Results (04_results_REVISED.Rmd) — Major Structural Surgery

**Files:**
- Modify: `paper/04_results_REVISED.Rmd`

This is the highest-impact task. Three operations:

### Step 1: Cut or drastically trim the "Model Overview" subsection (lines 110-116)

Current (lines 110-112): A 3-sentence subsection that re-summarizes DeSurv.

Replace with a single transitional sentence that goes directly to findings:
"We developed DeSurv, a survival-supervised deconvolution framework (Fig. \ref{fig:schema}; Methods), and evaluated it in controlled simulations and PDAC cohort analysis with external validation."

Remove the sentences about W/H gradient distinction (already in intro and methods) and the citation to lee1999learning/gaujoux2010flexible (cited elsewhere).

### Step 2: Strengthen prose in "Survival supervision clarifies NMF rank selection" (lines 119-127)

Line 127: Change "These results are consistent with the expectation that incorporating outcome information into the factorization itself improves recovery of the true factorization rank." →
"These results confirm that incorporating outcome information into factorization improves recovery of the true rank."

### Step 3: Strengthen prose in "DeSurv recovers prognostic gene programs" (lines 169-175)

Line 173: Change "This confirms that outcome-guided learning recovers prognostic programs more reliably when variance and prognosis diverge." → keep as is (already assertive with "confirms").

Line 175: Change "This indicates that the semi-supervised design does not impose spurious prognostic structure when none exists." → keep as is (already assertive).

### Step 4: Edit "Survival supervision reorganizes the learned factor structure in PDAC" (lines 187-197)

Line 195: Replace "A striking pattern emerged:" with "The pattern was clear:" (remove editorial flourish per Review 2).

Line 195: Replace "This apparent paradox — given the well-established prognostic role of basal-like biology in PDAC" with "This result — which might appear to contradict the prognostic role of basal-like biology in PDAC" (soften dramatic framing).

Line 189: The long sentence beginning "To enable a direct comparison..." (49 words) — split into two sentences:
"To enable a direct comparison of how survival supervision reorganizes factor structure at the same dimensionality, we fit standard NMF at DeSurv's BO-selected rank ($k = 3$). Results for standard NMF at its own BO-selected rank ($k = 7$) and at the elbow-selected rank ($k = 5$) are presented in the SI Appendix."

### Step 5: Edit "Survival-aligned programs generalize" (lines 345-351)

Line 349: Change "suggesting that variance-driven programs may capture cohort-specific variation that does not transfer as reliably" →
"consistent with variance-driven programs capturing cohort-specific variation that limits transferability"

### Step 6: COMPRESS the K-sensitivity section (lines 354-407) — MAJOR

This is the single highest-impact edit. The current section runs ~800 words (lines 354-407) with detailed per-cell p-values from the K×α grid.

Replace the ENTIRE prose content (lines 400-407, everything after the R code chunks) with approximately 4-5 sentences:

```
The 1-SE model selection rule chose $k = 3$ over the global BO maximum ($k = 7$; CV C-index 0.655 versus 0.646; margin 0.009). To assess whether this choice reflects a robust solution, we conducted an exhaustive sensitivity analysis across $K \in \{2, 3, 5, 7, 9\}$ and $\alpha \in \{0, 0.25, 0.35, 0.55, 0.75, 0.85, 0.95\}$ (35 combinations) with independent external validation (SI Appendix, Tables S1--S3). With the BO-selected gene-focusing parameter ($n_{\text{top}} = 270$), $K = 3$ achieved adjusted significance across a broad supervision range, while retaining all genes reduced validated combinations and eliminated $K = 3$ adjusted significance at the production supervision strength --- identifying gene focusing rather than regularization as the critical hyperparameter enabling generalization. The global optimum $k = 7$ did not achieve adjusted external validation with $n_{\text{top}} = 270$, while $k = 3$ validated robustly, providing empirical justification for the 1-SE rule's parsimony preference.
```

IMPORTANT: Keep ALL R code chunks (lines 356-398) untouched — they compute scalars used elsewhere. Only replace the prose that follows them.

### Step 7: Add synthesis sentence after the validation section

After the compressed K-sensitivity paragraph, add a brief synthesis that bridges back to the simulation framework:

"Together, the PDAC results are consistent with the intermediate regime between the primary and mixed simulation scenarios: survival supervision reorganizes factor structure, concentrates prognostic signal, and yields more transportable programs than the discover-then-evaluate paradigm."

### Step 8: Commit
```
git add paper/04_results_REVISED.Rmd
git commit -m "compress K-sensitivity to ~4 sentences, cut Model Overview, strengthen prose"
```

---

## Task 4: Revise Discussion (05_discussion_REVISED.Rmd) — Major Restructure

**Files:**
- Modify: `paper/05_discussion_REVISED.Rmd`

### Target structure (5 paragraphs, ~750 words):

**P1 (~100 words):** Key finding in biological terms + broader significance. Combines the current opening sentence with the current FINAL paragraph's broader implications.

**P2 (~150 words):** What DeSurv reveals about PDAC biology — the iCAF-classical co-program in D1, what joint tumor-stroma factor structure implies, connection to CAF subtype literature. Currently underdeveloped (Review 2 flagged this).

**P3 (~120 words):** When DeSurv is most useful and when it is not — simulation guidance, α=0 special case. Current P4-equivalent. Minor edits only.

**P4 (~120 words):** Semi-supervised design and parsimony mechanism — how supervision creates conditions for the 1-SE rule. Cut the K×α re-narration. Current P3-equivalent, trimmed.

**P5 (~250 words):** Limitations in 3 thematic groups + future directions. Consolidate current 7 limitations into: (a) statistical assumptions/scope, (b) generalizability, (c) biological interpretation of supervised factors. Close with multi-omics/spatial extension.

### Step 1: Write new P1 — biological significance + broader implications

Replace current P1 (lines 7) entirely. New P1:

"Aligning the discovery objective with the evaluation criterion — incorporating survival outcomes directly into NMF factorization rather than assessing them retrospectively — reorganizes the transcriptional landscape of PDAC in a clinically meaningful way. DeSurv concentrates prognostic signal into fewer, more interpretable gene programs that generalize across independent cohorts without retraining, achieving with three factors what standard NMF required seven to eight to approach. More broadly, the principle that outcome-guided dimensionality reduction targets different subspaces than variance-driven reduction extends beyond cancer genomics: sufficient dimension reduction theory [@cook2007fisher] and the information bottleneck framework [@tishby1999information] both predict that supervised compression retains outcome-relevant structure while discarding nuisance variation. DeSurv realizes this principle within the specific constraints of NMF deconvolution, where nonnegativity preserves biological interpretability and the factorization structure enables single-sample scoring by projection."

### Step 2: Write new P2 — PDAC biology (currently underdeveloped)

New P2:

"In PDAC, the reorganization produced by survival supervision is biologically informative. The survival-dominant factor (D1) integrates classical tumor programs with an iCAF-associated stromal signature, consistent with growing evidence that clinical outcomes depend on the joint configuration of tumor-intrinsic and cancer-associated fibroblast subtypes [@peng2024determination; @elyada2019cross]. This co-occurrence within a single factor suggests that the prognostically relevant biology in PDAC is not purely tumor-intrinsic but reflects tumor-stroma coupling that unsupervised methods distribute across multiple factors. Standard NMF, by contrast, allocated a dedicated factor to exocrine-compositional variation — a signal that reflects tumor purity rather than cancer biology [@rashid2020purity] — at the expense of resolving the microenvironmental programs that DeSurv separates. The methodological contribution lies not in the identity of these programs, which have been established through virtual microdissection [@moffitt2015virtual], experimental microdissection [@Maurer2019], and unsupervised deconvolution [@peng2019novo], but in how they are recovered: de novo, without pre-specified signatures, with direct quantification of their survival associations."

### Step 3: Write new P3 — when DeSurv helps and when it doesn't

Adapt current P4 (line 13) with minor edits. Keep the simulation guidance paragraph largely as is — it's well-written. Trim slightly.

### Step 4: Write new P4 — semi-supervised mechanism, cut K×α re-narration

Take current P3 (line 11) but CUT the second half that re-narrates the K×α sensitivity results (which now live in the compressed Results section and SI). Keep only the mechanistic explanation: supervision creates flat concordance plateau → 1-SE rule exploits it → parsimony is a property of the rule, not supervision per se.

### Step 5: Write new P5 — consolidated limitations + future directions

Consolidate the current 7 enumerated limitations (lines 15-17) into 3 thematic groups:

**(a) Statistical assumptions and scope (~80 words):** Cox PH assumption, alternative endpoint definitions, convergence gap between theory and implementation.

**(b) Generalizability (~60 words):** PDAC only, need cross-cancer benchmarking, computational cost of BO.

**(c) Biological interpretation (~80 words):** Supervised factors conditional on training cohort's treatment context, residual confounding possible, outcome data quality dependency.

Close with future directions (~30 words): multi-omics, spatial transcriptomics, alternative outcome models, single-cell annotation.

### Step 6: Commit
```
git add paper/05_discussion_REVISED.Rmd
git commit -m "restructure discussion: biology first, consolidate limitations, broader implications early"
```

---

## Task 5: Revise Methods (03_methods_REVISED.Rmd) — Targeted Edits

**Files:**
- Modify: `paper/03_methods_REVISED.Rmd`

### Step 1: Compress "The DeSurv Model" subsection (lines 11-23)

Cut the paragraph about optimization algorithm details (lines 22-23: "Optimization proceeds by alternating updates for $H$, $W$, and $\beta$, using multiplicative rules for $H$, multiplicative projected-gradient updates for $W$, and coordinate descent for $\beta$."). Replace with:
"Optimization alternates updates for $H$, $W$, and $\beta$ (SI Appendix). Although non-convex, these updates converge to a stationary point under mild conditions (SI Appendix)."

### Step 2: Remove redundancy in Hyperparameter Selection (lines 25-28)

The sentence "Because $\alpha$ is selected via cross-validated concordance, values large enough to overfit to cohort-specific survival patterns are penalized by poor out-of-sample performance, preventing the survival term from overwhelming the factorization." — this is stated identically in the Introduction. Cut it from Methods (keep in intro where it motivates the design).

### Step 3: Add simulation parameter values (lines 31-34)

The current simulation section says "marker genes were simulated to load strongly on a single factor and weakly on others" without parameter values. Add after the description of the three gene classes:

"Specifically, marker gene loadings were drawn from $\text{Gamma}(\text{shape}=10, \text{rate}=1)$ for the assigned factor and $\text{Gamma}(0.5, 1)$ for others; background gene loadings from $\text{Gamma}(5, 1)$ across all factors; and noise gene loadings from $\text{Gamma}(0.1, 1)$. Sample loadings $H$ were drawn from $\text{Gamma}(3, 1)$. Cox regression coefficients were $\beta = (1, 1, 0)^\top$ in the primary scenario, $\beta = 0$ in the null scenario, and $\beta = (1, 0.5, 0)^\top$ with shared background contribution in the mixed scenario."

NOTE: These values need to be verified against the actual simulation code. If the exact values differ, the student should correct them during rendering review. Add a comment:
`<!-- TODO: Verify gamma parameters against DeSurv simulation code -->`

### Step 4: Commit
```
git add paper/03_methods_REVISED.Rmd
git commit -m "compress methods model section, add simulation parameters, remove redundancy"
```

---

## Task 6: PNAS Formatting Pass — All Files

**Files:**
- Modify: all four child Rmd files + paper.Rmd

### Step 1: P not p for p-values

Search all files for lowercase "p =" and "$p" in the context of p-values and change to uppercase "$P".

IMPORTANT: Only change p-values, not the variable $p$ used for number of genes. The gene-count variable $p = 3{,}000$ should remain lowercase.

Specific patterns to change:
- `$p = ` when followed by a p-value → `$P = `
- `$P <` — already correct if present
- `p-value` in prose → `P value` (PNAS style, no hyphen)
- `$p$` when referring to significance → `$P$`

BUT KEEP lowercase $p$ when it refers to number of genes (e.g., "$p = 3{,}000$ genes").

### Step 2: Hyphenation

Search and fix:
- "survival driven" → "survival-driven" (compound modifier before noun)
- "variance dominant" → "variance-dominant"
- "outcome neutral" → "outcome-neutral"
- "survival supervised" → "survival-supervised"
- "semi supervised" → "semi-supervised" (verify consistent)

### Step 3: "Fig." not "Figure" in text references

PNAS convention uses "Fig. 1" not "Figure 1" in body text. The LaTeX \ref{} calls produce "Figure" by default from the PNAS class, so this may be handled automatically. Check the rendered PDF — in the PDF I read, references appear as "Fig. 1", "Fig. 2" etc., so this may already be correct. If any prose says "Figure" spelled out, change to "Fig."

### Step 4: Commit
```
git add paper/
git commit -m "PNAS formatting: P values, hyphenation, figure references"
```

---

## Task 7: Reference Formatting (references_30102025.bib)

**Files:**
- Modify: `paper/references_30102025.bib`

### Step 1: Fix journal name capitalization

Check that journal names use proper capitalization. Known issues from Review 2:
- Line 206: `journal={nature}` → `journal={Nature}` (for lee1999learning)
- Line 217: `journal={Proceedings of the national academy of sciences}` → `journal={Proceedings of the National Academy of Sciences}` (for brunet2004metagenes)

### Step 2: Verify Le Goff et al. publication status

Reference 22 (le2025survnmf) is listed as a PhD thesis. Check whether it has been formally published. If not, keep as is but the authors should note this to PNAS editors. Add a comment in the .bib:
`% NOTE: Verify if Le Goff et al. has been published in a journal by submission time`

### Step 3: Verify Peng et al. 2024 publication status

Reference 11 (peng2024determination) is a bioRxiv preprint. Check if published. Add similar comment.

### Step 4: Commit
```
git add paper/references_30102025.bib
git commit -m "fix reference formatting: journal capitalization, verify preprint status"
```

---

## Task 8: Final Review Pass — Prose Quality

**Files:**
- Modify: all four child Rmd files

### Step 1: Hedging reduction pass

Systematically find and replace hedging language per Review 2's table:

| Find | Replace |
|------|---------|
| "suggesting that" (when result is solid) | "demonstrating that" or "indicating that" |
| "consistent with the expectation that" | "confirms that" or "consistent with" (shorter) |
| "a striking pattern emerged" | Already addressed in Task 3 |
| "This apparent paradox" | Already addressed in Task 3 |

Apply judgment: some hedging is appropriate for genuinely uncertain claims (e.g., cross-cancer generalization). Only strengthen claims that are well-supported by the presented evidence.

### Step 2: Vary repeated phrasing

The phrase "variance and prognosis diverge" appears 5+ times. Replace 2-3 instances with alternatives:
- "the misalignment between variance and outcome relevance"
- "when the highest-variance programs are not the most prognostic"
- "the gap between reconstruction-optimal and survival-optimal subspaces"

### Step 3: Break overly long sentences

Any sentence >45 words should be examined for splitting. Key targets identified in reviews:
- The 49-word comparison sentence in Results (already split in Task 3, Step 4)
- Check Discussion for similar cases after restructuring

### Step 4: Add per-cohort KM reference

In the validation section of Results, after the pooled HR result, add parenthetical: "(per-cohort KM curves in SI Appendix, Fig. S5C--F)"

### Step 5: Commit
```
git add paper/
git commit -m "prose quality: reduce hedging, vary repeated phrases, break long sentences"
```

---

## Task 9: Supplement Cross-References — Verify Consistency

**Files:**
- Review (read-only): `paper/supplement.Rmd`
- Modify: `paper/04_results_REVISED.Rmd` (SI references in prose)

### Step 1: Verify SI figure numbering in main text matches supplement

Cross-check the SI figure mapping comment at the top of 04_results_REVISED.Rmd (lines 5-15) against actual supplement figure order. Ensure all "SI Appendix, Fig. S7" references point to the correct figure.

### Step 2: Note for student — supplement "Figure" → "Fig. S" prefix

The supplement currently uses "Figure 1", "Figure 2" etc. PNAS requires "Fig. S1", "Fig. S2". Flag this for the student to fix in supplement.Rmd (we should not modify supplement.Rmd's R code chunks since we can't render).

Add a comment at the top of supplement.Rmd noting this convention:
`<!-- PNAS CONVENTION: All figures should be labeled "Fig. S1", "Fig. S2", etc. not "Figure 1" -->`

### Step 3: Commit
```
git add paper/
git commit -m "verify SI cross-references, add supplement figure naming note"
```

---

## Execution Order and Dependencies

```
Task 1 (paper.Rmd: title/abstract/sig) ──────────┐
Task 2 (Introduction) ────────────────────────────┤
Task 3 (Results — major surgery) ─────────────────┤── All independent, can be
Task 4 (Discussion — major restructure) ──────────┤   done in any order
Task 5 (Methods — targeted edits) ────────────────┤
Task 7 (References .bib) ─────────────────────────┘
                                                   │
                                                   ▼
Task 6 (PNAS formatting pass — all files) ─────────┐── Depends on Tasks 1-5
Task 8 (Prose quality pass — all files) ───────────┤   (need final text first)
Task 9 (SI cross-references) ─────────────────────┘
```

Tasks 1-5 and 7 are independent and can be executed in any order.
Tasks 6, 8, and 9 should run after all content edits are complete.

---

## Word Budget Estimate

| Section | Current (~words) | Target (~words) | Change |
|---------|-----------------|----------------|--------|
| Abstract | ~230 | ~240 | +10 (add HR) |
| Introduction | ~850 | ~700 | -150 |
| Results: Model Overview | ~60 | ~25 | -35 |
| Results: Rank selection | ~250 | ~230 | -20 |
| Results: Simulations | ~350 | ~340 | -10 |
| Results: PDAC biology | ~650 | ~620 | -30 |
| Results: Validation | ~250 | ~260 | +10 |
| Results: K-sensitivity | ~800 | ~120 | **-680** |
| Discussion | ~1200 | ~750 | **-450** |
| Methods | ~550 | ~550 | 0 |
| **Total body** | **~5190** | **~3835** | **-1355** |

This brings the paper from ~5200 words to ~3800, well within PNAS's ~4500 word limit.

---

## Post-Edit Checklist (for student rendering)

After push, the student should verify:
1. [ ] Paper renders without errors from targets store
2. [ ] All `tar_load`/`tar_read` calls still resolve correctly
3. [ ] Figure numbering (Fig. 1-4) unchanged
4. [ ] SI figure cross-references (S1-S9) match supplement
5. [ ] Supplement figure labels changed to "Fig. S1" etc.
6. [ ] Reference numbers in rendered PDF match expectations
7. [ ] Word count of body text (excl. Methods, refs, figure legends) is within PNAS limits
8. [ ] Simulation gamma parameters verified against DeSurv code
9. [ ] Le Goff et al. and Peng et al. publication status checked
