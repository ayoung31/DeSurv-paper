# DeSurv Manuscript: Revised Introduction and Suggested Text for PNAS

**Purpose:** Drop-in replacement language for the Introduction and related sections, structured for PNAS audience. Builds a logical case that survival supervision should help, so simulation and real data results read as confirmation of testable predictions.

**Three arguments threaded throughout:**
1. The state of the art is a two-step approach (unsupervised NMF, then post-hoc survival evaluation)
2. Biological and statistical reasons why incorporating survival at training time should yield more generalizable factors
3. Survival information also simplifies the problem of tuning parameter selection in NMF

**Structure:** The Introduction makes predictions; Results validate them.

| Prediction (Introduction) | Validation (Results) |
|---------------------------|----------------------|
| Variance-dominant factors need not be prognostically relevant | Fig 4C: exocrine factor explains most variance but least survival signal |
| Supervision should concentrate survival signal into fewer factors | Fig 4A-B: DeSurv factors are more prognostically concentrated than NMF |
| Supervision should improve cross-cohort generalization | Fig 5: DeSurv factors show consistent HR across independent PDAC cohorts |
| Supervision should resolve rank selection ambiguity | Fig 2D: BO landscape shows clear optimum; contrast with disagreeing heuristics in Fig 2A-C |
| Survival-aligned programs may transfer across cancers | Fig 6B: PDAC-trained DeSurv factor retains prognostic signal in bladder cancer |

---

## A. Revised Introduction (5-Paragraph Structure)

### Paragraph 1: The biological problem

> Tumor transcriptomes are mixtures. In bulk RNA-sequencing, the observed expression profile of each sample reflects contributions from malignant cells, cancer-associated fibroblasts, immune infiltrates, and other microenvironmental components. Deconvolution methods based on nonnegative matrix factorization (NMF) have been instrumental in resolving this mixture into additive, interpretable gene programs that often correspond to recognizable cell types or transcriptional states [Lee & Seung 1999; Moffitt et al. 2015; Peng et al. 2019]. However, the programs that explain the most transcriptomic variance are not necessarily the programs that drive clinical outcomes. In pancreatic ductal adenocarcinoma (PDAC), for example, exocrine and tissue-composition signals dominate expression variance but contribute little to survival prediction, while lower-variance programs reflecting basal-like or activated stromal biology carry the strongest prognostic associations [Moffitt et al. 2015; Rashid et al. 2020; Peng et al. 2024]. This disconnect between variance and prognosis is not unique to PDAC: tumor purity alone accounts for a large fraction of observed expression variation across cancer types, and nearly half of genes used for subtyping in glioblastoma and lung adenocarcinoma are correlated with purity rather than biology [Aran et al. 2015]. When unsupervised methods preferentially capture these variance-dominant signals, the resulting programs require extensive downstream filtering to identify the subset with clinical relevance --- a process that is ad hoc, cohort-specific, and difficult to reproduce.

**Citable evidence:**
- Aran et al. 2015 (Nature Communications): 47.1% of GBM subtyping genes and 45.4% of LUAD subtyping genes correlated with tumor purity
- Bair & Tibshirani 2006 (JASA): "there is no guarantee that the resulting principal components will be associated with survival" when using all genes in unsupervised PCA
- Moffitt et al. 2015 (Nature Genetics): virtual microdissection of PDAC into tumor and stroma compartments
- Rashid et al. 2020: purity effects on PDAC subtyping
- Peng et al. 2019 (Nature Communications): DECODER framework for unsupervised compartment deconvolution

### Paragraph 2: Why the two-step approach falls short (tradeoff framing)

> The standard approach to molecular subtyping proceeds in two steps: first, discover latent programs through unsupervised factorization; then, evaluate their clinical relevance through retrospective survival analysis or downstream clustering [Brunet et al. 2004; Bailey et al. 2016; Peng et al. 2019]. This paradigm has yielded important biological insights, including the identification of basal-like and classical tumor programs in PDAC [Collisson et al. 2011; Moffitt et al. 2015] and compartment-specific deconvolution across 33 cancer types [Peng et al. 2019]. However, the two-step design creates a structural misalignment: the objective optimized during discovery (reconstruction error) differs from the criterion used during evaluation (survival association). Unsupervised NMF allocates representational capacity to the directions of greatest variance, which may be dominated by tissue composition, sample processing, or other outcome-neutral sources. Programs that are prognostically relevant but explain modest variance can be diluted across multiple factors or missed entirely. As a consequence, PDAC subtyping efforts spanning a decade produced between two and six proposed subtypes, with classification accuracy of individual signatures ranging from 35% to 90% across cohorts, and random gene sets sometimes performing comparably to published signatures [Schwarzova et al.; Rashid et al. 2020]. The field ultimately converged on a robust basal/classical dichotomy --- but only after extensive retrospective evaluation [Collisson et al. 2011; Moffitt et al. 2015; Bailey et al. 2016; TCGA 2017; consensus classifier 2025]. An alternative is to incorporate outcome information during discovery itself, so that the factorization is guided toward clinically relevant axes from the outset.

**Citable evidence:**
- Schwarzova et al.: PDAC signature cross-validation study showing random gene sets performing comparably to published signatures
- PDAC consensus molecular classifier (Genome Medicine, 2025): after ~14 years, field converged on two tumor classes (Consensus Classical / Consensus Non-classical) and two stroma classes
- survClust (Arora et al.): unsupervised k-means achieved 67.5% accuracy identifying true risk groups vs. 97.15% for survival-weighted clustering; liver cancer unsupervised subtypes p=0.42, survClust p<0.001

**Important note on tone:** This paragraph does NOT dismiss the two-step approach. It acknowledges its successes (basal/classical discovery, DECODER compartments) and frames the issue as a structural misalignment between the optimization objective and the evaluation criterion. This is a tradeoff, not a failure.

### Paragraph 3: Key insight and the case for supervision (1-2 sentences + supporting argument)

> There are both statistical and biological reasons to expect that incorporating survival information during factorization --- rather than only afterward --- can yield more generalizable gene programs. Statistically, sufficient dimension reduction theory establishes that response-guided subspace estimation targets the directions most relevant to the outcome, whereas variance-maximizing projections target different directions that may be orthogonal to the response [Cook & Forzani 2008; Bair & Tibshirani 2006]. When the outcome depends on low-variance features --- as is common when prognostically relevant tumor programs are obscured by high-variance microenvironmental or compositional signals --- unsupervised methods systematically miss these features, while outcome-supervised methods can recover them. Empirically, supervised principal components using only genes pre-filtered by Cox association outperformed standard PCA for survival prediction [Bair & Tibshirani 2006], and survival-weighted clustering identified prognostically distinct subtypes that unsupervised methods missed entirely [Arora et al.]. Biologically, supervision concentrates representational capacity on the programs that distinguish patients with different outcomes, effectively filtering out variance-dominant but prognostically neutral signals (e.g., exocrine content, normal tissue contamination) during learning rather than in a separate post-hoc step. This reorganization has a second benefit: by reducing the effective number of outcome-relevant factors, supervision can simplify the notoriously difficult problem of rank selection in NMF. Standard rank heuristics --- cophenetic correlation, silhouette width, dispersion --- frequently disagree and can overfit [Frigyesi & Hoglund 2008; Brunet et al. 2004], leaving the analyst without a principled criterion. When the objective includes a survival term, rank selection becomes a well-defined optimization problem: choose the rank (and other hyperparameters) that maximize cross-validated concordance with patient outcomes, providing a clinically meaningful model selection criterion.

**Citable evidence:**
- Cook & Forzani 2008 (Statistical Science): principal fitted components --- response-guided dimension reduction targets the central subspace
- Bair & Tibshirani 2006 (JASA): supervised PCA with 17 Cox-filtered genes outperformed all-gene PCA (R^2 = 0.113 vs 0.08); 5-gene supervised set outperformed 70-gene van't Veer predictor
- Arora et al. (survClust): survival-weighted clustering outperformed unsupervised in simulation (97.15% vs 67.5%) and real data (liver cancer survival logrank 21.56 vs 1.71)
- Frigyesi & Hoglund 2008 (Cancer Informatics): cophenetic correlation can overfit; proposed two-step rank determination

**Why this paragraph matters:** It sets up three testable predictions that the Results will validate:
1. DeSurv should concentrate survival signal into fewer factors (validated by Fig 4A-C)
2. DeSurv factors should generalize better across cohorts (validated by Fig 5)
3. BO-based rank selection should resolve the disagreement among heuristics (validated by Fig 2)

### Paragraph 4: What DeSurv does (technical approach)

> Here we present DeSurv, a survival-supervised deconvolution framework that integrates NMF with Cox proportional hazards modeling. The key architectural choice is where survival supervision enters the factorization. In DeSurv, factor scores are defined as linear projections of the expression data onto the gene program matrix: Z = W^T X. The Cox partial likelihood is then a function of W and the associated regression coefficients beta, so that the survival gradient acts directly on the gene programs. Sample-level loadings H are updated solely through the reconstruction objective and receive no survival gradient. This asymmetry is by design: it ensures that gene programs are shaped by both reconstruction fidelity and outcome association, while sample loadings remain unconstrained, preserving the biological interpretability of the deconvolution as a mixture model. Because the learned gene programs W are portable, new samples can be scored by simple projection (Z_new = W^T X_new) without requiring their survival data --- a property not shared by methods that supervise through H [Huang et al. 2020; Le Goff et al. 2025]. To address model selection, DeSurv employs Bayesian optimization over cross-validated concordance to jointly tune the factorization rank k, the supervision strength alpha, and regularization parameters, replacing heuristic rank criteria with a single, outcome-driven objective.

**Technical notes for reviewers:**
- The gradient of the Cox partial likelihood with respect to H is identically zero (because ell depends on W^T X, not on H). This is the mathematical basis for the W-vs-H supervision distinction.
- Model selection (k, alpha, lambda, nu, ngene, ntop) is performed entirely within training data using nested cross-validation. Validation cohorts are held out from all tuning. Reported external hazard ratios reflect purely out-of-sample generalization.

### Paragraph 5: What we show (setup for Results)

> We evaluate DeSurv in three settings that test the predictions made above. First, in simulations with known latent structure and survival associations, we show that DeSurv recovers the true factorization rank and the identity of lethal programs more reliably than unsupervised NMF, and that its advantage increases when prognostic programs explain low variance relative to outcome-neutral programs. Second, in PDAC, we demonstrate that survival supervision reorganizes the latent transcriptomic landscape: DeSurv suppresses variance-dominant but prognostically neutral signals (e.g., exocrine content) and concentrates survival association into a smaller set of biologically interpretable factors aligned with known tumor and microenvironmental programs. These survival-aligned factors generalize across independent PDAC cohorts with consistent hazard ratios and clearer survival separation than their unsupervised counterparts. Third, we show that a PDAC-trained DeSurv factor retains prognostic signal when projected into bladder cancer samples, consistent with prior reports that basal-like transcriptional structure generalizes across epithelial cancers [Damrauer et al. 2014; Hoadley et al. 2018].

---

## B. Related Work Paragraph

*(Place at the end of the Introduction or beginning of the Discussion)*

> DeSurv builds on a substantial body of prior work in transcriptomic deconvolution and molecular subtyping. Virtual microdissection first demonstrated that bulk PDAC expression can be decomposed into tumor- and stroma-specific compartments with distinct prognostic associations, including basal-like and classical tumor programs and activated versus normal stroma [Moffitt et al. 2015]. The DECODER framework formalized unsupervised NMF-based deconvolution and applied it across 33 TCGA cancer types, identifying seven PDAC compartments and showing that derived compartment ratios are associated with survival [Peng et al. 2019]. Two recent studies have proposed survival-aware NMF formulations: CoxNMF [Huang et al. 2020] and SurvNMF [Le Goff et al. 2025], both of which link the Cox objective to sample-level factor scores. Neither has been published in a peer-reviewed venue as of this writing, and neither provides a formal framework for hyperparameter selection or external validation. The goal of DeSurv is not to redefine established PDAC biology but to evaluate whether outcome-guided learning at the gene-program level can improve the stability, portability, and prognostic alignment of latent structure relative to variance-driven factorization. Accordingly, the biological programs recovered by DeSurv overlap with known tumor and microenvironmental states; the methodological contribution lies in how these programs are learned, selected, and validated.

---

## C. Setup-Payoff Connections

The Introduction makes specific, testable predictions. Below is how each prediction maps to Results content, with suggested transitional language.

### Prediction 1: Variance != Prognosis

**Intro setup:** "the programs that explain the most transcriptomic variance are not necessarily the programs that drive clinical outcomes"

**Results payoff (Fig 4C):** Suggested transition:
> "To directly test whether variance and prognostic relevance are aligned, we compared the fraction of expression variance explained by each factor against its survival association. As predicted, the dominant axis of transcriptomic variation --- corresponding to exocrine and tissue-composition signals --- showed negligible survival association, while lower-variance factors reflecting basal-like and activated stromal biology carried the strongest prognostic signal (Fig 4C)."

### Prediction 2: Supervision should concentrate survival signal

**Intro setup:** "supervision concentrates representational capacity on the programs that distinguish patients with different outcomes"

**Results payoff (Fig 4A-B):** Suggested transition:
> "Consistent with this expectation, DeSurv concentrated survival association into fewer factors than unsupervised NMF. In the gene overlap heatmaps, DeSurv factors showed clearer separation between prognostically relevant and irrelevant programs, whereas NMF distributed survival signal more diffusely across factors (Fig 4A-B)."

### Prediction 3: Supervision should improve cross-cohort generalization

**Intro setup:** "response-guided subspace estimation targets the directions most relevant to the outcome"

**Results payoff (Fig 5):** Suggested transition:
> "To assess whether survival-aligned programs generalize beyond the training cohort, we projected the learned gene programs onto independent PDAC datasets. The most prognostic DeSurv factor showed consistent hazard ratios and clearer survival separation than the corresponding NMF factor across all validation cohorts (Fig 5A-C), supporting the prediction that outcome-guided learning improves portability."

### Prediction 4: Supervision should resolve rank selection

**Intro setup:** "rank selection becomes a well-defined optimization problem: choose the rank that maximizes cross-validated concordance"

**Results payoff (Fig 2):** Suggested transition:
> "Standard rank selection heuristics --- cophenetic correlation, silhouette width --- provided ambiguous or conflicting guidance for the PDAC dataset (Fig 2A-C). In contrast, Bayesian optimization over cross-validated concordance identified a clear optimum at k = [value], with the survival-driven criterion resolving the disagreement among unsupervised metrics (Fig 2D)."

### Prediction 5: Cross-cancer transfer

**Intro setup:** "survival-aligned programs learned in PDAC recapitulate known basal-like transcriptional structure"

**Results payoff (Fig 6B):** Suggested transition:
> "Finally, we tested whether DeSurv programs trained on PDAC retain prognostic relevance in a different cancer type. A PDAC-trained DeSurv factor, when projected into bladder cancer samples, stratified patients by survival (Fig 6B), consistent with prior evidence that basal-like transcriptional programs generalize across epithelial cancers [Damrauer et al. 2014; Hoadley et al. 2018]."

---

## D. Targeted Claim Fixes

### Fix 1: Discovery-validation framing

**Current (paragraph 2, line 1):**
> "Separating discovery from validation can risk overfitting and limits biological and clinical generalizability."

**Replace with:**
> "The standard two-step approach --- unsupervised discovery followed by retrospective survival evaluation --- has produced foundational insights, but creates a structural misalignment between the optimization objective (reconstruction error) and the evaluation criterion (survival association). When variance-dominant signals are not prognostically relevant, this misalignment can cause the factorization to allocate capacity to outcome-neutral programs at the expense of clinically meaningful ones."

**Why:** Frames as a structural misalignment, not a criticism. Acknowledges the successes of the two-step approach.

### Fix 2: Prior survival-NMF methods

**Current (paragraph 4):**
> "both integrate the survival objective through the sample-specific loadings rather than the gene-level programs"

**Replace with:**
> "prior survival-aware NMF formulations [Huang et al. 2020; Le Goff et al. 2025] link the Cox objective to sample-level factor scores, with gene programs updated indirectly through their coupling to reconstruction error. In contrast, DeSurv defines the survival objective directly in terms of gene-program projections (Z = W^T X), so that Cox gradients act explicitly on the gene program matrix W. This architectural difference has a practical consequence: because the learned W is outcome-aligned, new samples can be scored by simple projection without requiring survival data, enabling direct external validation."

**Why:** Explains not just the mathematical difference but the practical consequence (portability).

### Fix 3: Cross-cancer claim

**Current:**
> "demonstrating cross-cancer transfer of prognostic programs"

**Replace with:**
> "demonstrating that survival-aligned programs learned in PDAC recapitulate known basal-like transcriptional structure and retain prognostic relevance when projected into bladder cancer, consistent with prior evidence of shared programs across epithelial cancers"

**Why:** Tempers novelty; acknowledges existing knowledge of cross-cancer basal-like programs.

### Fix 4: "Quintessential challenge" language

**Current:**
> "addressing the quintessential challenge of rank determination in matrix factorization"

**Replace with:**
> "replacing heuristic rank selection criteria with a single, outcome-driven objective optimized via Bayesian optimization"

**Why:** More precise; emphasizes the specific contribution (outcome-driven model selection) over flowery language.

---

## E. Revised Significance Statement

> Molecular subtyping of tumors relies on unsupervised discovery of gene programs, followed by retrospective evaluation of their clinical relevance. This two-step paradigm prioritizes expression variance, which can diverge from prognostic importance. We introduce DeSurv, a framework that incorporates patient survival directly into gene program discovery via nonnegative matrix factorization coupled with Cox regression. In DeSurv, survival supervision acts on gene programs rather than sample loadings, yielding signatures that are intrinsically outcome-aligned and portable to new cohorts by simple projection. Applied to pancreatic cancer, DeSurv concentrates survival signal into interpretable programs reflecting tumor and microenvironmental biology, suppresses variance-dominant but prognostically neutral structure, and generalizes across independent cohorts and cancer types.

---

## F. Key References for the Revised Introduction

| Reference | Key finding for DeSurv framing |
|-----------|-------------------------------|
| Bair & Tibshirani 2006 (JASA) | Supervised PCA outperforms unsupervised; "no guarantee principal components will be associated with survival" |
| Cook & Forzani 2008 (Stat Sci) | Sufficient dimension reduction: response-guided subspaces target outcome-relevant directions |
| Aran et al. 2015 (Nat Commun) | Tumor purity dominates expression variance; ~47% of subtyping genes correlate with purity |
| Moffitt et al. 2015 (Nat Genet) | Virtual microdissection of PDAC; basal/classical tumor + activated/normal stroma |
| Peng et al. 2019 (Nat Commun) | DECODER: unsupervised NMF deconvolution across 33 cancers; 7 PDAC compartments |
| Rashid et al. 2020 | Purity effects on PDAC subtyping; convergence on basal/classical dichotomy |
| Collisson et al. 2011 | Original 3-subtype PDAC classification |
| Bailey et al. 2016 (Nature) | 4-subtype PDAC classification |
| Frigyesi & Hoglund 2008 (Cancer Informatics) | Cophenetic correlation can overfit; rank selection instability |
| Arora et al. (survClust) | Survival-weighted clustering outperforms unsupervised (97% vs 68% accuracy in simulation) |
| Huang et al. 2020 (arXiv) | CoxNMF: survival-NMF through H; unpublished, cannot score new samples |
| Le Goff et al. 2025 (PhD thesis) | SurvNMF: survival-NMF through H; unpublished |
| Damrauer et al. 2014 | Basal-like subtype in bladder cancer parallels breast/PDAC |
| Hoadley et al. 2018 (Cell) | Pan-cancer atlas: shared transcriptional programs across epithelial cancers |
| Consensus PDAC Classifier 2025 (Genome Medicine) | After ~14 years, field converges on 2 tumor + 2 stroma classes |
| Bair & Tibshirani 2004 (PLoS Biology) | DLBCL: unsupervised clustering log-rank p=0.416; semi-supervised p=0.0001. Cleanest empirical precedent. |
| Cook 2007 (Fisher Lecture, Stat Sci) | "The principal component subspace and the central subspace are, in general, different objects." More quotable than Cook & Forzani 2008. |
| CoxNTF 2025 (arXiv:2506.06411) | Unsupervised NMF latent representations "degrade survival prediction in most datasets because NMF does not account for survival information." Very recent. |
| Kawaguchi et al. 2023 (ICML) | IB generalization bound: scales with degree of information compression, not number of parameters. Formal support for supervision-as-regularization. |
| Poirion et al. 2021 (Genome Medicine) | DeepProg: semi-supervised multi-omics C-index 0.73-0.80 across HCC cohorts; unsupervised SNF significant in only 13/32 cancers. |
| Li 1991 (JASA) | Sliced inverse regression: foundational SDR method; E(X\|Y) lies in the effective dimension reduction space |
| Ng & Jordan 2001 (NIPS) | Generative classifiers converge in O(log n) samples but discriminative achieve lower asymptotic error. Favors supervised factorization when labels available. |

---

## G. Usage Notes

1. **Paragraph 1** can be shortened for space by cutting the Aran et al. purity statistics; the core claim (variance != prognosis) stands without them
2. **Paragraph 3** is the most novel and critical section --- it provides the theoretical case that reviewers will scrutinize. The Cook/Bair/Tibshirani citations give it statistical grounding beyond hand-waving
3. **Section C (Setup-Payoff)** provides transitional language for the Results section; each transition explicitly connects back to an introduction prediction
4. **The Related Work paragraph (Section B)** is essential for addressing prior art concerns and should appear either at the end of the Introduction or the start of the Discussion
5. The significance statement (Section E) names the specific biological outcome (tumor/microenvironmental programs in PDAC) rather than speaking generically about "gene programs"

---

## H. The Semi-Supervised Argument: Why DeSurv Needs Both Terms

DeSurv is a **semi-supervised** method, not a fully supervised one. The alpha parameter explicitly controls the balance between unsupervised reconstruction (alpha=0) and supervised Cox association (alpha=1). This framing connects to established theory and strengthens the paper's positioning.

### Why not purely unsupervised (alpha = 0)?

Unsupervised NMF optimizes reconstruction error, allocating factors to explain variance. When the outcome depends on low-variance features --- as demonstrated by the exocrine-factor result in Fig 4C --- the factorization systematically underweights clinically meaningful programs. Survival labels steer the factorization toward the outcome-relevant subspace.

This parallels the **compatibility assumption** in semi-supervised learning (Chapelle et al. 2006): supervision helps when the generative structure of the data (expression patterns) and the discriminative structure (survival associations) share an underlying manifold but differ in their dominant axes. In PDAC, basal-like and activated stromal programs exist in the data but are obscured by higher-variance exocrine signals --- exactly the regime where outcome guidance is most beneficial.

### Why not purely supervised (alpha = 1)?

Without the reconstruction term, the model degenerates into penalized Cox regression on arbitrary linear projections. There is no constraint forcing W to represent real expression patterns. The reconstruction term ensures that learned gene programs correspond to actual transcriptomic programs, not statistical artifacts that happen to predict survival.

This is the **biological interpretability guarantee**: DeSurv discovers programs that both explain expression data and predict outcomes. Neither term alone achieves this. Dropping the unsupervised term sacrifices biological interpretability; dropping the supervised term sacrifices outcome alignment.

### The semi-supervised sweet spot

NMF has many equally good solutions --- the reconstruction error surface has many local optima with similar loss values (rotational/permutation ambiguity). Survival labels break this degeneracy by selecting, from among the comparably reconstructive solutions, the one most aligned with clinical outcomes. This is analogous to how semi-supervised NMF uses class labels to constrain the factor space [Liu et al. 2012; Cai et al. 2011].

From an information-theoretic perspective, supervised compression retains only outcome-relevant variation, discarding noise that unsupervised methods would preserve [Achille & Soatto 2018; Tishby et al. 1999]. The alpha parameter controls the compression strength: low alpha preserves more variance (biological completeness); high alpha compresses more aggressively (outcome focus). DeSurv's Bayesian optimization finds the alpha that maximizes cross-validated concordance, letting the data determine the optimal balance.

### Suggested text for the paper

Add to Introduction paragraph 4 (after describing the W-vs-H distinction):

> "DeSurv operates as a semi-supervised method: the reconstruction term (weight 1-alpha) ensures that learned programs correspond to genuine transcriptomic patterns, while the survival term (weight alpha) steers the factorization toward the decomposition --- among the many that explain the data comparably well --- that is most aligned with patient outcomes."

Add to Discussion (in the limitations or broader context paragraph):

> "The semi-supervised design reflects a deliberate tradeoff. Pure unsupervised factorization maximizes biological completeness but may miss outcome-relevant programs that explain modest variance. Pure supervised projection maximizes discrimination but may lose biological interpretability. By jointly optimizing both objectives, DeSurv targets the intermediate regime where programs are both biologically grounded and clinically informative. The alpha parameter, tuned via Bayesian optimization, lets the data determine the appropriate balance for a given cancer type and cohort size."

### Key references

| Reference | Argument |
|-----------|----------|
| Chapelle, Scholkopf, Zien 2006 (MIT Press) | Semi-supervised learning helps when p(x) and p(y\|x) share structure; the compatibility assumption |
| Liu et al. 2012 | Semi-supervised NMF: label constraints improve clustering over unsupervised NMF |
| Cai et al. 2011 | Graph-regularized NMF: incorporating label/structure information into NMF |
| Achille & Soatto 2018 (JMLR) | Information bottleneck: supervised compression retains minimal sufficient statistics |
| Tishby, Pereira, Bialek 1999 | Information bottleneck principle: optimal compression trades fidelity for relevance |

---

## I. Suggested Results Section Language

Each result should connect back to an introduction prediction, creating a **"prediction → finding"** structure that makes the paper read as hypothesis-driven rather than exploratory.

### Rank selection results (Fig 2)

> "Standard rank selection heuristics --- cophenetic correlation and silhouette width --- provided conflicting guidance for the PDAC dataset (Fig 2A-C). Cophenetic correlation suggested k = [X], while silhouette width favored k = [Y], illustrating the well-documented instability of these criteria [Frigyesi & Hoglund 2008]. In contrast, DeSurv's Bayesian optimization over cross-validated concordance identified a clear optimum at k = [Z] (Fig 2D), with the survival-driven criterion resolving the disagreement among unsupervised metrics. This supports the prediction that incorporating clinical outcomes provides a principled, clinically meaningful model selection criterion that sidesteps the ambiguity of purely data-driven heuristics."

### Simulation results (Fig 3)

> "To test the predictions from Section [Intro], we designed simulations with known ground-truth structure. DeSurv recovered the true factorization rank and the identity of lethal programs more reliably than unsupervised NMF, and this advantage was most pronounced when prognostic programs explained low variance relative to outcome-neutral programs (Fig 3A-B). Critically, under a null scenario with no survival signal, DeSurv did not outperform NMF, confirming that the observed advantages reflect genuine signal recovery rather than overfitting to noise. In the mixed scenario, where lethal programs overlapped partially with variance-dominant programs, DeSurv's advantage was present but attenuated --- consistent with the expectation that supervision helps most when variance and prognosis diverge."

**Note:** The null and mixed scenarios are essential. Without them, a reviewer can claim the results are cherry-picked. Frame the null result positively: it validates DeSurv's behavior under the conditions where supervision should NOT help.

### Variance vs. prognosis (Fig 4C)

> "To directly test the central prediction --- that the dominant axes of transcriptomic variation need not be the most prognostically relevant --- we compared the fraction of expression variance explained by each factor against its survival association (Fig 4C). The factor explaining the largest share of variance, corresponding to exocrine and tissue-composition signals, showed negligible survival association (C-index ~ 0.5). In contrast, factors reflecting basal-like and activated stromal biology, which individually explained less variance, carried the strongest prognostic signal. DeSurv reorganized this landscape by reallocating representational capacity: survival signal was concentrated into [N] factors rather than dispersed across the full factorization, while variance-dominant but prognostically neutral programs were deprioritized."

### External validation (Fig 5)

> "Consistent with the prediction from sufficient dimension reduction theory --- that outcome-guided subspaces should be more portable across datasets --- DeSurv factors showed more consistent hazard ratios across independent PDAC cohorts than the corresponding NMF factors (Fig 5A). The most prognostic DeSurv factor maintained significant survival separation in [dataset 1] (HR = [X], 95% CI: [Y-Z], p = [P]) and [dataset 2] (HR = [X], 95% CI: [Y-Z], p = [P]), whereas the best NMF factor showed attenuated or nonsignificant associations in the same cohorts (Fig 5B-C). This supports the hypothesis that outcome-aligned programs capture biology that generalizes, whereas variance-driven programs may capture cohort-specific noise."

### Cross-cancer transfer (Fig 6)

> "Finally, we tested whether DeSurv programs trained in one cancer type retain prognostic relevance in another. A PDAC-trained DeSurv factor, projected onto bladder cancer samples, stratified patients by survival (Fig 6B; HR = [X], p = [P]). This finding is consistent with prior evidence that basal-like transcriptional programs generalize across epithelial cancers [Damrauer et al. 2014; Hoadley et al. 2018] and suggests that survival supervision captures cross-cancer biology that unsupervised methods may distribute across multiple outcome-neutral factors."

---

## J. Suggested Discussion Language

### Paragraph 1: Summary of contributions

> "We have shown that incorporating survival information during NMF-based deconvolution reorganizes the latent transcriptional landscape, concentrating prognostic signal into fewer, more interpretable gene programs while suppressing variance-dominant but outcome-neutral structure. This reorganization is not a prediction improvement in the traditional sense --- DeSurv is not optimized primarily for maximal discrimination --- but rather a structural reallocation of representational capacity toward clinically relevant biology."

### Paragraph 2: PDAC biological findings

> "In PDAC, DeSurv recapitulated known tumor and microenvironmental programs, including the basal-like/classical distinction and activated versus normal stroma. The methodological contribution lies not in the identity of these programs --- which have been established through virtual microdissection [Moffitt et al. 2015], experimental microdissection [Maurer et al. 2019], and unsupervised deconvolution [Peng et al. 2019] --- but in how they are recovered. DeSurv discovers these programs de novo, without requiring pre-specified signatures or post-hoc filtering, and directly quantifies their survival associations during the factorization."

### Paragraph 3: The semi-supervised tradeoff and when supervision helps vs. hurts

> "The semi-supervised design reflects a deliberate balance between biological completeness and clinical relevance. The reconstruction term (weight 1-alpha) ensures that discovered programs explain real expression patterns; the survival term (weight alpha) steers the factorization toward the decomposition most aligned with outcomes. Our simulations provide guidance on when this tradeoff is most beneficial: DeSurv's advantage is greatest when prognostic programs explain modest variance relative to outcome-neutral signals, and vanishes under null conditions where no survival signal exists. DeSurv is therefore not intended to replace unsupervised NMF in all settings. When the goal is exploratory discovery without clinical endpoints, when survival data is sparse or unreliable, or when the application requires comprehensive biological characterization rather than prognostic stratification, unsupervised methods remain appropriate."

### Paragraph 4: Limitations

> "Several limitations should be noted. First, DeSurv assumes Cox proportional hazards, which may not hold in settings such as immunotherapy response where treatment effects are delayed; extensions to more flexible survival models (e.g., accelerated failure time) are a natural direction. Second, the computational cost of Bayesian optimization over cross-validated concordance exceeds that of standard NMF and may limit scalability to very large cohorts or rapid iterative analyses. Third, while the cross-cancer transfer result is encouraging, we tested only one transfer pair (PDAC to bladder), and broader benchmarking across cancer types is needed to establish generality. Fourth, DeSurv's convergence guarantee (Theorem 1 in the supplement) applies to an idealized algorithm; the implementation includes practical modifications (backtracking, gradient clamping) that depart from the theoretical analysis, though empirically these do not affect convergence behavior. Finally, as with any supervised method, DeSurv's value depends on the quality of the outcome data: heavily censored or misannotated survival information will degrade, not improve, the learned programs."

### Paragraph 5: Broader implications

> "More broadly, the principle instantiated by DeSurv --- that outcome-guided dimensionality reduction targets different subspaces than variance-driven reduction --- extends beyond cancer genomics. Sufficient dimension reduction theory [Cook & Forzani 2008] and the information bottleneck framework [Tishby et al. 1999] both predict that supervised compression retains outcome-relevant structure while discarding nuisance variation. DeSurv realizes this principle within the specific constraints of NMF deconvolution, where nonnegativity preserves biological interpretability and the factorization structure enables single-sample scoring. Extending this framework to other high-dimensional survival settings --- multi-omics integration, electronic health records, spatial transcriptomics --- is a natural next step."

---

## K. Accessibility for the General PNAS Audience

PNAS readership spans all scientific disciplines. The paper must be understood by an educated scientist outside computational genomics. Below are specific strategies.

### 1. Lead every section with the intuition, not the formalism

The revised Introduction (Section A) does this: "Tumor transcriptomes are mixtures" → "the programs that explain the most variance are not necessarily the programs that drive outcomes." No equations in the first three paragraphs.

### 2. Use a concrete analogy (optional, for Significance Statement or Introduction)

> "This is analogous to searching for disease risk factors in population data: the most variable trait in a population (e.g., height) is often not the one that predicts the outcome of interest (e.g., a rare genetic variant with large effect). Variance-driven methods find height; outcome-guided methods find the rare variant."

### 3. Frame the contribution as a general principle

The general insight, stated plainly:

> **"In high-dimensional data where outcome-relevant features explain modest variance, incorporating outcome information during dimensionality reduction yields more portable, more interpretable representations than unsupervised methods."**

This applies beyond cancer. Name this in one sentence in the Discussion to signal broader relevance to the PNAS audience.

### 4. Minimize jargon in key paragraphs

| Technical term | PNAS-friendly alternative |
|---------------|--------------------------|
| nonnegative matrix factorization | a deconvolution method that separates mixed signals into additive components |
| cophenetic correlation | standard model selection metrics |
| Cox partial likelihood | a standard survival model |
| factor scores Z = W^T X | how strongly each patient's tumor expresses each gene program |
| the gradient with respect to H is zero | survival information shapes the gene programs directly but does not constrain how individual tumors are decomposed |
| block coordinate descent | an alternating optimization procedure |
| Bayesian optimization | a systematic hyperparameter search strategy |
| concordance index | a measure of how well the model ranks patients by risk |

### 5. Make figures self-explanatory

For a general audience, three panels matter most:

- **Fig 4C** (variance vs. survival): The core insight. Label axes as "Expression variance explained (%)" vs. "Survival association (concordance index)." Add factor annotations directly on the plot ("Exocrine," "Basal-like," "Activated TME").

- **Fig 1** (model schematic): Add a callout: "Survival gradients act on W (gene programs), not H (sample loadings)." This single annotation conveys the key technical innovation visually.

- **Fig 5** (forest plot + KM curves): Add HR, 95% CI, and p-values directly on each KM panel. A general reader should be able to assess the result without reading the text.

### 6. Use the Significance Statement as the primary hook

The revised Significance Statement (Section E) is written at a level understandable by an undergraduate-educated scientist outside the field. It names the specific biology (PDAC, tumor/microenvironmental programs) rather than speaking generically. This is the first text a PNAS editor reads and is the primary determinant of whether the paper is sent for review.

### 7. One-sentence summaries for each results subsection

At the start of each results subsection, include a bolded one-sentence summary:

- **Rank selection:** "Standard metrics disagree; survival-driven optimization provides a clear, principled choice."
- **Simulations:** "DeSurv recovers ground-truth programs more reliably, but not under null conditions."
- **Biological structure:** "The dominant source of expression variance in PDAC --- exocrine content --- is not the dominant source of prognostic signal."
- **External validation:** "Survival-aligned programs generalize across independent PDAC cohorts."
- **Cross-cancer transfer:** "A PDAC-trained program retains prognostic relevance in bladder cancer."

---

## L. Does This Merit an Updated Narrative Arc?

The suggestions above represent a significant reframing of the paper's narrative. The existing NARRATIVE_ARC.md describes the current structure; the proposed changes would restructure the Introduction from a methods-first framing to a biology-first, prediction-validation framing. Key differences:

| Aspect | Current narrative | Proposed narrative |
|--------|-------------------|-------------------|
| **Opening** | Molecular subtyping + deconvolution challenge | Variance != prognosis as the biological problem |
| **Motivation** | Methods gap (no survival-supervised NMF) | Structural misalignment between optimization and evaluation objectives |
| **Theoretical grounding** | None (informal argument) | Cook/Bair/Tibshirani: sufficient dimension reduction theory |
| **Prior art framing** | Criticism of two-step approach | Tradeoff framing; acknowledgment of successes |
| **Results structure** | Sequential presentation | Prediction-validation pairs (intro predicts, results confirm) |
| **Discussion** | 3 short paragraphs | 5 paragraphs: contributions, biology, tradeoff, limitations, broader implications |
| **Semi-supervised framing** | Not present | Explicit: alpha controls the unsupervised-supervised balance |
| **PDAC subtyping context** | Cited but not discussed | Timeline of convergence motivates the need for supervision |

**Recommendation:** Yes, NARRATIVE_ARC.md should be updated to reflect the proposed prediction-validation structure. The updated arc would follow:

1. **Act I (Introduction):** The biological problem → why existing approaches create misalignment → the case for supervision → what DeSurv does → what we will show (predictions)
2. **Act II (Results):** Each finding validates a specific introduction prediction
3. **Act III (Discussion):** Contributions → biological context → semi-supervised tradeoff → limitations → broader implications

This structure transforms the paper from "here is a method and what it finds" to "here is why supervision should help, and here is the evidence that it does" --- the hypothesis-driven framing that PNAS expects.

---

---

## H. Methods Gap Additions (Issues #12, #13, #14)

The following draft text addresses three methods gaps identified during manuscript review. Each gap involves an analysis that appears in the results with no methods support in either the main text or SI Appendix.

### Gap 4: Factor selection criterion for external validation (Issue #12)

**Location:** Add to SI Appendix, "Survival Analysis" subsection, after the paragraph on median splits (line ~827).

> **Factor selection for external validation.** To compare the prognostic utility of individual factors across methods, we selected a single factor per method for the forest plot and Kaplan--Meier analyses (Fig. \ref{fig:extval}). For each method, we fit a univariate Cox proportional hazards model on the training data for each factor score $Z_j = W_j^\top X_{\text{train}}$ ($j = 1,\dots,k$) and computed the change in partial log-likelihood relative to a null (intercept-only) model:
> $$\Delta \ell_j = \ell(\hat{\beta}_j; Z_j) - \ell_0.$$
> The factor with the largest $\Delta \ell_j$ was selected as the most prognostically informative factor in the training data. This factor was then used to compute scores in each external validation cohort by projection ($Z_j = W_j^\top X_{\text{new}}$), and its survival associations were evaluated using hazard ratios, Kaplan--Meier curves, and log-rank statistics as described above. The same selection criterion was applied to both DeSurv and standard NMF factors to ensure a fair comparison. Importantly, factor selection was performed entirely on training data; validation cohort survival outcomes were used only for evaluation.

### Gap 5: Variance vs. survival quantification (Issue #13)

**Location:** Add as new SI Appendix subsection "Factor-level variance and survival decomposition" after "Survival Analysis."

> **Factor-level variance and survival decomposition.** To contrast the transcriptomic importance and prognostic relevance of individual factors (Fig. \ref{fig:bio}C, Fig. \ref{fig:bladder}A), we defined two complementary quantities for each factor $j \in \{1,\dots,k\}$.
>
> *Fraction of variance explained.* Given the gene program matrix $W$ and sample loadings $H$ from a fitted model, the fraction of total expression variance attributable to factor $j$ was computed as
> $$\text{VarExpl}_j = \frac{\lVert W_j H_j^\top \rVert_F^2}{\lVert X \rVert_F^2},$$
> where $W_j$ denotes the $j$th column of $W$, $H_j^\top$ the $j$th row of $H$, and $\lVert \cdot \rVert_F$ the Frobenius norm. The product $W_j H_j^\top$ is the rank-one reconstruction of $X$ from factor $j$ alone, and its squared Frobenius norm divided by the total sum of squares of $X$ gives the proportion of expression variation captured by that factor. Note that because NMF components are not orthogonal, these fractions need not sum to one and may exceed one.
>
> *Survival contribution (change in partial log-likelihood).* For each factor $j$, we computed sample-level factor scores $Z_j = W_j^\top X$ in the training data and fit a univariate Cox proportional hazards model $\text{Surv}(y, \delta) \sim Z_j$. The survival contribution was quantified as
> $$\Delta \ell_j = \ell(\hat{\beta}_j; Z_j) - \ell_0,$$
> where $\ell(\hat{\beta}_j; Z_j)$ is the maximized partial log-likelihood from the fitted univariate Cox model and $\ell_0$ is the partial log-likelihood from a null (intercept-only) model. Larger values of $\Delta \ell_j$ indicate stronger univariate association between factor $j$ and survival. This quantity was computed separately for DeSurv and standard NMF factors on the same training data.

### Gap 6: Gene program correlation analysis (Issue #14)

**Location:** Add as new SI Appendix subsection "Comparison with established PDAC gene programs" before "PDAC Datasets."

> **Comparison with established PDAC gene programs.** To characterize the biological identity of learned factors (Fig. \ref{fig:bio}A--B), we compared factor-specific gene rankings with published PDAC gene program signatures. Reference gene lists were compiled from the following sources: classical and basal-like tumor subtypes from Moffitt et al. [@moffitt2015virtual] and Collisson et al. [@collisson2011subtypes]; immune, activated stromal, and normal stroma programs from Moffitt et al. [@moffitt2015virtual]; tumor microenvironment programs from Elyada et al. [@elyada2019cross]; pro-inflammatory (proCAF) and restraining (restCAF) cancer-associated fibroblast signatures from Peng et al. [@peng2024determination]; and transcriptional subtypes from Puleo et al. [@puleo2018stratification]. Gene lists from Bailey et al. [@bailey2016genomic], DECODER [@peng2019compartmentalized], MSigDB immune signatures, and PurIST classifier genes were excluded to focus on gene programs derived from expression-based subtyping rather than genomic classification or purity estimation.
>
> For each method (DeSurv or standard NMF), the top 50 genes per factor were selected by ranking the gene weights in each column of the fitted $W$ matrix. The analysis was restricted to genes present in both the top-50 list (across all factors) and at least one reference gene set. For each factor $j$ and reference gene set $k$, we computed the Spearman rank correlation between the continuous gene weight vector $W_j$ (restricted to the shared gene set) and a binary indicator vector $v_k$ encoding membership in gene set $k$ (1 if the gene belongs to the set, 0 otherwise). This measures whether genes belonging to a given reference program tend to have higher weights in the factor.
>
> Statistical significance was assessed using Spearman's rank correlation test, with $p$-values adjusted for multiple comparisons within each factor column using the Benjamini--Hochberg procedure. For display, only reference gene programs with at least one factor correlation exceeding $|r| > 0.2$ were retained. Asterisks in Fig. \ref{fig:bio}A--B indicate adjusted $p < 0.1$.

---

**Note on Gap 6:** The current code hardcodes the deCAF (proCAF/restCAF) gene lists directly in `R/figure_targets.R:657-659` rather than loading them from a reference file. Amber should verify that these gene lists match the published versions in Peng et al. 2024 and consider loading them from a shared data file for reproducibility. Additionally, the SCISSORS citation is missing from `references_30102025.bib` — if SCISSORS subtypes appear in the heatmap, a citation needs to be added.

---

*Updated: 2026-02-10*
*Complements: PNAS_SENIOR_EDITOR_REVIEW.md, PNAS_REVIEW.md, FIGURE_ANALYSIS.md, NARRATIVE_ARC.md*
*Informed by: literature searches on NMF generalization, PDAC subtyping history, sufficient dimension reduction theory, semi-supervised NMF, and supervised vs. unsupervised representation learning*
