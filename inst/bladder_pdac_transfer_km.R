##############################################################################
# inst/bladder_pdac_transfer_km.R
#
# Cross-cancer transfer KM figure: project the PDAC-trained DeSurv W matrix
# (tcgacptac) onto the FULL IMVigor bladder cancer cohort (NE-like and
# BCG-treated samples excluded, n ≈ 260) and plot Kaplan-Meier curves
# using the PDAC log-rank-optimal cutpoint.
#
# Addresses GitHub Issue #23:  the existing cross-cancer result in the
# manuscript used only the 20% held-out test split from the bladder-specific
# pipeline (n ≈ 51).  Because bladder samples were NEVER in the PDAC training
# set, all bladder samples are genuinely held-out for the cross-cancer
# projection, so the full cohort is the appropriate evaluation set.
#
# The preprocessing applied here mirrors what targets_common_pipeline.R does
# for PDAC external validation cohorts:
#   1. Apply samp_keeps filter (NE-like / BCG exclusions from data loader)
#   2. Filter to samples with valid survival data
#   3. Restrict expression to PDAC training gene set
#   4. Apply rank transformation per sample (method_trans_train = "rank")
#   5. Compute linear predictor with PDAC W, beta, ntop
#   6. Normalise with PDAC training lp_mean / lp_sd
#   7. Stratify using PDAC log-rank-optimal z-cutpoint
#
# Run from the DeSurv-paper repo root:
#   Rscript inst/bladder_pdac_transfer_km.R
##############################################################################

suppressMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(cowplot)
  pkgload::load_all("../DeSurv", quiet = TRUE)
})

source("R/cv_grid_helpers.R")     # compute_lp()
source("R/load_data_bladder_vig.R")  # load_data_bladder_vig_nobcg()

STORE  <- "store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main"
OUT    <- "figures/bladder_pdac_transfer_km.pdf"

# ── 1. Load PDAC targets from store ─────────────────────────────────────────
cat("Loading PDAC targets from store...\n")
pdac_fit    <- readRDS(file.path(STORE, "objects/tar_fit_desurv_tcgacptac"))
pdac_stats  <- readRDS(file.path(STORE, "objects/desurv_lp_stats_tcgacptac"))
pdac_ntop   <- readRDS(file.path(STORE, "objects/tar_ntop_value_tcgacptac"))
bo_bundle   <- readRDS(file.path(STORE, "objects/bo_bundle_selected_tcgacptac"))
X_train     <- bo_bundle$data_filtered$ex   # rank-transformed PDAC training matrix

cat(sprintf("  W matrix dims:       %d genes x %d factors\n",
            nrow(pdac_fit$W), ncol(pdac_fit$W)))
cat(sprintf("  ntop:                %d\n", pdac_ntop))
cat(sprintf("  z_cutpoint (logrank):%.4f  (will be applied to factor-1-gene LP)\n",
            pdac_stats$optimal_z_cutpoint))

# ── 2. Load full IMVigor bladder cohort (NE-like + BCG excluded) ─────────────
cat("\nLoading bladder cohort...\n")
raw_bl <- load_data_bladder_vig_nobcg(
  "data/original/IMVigor210_sampInfo_ex_clinical.rds"
)

# Apply samp_keeps (NE-like exclusion from load_data_bladder_vig,
# plus BCG exclusion from load_data_bladder_vig_nobcg)
idx1 = which(raw_bl$clinical$Tissue=="bladder")
keep_idx <- raw_bl$samp_keeps
keep_idx = intersect(keep_idx,idx1)
X_bl     <- raw_bl$ex[, keep_idx, drop = FALSE]
sinfo_bl <- raw_bl$sampInfo[keep_idx, , drop = FALSE]

cat(sprintf("  Samples after NE-like + BCG filters: %d\n", ncol(X_bl)))

# ── 3. Filter to samples with valid survival data ────────────────────────────
valid_samp <- which(
  is.finite(sinfo_bl$time) &
    sinfo_bl$time > 0 &
    !is.na(sinfo_bl$event) &
    sinfo_bl$event %in% c(0, 1)
)
X_bl     <- X_bl[, valid_samp, drop = FALSE]
sinfo_bl <- sinfo_bl[valid_samp, , drop = FALSE]

cat(sprintf("  After survival filter:               %d (%d events)\n",
            ncol(X_bl), sum(sinfo_bl$event)))

# ── 4. Restrict to PDAC training gene set ────────────────────────────────────
# Mirrors preprocess_validation_data(genes = pdac_genes, ...) in the pipeline.
pdac_genes   <- rownames(pdac_fit$W)
common_genes <- intersect(pdac_genes, rownames(X_bl))
cat(sprintf("  Common genes (PDAC training ∩ bladder): %d / %d PDAC genes\n",
            length(common_genes), length(pdac_genes)))

X_bl_sub <- X_bl[common_genes, , drop = FALSE]

# ── 5. Rank-transform per sample (method_trans_train = "rank") ───────────────
# Replicates DeSurv::preprocess_data with method_trans_train = "rank":
# each sample's gene expression values are replaced by their within-sample ranks.
# This is the same transformation applied to PDAC training and external
# validation datasets.
X_bl_rank <- apply(X_bl_sub, 2, rank, ties.method = "average")
dimnames(X_bl_rank) <- dimnames(X_bl_sub)

# ── 6. Project onto factor 1 using all common genes ─────────────────────────
# Use all 1966 common genes with only the factor 1 column of W and beta[1].
w1    <- pdac_fit$W[common_genes, , drop = FALSE]   # all common genes x 1
beta1 <- pdac_fit$beta
lp_bl <- compute_lp(w1,beta1,X_bl_rank,pdac_ntop)

# ── 7. Recompute training LP stats using the same factor 1 projection ────────


z_bl     <- (lp_bl - pdac_stats$lp_mean) / pdac_stats$lp_sd
group_bl <- factor(
  ifelse(z_bl > z_cut, "High", "Low"),
  levels = c("Low", "High")
)

cat(sprintf("  Low / High split:    %d / %d\n",
            sum(group_bl == "Low"), sum(group_bl == "High")))

# ── 8. Assemble survival data frame ─────────────────────────────────────────
df_bl <- data.frame(
  time  = sinfo_bl$time,
  event = as.integer(sinfo_bl$event),
  group = group_bl
)

# ── 9. Compute survival statistics ──────────────────────────────────────────
sfit    <- survfit(Surv(time, event) ~ group, data = df_bl)
cox_fit <- coxph(Surv(time, event) ~ group, data = df_bl)
cox_ci  <- summary(cox_fit)$conf.int
hr      <- cox_ci[1, "exp(coef)"]
ci_lo   <- cox_ci[1, "lower .95"]
ci_hi   <- cox_ci[1, "upper .95"]
lr_test <- survdiff(Surv(time, event) ~ group, data = df_bl)
p_val   <- 1 - pchisq(lr_test$chisq, df = 1)

cat(sprintf("\nResults:\n"))
cat(sprintf("  n = %d  (%d events)\n", nrow(df_bl), sum(df_bl$event)))
cat(sprintf("  HR (High vs Low) = %.2f  (95%% CI %.2f-%.2f)\n", hr, ci_lo, ci_hi))
cat(sprintf("  Log-rank p = %.4f\n", p_val))

# ── 10. Build KM plot ────────────────────────────────────────────────────────
p_label  <- if (p_val < 0.001) "Log-rank p < 0.001" else
              sprintf("Log-rank p = %.3f", p_val)
hr_label <- sprintf(
  "HR (High vs Low) = %.2f\n(95%% CI %.2f-%.2f)\n%s",
  hr, ci_lo, ci_hi, p_label
)

x_max  <- max(df_bl$time, na.rm = TRUE)
breaks <- if (x_max < 30) 5 else 25

splot <- ggsurvplot(
  sfit,
  data              = df_bl,
  risk.table        = TRUE,
  xlab              = "Time (months)",
  palette           = c("violetred2", "turquoise4"),
  break.time.by     = breaks,
  legend.labs       = c("Low", "High"),
  risk.table.y.text = TRUE,
  fontsize          = 2.5,
  censor.size       = 2,
  font.legend       = 8,
  font.tickslab     = 8,
  font.x            = 10,
  font.y            = 10,
  tables.theme      = theme_classic(base_size = 10),
  title = sprintf(
    "PDAC -> Bladder transfer  (n = %d, %d events)",
    nrow(df_bl), sum(df_bl$event)
  )
)

splot$plot <- splot$plot +
  annotate(
    "text",
    x     = x_max * 0.98,
    y     = 0.85,
    hjust = 1,
    size  = 2.5,
    label = hr_label
  )

# ── 11. Save figure ──────────────────────────────────────────────────────────
fig <- plot_grid(splot$plot, splot$table,
                 nrow = 2, rel_heights = c(3, 1),
                 align = "v", axis = "lr")
ggsave(OUT, fig, width = 6, height = 5)
cat(sprintf("\nFigure saved to: %s\n", OUT))
