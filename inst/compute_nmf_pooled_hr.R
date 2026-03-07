#!/usr/bin/env Rscript
# compute_nmf_pooled_hr.R
#
# Computes the pooled continuous HR for standard NMF (at DeSurv's k)
# in the same way as the DeSurv pooled HR in 04_results_REVISED.Rmd.
#
# Requires: targets store with val_latent_std_desurvk_tcgacptac
#
# Usage:
#   Rscript inst/compute_nmf_pooled_hr.R [store_path]
#
# Output: results/nmf_pooled_hr.rds (list with hr, ci_lo, ci_hi, p, n, n_events, n_cohorts)

library(targets)
library(survival)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  Sys.setenv(TAR_STORE = args[1])
}

tar_load(val_latent_std_desurvk_tcgacptac)

pool_dfs <- lapply(val_latent_std_desurvk_tcgacptac, function(entry) {
  surv <- entry$survival
  if (is.null(surv) || !all(c("time", "event") %in% names(surv))) return(NULL)
  data.frame(
    time       = surv$time,
    event      = surv$event,
    dataset    = entry$dataset,
    risk_score = entry$risk_score
  )
})

pool_df <- do.call(rbind, Filter(Negate(is.null), pool_dfs))
pool_df <- pool_df[is.finite(pool_df$time) & is.finite(pool_df$risk_score), ]

cox_fit <- coxph(Surv(time, event) ~ risk_score + strata(dataset), data = pool_df)
cox_sum <- summary(cox_fit)

result <- list(
  hr        = round(cox_sum$coefficients[, "exp(coef)"], 2),
  ci_lo     = round(cox_sum$conf.int[, "lower .95"], 2),
  ci_hi     = round(cox_sum$conf.int[, "upper .95"], 2),
  p         = signif(cox_sum$coefficients[, "Pr(>|z|)"], 2),
  n         = nrow(pool_df),
  n_events  = sum(pool_df$event, na.rm = TRUE),
  n_cohorts = length(unique(pool_df$dataset))
)

dir.create("results", showWarnings = FALSE)
saveRDS(result, "results/nmf_pooled_hr.rds")

cat("NMF pooled continuous HR results:\n")
cat(sprintf("  HR = %.2f (95%% CI %.2f--%.2f), P = %s\n",
            result$hr, result$ci_lo, result$ci_hi, result$p))
cat(sprintf("  n = %d, events = %d, cohorts = %d\n",
            result$n, result$n_events, result$n_cohorts))
