#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survival)
})

resolve_repo_root <- function() {
  if (dir.exists("R/simulation_functions")) {
    return(getwd())
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  path <- args[startsWith(args, file_arg)]
  if (length(path)) {
    script_path <- sub(file_arg, "", path[[1]])
    root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
    if (dir.exists(file.path(root, "R/simulation_functions"))) {
      return(root)
    }
  }
  stop("Cannot locate repo root; run from the project root or via Rscript inst/...", call. = FALSE)
}

repo_root <- resolve_repo_root()
setwd(repo_root)

sim_files <- list.files(
  "R/simulation_functions",
  pattern = "[.]R$",
  full.names = TRUE
)
for (file in sim_files) {
  source(file)
}

set.seed(101)
sim <- simulate_desurv_scenario(
  scenario = "R_mixed",
  seed = 101,
  survival_gene_n = 150,
  survival_marker_frac = 0.7 # 70% markers, 30% background
)

expr_vals <- as.vector(sim$X)
if (length(expr_vals) > 200000) {
  expr_vals <- sample(expr_vals, size = 200000, replace = FALSE)
}

if (interactive()) {
  hist(expr_vals,
       breaks = 60,
       main = "R_mixed expression distribution",
       xlab = "Expression",
       col = "#4C78A8",
       border = "white")
}

surv_obj <- survival::Surv(sim$time, sim$status)
km_fit <- survival::survfit(surv_obj ~ 1)

if (interactive()) {
  plot(km_fit,
       xlab = "Time",
       ylab = "Survival probability",
       main = "R_mixed Kaplan-Meier curve",
       col = "#F58518",
       lwd = 2)
  grid(col = "gray85")
}

gene_names <- rownames(sim$W)
marker_genes <- unique(unlist(sim$marker_sets))
background_genes <- sim$background
noise_genes <- sim$noise_genes

gene_type <- rep("other", length(gene_names))
names(gene_type) <- gene_names
gene_type[marker_genes] <- "marker"
gene_type[background_genes] <- "background"
gene_type[noise_genes] <- "noise"

mean_loading <- rowMeans(sim$W)
type_levels <- c("marker", "background", "noise", "other")
gene_type <- factor(gene_type, levels = type_levels)

if (interactive()) {
  boxplot(mean_loading ~ gene_type,
          col = c("#54A24B", "#E45756", "#72B7B2", "#B279A2"),
          main = "Mean W loadings by gene type",
          xlab = "Gene type",
          ylab = "Mean loading")
}

if (!is.null(sim$survival_gene_sets) && length(sim$survival_gene_sets)) {
  survival_sets <- sim$survival_gene_sets
  counts <- vapply(
    survival_sets,
    function(genes) {
      genes <- genes[!is.na(genes)]
      type_counts <- table(factor(gene_type[genes], levels = type_levels))
      as.numeric(type_counts)
    },
    numeric(length(type_levels))
  )
  rownames(counts) <- type_levels

  if (interactive()) {
    barplot(counts,
            beside = FALSE,
            col = c("#54A24B", "#E45756", "#72B7B2", "#B279A2"),
            legend.text = TRUE,
            args.legend = list(x = "topright", bty = "n"),
            main = "Survival gene composition by factor",
            xlab = "Factor",
            ylab = "Gene count")
  }

  resolve_gene_indices <- function(genes, W) {
    if (is.null(genes) || !length(genes)) {
      return(integer(0))
    }
    if (is.character(genes)) {
      idx <- match(genes, rownames(W))
    } else {
      idx <- genes
    }
    idx <- idx[!is.na(idx) & idx >= 1 & idx <= nrow(W)]
    unique(as.integer(idx))
  }

  make_wtilde <- function(W, gene_sets) {
    Wtilde <- matrix(0, nrow = nrow(W), ncol = ncol(W), dimnames = dimnames(W))
    for (k in seq_len(ncol(W))) {
      idx <- resolve_gene_indices(gene_sets[[k]], W)
      if (length(idx)) {
        Wtilde[idx, k] <- W[idx, k]
      }
    }
    Wtilde
  }

  score_from_wtilde <- function(X, Wtilde, scale_scores = TRUE) {
    scores <- crossprod(X, Wtilde)
    if (scale_scores) {
      m <- colMeans(scores)
      s <- apply(scores, 2, sd)
      s[s == 0] <- 1
      scores <- sweep(scores, 2, m, "-")
      scores <- sweep(scores, 2, s, "/")
    }
    scores
  }

  compute_cindex <- function(risk, time, status) {
    surv_obj <- survival::Surv(time, status)
    cc <- tryCatch(
      survival::concordance(surv_obj ~ risk),
      error = function(e) survival::survConcordance(surv_obj ~ risk)
    )
    as.numeric(cc$concordance)
  }

  marker_surv_sets <- lapply(seq_along(survival_sets), function(k) {
    intersect(survival_sets[[k]], sim$marker_sets[[k]])
  })
  background_surv_sets <- lapply(seq_along(survival_sets), function(k) {
    intersect(survival_sets[[k]], sim$background)
  })

  Wtilde_marker <- make_wtilde(sim$W, marker_surv_sets)
  Wtilde_background <- make_wtilde(sim$W, background_surv_sets)
  Wtilde_all <- make_wtilde(sim$W, survival_sets)

  scores_marker <- score_from_wtilde(sim$X, Wtilde_marker)
  scores_background <- score_from_wtilde(sim$X, Wtilde_background)
  scores_all <- score_from_wtilde(sim$X, Wtilde_all)

  # Use negative linear predictor so higher values imply longer survival.
  risk_marker <- -drop(scores_marker %*% sim$beta)
  risk_background <- -drop(scores_background %*% sim$beta)
  risk_all <- -drop(scores_all %*% sim$beta)

  cindex_marker <- compute_cindex(risk_marker, sim$time, sim$status)
  cindex_background <- compute_cindex(risk_background, sim$time, sim$status)
  cindex_all <- compute_cindex(risk_all, sim$time, sim$status)

  cindex_vals <- c(
    marker_only = cindex_marker,
    background_only = cindex_background,
    combined = cindex_all
  )
  print(cindex_vals)

  if (interactive()) {
    barplot(cindex_vals,
            col = c("#59A14F", "#E15759", "#4E79A7"),
            ylim = c(0, 1),
            main = "Survival impact by gene set",
            ylab = "Concordance (c-index)")
    abline(h = 0.5, lty = 2, col = "gray60")
  }
}
