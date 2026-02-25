##############################################################################
# cv_grid_validation_analysis_allgenes.R
#
# Training + external validation analysis of cv_grid DeSurv fits over
# K=2,3,5,7,9 x alpha=0,0.25,0.35,0.55,0.75,0.85,0.95 (35 combos).
# Fixed parameters: ntop=NULL (all genes), lambda=0.349, nu=0.056,
# 100 initializations.
#
# Run from the DeSurv-paper repo root:
#   Rscript inst/cv_grid_validation_analysis_allgenes.R
##############################################################################

suppressMessages({
  library(survival)
  library(data.table)
})

STORE <- "store_PKG_VERSION=20260107bugfix_GIT_BRANCH=main"
TARGET_NTOP   <- NULL   # all genes
TARGET_LAMBDA <- 0.349
TARGET_NU     <- 0.056

load("data/derv/cmbSubtypes_formatted.RData")

orig_bundle <- readRDS(file.path(STORE, "objects/val_run_bundle_tcgacptac"))
orig_fit <- orig_bundle$fit_desurv
orig_W <- orig_fit$W

bo_bundle <- readRDS(file.path(STORE, "objects/bo_bundle_selected_tcgacptac"))
X_train <- bo_bundle$data_filtered$ex
surv_train <- bo_bundle$data_filtered$sampInfo

val_data <- readRDS(file.path(STORE, "objects/data_val_filtered_surv_tcgacptac"))
cat("Validation datasets:", paste(names(val_data), collapse=", "), "\n")
total_val <- sum(sapply(val_data, function(x) ncol(x$ex)))
cat("Total validation n:", total_val, "\n\n")

tcga_sub  <- read.csv("data/original/TCGA_PAAD_subtype.csv")
cptac_sub <- read.csv("data/original/CPTAC_subtype.csv")
all_train_sub <- rbind(
  tcga_sub[,  c("ID", "PurIST", "DeCAF")],
  cptac_sub[, c("ID", "PurIST", "DeCAF")]
)
all_train_sub$DeCAF[all_train_sub$DeCAF == "permCAF"] <- "proCAF"

source("R/get_top_genes.R")
orig_tops <- get_top_genes(orig_W, 270)  # original production factor top genes

elyada_icaf     <- top_genes$Elyada_CAF$iCAF
decoder_immune  <- top_genes$DECODER$Immune
decoder_basaltumor <- top_genes$DECODER$BasalTumor

cat("Loading all cv_grid_fit branches from store...\n")
fit_files <- list.files(
  file.path(STORE, "objects"),
  pattern = "^cv_grid_fit_[0-9a-f]+$",
  full.names = TRUE
)
all_fits <- lapply(fit_files, readRDS)
cat(sprintf("  Loaded %d branches.\n", length(all_fits)))

find_fit_null_ntop <- function(all_fits, k, alpha, lambda, nu) {
  for (fe in all_fits) {
    if (is.null(fe) || is.null(fe$k)) next
    fe_ntop   <- if (is.null(fe$ntop))   NA_real_ else as.numeric(fe$ntop)
    fe_lambda <- if (is.null(fe$lambda)) NA_real_ else as.numeric(fe$lambda)
    fe_nu     <- if (is.null(fe$nu))     NA_real_ else as.numeric(fe$nu)
    ntop_ok   <- is.na(fe_ntop)   # ntop=NULL stored as NA
    lam_ok    <- !is.na(fe_lambda) && abs(fe_lambda - lambda) < 1e-9
    nu_ok     <- !is.na(fe_nu)     && abs(fe_nu  - nu)     < 1e-9
    if (fe$k == k && abs(fe$alpha - alpha) < 1e-6 && ntop_ok && lam_ok && nu_ok) {
      return(fe)
    }
  }
  stop(sprintf("No fit found for k=%d, alpha=%.2f, ntop=NULL, lambda=%.3f, nu=%.3f",
               k, alpha, lambda, nu))
}

GRID_K     <- c(2, 3, 5, 7, 9)
GRID_ALPHA <- c(0, 0.25, 0.35, 0.55, 0.75, 0.85, 0.95)

fit_specs <- list()
for (ki in GRID_K) {
  for (ai in GRID_ALPHA) {
    nm <- sprintf("K%d_a%s", ki, sub("\\.", "", sprintf("%.2f", ai)))
    fit_specs[[nm]] <- list(k = ki, alpha = ai)
  }
}

fits <- lapply(fit_specs, function(spec) {
  find_fit_null_ntop(all_fits, k = spec$k, alpha = spec$alpha,
                     lambda = TARGET_LAMBDA, nu = TARGET_NU)
})

for (nm in names(fits)) {
  exp <- fit_specs[[nm]]
  got <- fits[[nm]]
  if (got$k != exp$k || abs(got$alpha - exp$alpha) > 1e-6) {
    stop(sprintf("Mismatch for %s", nm))
  }
  cat(sprintf("  %s: k=%d alpha=%.2f ntop=%s lambda=%.3f nu=%.3f OK\n",
    nm, got$k, got$alpha,
    ifelse(is.null(got$ntop), "NULL", got$ntop),
    ifelse(is.null(got$lambda), NA, got$lambda),
    ifelse(is.null(got$nu), NA, got$nu)))
}

rm(all_fits); gc(verbose = FALSE)
cat(sprintf("All %d fits loaded.\n\n", length(fits)))

# ---- helpers ----
build_val_df <- function(W, val_data, ntop = NULL) {
  if (!is.null(ntop)) {
    tops <- get_top_genes(W, ntop)
    W_focused <- W
    for (j in seq_len(ncol(W))) {
      keep_genes <- tops$top_genes[, j]
      non_top <- setdiff(rownames(W), keep_genes)
      W_focused[non_top, j] <- 0
    }
  }
  pooled <- list()
  for (ds_name in names(val_data)) {
    vd <- val_data[[ds_name]]
    X_val <- vd$ex
    si <- vd$sampInfo
    if (!is.null(ntop)) {
      genes <- intersect(rownames(W_focused), rownames(X_val))
      H_val <- t(X_val[genes, ]) %*% W_focused[genes, ]
    } else {
      genes <- intersect(rownames(W), rownames(X_val))
      H_val <- t(X_val[genes, ]) %*% W[genes, ]
    }
    colnames(H_val) <- paste0("Factor", seq_len(ncol(H_val)))
    decaf <- as.character(si$DeCAF)
    decaf[decaf == "permCAF"] <- "proCAF"
    ds_df <- data.frame(
      sample_id = colnames(X_val), time = si$time, event = si$event,
      PurIST = as.character(si$PurIST), DeCAF = decaf,
      dataset = ds_name, H_val, stringsAsFactors = FALSE
    )
    keep <- complete.cases(ds_df[, c("time", "event")])
    pooled[[ds_name]] <- ds_df[keep, ]
  }
  pooled_df <- do.call(rbind, pooled)
  rownames(pooled_df) <- NULL
  k <- sum(grepl("^Factor", names(pooled_df)))
  for (j in seq_len(k)) {
    fn <- paste0("Factor", j)
    pooled_df[[fn]] <- as.numeric(scale(pooled_df[[fn]]))
  }
  pooled_df
}

build_train_df <- function(W, X_train, surv_train, subtypes, ntop = NULL) {
  if (!is.null(ntop)) {
    tops <- get_top_genes(W, ntop)
    W_focused <- W
    for (j in seq_len(ncol(W))) {
      keep_genes <- tops$top_genes[, j]
      non_top <- setdiff(rownames(W), keep_genes)
      W_focused[non_top, j] <- 0
    }
    genes <- intersect(rownames(W_focused), rownames(X_train))
    H_train <- t(X_train[genes, ]) %*% W_focused[genes, ]
  } else {
    genes <- intersect(rownames(W), rownames(X_train))
    H_train <- t(X_train[genes, ]) %*% W[genes, ]
  }
  colnames(H_train) <- paste0("Factor", seq_len(ncol(H_train)))
  train_ids <- colnames(X_train)
  idx <- match(train_ids, subtypes$ID)
  train_df <- data.frame(
    sample_id = train_ids, time = surv_train$time, event = surv_train$event,
    PurIST = as.character(subtypes$PurIST[idx]),
    DeCAF = as.character(subtypes$DeCAF[idx]),
    H_train, stringsAsFactors = FALSE
  )
  k <- sum(grepl("^Factor", names(train_df)))
  for (j in seq_len(k)) {
    fn <- paste0("Factor", j)
    train_df[[fn]] <- as.numeric(scale(train_df[[fn]]))
  }
  train_df
}

# ---- storage matrices ----
hcor_matrix      <- matrix(NA_real_, nrow=length(GRID_K), ncol=length(GRID_ALPHA),
                            dimnames=list(paste0("K=",GRID_K), paste0("a=",sprintf("%.2f",GRID_ALPHA))))
val_p_matrix     <- hcor_matrix
val_p_adj_matrix <- hcor_matrix

master_rows <- list()

for (nm in names(fits)) {
  fit   <- fits[[nm]]
  W     <- fit$fit$W
  k     <- fit$k
  alpha <- fit$alpha
  ri <- match(k, GRID_K)
  ci <- match(round(alpha,2), round(GRID_ALPHA,2))

  cat(sprintf("\nProcessing %s (K=%d, alpha=%.2f)...\n", nm, k, alpha))

  # H-score correlation: use ntop=270 focused projection for ORIGINAL factor,
  # but use ALL genes (ntop=NULL) for the grid fit's W matrix
  # This measures how well the all-genes projection recovers the original factor
  W_orig_f <- orig_W
  for (j in seq_len(ncol(orig_W))) {
    keep_genes <- orig_tops$top_genes[, j]
    non_top <- setdiff(rownames(orig_W), keep_genes)
    W_orig_f[non_top, j] <- 0
  }
  genes_t <- Reduce(intersect, list(rownames(W), rownames(W_orig_f), rownames(X_train)))
  H_orig_f <- t(X_train[genes_t, ]) %*% W_orig_f[genes_t, ]
  # All-genes projection for grid fit
  H_new_all <- t(X_train[genes_t, ]) %*% W[genes_t, ]
  cors <- cor(H_orig_f[,1], H_new_all, method="spearman")
  best_f <- which.max(cors)
  hcor_val <- max(cors)
  hcor_matrix[ri, ci] <- hcor_val

  # Gene overlap: top 270 genes of best factor vs original top 270
  tops_new <- get_top_genes(W, 270)
  gene_ov <- length(intersect(tops_new$top_genes[, best_f], orig_tops$top_genes[, 1]))

  elyada_u  <- intersect(elyada_icaf, rownames(W))
  decoder_u <- intersect(decoder_immune, rownames(W))
  icaf_genes <- tops_new$top_genes[, best_f]
  elyada_ov  <- length(intersect(icaf_genes, elyada_u))
  decoder_ov <- length(intersect(icaf_genes, decoder_u))

  icaf_fn <- paste0("Factor", best_f)

  # Build train/val DFs using ALL genes (ntop=NULL)
  train_df <- build_train_df(W, X_train, surv_train, all_train_sub, ntop=NULL)
  val_df   <- build_val_df(W, val_data, ntop=NULL)

  # Training Cox unadjusted
  cx_t <- coxph(as.formula(paste0("Surv(time,event)~", icaf_fn)), data=train_df)
  s_t  <- summary(cx_t)
  hr_t <- s_t$conf.int[1,1]
  p_t  <- s_t$coefficients[1,"Pr(>|z|)"]

  # Training Cox adjusted
  adj_cols  <- c("time","event","PurIST","DeCAF", paste0("Factor",seq_len(k)))
  train_adj <- train_df[complete.cases(train_df[,adj_cols]),]
  cx_t_adj  <- coxph(as.formula(paste0("Surv(time,event)~PurIST+DeCAF+",icaf_fn)), data=train_adj)
  s_t_adj   <- summary(cx_t_adj)
  row_t     <- grep(paste0("^",icaf_fn,"$"), rownames(s_t_adj$coefficients))
  hr_t_adj  <- exp(s_t_adj$coefficients[row_t,"coef"])
  p_t_adj   <- s_t_adj$coefficients[row_t,"Pr(>|z|)"]

  # Training LRT
  base_t  <- coxph(Surv(time,event)~PurIST+DeCAF, data=train_adj)
  lrt_t   <- anova(base_t, cx_t_adj)
  lrt_p_t <- lrt_t[["Pr(>|Chi|)"]][2]

  # Validation Cox unadjusted
  cx_v <- coxph(as.formula(paste0("Surv(time,event)~", icaf_fn, "+strata(dataset)")), data=val_df)
  s_v  <- summary(cx_v)
  hr_v <- s_v$conf.int[1,1]
  p_v  <- s_v$coefficients[1,"Pr(>|z|)"]
  val_p_matrix[ri, ci] <- p_v

  # Validation Cox adjusted
  val_adj  <- val_df[complete.cases(val_df[,c("time","event","PurIST","DeCAF")]),]
  cx_v_adj <- coxph(as.formula(paste0("Surv(time,event)~PurIST+DeCAF+",icaf_fn,"+strata(dataset)")), data=val_adj)
  s_v_adj  <- summary(cx_v_adj)
  row_v    <- grep(paste0("^",icaf_fn,"$"), rownames(s_v_adj$coefficients))
  hr_v_adj <- exp(s_v_adj$coefficients[row_v,"coef"])
  p_v_adj  <- s_v_adj$coefficients[row_v,"Pr(>|z|)"]
  val_p_adj_matrix[ri, ci] <- p_v_adj

  # Validation LRT
  base_v    <- coxph(Surv(time,event)~PurIST+DeCAF+strata(dataset), data=val_adj)
  cx_v_full <- coxph(as.formula(paste0("Surv(time,event)~PurIST+DeCAF+",icaf_fn,"+strata(dataset)")), data=val_adj)
  lrt_v     <- anova(base_v, cx_v_full)
  lrt_p_v   <- lrt_v[["Pr(>|Chi|)"]][2]

  # KM
  val_df$icaf_group <- ifelse(val_df[[icaf_fn]] > median(val_df[[icaf_fn]]), "High", "Low")
  sf_v  <- survfit(Surv(time,event)~icaf_group, data=val_df)
  med_v <- summary(sf_v)$table[,"median"]
  lr_v  <- survdiff(Surv(time,event)~icaf_group, data=val_df)
  lr_p_v <- 1 - pchisq(lr_v$chisq, 1)
  km_hi <- med_v["icaf_group=High"]
  km_lo <- med_v["icaf_group=Low"]

  master_rows[[nm]] <- list(
    nm=nm, k=k, alpha=alpha, best_f=best_f,
    hcor=hcor_val, gene_ov=gene_ov,
    elyada_ov=elyada_ov, elyada_n=length(elyada_u),
    decoder_ov=decoder_ov, decoder_n=length(decoder_u),
    hr_t=hr_t, p_t=p_t,
    hr_t_adj=hr_t_adj, p_t_adj=p_t_adj, lrt_p_t=lrt_p_t,
    hr_v=hr_v, p_v=p_v,
    hr_v_adj=hr_v_adj, p_v_adj=p_v_adj, lrt_p_v=lrt_p_v,
    km_hi=km_hi, km_lo=km_lo, km_p=lr_p_v
  )
}

cat("\n\n========== H_cor MATRIX (ntop=NULL, lambda=0.349, nu=0.056) ==========\n")
cat(sprintf("%-6s","K"))
for (ai in GRID_ALPHA) cat(sprintf(" %8s", sprintf("a=%.2f",ai)))
cat("\n")
for (i in seq_along(GRID_K)) {
  cat(sprintf("%-6d", GRID_K[i]))
  for (j in seq_along(GRID_ALPHA)) {
    v <- hcor_matrix[i,j]
    cat(sprintf(" %8.3f", v))
  }
  cat("\n")
}

cat("\n========== UNADJUSTED VALIDATION P-VALUE MATRIX ==========\n")
cat(sprintf("%-6s","K"))
for (ai in GRID_ALPHA) cat(sprintf(" %10s", sprintf("a=%.2f",ai)))
cat("\n")
for (i in seq_along(GRID_K)) {
  cat(sprintf("%-6d", GRID_K[i]))
  for (j in seq_along(GRID_ALPHA)) {
    v <- val_p_matrix[i,j]
    cat(sprintf(" %10.4e", v))
  }
  cat("\n")
}

cat("\n========== ADJUSTED VALIDATION P-VALUE MATRIX ==========\n")
cat(sprintf("%-6s","K"))
for (ai in GRID_ALPHA) cat(sprintf(" %10s", sprintf("a=%.2f",ai)))
cat("\n")
for (i in seq_along(GRID_K)) {
  cat(sprintf("%-6d", GRID_K[i]))
  for (j in seq_along(GRID_ALPHA)) {
    v <- val_p_adj_matrix[i,j]
    cat(sprintf(" %10.4e", v))
  }
  cat("\n")
}

cat("\n========== MASTER TABLE ==========\n")
hdr <- sprintf("%-14s %3s %5s %4s | %6s %6s | %6s/%2s %6s/%2s | TrainHR(un) TrainHR(adj) TrainLRT | ValHR(un) ValHR(adj) ValLRT | KM_Hi  KM_Lo  KM_p",
  "Fit","K","Alpha","F","H-cor","GeneOv","Elyada","N","DECODER","N")
cat(hdr, "\n")
cat(paste(rep("-",length(hdr)+20), collapse=""), "\n")

for (nm in names(master_rows)) {
  r <- master_rows[[nm]]
  cat(sprintf(
    "%-14s %3d %5.2f %4d | %6.3f %6d | %6d/%2d %7d/%2d | %11s %12s %8.2e | %9s %10s %6.2e | %5.1f %6.1f %8.4e\n",
    r$nm, r$k, r$alpha, r$best_f,
    r$hcor, r$gene_ov,
    r$elyada_ov, r$elyada_n, r$decoder_ov, r$decoder_n,
    sprintf("%.3f(p=%.4f)", r$hr_t, r$p_t),
    sprintf("%.3f(p=%.4f)", r$hr_t_adj, r$p_t_adj),
    r$lrt_p_t,
    sprintf("%.3f(p=%.4f)", r$hr_v, r$p_v),
    sprintf("%.3f(p=%.3f)", r$hr_v_adj, r$p_v_adj),
    r$lrt_p_v,
    r$km_hi, r$km_lo, r$km_p))
}

cat("\nDone.\n")
