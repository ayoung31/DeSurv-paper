# Detect + clean multi-mapped symbols in an expression matrix
# - expr: numeric matrix/data.frame with rownames as symbols (possibly multi-mapped/concatenated)
# - gene_info (optional): data.frame with columns `symbol` and `biotype` (e.g., "protein_coding","lncRNA","snoRNA",...)
# - duplicate_strategy: how to handle duplicate symbols after cleaning
#     "highest_variance" (default), "highest_iqr", "mean", or "sum"
clean_multimapped_symbols <- function(
    expr,
    gene_info = NULL,
    duplicate_strategy = c("highest_variance", "highest_iqr", "mean", "sum"),
    verbose = TRUE
){
  stopifnot(!is.null(rownames(expr)))
  duplicate_strategy <- match.arg(duplicate_strategy)
  
  # 1) Delimiter pattern: common separators seen in GEO/array annotations
  delim_re <- "(?:\\.{5}|///|\\||;|,|\\s+/\\s+)"  # '.....', '///', '|', ';', ',', ' / '
  splitter  <- function(x) {
    parts <- unlist(strsplit(x, delim_re, perl = TRUE))
    parts <- trimws(parts)
    parts[nzchar(parts)]
  }
  
  # 2) Make a fast biotype lookup if provided
  biotype_lookup <- NULL
  if (!is.null(gene_info)) {
    stopifnot(all(c("symbol","biotype") %in% names(gene_info)))
    biotype_lookup <- setNames(as.character(gene_info$biotype), as.character(gene_info$symbol))
  }
  
  # 3) Scoring rules to pick the "best" single symbol from a set
  score_symbol <- function(sym){
    score <- 0
    # Prefer HGNC-looking, clean uppercase alnum with optional - or .
    if (grepl("^[A-Z0-9]+([-.][A-Z0-9]+)*$", sym)) score <- score + 1 else score <- score - 1
    # Penalize likely placeholders/predicted/LOC ids and scaffold-like
    if (grepl("^(LOC|LINC|AC[0-9]|AL[0-9]|AP[0-9]|RP[0-9])", sym)) score <- score - 2
    # Penalize common pseudogene suffixes
    if (grepl("P[0-9]*$", sym)) score <- score - 2
    # Lightly penalize small RNAs if no biotype info available
    if (is.null(biotype_lookup)) {
      if (grepl("^SNOR|^MIR|^MT-", sym)) score <- score - 1
    } else {
      bt <- biotype_lookup[[sym]]
      if (!is.null(bt)) {
        # Biotype preferences (tune as desired)
        score <- score + switch(bt,
                                "protein_coding" = 5,
                                "lncRNA"         = 2,
                                "antisense"      = 1,
                                "processed_transcript" = 1,
                                "miRNA"          = -1,
                                "snoRNA"         = -2,
                                "snRNA"          = -2,
                                "rRNA"           = -2,
                                "misc_RNA"       = -1,
                                0
        )
      }
    }
    # Tie-breakers: prefer longer symbols slightly (helps pick canonical gene over tiny RNA tags)
    score <- score + 0.01 * nchar(sym)
    return(score)
  }
  
  # 4) For each rowname, choose the best single symbol
  rn <- rownames(expr)
  parts_list <- lapply(rn, splitter)
  is_multi   <- vapply(parts_list, function(p) length(p) > 1, logical(1))
  chosen     <- character(length(rn))
  
  for (i in seq_along(rn)) {
    syms <- parts_list[[i]]
    if (length(syms) == 0L) {
      chosen[i] <- rn[i]  # fallback
    } else if (length(syms) == 1L) {
      chosen[i] <- syms
    } else {
      sc <- vapply(syms, score_symbol, numeric(1))
      # pick highest score; if tie, first in list
      chosen[i] <- syms[which.max(sc)]
    }
  }
  
  if (verbose) {
    n_multi <- sum(is_multi)
    message(sprintf("Detected %d multi-mapped rownames out of %d.", n_multi, length(rn)))
    if (!is.null(biotype_lookup)) {
      hit <- sum(chosen %in% names(biotype_lookup))
      message(sprintf("Biotype available for %d/%d chosen symbols.", hit, length(chosen)))
    }
  }
  
  # 5) Apply chosen symbols, then resolve duplicates
  expr2 <- as.matrix(expr)
  rownames(expr2) <- chosen
  
  # Duplicate handling
  if (anyDuplicated(rownames(expr2))) {
    if (verbose) {
      message("Resolving duplicate symbols with strategy: ", duplicate_strategy)
    }
    # split by symbol
    sp <- split(seq_len(nrow(expr2)), rownames(expr2))
    keep_idx <- integer(0)
    
    if (duplicate_strategy %in% c("highest_variance","highest_iqr")) {
      # choose the single row with max variance or IQR
      for (sym in names(sp)) {
        idx <- sp[[sym]]
        if (length(idx) == 1L) {
          keep_idx <- c(keep_idx, idx)
        } else {
          sub <- expr2[idx, , drop = FALSE]
          if (duplicate_strategy == "highest_variance") {
            s <- apply(sub, 1, var, na.rm = TRUE)
          } else {
            s <- apply(sub, 1, IQR, na.rm = TRUE)
          }
          keep_idx <- c(keep_idx, idx[which.max(s)])
        }
      }
      expr2 <- expr2[keep_idx, , drop = FALSE]
    } else {
      # "mean" or "sum" → aggregate duplicates row-wise
      agg_fun <- if (duplicate_strategy == "mean") colMeans else function(m) colSums(m, na.rm = TRUE)
      out <- lapply(sp, function(idx){
        if (length(idx) == 1L) expr2[idx, , drop = FALSE]
        else {
          m <- expr2[idx, , drop = FALSE]
          t(as.matrix(agg_fun(m)))
        }
      })
      expr2 <- do.call(rbind, out)
      expr2 <- expr2[order(rownames(expr2)), , drop = FALSE]
    }
  }
  
  return(expr2)
}
