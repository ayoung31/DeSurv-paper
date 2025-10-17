library(tidyestimate)
library(HGNChelper)

tcga=load_data("TCGA_PAAD")
X = tcga$ex[rownames(tcga$ex) != "?",]
X = as.data.frame(X)
genes=rownames(X)
# X = gene_filter(tcga$ex,ngene=5000)
# scores=estimateScore(X,is_affymetrix=FALSE)  # or "affymetrix"
hgnc_fix <- HGNChelper::checkGeneSymbols(genes)
genes_fixed <- hgnc_fix$Suggested.Symbol
genes_fixed[is.na(genes_fixed)] <- hgnc_fix$x[is.na(genes_fixed)]  # fallback to original if no suggestion
expr=X
expr$symbol    = rownames(expr)
expr          = dplyr::as_tibble(expr) |>
  dplyr::group_by(.data$symbol) |>
  dplyr::summarise(dplyr::across(dplyr::everything(), median), .groups = "drop")
expr          = as.data.frame(expr[!is.na(expr$symbol), ])
rownames(expr) = expr$symbol
expr$symbol    = NULL
expr=expr[,tcga$samp_keeps]
scores= expr |>
  filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE, find_alias = TRUE) |> 
  estimate_score(is_affymetrix=TRUE) 
