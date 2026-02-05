#' Overrepresentation analysis (ORA) on top genes
#' @param results A results list from select_model() and get_top_genes()
#' @param gene_lists A named list of gene sets
#' @param dplot Logical, whether to produce dotplots
#' @param save Logical, whether to save dotplots to PNG
#' @return results with $labels and $ora_res_list added
ORA <- function(results, gene_lists=NULL, plot = TRUE, p.adj=.05, msigdbr=FALSE) {
  tops <- results$tops
  msigdbr = is.null(gene_lists) | msigdbr
  
  if(msigdbr){
    list_name="ora_res_list_dbr"
    table_name="labels_dbr"
  }else{
    list_name="ora_res_list"
    table_name="labels"
  }
  if(!is.null(results[[table_name]]) & !is.null(results[[list_name]])){
    if(plot){
      ORA_plots(results,list_name,msigdbr,p.adj)
    }
    
    return(results)
  }
  
  # Prepare TERM2GENE data frame
  if(!msigdbr){
    gene_lists_flat <- stack(unlist(gene_lists, recursive = FALSE))
    colnames(gene_lists_flat) <- c("gene", "term")
    TERM2GENE <- gene_lists_flat[, c("term", "gene")]
  }else{
    library(msigdbr)
    msigdbr=TRUE
    
    # 1. Get MSigDB gene sets, for example Hallmark (H) gene sets for human
    cancer_sets <- msigdbr(species = "Homo sapiens", collection="C4")
    onco_sets = msigdbr(species = "Homo sapiens", collection="C6")
    immune_sets = msigdbr(species = "Homo sapiens", collection="C7",subcollection = "IMMUNESIGDB")
    hallmark_sets = msigdbr(species = "Homo sapiens", collection="H")
    c2_sets =  msigdbr(species = "Homo sapiens", collection="C2")
    c5_sets =  msigdbr(species = "Homo sapiens", collection="C5")
    c8_sets =  msigdbr(species = "Homo sapiens", collection="C8")
    # 2. Prepare TERM2GENE from msigdbr output
    TERM2GENE <- rbind(
      cancer_sets[, c("gs_name", "gene_symbol")],
      onco_sets[, c("gs_name", "gene_symbol")],
      immune_sets[, c("gs_name", "gene_symbol")],
      hallmark_sets[, c("gs_name", "gene_symbol")],
      c2_sets[, c("gs_name", "gene_symbol")],
      c5_sets[, c("gs_name", "gene_symbol")],
      c8_sets[, c("gs_name", "gene_symbol")]
    )
    
  }
  
  
  # Initialize outputs
  label_table <- data.frame(
    factor = seq_len(ncol(tops)),
    label = NA_character_,
    p.adjust = NA_real_,
    stringsAsFactors = FALSE
  )
  
  ora_res_list <- vector("list", length = ncol(tops))
  
  for (i in seq_len(ncol(tops))) {
    gene_set <- tops[[i]]
    
    
    ### change cutoff value here if type=tumor_strict
    ora_res <- clusterProfiler::enricher(
      gene = gene_set,
      TERM2GENE = TERM2GENE,
      minGSSize = 2,
      pvalueCutoff = p.adj
    )
    
    ora_res_list[[i]] <- ora_res
    ora_df = as.data.frame(ora_res)
    
    if (!is.null(ora_res) && nrow(ora_df) > 0) {
      label_table$label[i] <- ora_df$ID[1]
      label_table$p.adjust[i] <- ora_df$p.adjust[1]
    }
  }
  

  # Store into results
  results[[table_name]] <- label_table
  results[[list_name]] = ora_res_list

  return(results)
}
