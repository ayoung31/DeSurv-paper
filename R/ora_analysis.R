
ora = function(tops,universe,organism){
  if(!requireNamespace("clusterProfiler",quietly = TRUE)){
    stop("clusterProfiler must be installed to load this package")
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db must be installed to use ORA analysis")
  }
  universeKEGG <- suppressMessages(
    clusterProfiler::bitr(universe,
                          fromType = "SYMBOL",
                          toType = c("ENTREZID"),
                          OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  )
  universeKEGG <- unique(universeKEGG$ENTREZID)
  ego = list()
  enrich_KEGG = list()
  for(i in 1:ncol(tops)){
    gene_list_enrich = tops[,i]
    
    ego[[i]] <- clusterProfiler::enrichGO(
      gene = gene_list_enrich,
      OrgDb = organism,
      keyType = 'SYMBOL', 
      ont = "ALL",
      universe=universe,
      pAdjustMethod = "BH", 
      pvalueCutoff  = 1, 
      qvalueCutoff  = 1, 
      readable      = FALSE) 
    
    kegg_organism = "hsa"  #database what to use
    gene_list_enrichKEGG <- suppressMessages(
      clusterProfiler::bitr(
        gene_list_enrich,
        fromType = "SYMBOL",
        toType = c("ENTREZID"),
        OrgDb = org.Hs.eg.db::org.Hs.eg.db
      )
    )
    
    ek <- suppressMessages(
      clusterProfiler::enrichKEGG(
        gene = gene_list_enrichKEGG$ENTREZID,
        organism = kegg_organism,
        universe = universeKEGG,
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        keyType = "kegg"
      )
    )
    
    # Filter first, then convert to readable symbols
    # ek <- ek %>% dplyr::filter(p.adjust < 0.05, qvalue < 0.2)
    ek <- clusterProfiler::setReadable(ek, OrgDb = organism, keyType = "ENTREZID")
    
    
    enrich_KEGG[[i]] <- ek
    
  }
  
  return(list(enrich_GO = ego, enrich_KEGG=enrich_KEGG))
}
