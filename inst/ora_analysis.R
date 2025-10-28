tar_load(tops_best_50_1000)
tops=tops_best_50_1000$top_genes
organism = "org.Hs.eg.db"  #database what to use
tar_load(data_filtered_1000)
universe = rownames(data_filtered_1000$ex)
library(organism, character.only = TRUE)
library(clusterProfiler)
library(xlsx)
library(ggplot2)
library(enrichplot)

i=2
gene_list_enrich = tops[,i]
ego <- enrichGO(gene = gene_list_enrich,
                OrgDb = organism,
                keyType = 'SYMBOL', 
                ont = "ALL",
                universe=universe,
                pAdjustMethod = "BH", 
                pvalueCutoff  = 0.05, 
                qvalueCutoff  = 0.2, 
                readable      = FALSE) 
#setwd(dir_data_intersection)
saveRDS(ego, file = paste0("paper/figures/ORA_GO_",i,".Rds"))
write.xlsx(ego, paste0("paper/figures/ORA_analysis_",i,".xls"),
           sheetName = "ORA_GO", append = TRUE)

# 3. ORA_KEGG-------------------------
kegg_organism = "hsa"  #database what to use
gene_list_enrichKEGG <- bitr(gene_list_enrich, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
# BiocManager::install(kegg_organism, character.only = TRUE)
enrich_KEGG <- enrichKEGG(gene = gene_list_enrichKEGG$ENTREZID,
                          organism = kegg_organism,
                          minGSSize = 3,
                          maxGSSize = 800,
                          pvalueCutoff = 0.05,
                          keyType = "kegg")
enrich_KEGG <- setReadable(enrich_KEGG, 'org.Hs.eg.db', 'ENTREZID')%>%
  filter(p.adjust < 0.05,qvalue<0.2)
#setwd(dir_data_intersection)
saveRDS(enrich_KEGG, file = paste0("paper/figures/ORA_KEGG_",i,".Rds"))
write.xlsx(enrich_KEGG, paste0("paper/figures/ORA_analysis_",i,".xls"),
           sheetName = "ORA_KEGG", append = TRUE)


## 4 plot------------------------------------------------------------------
i=3
file=paste0("paper/figures/ORA_GO_FIGURES_",i,".pdf")
ORA_GO <- readRDS(paste0("paper/figures/ORA_GO_",i,".Rds"))
ORA_KEGG <- readRDS(paste0("paper/figures/ORA_KEGG_",i,".Rds"))

if(dim(ORA_GO@result)[1]>0){
  p1 <- dotplot(ORA_GO, showCategory=10) +
    ggtitle(paste0("ORA Analysis Factor ",i))+
    theme(plot.title = element_text(hjust = 0.5,size = 20))
  
  p2 <- emapplot(pairwise_termsim(ORA_GO), showCategory=40)+
    ggtitle(paste0("ORA_GO_erichment map (Top40)"))+
    theme(plot.title = element_text(hjust = 0.5,size = 20))
  
  if(dim(ORA_GO@result)[1]>1){
    p3 <- treeplot(pairwise_termsim(ORA_GO),showCategory=40,nCluster=2)+
      ggtitle(paste0("ORA_GO_treeplot (Top40)"))+
      theme(plot.title = element_text(hjust = 0.5,size = 20))
    plot_list = list(p1,p2,p3)
  }else{
    plot_list = list(p1,p2)
  }
  
  pdf(file)
  for (j in 1:length(plot_list)) {
    print(plot_list[[j]])
  }
  dev.off()
}

ggsave("paper/figures/ORA_factor3.png",p1,height=7,width=6)


file=paste0("ORA_KEGG_FIGURES_",i,".pdf")

if(dim(ORA_KEGG@result)[1]>0){
  p1 <- dotplot(ORA_KEGG, showCategory=20) +
    ggtitle(paste0("ORA_KEGG_dotplot (Top20)"))+
    theme(plot.title = element_text(hjust = 0.5,size = 20))
  p2 <- emapplot(pairwise_termsim(ORA_KEGG), showCategory=40)+
    ggtitle(paste0("ORA_KEGG_erichment map (Top40)"))+
    theme(plot.title = element_text(hjust = 0.5,size = 20))
  if(dim(ORA_KEGG@result)[1]>1){
    p3 <- treeplot(pairwise_termsim(ORA_KEGG),showCategory=40,nCluster=2)+
      ggtitle(paste0("ORA_KEGG_treeplot (Top40)"))+
      theme(plot.title = element_text(hjust = 0.5,size = 20))
    plot_list=list(p1,p2,p3)
  }else{
    plot_list=list(p1,p2)
  }
  
  pdf(file)
  for (j in 1:length(plot_list)) {
    print(plot_list[[j]])
  }
  dev.off()
}
