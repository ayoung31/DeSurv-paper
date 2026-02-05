library(tidyestimate)
library(HGNChelper)
library(targets)
library(dplyr)
library(ggplot2)

tar_load(data_val)
tar_load(aligned_clusters_desurv_50_2000)
clus=aligned_clusters_desurv_50_2000

sinfos=list()
for(i in 1:5){
  sinfos[[i]]=clus[[i]]$data$sampInfo
}
sampInfo=bind_rows(sinfos)

X = data_val[[1]]$ex

scores = lapply(1:5,function(i){
  dat=data_val[[i]]
  X=dat$ex
  X = X[rownames(X) != "?",]
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
  temp=clus[[i]]$data$sampInfo %>% left_join(scores)
  temp
  
},data_val,clus)

X = X[rownames(X) != "?",]
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
expr=expr[,data_val[[1]]$samp_keeps]
scores= expr |>
  filter_common_genes(id = "hgnc_symbol", tell_missing = FALSE, find_alias = TRUE) |> 
  estimate_score(is_affymetrix=TRUE) %>%
  mutate(sampID=sample)

temp=clus[[1]]$data$sampInfo %>% left_join(scores)

ggplot(temp,aes(x=samp_cluster,y=stromal, fill=samp_cluster, group=samp_cluster))+
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black", fill = "white") +
  labs(y = "Tumor Purity (ESTIMATE)", x = "Cluster") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

kruskal.test(purity ~ samp_cluster, data = temp)
pairwise.wilcox.test(temp$purity, temp$samp_cluster, p.adjust.method = "BH")
