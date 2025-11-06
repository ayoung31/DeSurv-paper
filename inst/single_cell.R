library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(targets)
library(VAM)
library(viridis)
library(khroma)
library(dichromat)
tar_load_globals()
tar_load(fit_consensus_2000)
W=fit_consensus_2000$W
tar_load(data_filtered_2000)
X = data_filtered_2000$ex
genes=intersect(rownames(W),rownames(X))
Z=t(X[genes,])%*%W[genes,]
dat=data.frame(scale(Z))
dat$y=data_filtered_2000$sampInfo$time
dat$d=data_filtered_2000$sampInfo$event
coxph(Surv(y,d)~X1,data=dat)

sc = readRDS("data/original/Elyada.Rds")
sc_caf = readRDS("data/original/Elyada_caf.Rds")
sc_tum=readRDS("data/original/Elyada_PDAC.Rds")

sc <- RunUMAP(sc, dims = 1:20) 
sc_caf= RunUMAP(sc_caf,dims=1:20)
sc_tum=RunUMAP(sc_tum,dims=1:20)

sc_save=sc
if ("SCT" %in% names(sc_save@assays)) {
  # remove stored SCTModel and misc stuff
  sc_save@assays$SCT@SCTModel.list <- list()
}
saveRDS(sc_save,file="data/derv/Elyada_umap.Rds",compress = "xz")

sc_save=sc_caf
if ("SCT" %in% names(sc_save@assays)) {
  # remove stored SCTModel and misc stuff
  sc_save@assays$SCT@SCTModel.list <- list()
}
saveRDS(sc_save,file="data/derv/Elyada_caf_umap.Rds",compress = xz)

sc_save=sc_tum
if ("SCT" %in% names(sc_save@assays)) {
  # remove stored SCTModel and misc stuff
  sc_save@assays$SCT@SCTModel.list <- list()
}
saveRDS(sc_save,file="data/derv/Elyada_PDAC_umap.Rds",compress = xz)

sc_all=sc

tar_load(fit_consensus_2000)
W=fit_consensus_2000$W
factor_genes = rownames(W)
sc_genes = rownames(sc)
int_genes = intersect(factor_genes,sc_genes)
W=W[int_genes,]
tops=get_top_genes(W,50)
create_table(tops$top_genes,top_genes,"DECODER",colors)
# tar_load(tops_desurv_50_2000)
# tops=tops_desurv_50_2000


desurv_genesets = as.list(tops$top_genes)
factor_names <- paste0("DeSurv_F", seq_along(desurv_genesets))
names(desurv_genesets) <- factor_names


sc <- AddModuleScore(
  object   = sc,
  features = desurv_genesets,
  name     = "DeSurv"
)

head(sc@meta.data[, grepl("^DeSurv", colnames(sc@meta.data)), drop = FALSE])

features_to_plot <- c("DeSurv1", "DeSurv2", "DeSurv3")

p_list <- lapply(features_to_plot, function(feat) {
  FeaturePlot(
    sc,
    features  = feat,
    reduction = "umap",
    pt.size   = 0.25,
    max.cutoff = "q98"    # trims extreme cells for nicer colors
  ) + ggtitle(feat)
})

## arrange in a row
p_fig2c <- wrap_plots(p_list, ncol = length(p_list))
p_fig2c


#### vam score
DefaultAssay(sc) <- "RNA"   # or "SCT" if that’s what you used

gene_ids <- rownames(sc)    # these should match your gene symbols

gs_collection <- createGeneSetCollection(
  gene.ids            = gene_ids,
  gene.set.collection = desurv_genesets,
  min.size            = 5   # or whatever minimum size you want
  # max.size          = 200 # optional
)

sc <- vamForSeurat(
  seurat.data        = sc,
  gene.set.collection = gs_collection,
  center             = FALSE,        # standard VAM settings
  gamma              = TRUE,
  sample.cov         = FALSE,
  return.dist        = FALSE         # we only want the CDF scores
)
DefaultAssay(sc) <- "VAMcdf"
head(rownames(sc))            # should include "DeSurv1", "DeSurv2", ...
GetAssayData(sc)[1:3, 1:5]    # a few scores

# 1) choose factors to plot
features_to_plot <- c("DeSurv-F1", "DeSurv-F2", "DeSurv-F3")

p_list <- lapply(features_to_plot, function(feat) {
  fp=FeaturePlot(
    sc,
    features   = feat,       # now interpreted as a "feature" in VAM.cdf
    reduction  = "umap",
    pt.size    = 0.25,
    slot = "data",
    max.cutoff = "q95",
  ) + ggtitle(paste0(feat, " (VAM)"))+
    scale_color_gradientn(
      colours = viridis(256,option="D"),
      limits  = c(0, 1)
    )
  fp[[1]]
})

p_fig2c <- wrap_plots(p_list, ncol = length(p_list))
p_fig2c

scplot=DimPlot(sc_all, group.by = "label_broad", reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("All Cell-type clusters")

scplot_caf=DimPlot(sc_caf, group.by = "label", reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("CAF Cell-type clusters")

scplot_tum=DimPlot(sc_tum, group.by = "label_broad", reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("All Cell-type clusters")

temp=list(scplot[[1]],scplot_caf[[1]],scplot_tum[[1]])
plot_grid(plotlist=temp,ncol=3)

temp=list(scplot[[1]],p_list[[1]],p_list[[2]],p_list[[3]])
plot_grid(plotlist=temp,ncol=4)

sc@meta.data %>%
  group_by(label_fine) %>%
  summarise(across(starts_with("DeSurv"), mean))



ct_col <- "label_fine"
if (!ct_col %in% colnames(sc@meta.data)) {
  stop("I expected a column named 'cell_type' in sc@meta.data. Rename or change ct_col.")
}

avg_scores <- sc@meta.data %>%
  group_by(.data[[ct_col]]) %>%
  summarise(
    across(c("DeSurv1","DeSurv2","DeSurv3"), ~ mean(.x, na.rm = TRUE))
  )

mat <- as.matrix(avg_scores[, -1, drop = FALSE])
rownames(mat) <- avg_scores[[ct_col]]

mat_capped=mat
upper <- quantile(mat, 0.95, na.rm = TRUE)
lower <- quantile(mat, 0.01, na.rm = TRUE)
mat_capped[mat_capped > upper] <- upper
mat_capped[mat_capped < lower] <- lower

col_fun <- circlize::colorRamp2(
  c(min(mat_capped), 0, max(mat_capped)),
  c("navy", "white", "firebrick")
)

ht <- Heatmap(
  mat,
  name = "score",
  col  = col_fun,
  cluster_rows    = TRUE,
  cluster_columns = TRUE,
  row_names_side  = "left"
)

ht


top_genes$deCAF = list(proCAF=c("IGFL2","NOX4","VSNL1","BICD1","NPR3","ETV1","ITGA11","CNIH3","COL11A1"),
                       restCAF=c("CHRDL1","OGN","PI16","ANK2","ABCA8","TGFBR3","FBLN5","SCARA5","KIAA1217"))
temp=purrr::list_flatten(top_genes)
ref_sigs=temp[!startsWith(names(temp),"Bailey") & !grepl("peri",names(temp)) & !startsWith(names(temp),"DECODER")]
W=fit_consensus_2000$W
common_genes <- Reduce(intersect, list(rownames(W), unique(unlist(ref_sigs))))
W <- W[common_genes, , drop = FALSE]

# initialize result matrix
cor_mat <- matrix(NA, nrow = ncol(W), ncol = length(ref_sigs),
                  dimnames = list(colnames(W), names(ref_sigs)))

# compute correlations
for (j in seq_len(ncol(W))) {
  wj <- W[, j]
  for (k in seq_along(ref_sigs)) {
    vk <- as.numeric(common_genes %in% ref_sigs[[k]])
    cor_mat[j, k] <- cor(wj, vk, method = "spearman")
  }
}

keep=apply(t(cor_mat),1,function(x) !any(is.na(x)))
mat=t(cor_mat[,which(keep)])
pheatmap::pheatmap(pmin(mat,.2),cluster_cols = FALSE)
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap::pheatmap(mat[,c(1,3)],
                   cluster_cols = FALSE,
                   color = my_colors,
                   breaks = seq(-.4, .4, length.out = 101))



library(fgsea)
j=1
ranks <- sort(W[, j], decreasing = TRUE)
fgseaRes <- fgsea(pathways = ref_sigs, stats = ranks, nperm = 1000)


all_genes <- rownames(W)

# init matrix: factors x signatures
nes_mat <- matrix(NA_real_,
                  nrow = ncol(W),
                  ncol = length(ref_sigs),
                  dimnames = list(colnames(W), names(ref_sigs)))

for (j in seq_len(ncol(W))) {
  # 1) get weights for this factor
  wj <- W[, j]
  
  # 2) make a named numeric vector, sorted desc
  ranks <- sort(wj, decreasing = TRUE)
  
  # 3) run fgsea
  fg <- fgsea(
    pathways = ref_sigs,
    stats    = ranks,
    nperm    = 2000  # bump to 5k/10k for final paper
  )
  
  # 4) store NES in matrix
  nes_mat[j, fg$pathway] <- fg$NES
}


temp=t(nes_mat)[which(apply(t(nes_mat),1,function(x) !any(is.na(x)))),]
temp2=get_top_genes(temp,10)$top_diffs
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pheatmap(temp2,cluster_cols = FALSE,
         color = my_colors,
         breaks = seq(-1, 1, length.out = 101))
