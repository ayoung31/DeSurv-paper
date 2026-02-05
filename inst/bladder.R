mot = readRDS("data/original/IMmotion150_sampInfo_ex.rds")
vig = readRDS("data/original/IMVigor210_sampInfo_ex_clinical.rds")

save(mot,"data/original/")

hist(as.numeric(log2(as.matrix(mot$ex[,3:ncol(mot$ex)])+1)))
hist(as.numeric(log2(as.matrix(vig$ex)+1)))

# tpms not log transformed
tcga=readRDS("data/original/TCGA_PAAD.rds")
hist(as.numeric(log2(as.matrix(tcga$ex)+1)))

# already log transformed tpms
cptac=readRDS("data/original/CPTAC.rds")
hist(2^as.numeric(as.matrix(cptac$ex)))

# raw counts
dijk=readRDS("data/original/Dijk.rds")
hist(as.numeric(log2(as.matrix(dijk$ex)+1)))

# log transformed microarray
# note microarray data is already quantile normalized after log transform
moff = readRDS("data/original/Moffitt_GEO_array.rds")
hist(as.numeric(as.matrix(moff$ex)))


hist(as.numeric(as.matrix(dijk$ex)))
hist(as.numeric(as.matrix(tcga$ex)))


# 
# library(biomaRt)
# 
# # Connect to Ensembl
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# # Example: get gene length for a vector of Ensembl IDs
# genes <- tcga$featInfo$SYMBOL
# X=tcga$ex
# rownames(X) = genes

# genes = vig$featInfo
# X=vig$ex
# rownames(X) = genes


# df_tx <- getBM(
#   attributes = c("external_gene_name","ensembl_gene_id",
#                  "ensembl_transcript_id","transcript_length"),
#   filters = "external_gene_name",
#   values  = genes,
#   mart    = mart
# )
# 
# # Collapse transcript lengths to a per-gene length (max or median are common)
# df_gene <- aggregate(transcript_length ~ ensembl_gene_id + external_gene_name,
#                      data = df_tx, FUN = max)  # or median
# names(df_gene)[3] <- "gene_length_tx_based"
# head(df_gene)
# df_gene$symbol = df_gene$external_gene_name
# 
# gene_avg=rowMeans(X)
# 
# gene_avg = data.frame(avg_exp = gene_avg,symbol=names(gene_avg))
# 
# dat = gene_avg %>% left_join(df_gene)
# 
# plot(dat$avg_exp,dat$gene_length_tx_based)
# 
# cor(dat$avg_exp,dat$gene_length_tx_based,use="complete.obs",method="spearman")
