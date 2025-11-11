mot = readRDS("data/original/IMmotion150_sampInfo_ex.rds")
vig = readRDS("data/original/IMVigor210_sampInfo_ex_clinical.rds")
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
