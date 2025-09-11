library(NMF)
library(coxNMF)


purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)

raw = load_data("TCGA_PAAD")

data = preprocess_data(raw)
nmf_fit=NMF::nmf(as.matrix(data$ex),rank=2:12,method="lee",nrun=30, .options="v4")
plot(nmf_fit)
