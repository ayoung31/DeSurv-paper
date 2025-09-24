purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)

### load data
VAL_DATASETS = c("Dijk","Linehan","Moffitt_GEO_array",
                 "PACA_AU_array","PACA_AU_seq","Puleo_array")
VAL_PREFIX = paste0(VAL_DATASETS, collapse = ".")
data = readRDS(paste0("data/derv/",VAL_PREFIX,"_formatted.rds"))

data = readRDS(paste0("data/derv/","dijk","_formatted.rds"))


### load top genes
load("data/top_genes/selected_model.RData")
load("data/derv/cmbSubtypes_formatted.RData")

# specify the factors you want to cluster on
facs = c(1,2)
res = run_clustering(tops,data,top_genes,colors,facs=facs)


