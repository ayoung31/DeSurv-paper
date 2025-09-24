purrr::walk(list.files("R", full.names = TRUE, pattern = "[.]R$"), source)

### load data (PICK ONE to uncomment below)

## combined datasets
# VAL_DATASETS = c("Dijk","Linehan","Moffitt_GEO_array",
#                  "PACA_AU_array","PACA_AU_seq","Puleo_array")
# VAL_PREFIX = paste0(VAL_DATASETS, collapse = ".")
# data = readRDS(paste0("data/derv/",VAL_PREFIX,"_formatted.rds"))

## dijk
# data = readRDS(paste0("data/derv/","dijk","_formatted.rds"))

## moffitt
# data = readRDS(paste0("data/derv/","moff","_formatted.rds"))

## paca seq
# data = readRDS(paste0("data/derv/","paca_seq","_formatted.rds"))

## paca array
# data = readRDS(paste0("data/derv/","paca_array","_formatted.rds"))

## puleo
# data = readRDS(paste0("data/derv/","puleo","_formatted.rds"))

### load top genes
load("data/top_genes/selected_model.RData")
load("data/derv/cmbSubtypes_formatted.RData")

# specify the factors you want to cluster on
facs = c(1,2)
res = run_clustering(tops,data,top_genes,colors,facs=facs)


