purrr::walk(
  list.files("R/simulation_functions", full.names = TRUE, pattern = "[.]R$", recursive = TRUE),
  source
)

purrr::walk(
  list.files("R", full.names = TRUE, pattern = "[.]R$", recursive = TRUE),
  source
)

library(DeSurv)



data$surv$dataset = "sim"
data_clean=preprocess_data(
  data$counts,
  data$surv$time,
  data$surv$status,
  dataset=data$surv$dataset,
  ngene=300,
  method_trans_train = "none")
kept_genes = rownames(data_clean$ex)
pairs(t(data$X_mean)%*%data$W_true)
pairs(t(data_clean$ex)%*%data$W_true[kept_genes,])
pairs(t(data$counts[marker_genes,])%*%data$W_true[marker_genes,])


marker_genes = unlist(lapply(data$marker_info,function(x) x$genes))
length(intersect(marker_genes,kept_genes))


temp = data.frame(scale(t(data$X_mean)%*%data$W_true))
facs = colnames(temp)
temp$time = data$surv$time
temp$event = data$surv$status
cox_formula <- as.formula(
  paste0("Surv(", "time", ", ", "event", ") ~ ",
         paste(facs, collapse = " + "))
)
sfit = coxph(cox_formula,data=temp)
summary(sfit)

sfit2 = coxph(Surv(time,event)~prog1,data=temp)
summary(sfit2)

sfit3 = coxph(Surv(time,event)~prog2+prog3+prog4,data=temp)
summary(sfit3)

#check that we can recover marker genes from true W
tops=get_top_genes(data$W_true,100)
length(intersect(unlist(tops$top_genes),marker_genes))
