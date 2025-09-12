tar_load(CV_metrics_full_file)

metrics = list()
for(i in 1:length(CV_metrics_full_file)){
  print(i)
  metrics[[i]] = readRDS(CV_metrics_full_file[[i]])
}

CV_metrics_full = dplyr::bind_rows(metrics)

CV_metrics = CV_metrics_full %>%
  group_by(k,alpha,lambda,eta,lambdaW,lambdaH) %>%
  summarise(bic_mean = mean(bic,na.rm=TRUE),
            bic_sd = sd(bic,na.rm=TRUE),
            nfold = sum(!is.na(bic))) %>%
  filter(nfold==5) %>%
  ungroup()

selected_params = CV_metrics %>%
  slice_min(order_by = bic_mean, n = 50, with_ties = FALSE) 

temp1 = CV_metrics_full %>% filter(fold==1) %>%
  slice_min(order_by = bic, n = 1, with_ties = FALSE) 

temp2 = CV_metrics_full %>% filter(fold==2) %>%
  slice_min(order_by = bic, n = 1, with_ties = FALSE) 

temp3 = CV_metrics_full %>% filter(fold==3) %>%
  slice_min(order_by = bic, n = 1, with_ties = FALSE) 

temp4 = CV_metrics_full %>% filter(fold==4) %>%
  slice_min(order_by = bic, n = 1, with_ties = FALSE) 

temp5 = CV_metrics_full %>% filter(fold==5) %>%
  slice_min(order_by = bic, n = 1, with_ties = FALSE) 

temp = temp1