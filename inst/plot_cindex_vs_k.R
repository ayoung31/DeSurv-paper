tar_load(desurv_bo_results_tcgacptacupper)
Xcand=as.data.frame(desurv_bo_results_tcgacptacupper$history)
Xcand1 = Xcand %>% filter(run_id==1) %>% dplyr::select(k_grid,alpha_grid,lambda_grid,nu_grid,ntop)
Xcand2 = Xcand %>% filter(run_id==2) %>% dplyr::select(k_grid,alpha_grid,lambda_grid,nu_grid,ntop)

km_fit=desurv_bo_results_tcgacptacupper$runs[[2]]$km_fit
# Xcand=km_fit@X

lower1=sapply(desurv_bo_results_tcgacptacupper$bounds[[1]],function(x)x$lower)
upper1=sapply(desurv_bo_results_tcgacptacupper$bounds[[1]],function(x)x$upper)

lower2=sapply(desurv_bo_results_tcgacptacupper$bounds[[2]],function(x)x$lower)
upper2=sapply(desurv_bo_results_tcgacptacupper$bounds[[2]],function(x)x$upper)

Xcand1_scale = as.data.frame(t(apply(Xcand1,1,function(x) (x-lower1)/(upper1-lower1))))
Xcand2_scale = as.data.frame(t(apply(Xcand2,1,function(x) (x-lower2)/(upper2-lower2))))

Xcand_scale = rbind(Xcand1_scale,Xcand2_scale)

Xcand1_lowk=Xcand1 %>% dplyr::filter(k_grid<=8)
Xcand1_lowk_scale = as.data.frame(t(apply(Xcand1_lowk,1,function(x) (x-lower1)/(upper1-lower1))))

pred_lowk = predict(desurv_bo_results_tcgacptacupper$runs[[1]]$km_fit,newdata=Xcand1_lowk_scale,type="UK")
Xcand1_lowk$mu = pred_lowk$mean
Xcand1_lowk$sd = pred_lowk$sd
Xcand1_lowk$upper = pred_lowk$upper95
Xcand1_lowk$lower = pred_lowk$lower95

pred=predict(km_fit,newdata=Xcand_scale,type="UK")

Xcand$mu=pred$mean
Xcand$sd=pred$sd
Xcand$lower = pred$lower95
Xcand$upper = pred$upper95

df1 = Xcand1_lowk %>% group_by(k_grid) %>% slice_max(mu, n = 1, with_ties = FALSE) %>%ungroup()
df2 = Xcand %>% filter(k_grid>=9)%>% group_by(k_grid) %>% slice_max(mu, n = 1, with_ties = FALSE) %>%
  ungroup() %>% dplyr::select(k_grid,alpha_grid,lambda_grid,nu_grid,ntop,mu,sd,lower,upper)

df = rbind(df1,df2)

ggplot(df, aes(x = k_grid, y = mu)) +
  geom_line() +
  geom_point() +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2
  ) +
  labs(
    x = "Number of factors (k)",
    y = "Cross-validated C-index"
  ) +
  theme_minimal()
