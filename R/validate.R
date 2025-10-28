validate = function(data,W,beta,sdZ,meanZ){
  genes = intersect(rownames(W),rownames(data$ex))
  W = W[genes,]
  data_filtered = preprocess_data(data,genes=genes)
  X = data_filtered$ex
  
  Z = t(X)%*%W
  lp = Z%*%beta
  Z_bin = apply(Z,2,function(x)x>median(x))
  dat = as.data.frame(Z_bin)
  dat$lp=lp
  dat$lp_bin = lp>median(lp)
  dat$y = data_filtered$sampInfo$time
  dat$delta = data_filtered$sampInfo$event
  dat$dataset = data_filtered$sampInfo$dataset
  # cvwrapr::getCindex(Surv(dat$y,dat$delta))
  fitph=coxph(Surv(y,delta)~lp+strata(dataset),data=dat)
  # print(summary(fitph))
  
  fitsurv1 = survfit(Surv(y,delta)~lp_bin,data=dat)
  names(fitsurv1$strata) = c("high","low")
  p1 = ggsurvplot(fitsurv1,pval=TRUE)+
          labs(title="Linear predictor")
  
  fitsurv2 = survfit(Surv(y,delta)~V2,data=dat)
  names(fitsurv2$strata) = c("high","low")
  p2 = ggsurvplot(fitsurv2,pval=TRUE)+
          labs(title="Factor 2")
  
  fitsurv3 = survfit(Surv(y,delta)~V3,data=dat)
  names(fitsurv3$strata) = c("high","low")
  p3 = ggsurvplot(fitsurv3,pval=TRUE)+
          labs(title="Factor 3")
  
  return(list(p1=p1$plot,p2=p2$plot,p3=p3$plot))
}
