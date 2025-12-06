dat=simulate_desurv_easy()
fit_true=coxph(Surv(time,status)~linpred_true,data=dat$surv)
summary(fit_true)
X = dat$counts
dat_new=preprocess_data(X,dat$surv$time,dat$surv$status,
                        dataset=rep("sim",length(dat$surv$time)),
                        ngene = 2000,
                        method_trans_train = "quant")
# Xrank=apply(X,2,rank)
genes=intersect(rownames(dat$W_true),rownames(dat_new$ex))
# ex_rank = apply(t(dat_new$ex),2,rank)
Z=scale(t(dat_new$ex[genes,]) %*%dat$W_true[genes,])
# Z=scale(t(X)%*%dat$W_true)
dat$surv[,colnames(Z)] = as.data.frame(Z)
sfit=coxph(Surv(time,status)~prog1+prog2+prog3+prog4,data=dat$surv)
summary(sfit)
summary(dat$surv$linpred_true)
sd(dat$surv$linpred_true)
quantile(dat$surv$linpred_true, c(.05,.95))


fits=list()
i=1
for(a in seq(0,.95,by=.05)){
fits[[i]]=glmnet::cv.glmnet(t(dat_new$ex),
                          Surv(dat_new$sampInfo$time,dat_new$sampInfo$event),
                          family="cox",
                          alpha=a,
                          type.measure="C")
i=i+1
}

which.max(sapply(fits,function(x){
  max(x$cvm)
}))

beta=coef(fits[[8]]$glmnet.fit,s=fits[[8]]$lambda.min)

linpred=t(dat_new$ex) %*% beta
dat_glm = data.frame(linpred=linpred)
dat_glm$time=dat_new$sampInfo$time
dat_glm$event = dat_new$sampInfo$event

fit_glm=coxph(Surv(time,event)~X1,data=dat_glm)
summary(fit_glm)
