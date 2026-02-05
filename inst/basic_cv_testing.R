
library(devtools)
library(targets)
tar_load(data_filtered)
load_all('../DeSurv')

source("R/load_data.R")
source("R/load_data_internal.R")
source("R/set_folds.R")

data = load_data("TCGA_PAAD")

data_folds = set_folds(data,nfold = 5,method_trans_train="rank")

dat = data.frame(matrix(ncol=7,nrow=0))
colnames(dat) = c("alpha","fold","ctrain","ctest","sloss","nloss","pen")
for(a in seq(0,.3,by=.05)){
  print(a)
  for(f in 1:5){
    print(f)
    data_train = data_folds$data_train[[f]]
    data_test = data_folds$data_test[[f]]
    
    fit = desurv_fit(X=data_train$ex,y=data_train$sampInfo$time,d=data_train$sampInfo$event,
                     k=3, alpha=a, lambda=.001, nu=.1,lambdaW=0,lambdaH=0,ninit=10,
                     verbose = FALSE,maxit=500,
                     max_iter_beta = 1000,tol=1e-5)
    ctrain = fit$cindex
    
    ctest = cvwrapr::getCindex(t(data_test$ex)%*%fit$W%*%fit$beta,
                       survival::Surv(data_test$sampInfo$time,data_test$sampInfo$event))
    dat[nrow(dat)+1,] = c(a,f,ctrain,ctest,fit$loss$surv_loss,fit$loss$nmf_loss,fit$loss$penalty_beta)
    
  }
}


library(ggplot2)
library(dplyr)
dat_mean=dat %>% group_by(alpha) %>% summarise(mctr=mean(ctrain),mcte=mean(ctest),sctr=sd(ctrain),scte=sd(ctest))
ggplot(dat_old, aes(x=alpha,y=ctest,color=as.factor(fold)))+
  geom_point()+
  geom_line()
ggplot(dat_mean, aes(x=alpha,y=mcte))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=mcte-scte,ymax=mcte+scte))

plot(fit$Wlossit)
plot(fit$lossit)
plot(fit$slossit)
plot(fit$nlossit)


dat = data.frame(slossit = c(fit$slossit,fit$slossitW,fit$slossitH),
           nlossit = c(fit$nlossit,fit$nlossitW,fit$nlossitH),
           lossit = c(fit$lossit,fit$lossitW,fit$lossitH),
           iter = rep(1:length(fit$lossit),3),
           update = rep(c("beta","W","H"),each=length(fit$lossit)))

library(ggplot2)
library(dplyr)

dat2=dat %>% filter(iter < 100 & iter>80)

ggplot(dat2)+
  geom_point(aes(x=iter,y=slossit,color=update))
