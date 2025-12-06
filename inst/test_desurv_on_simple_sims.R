# library(DeSurv)
library(glmnet)
library(NMF)
library(survival)

purrr::walk(
  list.files("../R/simulation_functions", full.names = TRUE, pattern = "[.]R$", recursive = TRUE),
  source
)

purrr::walk(
  list.files("../R", full.names = TRUE, pattern = "[.]R$", recursive = TRUE),
  source
)

set.seed(123)
data=simulate_desurv_easy(
  G=1000,
  N=100,
  low_mean =  1,
  high_mean = 3,
  lethal_effect =4,
  big_prog_multiplier = 1,
  lib_sdlog = .2,
  h_shape_bg = 3,
  h_shape_noise = .5,
  correlated_pairs = list()
)

Z = t(data$counts)%*%data$W_true
psfit=cv.glmnet(Z,Surv(data$surv$time,data$surv$status),family="cox",
                type.measure = "C",alpha=.1)
plot(psfit)
coef(psfit$glmnet.fit,s=psfit$lambda.min)*apply(Z,2,sd)
max(psfit$cvm)

#subset to nonzero rows
keeps=which(rowSums(data$counts)>0)
x=data$counts[keeps,]
nfit=nmf(x,4)
Z = t(data$counts)%*%nfit@fit@W
psnfit = cv.glmnet(Z,Surv(data$surv$time,data$surv$status),family="cox",
                   type.measure = "C",alpha=.9)
plot(psnfit)
coef(psnfit$glmnet.fit,s=psnfit$lambda.min)*apply(Z,2,sd)
max(psnfit$cvm)

cor(nfit@fit@W,data$W_true[keeps,])

dfit0 = desurv_fit(X=log2(data$counts),y=data$surv$time,d=data$surv$status,k=2,
                  alpha=0,lambda=.3,nu=.99,lambdaW=0,lambdaH=0,ninit=1,maxit=3000,tol=1e-7)

dfit0$cindex
Z0 = t(data$counts)%*%dfit0$W
dfit0$beta*apply(Z0,2,sd)

cor(dfit0$W,data$W_true)
intersect(get_top_genes(dfit0$W,100)$top_genes[,1],
          data$marker_info[[1]]$genes)

dfit = desurv_fit(X=log2(data$counts),y=data$surv$time,d=data$surv$status,k=2,
           alpha=.7,lambda=.3,nu=.99,lambdaW=0,lambdaH=0,ninit=1,maxit=3000,tol=1e-7)

dfit$cindex
Z = t(data$counts)%*%dfit$W
dfit$beta*apply(Z,2,sd)

cor(dfit$W,data$W_true)
intersect(get_top_genes(dfit$W,100)$top_genes[,1],
data$marker_info[[1]]$genes)
data$marker_info[[1]]