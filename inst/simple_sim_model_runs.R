library(NMF)
purrr::walk(
  list.files("R/simulation_functions", full.names = TRUE, pattern = "[.]R$", recursive = TRUE),
  source
)

data=simulate_desurv_scenario(scenario = "R0",seed=14537)
markers = data$marker_sets
split=split_simulation_samples(sim_dataset[[1]],.7,10000)
dat=prepare_simulation_data(sim_dataset[[1]],sample_ids = split$train_ids,ngene=NULL,transform_method = "none")
#oracle
temp=data.frame(time=data$time,status=data$status,linpred=data$linpred)
sfit = coxph(Surv(time,status)~linpred,data=temp)
summary(sfit)



#alpha=0 BO
bounds=desurv_bo_default_bounds()
bounds$n_starts=NULL
bounds$nfolds=NULL
bounds$ngene=NULL
bounds$ntop=NULL
bounds$tol=NULL
bounds$maxit=NULL
bounds$lambdaW_grid=NULL
bounds$lambdaH_grid=NULL
bounds$k_grid$upper = 8
bounds$lambda_grid$lower=1e-2
bounds$lambda_grid$upper=1e2
bounds$alpha_grid=NULL
bounds
bo0 = desurv_bayesopt(X=data$X,y=data$time,d=data$status,
                      bo_bounds = bounds,
                      bo_fixed= list(n_starts=1,
                                     nfolds=5,
                                     ngene=5000,
                                     ntop=NULL,
                                     tol=1e-6,
                                     maxit=600,
                                     lambdaW_grid=0,
                                     lambdaH_grid=0,
                                     alpha_grid=0),
                      preprocess = FALSE)

#alpha=0 fit
dfit0 = desurv_fit(X=data$X,y=data$time,d=data$status,
                   k=bo0$best$params['k_grid'],
                   alpha=0,
                   lambda=bo0$best$params['lambda_grid'],
                   nu=bo0$best$params['nu_grid'],
                   lambdaW=0,lambdaH=0,
                   ninit=20,maxit=5000,tol=1e-6)
plot(dfit0$lossit)
plot(dfit0$nlossit)
plot(dfit0$slossit)

tops0=get_top_genes(dfit0$W,150)$top_genes
intersect(unlist(tops0),unlist(markers))
precision_recall(list(unlist(tops0[,4])),list(unlist(markers[[1]])))

dfit0$cindex

dfit0$beta*apply(t(data$X)%*%dfit0$W,2,sd)

recon0 = norm(data$X - dfit0$W %*% dfit0$H, "F")


#alpha>0 bo
bounds=desurv_bo_default_bounds()
bounds$n_starts=NULL
bounds$nfolds=NULL
bounds$ngene=NULL
bounds$ntop=NULL
bounds$tol=NULL
bounds$maxit=NULL
bounds$lambdaW_grid=NULL
bounds$lambdaH_grid=NULL
bounds$k_grid$upper = 8
bounds$lambda_grid$lower=1e-2
bounds$lambda_grid$upper=1e2
bounds$alpha_grid$upper=1
bounds

bo = desurv_bayesopt(X=data$X,y=data$time,d=data$status,
                      bo_bounds = bounds,
                      bo_fixed= list(n_starts=1,
                                     nfolds=5,
                                     ngene=5000,
                                     ntop=NULL,
                                     tol=1e-6,
                                     maxit=600,
                                     lambdaW_grid=0,
                                     lambdaH_grid=0),
                      preprocess = FALSE)


#alpha>0 fit
dfit = desurv_fit(X=data$X,y=data$time,d=data$status,
                  k=bo$best$params['k_grid'],
                  alpha=bo$best$params['alpha_grid'],
                  lambda=bo$best$params['lambda_grid'],
                  nu=bo$best$params['nu_grid'],
                  lambdaW=0,lambdaH=0,
                  ninit = 20, maxit=5000,tol=1e-6)
plot(dfit$lossit)
plot(dfit$nlossit)
plot(dfit$slossit)

dfit$cindex

tops=get_top_genes(dfit$W,150)$top_genes
intersect(unlist(tops),unlist(markers))
dfit$beta*apply(t(data$X)%*%dfit$W,2,sd)
precision_recall(list(unlist(tops[,2])),list(unlist(markers[[1]])))

dfit$beta*apply(t(data$X)%*%dfit$W,2,sd)

recon = norm(data$X - dfit$W %*% dfit$H, "F")


## specific alpha>0 fit
ids=sample(1:ncol(data$X),140,replace = FALSE)
data$X = data$X[,ids]
data$time = data$time[ids]
data$status=data$status[ids]
c=numeric()
i=1
for(a in seq(0,.95,by=.05)){
  dfit = desurv_fit(X=data$X,y=data$time,d=data$status,
                    k=3,
                    alpha=.6,
                    lambda=.1,
                    nu=.3,
                    lambdaW=0,lambdaH=0,
                    ninit = 1, maxit=2000,tol=1e-5)
  c[i] = dfit$cindex
  i=i+1
}

plot(seq(0,.95,by=.05),c)

plot(dfit$lossit)
plot(dfit$nlossit)
plot(dfit$slossit)

dfit$cindex

tops=get_top_genes(dfit$W,150)$top_genes
intersect(unlist(tops),unlist(markers))
dfit$beta*apply(t(data$X)%*%dfit$W,2,sd)
precision_recall(list(unlist(tops[,2])),list(unlist(markers[[1]])))

dfit$beta*apply(t(data$X)%*%dfit$W,2,sd)

recon = norm(data$X - dfit$W %*% dfit$H, "F")
