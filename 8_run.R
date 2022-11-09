##### use this file to run a simulation

library(bayesm)
library(MASS)
library(EnvStats)
library(DIRECT)
library(nnet)
library(LaplacesDemon)
library(parallel)
library(doParallel)
library(foreach)
library(R.utils)
library(igraph)
library(dplyr)
library(tcltk)

library(ggplot2)
library(corrplot)
library(mclust)
library(cowplot)
setwd("D:\\5G\\VEM_code")

source("7_baseline.R")
#modify here to change a simulation setting
## setting:
# Q: number of responses
# N: number of nodes
# P: number of covariates
# K:****parameter**** number of states K>=2
# T_: number of time points

MINLOG = log(.Machine$double.xmin)
MAXLOG = log(.Machine$double.xmax)


set.seed(1)
Q = 2
N = 10
K = 2
P = 2
T_ = 100
T_ = as.integer(T_/0.8)
spar = 0.3
rho = 0.5
stepsize = 0.01
g_n = 5
gc = 1:N
Z = matrix(0,N,N)
for(i in 1:g_n)
{
  nlist = sample(gc,N/g_n,replace = F)
  Z[nlist,nlist] = rbinom((N/g_n)^2,1,prob = spar) * (1-diag(N/g_n))
  gc = gc[-nlist]
}
Z[which(rowSums(Z)==0),sample(1:N,1)] = 1
Z[which(rowSums(Z)==0),sample(1:N,1)] = 1
diag(Z) = 0
#Z = matrix(rbinom(N*N,1,prob = spar),N,N) * (1-diag(N))
BOUND=pmax(1/max(rowSums(Z)),1/max(colSums(Z)))

#THETA = generate_theta(N,K,Q,P,Z,T)
KAPPA = rWishart(1,Q,diag(rep(0.025,Q)))[,,1]

# parameters setting
# THETA is a list of real parameters
set.seed(829)
THETA = list()
THETA = generate_theta(N,K,Q,P,Z,is_network=TRUE,seed=123)
THETA$var_pi = c(0.7,0.3)
THETA$var_A = matrix(c(0.7,0.2,0.3,0.8),K,K)
THETA$delta = matrix(c(5,-5,2,-2),nrow=K,ncol=Q)
THETA$beta = matrix(rnorm(Q*P,0,0.25),P,Q)
THETA$sigma = rWishart(1,Q,diag(rep(0.0025,Q)))[,,1]
THETA$W = matrix(rnormTrunc(N*N,0,0.5,0,1),N,N) * Z
THETA$W = t(apply(THETA$W,1,function(x){
  if(sum(x)>BOUND)
  {
    BOUND*x/sum(x)
  }
  else
    x
}))
#THETA$W = diag(rnormTrunc(N,0,2/30,0,1)) %*% Z

#generate covariates
# x21 = apply(array(1:T_),1,function(t) abs(3*sin(7*t)))
# x1 = apply(array(1:N),1,function(n) ifelse(runif(T_)<0.1,1,0))  #T*N
# x2 = apply(array(1:N),1,function(n) rnorm(T_,0,0.025)) #T*N
# x2 = apply(array(1:N),1,function(n) rnorm(T_,0,0.025))
# xt = array(0,dim=c(N,T_,P))
# xt[,,1] = scale(x1)
# xt[,,2] = scale(x21+x2)

source("app_all.R")
source("1_generate_sim_case.R")
source("2_initialize_function.R")

#----------------include the functions-----------------

#######whole run
#outer: number of case
#rep_res: estimation result
#pred: prediction result
outer = 30
mod_num = 6
nhmm_theta = data.frame()
#nhmm_mae = data.frame()
#fit_rmse = data.frame()
#fit_ari = data.frame()
pred_rmse = array(0,dim=c(Q,mod_num,outer))
time_mark = array(0,dim=c(mod_num,outer))
yfit = list()
pb = tkProgressBar("Simulation"," finished %",0,100)
for(gg in 1:outer)
{
  progress_info = sprintf("finished %d%%", round((gg-1)/outer*100))
  setTkProgressBar(pb,round((gg-1)*100/outer),sprintf("Progress (%s)",progress_info),progress_info)
  prmse = c()
  ptime = c()
  #pari = c()
  #-------------run this code to generate a simulation case-----------------------#
  K=2
  theta = generate_theta(N,K,Q,P,Z,seed = gg^2)
  # theta = THETA
  train_test_ratio = 1
  data = generate_data_set(q_=Q,t_=as.integer(T_/0.8),
                           n_=N,k_=K,p_=P,pars=THETA,seed=gg^2)
  
  sep = as.integer(T_)
  train_set = list(obs=data$obs[,,(1:sep)],cov = data$cov[,(1:sep),])
  test_set = list(obs=data$obs[,,((sep+1):T_)],cov = data$cov[,((sep+1):T_),])

  ypred = data.frame(as.vector(aperm(test_set$obs,c(3,2,1))))
  #----------------------------estimate the parameters---------------------------#
  ##parameter of VEM algorithm
  # G:generation of VEM
  # precision: convergence condition
  # train
  G=20
  precision = 1e-4
  cl.cores = detectCores()-1
  cl = makeCluster(cl.cores)
  clusterEvalQ(cl,library(DIRECT))
  clusterEvalQ(cl,library(MASS))
  clusterEvalQ(cl,library(foreach))
  clusterEvalQ(cl,library(iterators))
  clusterEvalQ(cl,library(nnet))
  clusterExport(cl,varlist = ls(environment(VEM_unify)),envir = environment())
  clusterExport(cl,varlist = ls(environment(app_po)),envir = environment())
  clusterExport(cl,varlist = ls(environment(app_update_tau)),envir = environment())
  #clusterExport(cl,varlist = ls(environment(gradw)),envir = environment())
  K=2
  #cat("nhmm\n")
  p1=proc.time()
  res_app_vem = VEM_unify(G,precision,train_set$obs,train_set$cov,N,is_network = T)
  #nhmm_theta = rbind(nhmm_theta,unlist(res_app_vem$theta_est))
  #nhmm_mae=rbind(nhmm_mae,unlist(mapply(function(x,y) abs(x-y),theta,res_app_vem$theta_est)))
  #pari = c(pari,adjustedRandIndex(res_app_vem$infs$lt,data$ltn[,1:sep]))
  ptime = c(ptime,(proc.time()-p1)[3])
  #cat("NHMM,",proc.time()-p1,"\n")
  #yfitt = fitting_nhmm(train_set$cov,res_app_vem$infs$lt,res_app_vem$theta_est,T)
  #pyfit = cbind(pyfit,as.vector(aperm(yfitt[,,-1],c(3,2,1))))
  #prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((yfitt[,q,]-train_set$obs[,q,])^2)))))
  # cat(prmse)

  progress_info = sprintf("finished %d%%", round((gg-1+1/mod_num)/outer*100))
  setTkProgressBar(pb,round((gg-1+2/mod_num)/outer*100),sprintf("Progress (%s)",progress_info),progress_info)

  #------------------------------one step predict-----------------------------#
  s0 = res_app_vem$inf$lt[,sep]
  clusterExport(cl,varlist = ls(environment(app_predict_nhmm)),envir = environment())
  py = app_predict_nhmm(s0,test_set$obs,test_set$cov,res_app_vem$theta_est,as.integer(N/2),is_network = T)
  ypred = cbind(ypred,as.vector(aperm(py$y,c(3,2,1))))
  #nhmm_predict[[gg]] = py
  pre_prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((py$y[,q,]-test_set$obs[,q,])^2)))))
  #pari = cbind(pari,adjustedRandIndex(py$s,data$ltn[,(sep+1):T_]))

  #--------------------baseline1: naive model, MLE
  cat("naive\n")
  K=1
  p1=proc.time()
  res_bsl = MLE_naive(G,precision,train_set$obs,train_set$cov)
  ptime = c(ptime,(proc.time()-p1)[3])
  yfitt = fitting_naive(train_set$cov,res_bsl$theta_est)
  pyfit = cbind(pyfit,as.vector(aperm(yfitt[,,-1],c(3,2,1))))
  prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((yfitt[,q,]-train_set$obs[,q,])^2)))))
  # py_bsl = predict_naive(test_set$obs,test_set$cov,res_bsl$theta_est)
  # prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((py_bsl[,q,]-test_set$obs[,q,])^2)))))
  progress_info = sprintf("finished %d%%", round((gg-1+2/mod_num)/outer*100))
  setTkProgressBar(pb,round((gg-1+2/mod_num)/outer*100),sprintf("Progress (%s)",progress_info),progress_info)

  # 
  # #--------------------baseline2: feHMM
  K=2
  #N*
  #cat("fehmm\n")
  p1=proc.time()
  res_bsl = VEM_unify(G,precision,train_set$obs,train_set$cov,N,is_network = F)
  ptime = c(ptime,(proc.time()-p1)[3])
  yfitt = fitting_nhmm(train_set$cov,res_bsl$infs$lt,res_bsl$theta_est,F)
  pyfit = cbind(pyfit,as.vector(aperm(yfitt[,,-1],c(3,2,1))))
  prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((yfitt[,q,]-train_set$obs[,q,])^2)))))
  #pari = c(pari,adjustedRandIndex(res_bsl$infs$lt,data$ltn[,1:sep]))
  # s0 = res_bsl$infs$lt[,sep]
  # py_bsl = app_predict_nhmm(s0,test_set$obs,test_set$cov,res_bsl$theta_est,N-1,is_network = F)
  # prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((py_bsl$y[,q,]-test_set$obs[,q,])^2)))))
  # pari = cbind(pari,adjustedRandIndex(py_bsl$s,data$ltn[,(sep+1):T_]))
  progress_info = sprintf("finished %d%%", round((gg-1+3/mod_num)/outer*100))
  setTkProgressBar(pb,round((gg-1+3/mod_num)/outer*100),sprintf("Progress (%s)",progress_info),progress_info)

  # 
  # #--------------------baseline3: mixVAR
  #cat("mixVAR\n")
  p1=proc.time()
  res_bsl = estimate_MVAR(train_set$obs,train_set$cov,K)
  ptime = c(ptime,(proc.time()-p1)[3])
  yfitt = predict_MVAR(train_set$obs,train_set$cov,res_bsl,1)
  pyfit = cbind(pyfit,as.vector(aperm(yfitt,c(3,2,1))))
  prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((yfitt[,q,]-train_set$obs[,q,2:sep])^2)))))

  progress_info = sprintf("finished %d%%", round((gg-1+4/mod_num)/outer*100))
  setTkProgressBar(pb,round((gg-1+4/mod_num)/outer*100),sprintf("Progress (%s)",progress_info),progress_info)

  # mvar_obs = data$obs[,,sep:T_]
  # mvar_cov = xt[,sep:T_,]
  # py_bsl = predict_MVAR(mvar_obs,mvar_cov,res_bsl,1)
  # prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((py_bsl[,q,]-test_set$obs[,q,])^2)))))

  # #--------------------baseline4: NAR
  #cat("NAR\n")
  p1=proc.time()
  W_nar = t(apply(Z,1,function(x) x/sum(x)))
  # ypred = array(0,dim=c(N,Q,dim(test_set$obs)[3]))
  yfitt = array(0,dim=c(N,Q,dim(train_set$obs)[3]-1))
  for(q in 1:Q)
  {
    Ymat = train_set$obs[,q,]
    theta_nar= betaOLS(Ymat,W_nar,train_set$cov)
    yfitt[,q,] = pred_nar(train_set$obs[,q,],train_set$cov[,2:sep,],theta_nar,W_nar)
    # ypred[,q,] = pred_nar(nar_obs,test_set$cov,theta_nar,W_nar)
  }
  ptime = c(ptime,(proc.time()-p1)[3])
  pyfit = cbind(pyfit,as.vector(aperm(yfitt,c(3,2,1))))
  prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((yfitt[,q,]-train_set$obs[,q,2:sep])^2)))))
  # prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((ypred[,q,]-test_set$obs[,q,])^2)))))
  progress_info = sprintf("finished %d%%", round((gg-1+5/6)/outer*100))
  setTkProgressBar(pb,round((gg-1+5/6)/outer*100),sprintf("Progress (%s)",progress_info),progress_info)

  # 
  # #--------------------baseline5: CHMM(py)
  #cat("GChmm\n")
  p1=proc.time()
  res_bsl = estimate_CHMM(train_set$obs,train_set$cov,K,Z)
  ptime = c(ptime,(proc.time()-p1)[3])
  yfitt = fitting_chmm(train_set$cov,res_bsl$status[,,1],res_bsl$theta)
  pyfit = cbind(pyfit,as.vector(aperm(yfitt[,,-1],c(3,2,1))))

  bsl_s = t(res_bsl$status[,,1])
  prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((yfitt[,q,]-train_set$obs[,q,])^2)))))
  #pari = c(pari,adjustedRandIndex(bsl_s,data$ltn[,1:sep]))

  # s0 = res_bsl$status[sep,,1]
  # py_bsl = predict_CHMM(s0,test_set$obs,test_set$cov,res_bsl$theta)
  # prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((py_bsl$y[,q,]-test_set$obs[,q,])^2)))))
  # py_bsl$s = sum(1:K) - py_bsl$s
  # pari = cbind(pari,adjustedRandIndex(py_bsl$s,data$ltn[,(sep+1):T_]))
  yfit[[gg]] = pyfit
  fit_rmse = rbind(fit_rmse,prmse)
  time_mark[,gg] = ptime
  #fit_ari = rbind(fit_ari,pari)
  stopCluster(cl)
  progress_info = sprintf("finished %d%%", round(gg*100/outer))
  setTkProgressBar(pb,round(gg*100/outer),sprintf("Progress (%s)",progress_info),progress_info)
}
close(pb)
filename=paste("yfit_",as.character(N),".Rdata",sep="")
save(yfit,file=filename)
rmse_model_name = c("NHMM","Naive","FEHMM","mVAR","NAR","GCHMM")
ari_model_name = c("NHMM","FEHMM","GCHMM")
rmse_name=c()

for(a in rmse_model_name)
{
  for(q in 1:Q)
    rmse_name =c(rmse_name,paste(a,"_",as.character(q),sep=""))
}

names(fit_rmse) = rmse_name
#names(fit_ari) = ari_model_name

#table_est(nhmm_mae)
table_fit(fit_rmse,rmse_model_name)


#visual
G = graph_from_adjacency_matrix(Z,mode='directed')
plot(sub1,layout=layout.gem(sub1),edge.arrow.size=0,vertex.color="#6E9ECE",
     vertex.frame.color="white",vertex.label="",
     vertex.size=5,edge.width=0.2,edge.color="black")

corrplot(Z,method = "square",tl.col = "black",is.corr = F,col= c("white","#6E9ECE"))

corrplot(cor(t(as.data.frame(array(yobs,c(N,Q*T_))))),method="color",tl.col="black",tl.cex=0.5)

draw_obs_dens_st(data$obs)

node = sample(1:N,1)
nodelist=sort(c(node,which(Z[node,]>0)[1:3]),decreasing = F)
draw_series(yobs,nodelist)

draw_rmse(fit_rmse,rmse_model_name)

dffit = yfit[[which.min(fit_rmse[,1])]]
names(dffit) = c("Real",rmse_model_name)
draw_comp(pyfit,names(dffit),N,T_-1,Q)

rmse_model_name=c("NHMM","Naive","FEHMM","NAR","GCHMM","mVAR")
res_name = c("Accessibility","Downlin_Latency","Downlin_Throughput","Drop_Call_Rate","Uplink_Throughput")

fit_rmse = data.frame()
for(i in 2:7)
{
  rmse_n = (pyfit[,i] - pyfit[,1])^2
  rmse_n = array(rmse_n,c(T_,Q,N))
  rmse_n = cbind(as.vector(sqrt(apply(rmse_n,2,mean))),rep(rmse_model_name[(i-1)],Q),res_name)
  fit_rmse = rbind(fit_rmse,rmse_n)
}
names(fit_rmse) = c("RMSE","Model","Response")
fit_rmse$Model = factor(fit_rmse$Model,levels = rmse_model_name)
fit_rmse$Response = as.factor(fit_rmse$Response)

draw_rmse(fit_rmse,rmse_model_name)


#---------------------pred------------------------------------------

outer = 1
mod_num = 6

pred_rmse = array(0,dim=c(Q,mod_num,outer))
time_mark = array(0,dim=c(mod_num,outer))
yp = list()

pb = tkProgressBar("Simulation"," finished %",0,100)
for(gg in 1:outer)
{
  progress_info = sprintf("finished %d%%", round((gg-1)/outer*100))
  setTkProgressBar(pb,round((gg-1)*100/outer),sprintf("Progress (%s)",progress_info),progress_info)
  
  ptime = c()
  m = 1

  #-------------run this code to generate a simulation case-----------------------#
  K=2
  theta = generate_theta(N,K,Q,P,Z,seed = gg^2)
  train_test_ratio = 1
  data = generate_data_set(q_=Q,t_=as.integer(T_),
                           n_=N,k_=K,p_=P,pars=THETA,seed=gg^2)
  
  sep = as.integer(T_*0.8)
  train_set = list(obs=data$obs[,,(1:sep)],cov = data$cov[,(1:sep),])
  test_set = list(obs=data$obs[,,((sep+1):T_)],cov = data$cov[,((sep+1):T_),])
  
  ypred = data.frame(as.vector(aperm(test_set$obs,c(3,2,1))))
  #----------------------------estimate the parameters---------------------------#
  ##parameter of VEM algorithm
  # G:generation of VEM
  # precision: convergence condition
  # train
  G=20
  precision = 1e-4
  cl.cores = detectCores()-1
  cl = makeCluster(cl.cores)
  clusterEvalQ(cl,library(DIRECT))
  clusterEvalQ(cl,library(MASS))
  clusterEvalQ(cl,library(foreach))
  clusterEvalQ(cl,library(iterators))
  clusterEvalQ(cl,library(nnet))
  clusterExport(cl,varlist = ls(environment(VEM_unify)),envir = environment())
  clusterExport(cl,varlist = ls(environment(app_po)),envir = environment())
  clusterExport(cl,varlist = ls(environment(app_update_tau)),envir = environment())
  #clusterExport(cl,varlist = ls(environment(gradw)),envir = environment())
  K=2
  #cat("nhmm\n")
  p1=proc.time()
  res_app_vem = VEM_unify(G,precision,train_set$obs,train_set$cov,N,is_network = T)

  ptime = c(ptime,(proc.time()-p1)[3])
  
  progress_info = sprintf("finished %d%%", round((gg-1+1/mod_num)/outer*100))
  setTkProgressBar(pb,round((gg-1+2/mod_num)/outer*100),sprintf("Progress (%s)",progress_info),progress_info)
  
  #------------------------------one step predict-----------------------------#
  s0 = res_app_vem$inf$lt[,sep]
  clusterExport(cl,varlist = ls(environment(app_predict_nhmm)),envir = environment())
  py = app_predict_nhmm(s0,test_set$obs,test_set$cov,res_app_vem$theta_est,as.integer(N/2),is_network = T)
  ypred = cbind(ypred,as.vector(aperm(py$y,c(3,2,1))))
  
  pred_rmse[,m,gg] = sapply(1:Q, function(q) sqrt(mean((py$y[,q,]-test_set$obs[,q,])^2)))

  
  #--------------------baseline1: naive model, MLE
  cat("naive\n")
  m = m + 1
  K=1
  p1=proc.time()
  res_bsl = MLE_naive(G,precision,train_set$obs,train_set$cov)
  ptime = c(ptime,(proc.time()-p1)[3])
  py_bsl = predict_naive(test_set$obs,test_set$cov,res_bsl$theta_est)
  ypred = cbind(ypred,as.vector(aperm(py_bsl,c(3,2,1))))
  pred_rmse[,m,gg] = sapply(1:Q, function(q) sqrt(mean((py_bsl[,q,]-test_set$obs[,q,])^2)))
  
  # 
  # prmse = cbind(prmse,t(sapply(1:Q, function(q) sqrt(mean((py_bsl[,q,]-test_set$obs[,q,])^2)))))
  progress_info = sprintf("finished %d%%", round((gg-1+2/mod_num)/outer*100))
  setTkProgressBar(pb,round((gg-1+2/mod_num)/outer*100),sprintf("Progress (%s)",progress_info),progress_info)
  
  # 
  # #--------------------baseline2: feHMM
  K=2
  m = m + 1
  #N*
  #cat("fehmm\n")
  p1=proc.time()
  res_bsl = VEM_unify(G,precision,train_set$obs,train_set$cov,N,is_network = F)
  ptime = c(ptime,(proc.time()-p1)[3])

  s0 = res_bsl$infs$lt[,sep]
  py = app_predict_nhmm(s0,test_set$obs,test_set$cov,res_bsl$theta_est,N-1,is_network = F)
  pred_rmse[,m,gg] = sapply(1:Q, function(q) sqrt(mean((py$y[,q,]-test_set$obs[,q,])^2)))
  ypred = cbind(ypred,as.vector(aperm(py$y,c(3,2,1))))

  progress_info = sprintf("finished %d%%", round((gg-1+3/mod_num)/outer*100))
  setTkProgressBar(pb,round((gg-1+3/mod_num)/outer*100),sprintf("Progress (%s)",progress_info),progress_info)
  
  # 
  # #--------------------baseline3: mixVAR
  #cat("mixVAR\n")
  m = m + 1
  p1=proc.time()
  res_bsl = estimate_MVAR(train_set$obs,train_set$cov,K)
  ptime = c(ptime,(proc.time()-p1)[3])
  
  mvar_obs = data$obs[,,sep:T_]
  mvar_cov = data$cov[,sep:T_,]
  py_bsl = predict_MVAR(mvar_obs,mvar_cov,res_bsl,1)
  pred_rmse[,m,gg] = sapply(1:Q, function(q) sqrt(mean((py_bsl[,q,]-test_set$obs[,q,])^2)))
  ypred = cbind(ypred,as.vector(aperm(py_bsl,c(3,2,1))))
  
  progress_info = sprintf("finished %d%%", round((gg-1+4/mod_num)/outer*100))
  setTkProgressBar(pb,round((gg-1+4/mod_num)/outer*100),sprintf("Progress (%s)",progress_info),progress_info)
  

  
  # #--------------------baseline4: NAR
  #cat("NAR\n")
  m = m + 1
  p1=proc.time()
  W_nar = t(apply(Z,1,function(x) x/sum(x)))
  py_bsl = array(0,dim=c(N,Q,dim(test_set$obs)[3]))
  # yfitt = array(0,dim=c(N,Q,dim(train_set$obs)[3]-1))
  for(q in 1:Q)
  {
    nar_obs = data$obs[,q,sep:(T_-1)]
    Ymat = train_set$obs[,q,]
    theta_nar= betaOLS(Ymat,W_nar,train_set$cov)
    py_bsl[,q,] = pred_nar(nar_obs,test_set$cov,theta_nar,W_nar)
  }
  ptime = c(ptime,(proc.time()-p1)[3])
  pred_rmse[,m,gg] = sapply(1:Q, function(q) sqrt(mean((py_bsl[,q,]-test_set$obs[,q,])^2)))
  ypred = cbind(ypred,as.vector(aperm(py_bsl,c(3,2,1))))
  
  progress_info = sprintf("finished %d%%", round((gg-1+5/6)/outer*100))
  setTkProgressBar(pb,round((gg-1+5/6)/outer*100),sprintf("Progress (%s)",progress_info),progress_info)
  
  # 
  # #--------------------baseline5: CHMM(py)
  #cat("GChmm\n")
  m = m + 1
  p1=proc.time()
  res_bsl = estimate_CHMM(train_set$obs,train_set$cov,K,Z)
  ptime = c(ptime,(proc.time()-p1)[3])
  
  s0 = res_bsl$status[sep,,1]
  py = predict_CHMM(s0,test_set$obs,test_set$cov,res_bsl$theta)
  pred_rmse[,m,gg] = sapply(1:Q, function(q) sqrt(mean((py$y[,q,]-test_set$obs[,q,])^2)))
  ypred = cbind(ypred,as.vector(aperm(py$y,c(3,2,1))))

  
  yp[[gg]] = ypred
  time_mark[,gg] = ptime

  
  stopCluster(cl)
  progress_info = sprintf("finished %d%%", round(gg*100/outer))
  setTkProgressBar(pb,round(gg*100/outer),sprintf("Progress (%s)",progress_info),progress_info)
}
close(pb)
filename=paste("ypred_",as.character(N),".Rdata",sep="")
save(yp,file=filename)
model_name = c("NHMM","Naive","FEHMM","mVAR","NAR","GCHMM")

mean_rmse = apply(pred_rmse, c(1,2), mean)
names(mean_rmse)[2] = model_name
print(mean_rmse)
names(time_mark)[2] = model_name
print(time_mark)

