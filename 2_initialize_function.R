## initialize function
# ***initialize tau and theta

#correct taul
correct_taul<-function(taul)
{
  taul = ifelse(taul<.Machine$double.xmin,.Machine$double.xmin,taul)
  taul_sum = apply(taul,1,sum)
  taul = taul/taul_sum
  return(taul)
}

#correct taut
correct_taut<-function(taut)
{
  taut = ifelse((taut<.Machine$double.xmin),.Machine$double.xmin,taut)

  k_ = dim(taut)[1]
  n_ = dim(taut)[3]
  t_ = dim(taut)[4]
  taut = array(apply(taut,c(3,4),function(x){
    t(apply(x,1,function(y) y/sum(y)))}),dim=c(k_,k_,n_,t_))
  return(taut)
}

#correct taum
correct_taum<-function(taum)
{
  taum = ifelse((taum<.Machine$double.xmin),.Machine$double.xmin,taum)
  for(t in 1:dim(taum)[1])
  {
    for(n in 1:dim(taum)[2])
    {
      taum[t,n,] = taum[t,n,]/sum(taum[t,n,])
    }
  }
  return(taum)
}

# given tau update taum
update_taum<-function(taul,taut)
{
  n_=dim(taul)[1]
  k_=dim(taul)[2]
  t_=dim(taut)[4]
  taum = array(0,dim=c(t_,n_,k_))
  taum[1,,] = taul
  for(t in 2:t_)
  {
    for(n in 1:n_)
    {
      for(j in 1:k_)
      {
        taum[t,n,j] = sum(taum[(t-1),n,]*taut[,j,n,t])
      }
    }
  }
  taum = correct_taum(taum)
  return(taum)
}

# init tau
init_tau<-function(n_,k_,t_)
{
  taul = array(1/k_,dim=c(n_,k_))
  taut = array(1/k_,dim=c(k_,k_,n_,t_))
  taum = update_taum(taul,taut)
  res = list(taul = taul,taut=taut,taum=taum)
  return(res)
}

#init theta
init_theta<-function(n_,k_,q_,p_,z,is_network=TRUE)
{
  var_pi = rep(1/k_,k_)
  var_A = matrix(1/k_,k_,k_)
  delta = matrix(rnormTrunc(k_*q_,2,min=0),k_,q_)
  delta = apply(delta,2,sort)
  delta = matrix(delta,k_,q_)
  sigma = diag(rep(0.0025,q_))
  if(!is.null(p_))
    beta = matrix(rep(1,p_*q_),p_,q_)
  else
    beta = rnorm(0)
  if(is_network)
  {
    W = matrix(rnormTrunc(n_*n_,0,0.5,0,1),n_,n_) * z
    W = t(apply(W,1,function(x){
      if(sum(x)>BOUND)
      {
        BOUND*x/sum(x)
      }
      else
        x
    }))
    res=list(var_pi=var_pi, var_A=var_A, delta=delta, beta=beta,sigma=sigma, W=W)
  }
  else
  {
    fe = matrix(0,n_,q_)
    W = 0*Z
    res=list(var_pi=var_pi, var_A=var_A, delta=delta, beta=beta,sigma=sigma, W=W, fe=fe)
  }

  #W = matrix(rnorm(n_*n_,0,0.01),n_,n_) * z

  #W = THETA$W
  return(res)
}



#parallel density(O(n))
po<-function(yt,st,covt,var_theta,log=TRUE)
{
  # yt: observations, NQ
  # st: all combination of St, ns*N*K
  # covt: NP
  # theta

  if(is.null(dim(yt)))
  {
    q_ = length(yt)
    n_=1
    p_ = length(covt)
    yt = matrix(yt,n_,q_)
    covt = matrix(covt,n_,p_)
  }

  q_ = ncol(yt)
  n_ = nrow(yt)

  yt = as.vector(yt)
  covt = kronecker(diag(q_),covt)

  xb = covt %*% as.vector(var_theta$beta)
  revD = ginv(diag(n_) - var_theta$W)

  bigsgm = kronecker(var_theta$sigma,diag(n_))


  logdens = apply(array(1:nrow(st)),1,function(x)
  {
   rt = kronecker(diag(q_),(revD %*% st[x,,])) %*% as.vector(var_theta$delta)
   y = as.vector(yt - xb - rt)
   #cat(x,apply(st[x,,],1,function(z) which.is.max(z)),dMVNorm(y,mean=rep(0,n_*q_),sigma = bigsgm,log=TRUE),"\n")
   dMVNorm(y,mean=rep(0,n_*q_),sigma = bigsgm,log=T)
  }
    )
  if(!log)
    return(exp(logdens))
  else
    return(logdens)
  #result: N(yt|xt,var_theta,all_S_combination) ns*1
}

#log_complete likelihood
elbo_log<-function(obs,covar,tau,theta,st_num,log_dens)
{
  n_=dim(obs)[1]
  q_=dim(obs)[2]
  t_=dim(obs)[3]
  taul = tau$taul
  taut = tau$taut
  taum = tau$taum
  #term pi
  loglik_pi = sum(taul * (rep(log(theta$var_pi),each=n_) - log(taul)))
  #print(loglik_pi)
  #term trans
  loglik_tr = 0.0
  taut_log = taut * (rep(log(theta$var_A),(t_*n_)) - log(taut))

  #term response

  taum_log_sum = array(0.0,dim=c(t_,nrow(st_num)))
  for(n in 1:n_)
  {
    taum_log_sum = taum_log_sum + log(taum[,n,st_num[,n]])
    loglik_tr = loglik_tr + sum(apply(array(2:t_),1,function(t) tr(diag(taum[(t-1),n,])%*%taut_log[,,n,t])))
  }
  #print(loglik_tr)
  #taum_log_sum=ifelse(taum_log_sum<.Machine$double.min.exp,.Machine$double.min.exp,taum_log_sum)
  taum_prod = exp(taum_log_sum)

  #print(t(taum_prod)*log_dens)
  loglik_rp = sum(t(taum_prod) * log_dens)
  #print(loglik_rp)
  loglik = loglik_pi + loglik_tr + loglik_rp

  return(loglik)
}

