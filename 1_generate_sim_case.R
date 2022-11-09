
trans<-function(num,n,k)
{
  num = num -1
  vec = c()
  while(num!=0)
  {
    vec = c(num%%k,vec)
    num = floor(num/k)
  }
  if(length(vec)<n)
  {
    sup = rep(0,(n-length(vec)))
    vec = c(sup,vec)
  }
  return(vec+1)
}

revtrans<-function(vec,n,k)
{
  vec = vec - 1
  a = k^(seq(n-1,0,-1))
  num = sum(vec*a)
  return(num+1)
}

gen_tr<-function(n_,k_)
{
  if(k_==1)
  {
    st_n = array(1,dim=c(1,n_))
    st01 = array(1,dim=c(1,n_,k_))
    res = list(A = 1.0, pi=1.0, st = st_n, st_01=st01)
    return(res)
  }

  stlist = list()
  for(n in 1:n_)
  {
    stlist[[n]] = 1:k_
  }
  st_n = as.matrix(expand.grid(stlist)[,n_:1])


  st01 = array(class.ind(st_n),dim=c(nrow(st_n),n_,k_))
  ns = nrow(st_n)
  pi = rep(1/ns,ns)
  A = matrix(1/ns,ns,ns)
  res = list(A = A, pi=pi, st = st_n, st_01=st01)
  return(res)
}
#mean of lambda:Q*K

#generate simulation data
# obs: observation Y:N*Q*T
# lt: latent state S:N*M*T

generate_data_set<-function(q_=Q,t_=T_,n_=N,k_=K,p_=P,pars=THETA,seed=0) #generate data set: given parameters and covariates return latent state sequence and data
{
  if(seed!=0)
    set.seed(seed)
  xb = array(0,dim=c(n_,q_,t_))

  if(p_==0)
  {
    covs = NULL
  }
  if(p_>0)
  {
    covs = array(0,dim=c(n_,t_,p_))
    for(p in 1:P)
    {
      sd_rand = runif(1)
      xp = rnorm((n_*t_),sd=sd_rand)
      type = sample(1:3,1)
      if(type==1)
      {
        covs[,,p] = array(scale(xp),dim=c(n_,t_))
      }
      else
      {
        if(type==2)
        {
          xpt = abs(3*sin(7*(1:t_)))
          covs[,,p] = array(scale(xp+xpt),dim=c(n_,t_))
        }
        if(type==3)
        {
          xpt = rnorm(n_)
          covs[,,p] = array(scale(rep(xpt,each=t_)+xp),dim=c(n_,t_))
        }
        if(type==4)
        {
          covs[,,p] = array(scale(ifelse(xp>0.7,1,0)),dim=c(n_,t_))
        }
      }
    }
    xb = array(apply(covs,2,function(x) x%*%pars$beta),dim=c(n_,q_,t_))
  }


  obs = array(0,dim=c(n_,q_,t_))
  ltn = array(0,dim=c(n_,t_))

  revD = ginv(diag(n_) - pars$W)
  bigsgm = kronecker(pars$sigma,diag(n_))
  ## generate latent state sequence and observations

  # t=1
  temp = sample(1:k_,n_,replace = T,prob = pars$var_pi)
  ltn[,1] = temp

  rt = revD%*%pars$delta[ltn[,1],]
  ep = as.vector(rMVNorm(1,mean=rep(0,n_*q_),sigma = bigsgm))

  obs[,,1] = as.vector(xb[,,1] + rt) + ep
  #cat(xb,"+\n",rt,"+\n",ep,"+\n",obs[,,1],"\n")

  # t>1

  for(t in 2:t_)
  {
    temp = ltn[,(t-1)]
    newtemp = sapply(temp,function(i) sample(1:k_,1,replace = T,prob=pars$var_A[i,]))
    ltn[,t] = newtemp

    rt = revD%*%pars$delta[ltn[,t],]
    ep = as.vector(rMVNorm(1,mean=rep(0,n_*q_),sigma = bigsgm))
    obs[,,t] = as.vector(xb[,,t] + rt) + ep

  }

  res=list(obs=obs,cov=covs,ltn=ltn)

  return(res)
}


generate_theta<-function(n_,k_,q_,p_,z,is_network=TRUE,seed=0)
{
  if(seed!=0)
    set.seed(seed)
  if(k_>1)
  {
    rand_num = sample(1:10,k_)
    var_pi = rDirichlet(1,rand_num)
    var_A = rDirichlet(k_,rand_num)
    diag(var_A) = diag(var_A) + 5
    var_A = var_A/rowSums(var_A)
    delta = rnorm(q_)

    deltai = t(sapply(1:(k_-1),function(x) runif(q_,(x-1)*5+1,x*5)))
    delta = rbind(delta, rep(delta,each=k_-1)+deltai)
    
    delta = matrix(delta,k_,q_)
  }
  if(k_==1)
  {
    var_pi = 1.0
    var_A = matrix(1,1,1)
    delta = matrix(rnorm(q_),k_,q_)
  }

  rand_num = runif(q_,0,0.5)
  sigma = rWishart(1,q_,diag(rand_num))[,,1]

    
  beta = matrix(rnorm(q_*p_),p_,q_)

  if(is_network)
  {
    rand_num = runif(1,0,0.5)
    W = matrix(rnormTrunc(n_*n_,rand_num,0.5,0,1),n_,n_) * z
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
    fe = matrix(rnorm(n_*q_),n_,q_)
    W = 0*z
    res=list(var_pi=var_pi, var_A=var_A, delta=delta, beta=beta,sigma=sigma, W=W, fe=fe)
  }

  #W = matrix(rnorm(n_*n_,0,0.01),n_,n_) * z

  #W = THETA$W
  return(res)
}
