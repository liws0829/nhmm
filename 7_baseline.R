#baseline

##1. MOHMM
## y_it = F_i + x_it*B + lambda_sti + epsilon_it (N*Q)
## set W=0


##2. Naive model
##set K=1
MLE_naive<-function(iteration,presision,obs,covar)
{
  n_ = dim(obs)[1]
  q_=dim(obs)[2]
  t_=dim(obs)[3]
  p_=dim(covar)[3]
  k_ = 1
  var_theta = init_theta(n_,k_,q_,p_,Z,T)
  old_theta = var_theta
  loglik= 0.0
  for(g in 1:iteration)
  {
    # beta
    revD = ginv(diag(n_) - var_theta$W)
    lbd = as.vector(rowSums(revD)) %*% var_theta$delta
    ut = obs - array(rep(lbd,t_),dim(obs))

    if(!is.null(p_))
    {
      term1 = kronecker(ginv(var_theta$sigma),matrix(rowSums(apply(covar,2,crossprod)),p_,p_))  #pq*pq
      term2 = as.vector(matrix(rowSums(sapply(1:t_,function(t) crossprod(covar[,t,],ut[,,t]))),p_,q_)%*%ginv(var_theta$sigma)) #pq*1
      newbeta=matrix(solve(term1,term2),p_,q_)
      # cat("beta\n")
      # print(mean(abs(newbeta-THETA$beta)))
    }
    else
      newbeta = rnorm(0)


    #delta
    if(!is.null(p_))
      xb = array(apply(covar,2,function(x) x%*%newbeta),dim=c(n_,q_,t_))
    else
      xb = array(0,c(n_,q_,t_))
    ut = obs - xb

    rt = as.vector(rowSums(revD))
    term1 = t(as.vector(apply(ut,c(1,2),sum))) %*%  kronecker(ginv(var_theta$sigma),rt)
    term2 = t_*kronecker(ginv(var_theta$sigma),crossprod(rt))

    newdelta = matrix(solve(t(term2),t(term1)),1,q_)

    #W
    ta = apply(ut,c(1,2),sum)
    tb = matrix(rep(newdelta,each=n_),n_,q_)
    baterm = matrix(0,n_*n_,n_*n_)
    cterm = as.vector(rep(0,n_*n_))

    for(q1 in 1:q_)
    {
      for(q2 in 1:q_)
      {
        b_ = 0.5 * (tb[,q1]%*%t(tb[,q2]) + tb[,q2]%*%t(tb[,q1]))
        baterm = baterm + kronecker(b_,ginv(diag(rep(var_theta$sigma[q1,q2],n_))))
        cterm = cterm + as.vector(ginv(diag(rep(var_theta$sigma[q1,q2],n_)))%*%ta[,q1]%*%t(tb[,q2]))
      }
    }

    vecD = matrix(1/t_*ginv(baterm)%*%cterm,n_,n_)
    newW = (diag(n_) - ginv(vecD)) * Z
    newW = ifelse(newW<0,0,newW)
    newW = ifelse(newW>1,1,newW)

    #sigma
    revD = ginv(diag(n_)-newW)
    lbd =  as.vector(rowSums(revD)) %*% newdelta
    et = ut - array(rep(lbd,t_),dim=c(n_,q_,t_))

    temp = matrix(rowSums(apply(et,3,function(x) crossprod(t(x)))),n_*q_,n_*q_)
    newsigma = matrix(0,q_,q_)
    for(q1 in 1:q_)
    {
      for(q2 in 1:q_)
      {
        newsigma[q1,q2] = tr(temp[((q1-1)*n_+1):(q1*n_),][,((q2-1)*n_+1):(q2*n_)])
      }
    }
    newsigma = newsigma/(n_*t_)
    newsigma = pmin(newsigma,KAPPA)

    var_theta = list(var_pi=1.0, var_A=1.0, delta=newdelta, beta=newbeta,sigma=newsigma, W=newW)
    convg = max(unlist(mapply(function(x,y) abs(x-y),old_theta,var_theta)))
    loglik = sum(sapply(1:t_,function(t) stpo(obs[,,t],NULL,covar[,t,],var_theta,log=TRUE,is_network=T)))
    old_theta = var_theta

    if(g>1 && convg<=presision) break()
  }
  res=list(theta_est= var_theta,loglik=loglik)
  return(res)
}
#ynt = xnt*b+revD*delta+ep
predict_naive<-function(obs,covar,var_theta)
{
  n_ = dim(obs)[1]
  q_ = dim(obs)[2]
  t_ = dim(obs)[3]
  predict_y = array(0,dim=c(n_,q_,t_))

  revD = ginv(diag(n_)-var_theta$W)

  for(t in 1:t_)
  {
    #given inference state of t, and covar of t, predict y_t+1
    xb = covar[,t,]%*%var_theta$beta
    rt = revD %*% as.vector(rep(1,n_))%*%matrix(var_theta$delta,ncol=q_)
    predict_y[,,t] = xb + rt
  }
  res = predict_y
  return(res)
}
fitting_naive<-function(covar,var_theta)
{
  n_ = dim(covar)[1]
  t_ = dim(covar)[2]
  q_ = ncol(var_theta$beta)
  revD = ginv(diag(n_)-var_theta$W)
  lbd = matrix(rep(var_theta$delta,each=n_),n_,q_)

  xb = array(apply(covar,2,function(x) x%*%var_theta$beta),c(n_,q_,t_))
  rt = revD%*%lbd

  yfit = xb + rep(rt,t_)
  return(yfit)
}

##3. GMVAR(1,2)
source("modi_mAR.R")
estimate_MVAR<-function(obs,covar,k_)
{
  n_=dim(obs)[1]
  q_=dim(obs)[2]
  t_=dim(obs)[3]
  p_=dim(covar)[3]
  AR = list()
  for(k in 1:k_)
  {
    AR[[k]] = array(kronecker(diag(q_),diag(n_)),dim=c(n_*q_,n_*q_,1))
  }
  prob <- rep(1/k_,k_)
  shift <- matrix(0,n_*q_,k_)
  vcov <- array(rep(kronecker(diag(q_),diag(n_)),k_),dim=c(n_*q_,n_*q_,k_))
  mod <- new("MixVARGaussian", prob = prob, vcov = vcov, arcoef = AR, shift = shift)
  MVAR_fit <- mixVARreg(obs,covar,mod,niter=10)
  return(MVAR_fit)
}
predict_MVAR<-function(obs,covar,var_theta,p)
{
  #t_ = (sep - p) : T_
  n_=dim(obs)[1]
  q_=dim(obs)[2]
  t_=dim(obs)[3]
  p_=dim(covar)[3]
  xb = apply(covar,2,function(x) x%*%var_theta$reg) #NQ*(test+p)
  ut = t(array(obs,dim=c(n_*q_,t_)) - xb ) #(test+p)*NQ
  pk <- var_theta$mixVARmodel@order
  g <- length(pk)

  ypred <- matrix(0,t_- p,n_*q_)
  for(k in 1:g){
    for(t in 1:(t_ - p)){
      for(u in 1:pk[k])
      {
        ypred[t,] <- ypred[t,] + var_theta$mixVARmodel@prob[k] *
          (var_theta$mixVARmodel@arcoef@a[[k]][,,u] %*% ut[(t+p-u),])
      }
      ypred[t,] <- ypred[t,] + xb[,(t+p)]
    }
  }
  ypred = array(t(ypred),dim=c(n_,q_,t_-p))
  return(ypred)
}

##4. GCHMM
source("modi_CHMM.R")
estimate_CHMM<-function(obs,covar,k_,adj)
{
  n_ = dim(obs)[1]
  q_ = dim(obs)[2]
  p_ = dim(covar)[3]
  t_ = dim(obs)[3]
  var_pi = rep(0,k_)
  var_A = matrix(0,k_,k_)
  beta = matrix(0,p_,q_)
  delta = matrix(0,k_,q_)
  sigma = c()
  omega = c()
  scor = matrix(apply(sapply(1:q_, function(q) cor(t(obs[,q,]))),1,mean),n_,n_)*adj
  ltn = array(0,dim=c(t_,n_,q_))
  for(q in 1:q_)
  {
    Xmat = t(obs[,q,])
    chmm = coupledHMM(Xmat,covar,nb.states = k_, S=scor, itmax = 20)
    res_chmm = chmm$model
    var_pi = var_pi + res_chmm$initPr
    var_A = var_A + res_chmm$transPr
    beta[,q] = res_chmm$esBeta
    delta[,q] = res_chmm$esAvg
    sigma = c(sigma , unique(res_chmm$esVar))
    omega = c(omega,chmm$omega)
    ltn[,,q] = chmm$status
  }
  omega = mean(omega)
  var_pi = var_pi/q_
  var_pi = var_pi/sum(var_pi)
  var_A = var_A/q_
  var_A = t(apply(var_A,1,function(x) x/sum(x)))
  sigma = diag(sigma)
  theta = list(var_pi = var_pi, var_A = var_A, beta = beta,
             delta = delta, sigma = sigma, scor= scor, omega = omega)
  res = list(theta=theta,status = ltn)
  return(res)
}
predict_CHMM<-function(s0,obs,covar,var_theta)
{

  # s0: n*p
  #var_theta: res_chmm$theta
  k_ = length(var_theta$var_pi)
  n_ = dim(obs)[1]
  q_ = dim(obs)[2]
  t_ = dim(obs)[3]
  p = length(s0)/n_

  st = gen_tr(n_,k_)$st
  sl = apply(st,1,function(x) sum(var_theta$scor[lower.tri(var_theta$scor)][which(dist(x)>0)]))
  logw = sl*log(var_theta$omega)
  xb = array(apply(covar,2,function(x) x%*%var_theta$beta),dim=c(n_,q_,t_)) #NQ*T
  py = array(0,dim=c(n_,q_,t_))
  inf_s = array(0,dim=c(n_,t_))
  #wl
  #reshape

  #dens
  logdens = c()
  for(n in 1:n_)
  {
    tempq = sapply(1:q_,function(q) Emis.Gauss(obs[n,q,]-xb[n,q,],var_theta$delta[,q],rep(diag(var_theta$sigma)[q],k_)))
    tempq = apply(tempq,1,sum)
    logdens = c(logdens,tempq)
  }
  logdens = array(logdens,dim=c(t_,k_,n_))

  for(t in 1:t_)
  {
    if(t<=p)
    {
      if(!is.null(dim(s0)))
        s = s0[,t]
      else
      {
        if(length(s0)==n_)
          s = s0
        if(length(s0)!=n_)
          stop("s0 input must be n*t")
      }
    }
    else
    {
      s = inf_s[,t-p]
    }
    #predict
    temp = logw + apply(st,1,function(x) tr(p*log(var_theta$var_A[s,x])))  #NS
    temp = exp(temp - max(temp))
    temp = temp/sum(temp)   #NS
    dt = apply(rep(rep(temp,each=n_),q_) * array(var_theta$delta[t(st),],dim=c(n_,nrow(st),q_)),c(1,3),sum)

    py[,,t] = xb[,,t] + dt
    #inference
    for(i in 1:n_)
    {
      s[i] = which.is.max(sapply(1:k_,function(k) sum(temp[which(st[,i]==k)]) * logdens[t,k,i]))
    }
    inf_s[,t] = s
  }
  res = list(y=py,s=inf_s)
  return(res)
}
fitting_chmm<-function(covar,infs,var_theta)
{
  #infs:t*n
  n_ = dim(covar)[1]
  t_ = dim(covar)[2]
  q_ = ncol(var_theta$beta)
  xb = array(apply(covar,2,function(x) x%*%var_theta$beta),c(n_,q_,t_))
  rt = array(apply(infs,1,function(x) var_theta$delta[x,]),c(n_,q_,t_))
  yfit = xb + rt
  return(yfit)
}
##5. NAR
source("NAR_helper.R")
pred_nar<-function(obs,covar,theta_nar,w_nar)
{
  #time:[sep+1:(T_-1)]
  #pred:[sep+2:T_]
  #obs:N*pt
  Ymat = w_nar %*% obs
  #X:(N*pt-1)*(P+3)
  X = cbind(rep(1, nrow(obs) * (ncol(obs)-1)), as.vector(Ymat[, -ncol(Ymat)]),
            as.vector(Ymat[, -ncol(Ymat)]), apply(covar[,-ncol(Ymat),],3,rbind))
  Ypred = matrix(X %*% theta_nar$theta,nrow=N)
  return(Ypred)
}


