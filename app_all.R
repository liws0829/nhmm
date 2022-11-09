# require
stpo<-function(yt,st,covt,var_theta,log=TRUE,is_network)
{
  # yt: observations, NQ
  # st: all combination of St, ns*N
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
  if(!is.null(st))
    st = matrix(st,ncol=n_)

  if(!is_network)
    yt = yt - var_theta$fe

  # yt = as.vector(yt)
  # covt = kronecker(diag(q_),covt)

  if(!is.null(covt))
    xb = covt %*% var_theta$beta
  else
    xb = 0
  revD = ginv(diag(n_) - var_theta$W)
  bigsgm = kronecker(var_theta$sigma,diag(n_))
  if(!isSymmetric(bigsgm))
    bigsgm[lower.tri(bigsgm)] = t(bigsgm)[lower.tri(bigsgm)]

  if(K==1)
  {
    rt = revD %*% matrix(rep(var_theta$delta,each=n_),nrow=n_)
    # rt = kronecker(diag(q_),rowSums(revD)) %*% as.vector(var_theta$delta)
    y = as.vector(yt-xb-rt)
    logdens = dMVNorm(y,mean=rep(0,n_*q_),sigma = bigsgm,log=T)
    return(logdens)
  }

  if(K>1)
  {
    if(!is.null(nrow(st)))
      ns = nrow(st)
    else
    {
      if(length(ns)==n_)
        ns = 1
    }

    logdens = apply(array(1:ns),1,function(x)
    {
      rt = revD %*% var_theta$delta[st[x,],]
      # rt = kronecker(diag(q_),revD) %*% as.vector(var_theta$delta[st[x,],])
      y = as.vector(yt - xb - rt)
      #cat(x,apply(st[x,,],1,function(z) which.is.max(z)),dMVNorm(y,mean=rep(0,n_*q_),sigma = bigsgm,log=TRUE),"\n")
      dMVNorm(y,mean=rep(0,n_*q_),sigma = bigsgm,log=T)
    }
    )
    return(logdens)
  }

  #result: N(yt|xt,var_theta,all_S_combination) ns*1
}

app_po<-function(yt,covt,var_theta,limit,log=TRUE,is_network)
{
  k_ = K
  n_ = nrow(yt)
  q_ = ncol(yt)
  if(!is.null(covt))
    ut = yt - (covt %*% var_theta$beta) #N*Q
  else
    ut = yt

  if(k_==1)
  {
    st_new = rep(1,n_)
    logdens = stpo(yt,st_new,covt,var_theta,is_network=is_network)
    res = cbind(st_new,logdens)
    colnames(res) = c(1:n_,"logdens")
    return(as.data.frame(res))
  }

  if(!is_network)
    ut = ut - var_theta$fe

  if(limit>0)
  {
    bigW = diag((diag(n_)-var_theta$W)%*%t(diag(n_)-var_theta$W))
    mid = (diag(n_)-var_theta$W) %*% ut
    temp = sapply(1:n_,function(n)
      {
      sapply(1:k_,function(i){
        dMVNorm(as.vector(mid[n,]-var_theta$delta[i,]),rep(0,q_),sigma = bigW[n]*var_theta$sigma,log=T)
      })
    }) #K*N
    n_star = sort(order(apply(apply(temp,2,function(x) (max(x)-x)),2,max),decreasing=T)[1:limit])
    k_star = rep(0,n_)
    k_star[n_star]  = apply(temp,2,which.is.max)[n_star]
    if(limit<n_)
    {
      st_num = gen_tr((n_-limit),k_)$st
      st_new = matrix(0,nrow=nrow(st_num),ncol=N)
      st_new[,n_star] = rep(k_star[n_star],each=nrow(st_new))
      st_new[,which(colSums(st_new)==0)] = st_num
    }
    else
      st_new = matrix(k_star,nrow=1)
  }
  else
    st_new  = gen_tr(n_,k_)$st


  logdens = stpo(yt,st_new,covt,var_theta,is_network=is_network)
  res = cbind(st_new,logdens)
  colnames(res) = c(1:n_,"logdens")
  return(as.data.frame(res))
}

app_elbo_log<-function(tau,theta,log_dens)
{
  n_=dim(tau$taum)[2]
  t_=dim(tau$taum)[1]
  k_=dim(tau$taum)[3]
  taul = tau$taul
  taut = tau$taut
  taum = tau$taum
  #term pi
  loglik_pi = sum(taul * (rep(log(theta$var_pi),each=n_) - log(taul)))
  #print(loglik_pi)
  #term trans
  loglik_tr = 0.0
  taut_log = taut * (rep(log(theta$var_A),t_*n_) - log(taut))
  # taut_log = taut * (array(rep(log(theta$var_A),each=(t_*n_)),dim=c(t_,n_,K,K)) - log(taut))

  #term response
  #registerDoParallel(cl)
  loglik_rp = sum(log_dens$logdens * log_dens$gama)

  for(n in 1:n_)
  {
    #taum_log_sum = taum_log_sum + log(taum[,n,st_num[,n]])
    loglik_tr = loglik_tr + sum(apply(array(2:t_),1,function(t) tr(diag(taum[(t-1),n,])%*%taut_log[,,n,t])))
  }

  #taum_log_sum=ifelse(taum_log_sum<.Machine$double.min.exp,.Machine$double.min.exp,taum_log_sum)

  #print(loglik_rp)
  loglik = loglik_pi + loglik_tr + loglik_rp

  return(loglik)
}

app_update_tau<-function(log_dens,tau,var_pi,var_A)
{

  # log_dens:ns*T, pre calculate log-density
  # st_num: combination of all state
  # tau:old tau, list of taul,taut,taum

  t_ = dim(tau$taum)[1]
  n_ = dim(tau$taum)[2]
  k_ = K
  newtaul = array(0,dim=c(n_,k_))
  #newtaut = array(0,dim=c(t_,n_,k_,k_))

  #update taul(t=1)

  for(n in 1:n_)
  {
    #taul_prod
    tauln = sapply(1:k_,function(k){
      nsi = log_dens[(log_dens$t==1) & (log_dens[,n]==k),]
      if(nrow(nsi)==0)
        -Inf
      else
      {
        #nrow(nsi)
        temp = apply(nsi,1,function(x)
        {
          u = matrix(0,(n_-1),k_)
          u[,sort(unique(x[(1:n_)[-n]]))] = class.ind(x[(1:n_)[-n]])
          sum(log(tau$taul[-n,])*u)
        })
        temp = exp(temp-max(temp))
        sum((temp/sum(temp))*nsi['logdens']) + log(var_pi[k])
      }
    })
    tauln = exp(tauln - max(tauln))
    newtaul[n,] = tauln/sum(tauln)
  }

  alltaut = parApply(cl,array(2:t_),1,function(t)
  {
    sapply(1:n_,function(n)
    {
      log_taut = sapply(1:k_,function(k)
      {
        nsi = log_dens[(log_dens$t==t)&(log_dens[,n]==k),]
        if(nrow(nsi)==0)
          rep(-Inf,k_)
        else
        {
          temp = apply(nsi,1,function(x)
          {
            u = matrix(0,(n_-1),k_)
            u[,sort(unique(x[(1:n_)[-n]]))] = class.ind(x[(1:n_)[-n]])
            sum(log(tau$taum[t,-n,])*u)
          })
          temp = exp(temp-max(temp))
          sum((temp/sum(temp))*nsi['logdens']) + log(var_A[,k])
        }
      })
      log_taut = apply(log_taut,1,function(x) exp(x-max(x)))
      t(log_taut/colSums(log_taut))
    })
  })
  alltaut = cbind(rep(0,n_*k_*k_),alltaut)
  newtaut = array(alltaut,dim=c(k_,k_,n_,t_))


  newtaul = correct_taul(newtaul)
  newtaut = correct_taut(newtaut)


  # #print(apply(newtaul,1,sum))
  # #print(apply(newtaut,c(1,2,3),sum))
  #
  newtaum = update_taum(newtaul,newtaut)
  res=list(taul=newtaul,taut=newtaut,taum=newtaum)
  return(res)
}

app_inference_s<-function(var_tau)
{
  t_=dim(var_tau$taum)[1]
  n_=dim(var_tau$taum)[2]
  k_ = dim(var_tau$taum)[3]

  inf_lt = array(apply(var_tau$taum,1,function(x) {x = matrix(x,nrow=n_)
  apply(x,1,function(y) which.is.max(y))}),dim=c(n_,t_))   #N*T


  res = list(lt=inf_lt)
  return(res)
}

gradw<-function(varw,ut,rt,sig,rho,z)
{
  n_ = dim(ut)[1]
  q_ = dim(ut)[2]
  t_ = dim(ut)[3]

  revD = ginv(diag(n_)-varw)

  temp = matrix(0,n_,n_)
  for(q1 in 1:q_)
  {
    for(q2 in 1:q_)
    {
      temp = temp + matrix(rowSums(sapply(1:t_,function(t)
      {
        y = tcrossprod(rt[,q1,t],rt[,q2,t])
        b_ = 0.5*(y+t(y))
        a_ = tcrossprod(ut[,q1,t],rt[,q2,t])
        ginv(sig)[q1,q2]*( varw %*% b_ - a_)
      })),n_,n_)
    }
  }
  grad = t(revD) %*% temp %*% t(revD) + rho * varw
  return(grad*z)
}

VEM_unify<-function(iteration,presision,obs,covar,limit,is_network=TRUE)
{
  is_grad = F
  if(K<=1)
  {
    stop("K must larger than 1, if K==1, use MLE_navie() instead")
  }
  loglik = rep(0,iteration)

  n_ = dim(obs)[1]
  q_ = dim(obs)[2]
  t_ = dim(obs)[3]
  p_ = dim(covar)[3]

  #pb = txtProgressBar(style = 3)
  #start_time = Sys.time()
  #init
  var_tau = init_tau(n_,K,t_)
  var_theta = init_theta(n_,K,q_,p_,Z,is_network)

  kc = kmeans(array(aperm(obs,c(1,3,2)),c(n_*t_,q_)),K)
  var_theta$delta = kc$centers[order(kc$centers[,1]),]

  # var_theta = THETA

  registerDoParallel(cl)
  log_dens = foreach(t=1:t_,.combine = rbind) %dopar%
    {
      if(!is.null(p_))
        temp = app_po(obs[,,t],covar[,t,],var_theta,limit,log=T,is_network=is_network)
      else
        temp = app_po(obs[,,t],NULL,var_theta,limit,log=T,is_network=is_network)
      temp$t = t
      temp
    }
  log_dens$gama = apply(log_dens,1,function(x)
  {
    u = matrix(0,n_,K)
    u[,sort(unique(as.numeric(x[1:n_])))] = class.ind(as.numeric(x[1:n_]))
    sum(log(var_tau$taum[x['t'],,])*u)
  })

  log_dens$gama = (group_by(log_dens,t) %>%
                     dplyr::summarise(gama=exp(gama - max(gama)),.groups = 'drop')  %>%
                     group_by(t) %>%
                     dplyr::summarise(gama=gama/sum(gama),.groups = 'drop'))$gama

  gama = log_dens[log_dens$gama!=0,]
  temp = t(apply(gama,1,function(x)
  {
    u = matrix(0,n_,K)
    u[,sort(unique(x[1:n_]))] = class.ind(x[1:n_])
    c(u*x['gama'],x['t'])

  }))
  gama = array(sapply(1:t_,function(t) colSums(matrix(temp[temp[,'t']==t,-ncol(temp)],ncol=ncol(temp)-1))),c(n_,K,t_))  #NK*T
  loglik[1] = app_elbo_log(var_tau,var_theta,log_dens)

  epsilon=c()
  old_theta = var_theta
  # cat("1,loglik,",loglik[1],"\n")

  for(g in 2:iteration)
  {

    #VE
    var_tau = app_update_tau(log_dens,var_tau,var_theta$var_pi,var_theta$var_A)

    # temp = abs(apply(var_tau$taum,1,function(x) apply(matrix(x,nrow=n_),1,which.is.max)) - data$ltn[,1:sep])
    # print(length(which(colSums(temp)>0)))

    #M

    #pi
    newpi = apply(var_tau$taul,2,mean)
    # cat("pi\n")
    # print(mean(abs(newpi-THETA$var_pi)))

    #A
    newA = t(apply(rep(aperm(var_tau$taum[-t_,,],c(3,2,1)),each=K)
                   * aperm(var_tau$taut[,,,-1],c(2,1,3,4)),c(1,2),sum))
    denom = apply(var_tau$taum[-t_,,], 3, sum)
    newA = newA/denom
    newA = matrix(pmax(.Machine$double.xmin,newA),K,K)
    newA = newA/rowSums(newA)
    # cat("A\n")
    # print(mean(abs(newA-THETA$var_A)))


    #gama
    log_dens$gama = apply(log_dens,1,function(x)
    {
      u = matrix(0,n_,K)
      u[,sort(unique(as.numeric(x[1:n_])))] = class.ind(as.numeric(x[1:n_]))
      sum(log(var_tau$taum[as.numeric(x['t']),,])*u)
    })
    log_dens$gama = (group_by(log_dens,t) %>%
                       dplyr::summarise(gama=exp(gama - max(gama)),.groups = 'drop')  %>%
                       group_by(t) %>%
                       dplyr::summarise(gama=gama/sum(gama),.groups = 'drop'))$gama
    gama = log_dens[log_dens$gama!=0,]
    temp = t(apply(gama,1,function(x)
    {
      u = matrix(0,n_,K)
      u[,sort(unique(x[1:n_]))] = class.ind(x[1:n_])
      c(u*x['gama'],x['t'])

    }))
    gama = array(sapply(1:t_,function(t) colSums(matrix(temp[temp[,'t']==t,-ncol(temp)],ncol=ncol(temp)-1))),c(n_,K,t_))  #NK*T
    revD = ginv(diag(n_)-var_theta$W)
    # bgs = kronecker(ginv(var_theta$sigma),diag(n_))

    #beta
    if(!is.null(p_))
    {
      lbd = array(apply(gama,3,function(x) revD%*%x%*%var_theta$delta),c(n_,q_,t_))
      if(!is_network)
        ut = obs - rep(var_theta$fe,t_) - lbd
      else
        ut = obs - lbd
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
      xb = array(0,dim=c(n_,q_,t_))
    if(!is_network)
      ut = obs - rep(var_theta$fe,t_) - xb
    else
      ut = obs - xb
    rt = array(apply(gama,3,function(x) revD%*%x),c(n_,K,t_))
    term1 = kronecker(ginv(var_theta$sigma), matrix(rowSums(apply(rt,3,crossprod)),K,K))
    term2 = as.vector(matrix(rowSums(sapply(1:t_,function(t) crossprod(rt[,,t],ut[,,t]))),K,q_)%*%ginv(var_theta$sigma))
    newdelta = matrix(solve(term1-rho*diag(K*q_),term2),K,q_)
    # cat("delta\n")
    # print(mean(abs(newdelta-THETA$delta)))

    rt = array(apply(gama,3,function(x) x%*%newdelta),c(n_,q_,t_))
    #if not network, fe
    newfe = matrix(0,n_,q_)
    newW = matrix(0,n_,n_)
    if(!is_network)
    {
      term2 = apply(ut-rt,c(1,2),mean)
      newfe = matrix(term2,n_,q_)
      newW = matrix(0,n_,n_)
    }

    #else, W
    if(is_network)
    {
      if(is_grad==T)
      {
        batch_size = as.integer(0.2 * t_)
        tlist = ((1:floor(t_/batch_size))-1)*batch_size +1
        registerDoParallel(cl)
        MiniBat=foreach(ti=iter(tlist),.combine = "+") %dopar%
          {
            te = pmin(((ti + batch_size)-1),t_)
            f_ = gradw(var_theta$W,ut[,,ti:te],rt[,,ti:te],var_theta$sigma,rho,Z)
            f_
          }

        grad = MiniBat/length(tlist)* Z
        print(grad)

        newW = (var_theta$W + stepsize * grad)
        newW = ifelse(newW<0,0,newW)
        newW = ifelse(newW>1,1,newW)
        newW[which(rowSums(newW)>BOUND),] = BOUND *
          newW[which(rowSums(newW)>BOUND),] /
          rowSums(newW)[which(rowSums(newW)>BOUND)]
      }
      else
      {
        term1 = matrix(0,n_^2,n_^2)
        term2 = as.vector(rep(0,n_^2))
        for(q1 in 1:q_)
        {
          for(q2 in 1:q_)
          {
            temp = matrix(0.5* (tcrossprod(rt[,q1,],rt[,q2,]) + t(tcrossprod(rt[,q1,],rt[,q2,]))),n_,n_)
            term1 = term1 + kronecker(temp,diag(rep(ginv(var_theta$sigma)[q1,q2],n_)))
            term2 = term2 + as.vector(ginv(var_theta$sigma)[q1,q2] * tcrossprod(ut[,q1,],rt[,q2,]))
          }
        }
        revD = matrix(solve(term1,term2),n_,n_)

        newW = (diag(n_) - ginv(revD)) * Z
        newW = ifelse(newW>1,1,newW)
        newW = ifelse(newW<0,0,newW)
        newW[which(rowSums(newW)>BOUND),] = BOUND *
          newW[which(rowSums(newW)>BOUND),] /
          rowSums(newW)[which(rowSums(newW)>BOUND)]
      }
      # print(mean(abs(newW-THETA$W)))
    }

    # cat("W",newW,"\n")
    #sigma

    revD = ginv(diag(n_)-newW)
    lbd = array(apply(gama,3,function(x) revD%*%x%*%var_theta$delta),c(n_,q_,t_))
    et = obs - xb - lbd
    if(!is_network)
      et = et - rep(newfe,t_)

    temp = matrix(rowSums(apply(et,3,function(x) tcrossprod(as.vector(x)))),n_*q_,n_*q_)
    newsigma = matrix(0,q_,q_)
    for(q1 in 1:q_)
    {
      for(q2 in 1:q_)
      {
        newsigma[q1,q2] = tr(temp[((q1-1)*n_+1):(q1*n_),((q2-1)*n_+1):(q2*n_)])
      }
    }
    newsigma = newsigma/(n_*t_)
    if(all(newsigma>KAPPA))
      newsigma = KAPPA

    #cat("sigma",newsigma,"\n")
    # print(mean(abs(newsigma)))

    if(!is_network)
      var_theta = list(var_pi=newpi, var_A=newA, delta=newdelta, beta=newbeta,sigma=newsigma, W=newW, fe=newfe)
    else
      var_theta = list(var_pi=newpi, var_A=newA, delta=newdelta, beta=newbeta,sigma=newsigma, W=newW)

    # print(mapply(function(x,y) mean(abs(x-y)),var_theta,THETA))
    cvg = max(unlist(mapply(function(x,y) abs(x-y),old_theta,var_theta)))
    old_theta = var_theta

    registerDoParallel(cl)
    log_dens = foreach(t=1:t_,.combine = rbind) %dopar%
      {
        if(!is.null(p_))
          temp = app_po(obs[,,t],covar[,t,],var_theta,limit,log=T,is_network=is_network)
        else
          temp = app_po(obs[,,t],NULL,var_theta,limit,log=T,is_network=is_network)
        temp$t = t
        temp
      }
    loglik[g] = app_elbo_log(var_tau,var_theta,log_dens)
    ep = abs(loglik[g]-loglik[g-1])
    epsilon = c(epsilon, ep)

    # cat(g,",loglik,",loglik[g],"\n")
    #setTxtProgressBar(pb, g/iteration)
    if(g>1 && cvg <= presision) break()

  }

  #inference
  s_star = app_inference_s(var_tau)
  res=list(theta_est=var_theta, converge = loglik,  infs = s_star)
  #end_time = Sys.time()
  #close(pb)
  return(res)
}

app_predict_nhmm<-function(s0,obs,covar,var_theta,limit,is_network)
{
  # s0:n_*p

  n_ = dim(obs)[1]
  q_ = dim(obs)[2]
  t_ = dim(obs)[3]
  k_ = K
  p = length(s0)/n_

  predict_y = array(0,dim=c(n_,q_,t_))
  inf_s = matrix(0,n_,t_)
  revD = ginv(diag(n_)-var_theta$W)
  xb = array(apply(covar,2,function(x) x%*%var_theta$beta),dim=c(n_,q_,t_))

  registerDoParallel(cl)
  log_dens = foreach(t=1:t_,.combine = rbind) %dopar%
    {
      temp = app_po(obs[,,t],covar[,t,],var_theta,limit,log=T,is_network=is_network)
      temp$t = t
      temp
    }

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
      s = inf_s[,(t-p)]
    }
    temp_p = p*log(var_theta$var_A[s,])
    temp_p = apply(temp_p,1,function(x) exp(x-max(x)))
    temp_p = apply(temp_p,2,function(x) x/sum(x))
    temp = t(apply(temp_p,2,function(x) apply(x*var_theta$delta,2,sum)))  #N*Q

    if(!is_network)
    {
      predict_y[,,t] = xb[,,t] + var_theta$fe + temp
    }
    else
    {
      predict_y[,,t] = xb[,,t] + revD%*%temp
    }

    prob = apply(as.matrix(log_dens[log_dens$t==t,]),1,function(x){
      tr(log(var_theta$var_A[s,x[1:n_]])) + x['logdens']
    })
    s = as.matrix(log_dens[log_dens$t==t,1:n_][which.is.max(prob),])
    inf_s[,t] = s
  }

  res = list(y=predict_y,s=inf_s)
  return(res)
}

fitting_nhmm<-function(covar,infs,var_theta, is_network)
{
  #infs:N*T
  n_ = nrow(var_theta$W)
  q_ = ncol(var_theta$delta)
  t_ = ncol(infs)
  revD = ginv(diag(n_)-var_theta$W)

  xb = array(apply(covar,2,function(x) x%*%var_theta$beta),c(n_,q_,t_))
  rt = array(apply(infs,2,function(x){
    revD%*%var_theta$delta[x,]
  }),c(n_,q_,t_))

  yfit = xb + rt
  if(!is_network)
    yfit = yfit + rep(var_theta$fe,t_)

  return(yfit)
}



