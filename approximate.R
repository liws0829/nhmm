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
  
  if(!is_network)
    yt = yt - var_theta$fe
  
  yt = as.vector(yt)
  covt = kronecker(diag(q_),covt)
  
  xb = covt %*% as.vector(var_theta$beta)
  revD = ginv(diag(n_) - var_theta$W)
  bigsgm = kronecker(var_theta$sigma,diag(n_))
  
  if(K==1)
  {
    rt = kronecker(diag(q_),rowSums(revD)) %*% as.vector(var_theta$delta)
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
      rt = kronecker(diag(q_),revD) %*% as.vector(var_theta$delta[st[x,],])
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
  ut = yt - (covt %*% var_theta$beta) #N*Q
  
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
    test = apply(var_theta$delta,1,function(x) ut-rep(x,each=n_))  #NQ*K
    
    nnorm = apply(test,2,function(x) 
    {
      mid = matrix(x,nrow=n_)
      apply(mid,1,function(y) dMVNorm(as.vector(y),rep(0,q_),sigma = var_theta$sigma,log=TRUE))
    })#N*K
    
    minu = apply(nnorm,1,function(x) max(max(x)-x))
    
    n_star = which((minu/max(minu))>=0.1)
    if(length(n_star)>limit)
    {
      n_star = sort(order(minu)[1:limit])
    }

    
    k_star = apply(matrix(nnorm[n_star,],nrow =length(n_star) ),1,which.is.max)
    st_temp = gen_tr((n_-length(n_star)),k_)$st
  
    st_new = matrix(0,nrow(st_temp),n_)
    
    for(i in 1:n_)
    {
      if(i %in% n_star)
      {
        st_new[,i] = k_star[which(n_star==i)]
      }
      else
      {
        if(!is.null(dim(st_temp)))
        {
          st_new[,i] = st_temp[,1]
          st_temp = st_temp[,-1]
        }
        else
        {
          if(length(st_temp)>0)
          {

            st_new[,i] = st_temp
            st_temp = NULL
          }
        }
      }
        
    }
    
    # for(i in 1:length(n_star))
    # {
    #   #print(n_star[i])
    #   st_new = rep(k_star[i],nrow(st_temp))
    #   if(n_star[i]==1)
    #     st_temp = cbind(st_new,st_temp)
    #   else
    #   {
    #     if(n_star[i]>ncol(st_temp))
    #       st_temp = cbind(st_temp,st_new)
    #     else
    #     {
    #       st_temp = cbind(st_temp[,1:(n_star[i]-1)],st_new,st_temp[,(n_star[i]:ncol(st_temp))])
    #     }
    #   }
    #   
    #   #if((n_star[i]>1)&(n_star[i]<=ncol(st_temp)))
    # }
  }
  else
    st_new = gen_tr(n_,k_)$st
  
  # nslist = apply(st_temp,1,function(x) revtrans(x,n_,k_))
  #print(st_temp)
  logdens = stpo(yt,st_new,covt,var_theta,is_network=is_network)
  res = cbind(st_new,logdens)
  colnames(res) = c(1:n_,"logdens")
  return(as.data.frame(res))
}

app_elbo_log<-function(obs,covar,tau,theta,log_dens)
{
  n_=dim(obs)[1]
  q_=dim(obs)[2]
  t_=dim(obs)[3]
  k_=K
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
  loglik_rp = foreach(t=1:t_,.combine = "+")%do%
    {
      # st_num = t(sapply(log_dens[[t]]$nslist,function(x) trans(x,n_,k_)))
      st_num = as.matrix(log_dens[log_dens$t==t,1:n_])
      #print(st_num)
      taum_log_sum = apply(apply(array(1:n_),1,function(n) log(taum[t,n,st_num[,n]])),1,sum)
      sum(exp(taum_log_sum) * log_dens[log_dens$t==t,'logdens'])
    }

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
  # st_num = t(sapply(log_dens[[1]]$nslist,function(x) trans(x,n_,k_)))
  st_num = log_dens[log_dens$t==1,1:n_]
  for(n in 1:n_)
  {
    taul_prod = foreach(i=1:k_,.combine = cbind) %do%
      {
        nsi = which(st_num[,n]==i)
        if(length(nsi)==0)
          0
        else
        {
          temp_sum = foreach(n2=(1:n_)[-n],.combine = "+") %do%
            {
              log(tau$taul[n2,st_num[nsi,n2]])
            }
          a = exp(temp_sum-max(temp_sum))
          a = a/sum(a)
          a
        }
      }

    log_taul = sapply(1:k_,function(i) 
      {
        if(all(taul_prod[,i]==0))
          -Inf
        else
          sum(taul_prod[,i]*log_dens[(log_dens$t==1)&(log_dens[,n]==i),'logdens']) + log(var_pi[i])
      })
    
    log_taul = exp(log_taul-max(log_taul))
    log_taul = log_taul/sum(log_taul)
    newtaul[n,] = log_taul
  }

  #update taut(t>1)
  #registerDoParallel(cl)
  df_taut = foreach(t = 2:t_, .combine = cbind) %dopar%
  {
    # st_num = t(sapply(log_dens[[t]]$nslist,function(x) trans(x,n_,k_)))
    df = log_dens[log_dens$t==t,]
    st_num = df[,1:n_]
    ntemp = foreach(n=1:n_,.combine = cbind) %do%
    {
      temp_prod = foreach(i=1:k_,.combine = cbind) %do%
      {
        nsi = which(st_num[,n]==i)
        if(length(nsi)==0)
          0
        else
        {
          logsum = foreach(n2=iter((1:n_)[-n]),.combine = "+") %do%
            {
              log(tau$taum[t,n2,st_num[nsi,n2]])
            }
          a= exp(logsum - max(logsum))
          a=a/sum(a)
        }
      }
      log_taut = sapply(1:k_,function(j) 
      {
        if(all(temp_prod[,j]==0))
          rep(-Inf,k_)
        else
          sum(temp_prod[,j]*df[df[,n]==j,'logdens']) + log(var_A[,j])
      })
      temp = apply(log_taut,1,function(x) exp(x-max(x)))
      temp = t(apply(temp,2,function(x) x/sum(x)))
      temp
    }
    ntemp
  }
  #print(dim(df_taut))
  df_taut = cbind(matrix(0,k_,k_*n_),df_taut)
  #print(dim(df_taut))
  
  newtaut = array(df_taut,dim=c(k_,k_,n_,t_))
  
  newtaul = correct_taul(newtaul)
  newtaut = correct_taut(newtaut)


  # #print(apply(newtaul,1,sum))
  # #print(apply(newtaut,c(1,2,3),sum))
  #
  newtaum = update_taum(newtaul,newtaut)
  res=list(taul=newtaul,taut=newtaut,taum=newtaum)
  return(res)
}

#baseline: fixed effect
app_update_fe<-function(obs,covar,var_theta,var_gama)
{
  #var_gama: df with col[1:N,'t','gama']
  #beta:P*Q
  n_ = dim(obs)[1]
  q_ = dim(obs)[2]
  t_ = dim(obs)[3]
  k_ = K
  
  revD = ginv(diag(n_)-var_theta$W)
  
  registerDoParallel(cl)
  mid=foreach(t=1:t_,.combine = "+") %dopar%
    {
      xb = kronecker(diag(q_),matrix(covar[,t,],nrow=n_)) %*% as.vector(var_theta$beta)
      nslist = as.matrix(var_gama[(var_gama$t==t)&(var_gama$gama!=0),c(1:n_,'gama')])
      term2 = foreach(i=1:nrow(nslist),.combine = "+") %do%
        {
          nslist[i,'gama']*(as.vector(obs[,,t]) - xb - kronecker(diag(q_),revD)%*%as.vector(var_theta$delta[nslist[i,1:n_],]))
        }
      term2
    }#PQ*(PQ+1)

  newfe=matrix(mid/t_,n_,q_)
  return(newfe)
}

#max beta
app_update_B<-function(obs,covar,var_theta,var_gama,is_network)
{
  #var_gama: df with col[1:N,'t','gama']
  #beta:P*Q
  n_ = dim(obs)[1]
  q_ = dim(obs)[2]
  t_ = dim(obs)[3]
  p_ = dim(covar)[3]
  k_ = K

  if(!is_network)
    obs = obs - rep(var_theta$fe,t_)

  bigsgm = kronecker(ginv(var_theta$sigma),diag(n_))  #NQ*NQ
  revD = ginv(diag(n_)-var_theta$W)
  
  registerDoParallel(cl)
  mid=foreach(t=1:t_,.combine = "+") %dopar%
    {
      xtbar = kronecker(diag(q_),matrix(covar[,t,],nrow=n_))
      term1 = (t(xtbar) %*% bigsgm %*% xtbar)   #PQ*PQ
      nslist = as.matrix(var_gama[(var_gama$t==t)&(var_gama$gama!=0),c(1:n_,'gama')])
      term2 = foreach(i=1:nrow(nslist),.combine = "+") %do%
        {
          nslist[i,'gama']*t(xtbar)%*%bigsgm%*%(as.vector(obs[,,t]) - kronecker(diag(q_),revD)%*%as.vector(var_theta$delta[nslist[i,1:n_],]))
        }
      cbind(term1,term2)
    }#PQ*(PQ+1)
  term1 = as.matrix(mid[,-ncol(mid)])
  term2 = as.matrix(mid[,ncol(mid)])
  newbeta=matrix(ginv(term1)%*%term2,p_,q_)
  return(newbeta)
}

#max sigma
app_update_sigma<-function(obs,covar,var_theta,var_gama,is_network)
{
  n_ = dim(obs)[1]
  q_ = dim(obs)[2]
  t_ = dim(obs)[3]
  p_ = dim(covar)[3]
  k_ = K
  
  revD = ginv(diag(n_)-var_theta$W)
  #eti = array(0.0,dim=c(t_,ns,n_*q_))
  termphi = array(0,dim=c(q_,q_))
  #temp = matrix(0,n_*q_,n_*q_)
  
  if(!is_network)
    obs = obs - rep(var_theta$fe,t_)
  
  a_tilta = obs - array(apply(covar,2,function(x) x%*%var_theta$beta),dim=c(n_,q_,t_))
  
  registerDoParallel(cl)
  temp = foreach(t=1:t_,.combine = "+") %dopar%
    {
      yb=as.vector(a_tilta[,,t])
      nslist=as.matrix(var_gama[(var_gama$t==t)&(var_gama$gama!=0),c(1:n_,'gama')])
      mid = foreach(i=1:nrow(nslist),.combine = "+") %do%
        {
          eti = yb - kronecker(diag(q_),revD) %*% as.vector(var_theta$delta[nslist[i,1:n_],])
          nslist[i,'gama']*crossprod(t(eti))
        }
      mid
    }
  
  for(i in 1:q_)
  {
    for(j in 1:q_)
    {
      termphi[i,j] = tr(temp[((i-1)*n_+1):(i*n_),((j-1)*n_+1):(j*n_)])
    }
  }
  #print(termphi)
  newsigma = termphi / (n_*t_)
  if(all(newsigma>KAPPA))
    newsigma = KAPPA
  #newsigma = ifelse(newsigma>KAPPA,KAPPA,newsigma)
  #print(newsigma)
  return(newsigma)
}

#max lambda
app_update_lambda<-function(obs,covar,var_theta,var_gama,is_network)
{
  n_ = dim(obs)[1]
  q_ = dim(obs)[2]
  t_ = dim(obs)[3]
  p_ = dim(covar)[3]
  k_ = K
  
  if(!is_network)
    obs = obs - rep(var_theta$fe,t_)
  bigsgm = kronecker(ginv(var_theta$sigma),diag(n_))  #NQ*NQ
  #zt = matrix(0,ns,n_*q_)
  revD = ginv(diag(n_)-var_theta$W)
  registerDoParallel(cl)
  mid = foreach(t=1:t_,.combine = "+") %dopar%
    {
      xb = as.vector(obs[,,t])- kronecker(diag(q_),matrix(covar[,t,],nrow=n_)) %*% as.vector(var_theta$beta)
      nslist = as.matrix(var_gama[(var_gama$t==t)&(var_gama$gama!=0),c(1:n_,'gama')])
      mmid = foreach(i=1:nrow(nslist),.combine = "+") %do%
        {
          st = matrix(0,n_,k_)
          for(j in 1:k_)
          {
            st[which(nslist[i,1:n_]==j),j] = 1
          }
          rts = kronecker(diag(q_),revD%*% st)
          term1 = nslist[i,'gama']*t(rts)%*%bigsgm%*%rts    #QK*QK
          term2 = nslist[i,'gama']*t(rts)%*%bigsgm%*%(xb)  #QK*1
          cbind(term1,term2)
        }
      mmid
    }
  term1 = mid[,-ncol(mid)]
  term2 = mid[,ncol(mid)]
  newdelta = matrix(ginv(term1-rho*diag(q_*k_))%*%term2,k_,q_)

  return(newdelta)
}

#max network
app_update_W<-function(obs,covar,var_theta,var_gama,adj)
{
  n_ = dim(obs)[1]
  q_ = dim(obs)[2]
  t_ = dim(obs)[3]
  p_ = dim(covar)[3]
  k_ = K
  
  a_tilta = obs - array(apply(covar,2,function(x) x%*%var_theta$beta),dim=c(n_,q_,t_))
  
  baterm = matrix(0,n_*n_,n_*n_)
  cterm = as.vector(rep(0,n_*n_))
  
  
  for(t in 1:t_)
  {
    nslist = as.matrix(var_gama[(var_gama$t==t)&(var_gama$gama!=0),c(1:n_,'gama')])
    for(i in 1:nrow(nslist))
    {
      tempb = matrix(0,n_*n_,n_*n_)
      tempc = as.vector(rep(0,n_*n_))
      
      b_tilta = var_theta$delta[nslist[i,1:n_],]   #N*Q
      for(q1 in 1:q_)
      {
        for(q2 in q1:q_)
        {
          b_ = 0.5 * (b_tilta[,q1]%*%t(b_tilta[,q2]) + b_tilta[,q2]%*%t(b_tilta[,q1]))
          tempb = tempb + kronecker(b_,ginv(diag(rep(var_theta$sigma[q1,q2],n_))))
          tempc = tempc + as.vector(ginv(diag(rep(var_theta$sigma[q1,q2],n_)))%*%a_tilta[,q1,t]%*%t(b_tilta[,q2]))
        }
      }
      tempb = nslist[i,'gama']*tempb
      baterm = baterm + tempb
      tempc = nslist[i,'gama']*tempc
      cterm = cterm + tempc
    }
  }
  vecD = ginv(baterm)%*%cterm
  vecD = matrix(vecD,n_,n_)
  newW = (diag(n_) - ginv(vecD)) * adj
  newW = ifelse(newW<0,0,newW)
  newW = ifelse(newW>1,1,newW)

  return(newW)
}

app_M_step<-function(obs,covar,var_theta,var_tau,nslist,is_network)
{
  n_ = dim(obs)[1]
  t_ = dim(obs)[3]
  #nslist = lapply(log_dens,function(x) x[,1:n_]) list of length t_
  var_theta$var_pi = update_pi(var_tau)
  var_theta$var_A = update_A(var_tau)
  
  #list of length t_, each element is num of state and gama value
  if(K>1)
  {
    gama = foreach(t=1:t_,.combine=rbind) %do%
      {
        ns = nslist[nslist$t==t,c(1:n_,'t')]
        value = apply(ns[,1:n_],1,function(x)
          {
            foreach(j = 1:K,.combine = "+") %do%
            {
              sum(log(var_tau$taum[t,which(x==j),j]))
            }
        })
        value = exp(value-max(value))
        value = value/sum(value)
        ns$gama = value
        ns[!(ns$gama==0),]
      }
  }
  if(K==1)
  {
    gama = cbind(matrix(1,t_,n_),(1:t_),1)
    colnames(gama) = c(1:n_,"t","gama")
  }
  
  # print(sapply(1:t_, function(t) revtrans(gama[gama$t==t,][which.is.max(gama[gama$t==t,'gama']),1:n_],n_,K) - data$lti[t]))
  
  var_theta$beta = app_update_B(obs,covar,var_theta,gama,is_network)
  var_theta$delta = app_update_lambda(obs,covar,var_theta,gama,is_network)
  # var_theta$delta = THETA$delta
  if(!is_network)
  {
    var_theta$W = matrix(0,n_,n_)
    var_theta$fe = app_update_fe(obs,covar,var_theta,gama)
  }
  else
  {
    if(n_>1)
    {
      var_theta$W = app_update_W(obs,covar,var_theta,gama,Z)
    }
    else
      var_theta$W = 1.0
  }
  
  # var_theta$W = THETA$W
  
  # var_theta$beta = THETA$beta
  var_theta$sigma = app_update_sigma(obs,covar,var_theta,gama,is_network)
  # var_theta$sigma = THETA$sigma
  return(var_theta)
}

app_inference_s<-function(var_tau)
{
  t_=dim(var_tau$taum)[1]
  n_=dim(var_tau$taum)[2]
  k_ = dim(var_tau$taum)[3]
  
  inf_lt = array(apply(var_tau$taum,1,function(x) {x = matrix(x,nrow=n_)
  apply(x,1,function(y) which.is.max(y))}),dim=c(n_,t_))   #N*T
  
  inf_lti = apply(inf_lt,2,function(x) revtrans(x,n_,k_))
  
  res = list(lt=inf_lt,lti=inf_lti)
  return(res)
}

app_VEM_algo<-function(iteration,presision,obs,covar,limit,is_network=TRUE)
{
  loglik = rep(0,iteration)
  
  n_ = dim(obs)[1]
  q_ = dim(obs)[2]
  t_ = dim(obs)[3]
  p_ = dim(covar)[3]
  
  #init
  var_tau = init_tau(n_,K,t_)
  var_theta = init_theta(n_,K,q_,p_,Z,is_network)
  
  #var_theta = THETA
  clusterExport(cl,varlist = ls(environment(app_po)),envir = environment())
  clusterExport(cl,varlist = ls(environment(app_M_step)),envir = environment())
  clusterExport(cl,varlist = ls(environment(app_update_tau)),envir = environment())
  epsilon=c()
  old_theta = var_theta
  for(g in 1:iteration)
  {
    registerDoParallel(cl)
    log_dens = foreach(t=1:t_,.combine = rbind) %dopar%
      {
        temp = app_po(obs[,,t],covar[,t,],var_theta,limit,log=T,is_network=is_network)
        temp$t = t
        temp
      }
    #logdens: T*ns(t) x N+2
    #cat("logdens\n")
    #print(sapply((1:t_),function(t) revtrans(log_dens[log_dens$t==t,][which.is.max(log_dens[log_dens$t==t,'logdens']),1:n_],n_,K)- data$lti[t]) )
    loglik[g] = app_elbo_log(obs,covar,var_tau,var_theta,log_dens)
    
    cat(g,",loglik,",loglik[g],"\n")
    if(g>1)
    {
      ep = abs(loglik[g]-loglik[g-1])
      epsilon = c(epsilon, ep)
      if(loglik[g]<loglik[g-1])
      {
        cat("Likelihood decrease!\n")
      }
    }
    
    #VE
    if(K>1)
    {
      var_tau = app_update_tau(log_dens,init_tau(n_,K,t_),var_theta$var_pi,var_theta$var_A)
      # cat(g,",after VE,",app_elbo_log(obs,covar,var_tau,var_theta,log_dens),"\n")
      # print(sapply(1:t_, function(t) revtrans(apply(var_tau$taum[t,,],1,which.is.max),n_,K) - data$lti[t]))
    }
    # print(sum(abs(apply(var_tau$taum,1,function(x) apply(x,1,function(y) which.is.max(y))) - apply(data$lt[,,1:sep],3,function(x) apply(x,1,function(y) which.is.max(y))))))
    #M
    # nslist = lapply(log_dens,function(x) x[,1:n_])
    var_theta = app_M_step(obs,covar,var_theta,var_tau,log_dens,is_network)
    
    cvg = max(unlist(mapply(function(x,y) abs(x-y),old_theta,var_theta)))
    old_theta = var_theta
    
    if(g>1 && cvg <= presision) break()
    
  }
  
  #inference
  s_star = app_inference_s(var_tau)
  res=list(theta_est=var_theta, converge = loglik, tau=var_tau, infs = s_star)
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

