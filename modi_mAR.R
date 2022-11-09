library(mixAR)

fit_mixVARreg <- function(x, y, mixVARmodel, EMinit, ...){
  stop("There is currently no default for this funciton")
}
setGeneric("fit_mixVARreg")

setMethod("fit_mixVARreg", c(x = "matrix", y = "matrix", mixVARmodel = "MixVAR",
                            EMinit = "missing"),
          function(x, y, mixVARmodel, EMinit, ...){
            mixVARreg(x, y, mixVARmodel, ...)
          }
)

mixVARreg <- function(x, y, mixVARmodel, tol = 1e-6, niter = 200){
  
  #x:(N*Q)*T
  #y:N*T*P
  n_ <- dim(x)[1]
  q_ <- dim(x)[2]
  t_ <- dim(x)[3]
  p_ <- dim(y)[3]
  
  
  if(t_ != ncol(y)){ 
    stop("x and y must contain the same number of observations")
  } ## Need a first iteration to initialise all variables
  
  # reg_model <- lm(data.frame(x, y))
  temp <- foreach(t=1:t_,.combine = "+") %do%
  {
    term1 = crossprod(x[,,t],y[,t,])  #q*p
    term2 = crossprod(y[,t,])     #p*p
    rbind(term1,term2)    #(q+p)*p
  }
  
  
  #beta: (P*Q)
  reg_parameters <- beta_old <- t(temp[1:q_,] %*% ginv(temp[(q_+1):nrow(temp),]))
  
  
  #residuals: T*(NQ)
  reg_residuals <- t(array(x,dim=c(n_*q_,t_)) - apply(y,2,function(x) x%*%reg_parameters))
  
  MVAR_mod <- MVAR_old <- mixVARfit(reg_residuals,mixVARmodel, fix="shift")$model
 
  pk <- MVAR_mod@order
  g <- length(pk)
  p <- max(pk)
  ufit <- matrix(0,t_- p,n_*q_)
  for(k in 1:g){
    for(t in 1:(t_ - p)){
      for(u in 1:pk[k])
      {
        ufit[t,] <- ufit[t,] + MVAR_mod@prob[k] * 
          (MVAR_mod@arcoef@a[[k]][,,u] %*% reg_residuals[(t+p-u),])
      }
    }
  }
  xstar <- x ;
  xstar[,,(p+1):t_] <- x[,,(p+1):t_] - array(t(ufit),dim=c(n_,q_,t_-p))
  diff <- 1
  i <- 1
  ## @Davide: not urgent but beware - 
  ##    (b) 'diff' is the sum of quantities, which are potentially on different
  ##         scales, whi may cause trouble.
  ## 24/04/2020 changed this, now stops if all single differences are smaller than tol

  while(diff > tol){  ## From iteration 2, basically repeats the code above
    i <- i+1
    temp <- foreach(t=1:t_,.combine = "+") %do%
      {
        term1 = crossprod(xstar[,,t],y[,t,])  #q*p
        term2 = crossprod(y[,t,]) #p*p
        rbind(term1,term2)
      }#(q+p)*p
    
    
    #beta: (P*Q)
    reg_parameters <- beta_old <- t(temp[1:q_,] %*% ginv(temp[(q_+1):nrow(temp),]))
    #residuals: T*(NQ)
    reg_residuals <- t(array(x,dim=c(n_*q_,t_)) - apply(y,2,function(x) x%*%reg_parameters))
    MVAR_mod <- MVAR_old <- mixVARfit(reg_residuals,mixVARmodel, fix="shift")$model
    # MAR_mod <- MAR_old <- fit_mixAR(reg_residuals, mixARmodel, fix = "shift")$model
    pk <- MVAR_mod@order
    g <- length(pk)
    p <- max(pk)
    ufit <- matrix(0,t_- p,n_*q_)
    for(k in 1:g){
      for(t in 1:(t_ - p)){
        for(u in 1:pk[k])
        {
          ufit[t,] <- ufit[t,] + MVAR_mod@prob[k] * 
            (MVAR_mod@arcoef@a[[k]][,,u] %*% reg_residuals[(t+p-u),])
        }
      }
    }
    
    # diff <- sum(abs(reg_parameters / beta_old - 1)) +
    #                sum(abs(MAR_mod@prob / MAR_old@prob - 1)) +
    #                sum(abs(MAR_mod@scale / MAR_old@scale - 1)) +
    #                sum(abs(unlist(MAR_mod@arcoef@a) / unlist(MAR_old@arcoef@a) -1))
    diff <- c(abs(reg_parameters - beta_old),
              abs(MVAR_mod@prob - MVAR_old@prob),
              abs(MVAR_mod@vcov - MVAR_old@vcov),
              abs(unlist(MVAR_mod@arcoef@a) - unlist(MVAR_old@arcoef@a)))
    diff <- ifelse(all(diff < tol), 0, 1)
    beta_old <- reg_parameters
    MVAR_old <- MVAR_mod
    xstar <- x ;
    xstar[,,(p+1):t_] <- x[,,(p+1):t_] - array(t(ufit),dim=c(n_,q_,t_-p))
    if(i >= niter) 
      break
  }
  list(reg = reg_parameters, mixVARmodel = MVAR_mod, niter = i, convergence = diff)
}

cond_loglikV <- function(model, x, index){
  if(ncol(x) > nrow(x)) x <- t(x)
  mult <- ncol(x)
  n <- nrow(x)
  g <- length(model@prob)
  pk <- sapply(model@arcoef@a, function(x) dim(x)[3])
  p <- max(pk) 
  
  if(missing(index))
    index <- (p+1):nrow(x)
  
  yhat <- array(matrix(0, ncol=ncol(x), nrow=length(index)), 
                dim=c(length(index), ncol(x), g))
  dens <- matrix(0, ncol=g, nrow=length(index))
  
  for(k in 1:g){
    armat <- model@arcoef[k]
    
    for(i in index){
      yhat[(i + 1 - min(index)), ,k] <- model@shift[,k] + 
        rowSums(t(t(armat) * c( t(x[(i-1):(i-pk[k]),]) )))
      
      dens[i + 1 - min(index), ] <- model@prob[k] * 
        dMVNorm(x[i, ], yhat[(i + 1 - min(index)), , k], model@vcov[,,k], log = TRUE)
    }
    
  }
  sum(dens)
}

mixVARfit <- function(y, model, fix = FALSE, tol = 10^-6, verbose = FALSE){
  verbose <- verbose && interactive()
  
  prob <- model@prob
  Scale <- model@vcov
  arcoef <- model@arcoef@a
  shift <- model@shift
  one_or_zero <- ifelse(identical(fix, "shift"), FALSE, TRUE)
  if(ncol(y) > nrow(y)) y <- t(y)
  mult <- ncol(y)
  n <- nrow(y)
  g <- length(prob)
  pk <- sapply(arcoef, function(x) dim(x)[3])
  p <- max(pk) 
  e <- array(matrix(0, ncol=ncol(y), nrow=(n-p)), dim=c(n-p, ncol(y), g))
  diff <- 1
  count <- 0
  
  while(diff>tol){
    count <- count +1; if(count%%25 == 0) {
      loglik <- cond_loglikV(new("MixVARGaussian", prob=prob, vcov=Scale, arcoef=arcoef,
                                 shift=shift), y)
      if(verbose)
        cat("niter:", count, "\tvallogf:", loglik, "\n")
    }
    ##Calculate residuals. Each residual is a mult x 1 vector
    
    for(k in 1:g){
      armat <- matrix(arcoef[[k]], nrow=ncol(y), ncol=(ncol(y) * pk[k]))
      
      for(t in 1:(n-p)){
        e[t, ,k] <- y[t+p, ] - shift[,k] - 
          rowSums(t(t(armat) * c( t(y[(t+p-1):(t+p-pk[k]),]) )))
      }
      
    }
    dens <- matrix(ncol=g, nrow=(n-p))
    for(k in 1:g){
      d <- sqrt(det(Scale[,,k]))
      inv <- ginv(Scale[,,k])
      for(t in 1:(n-p)){
        dens[t, k] <- prob[k] / d * 
          exp(-0.5*e[t,,k] %*% inv %*% e[t,,k]) 
      }
    }
    dens <- pmax(dens,.Machine$double.xmin)
    dens <- pmin(dens,.Machine$double.xmax)
    #cat(sum(is.na(dens)),"\t")
    #cat(dens[1,],max(dens[1,]),dens[1,]-max(dens[1,]))
    dens <- t(apply(dens,1,function(x) x-max(x)))
    dens[which(rowSums(dens)==0),] = prob
    dens_sum <- ifelse(rowSums(dens)==0,0,1/pmax(rowSums(dens),.Machine$double.xmin))
    #cat(length(which(dens_sum==0)))
    tau <- dens*dens_sum
    prob_new <- colMeans(tau)
    # cat(sum(is.na(dens)),sum(is.na(tau)),sum(is.na(prob_new)),"\n")
    
    Scale_new <- Scale 
    for(k in 1:g){
      sig_sum <-matrix(0, nrow=mult, ncol=mult)
      for(t in 1:(n-p)){
        sig_sum <- sig_sum + tau[t, k] * e[t,,k] %*% t(e[t,,k])
      }
      Scale_new[,,k] <- sig_sum / sum(tau[,k])
      Scale_new[,,k] <- pmax(Scale[,,k],1e-20)
      if(!isSymmetric(Scale_new[,,k]))
        Scale_new[,,k][lower.tri(Scale_new[,,k])] <- t(Scale_new[,,k])[lower.tri(Scale_new[,,k])]
    }
    
    AR_new <- arcoef
    if(one_or_zero) shift_new <- shift
    for(k in 1:g){
      sum1 <- matrix(0,ncol=(mult * pk[k]+one_or_zero), nrow=(mult*pk[k]+one_or_zero))
      sum2 <- matrix(0,nrow=(mult * pk[k]+one_or_zero), ncol=mult)
      for(t in 1:(n-p)){
        x <- c(if(one_or_zero) 1, c(t(y[(t+p-1):(t+p-pk[k]),])))
        sum1 <- sum1 + tau[t,k] * x %*% t(x)
        sum2 <- sum2 + tau[t,k] * x %*% t(y[t+p,])
      }
      theta <- t(solve(sum1) %*% sum2)
      if(one_or_zero) {shift_new[,k] <- theta[,1] ; theta <- theta[,-1]}
      for(i in 1:pk[k]){
        AR_new[[k]][,,i] <- theta[ ,1:mult] ; theta <- theta[ ,-c(1:mult)]  
      }  
    }
    
    #cat(sum(is.na(prob_new)),sum(is.na(unlist(Scale_new))),sum(is.na(unlist(AR_new))),"\n")
    
    diff <- sum(abs(prob_new-prob)) + 
      sum(abs(unlist(Scale) - unlist(Scale_new))) +
      sum(abs(unlist(arcoef) - unlist(AR_new)))
    
    
    Scale  <- Scale_new
    arcoef <- AR_new
    if(one_or_zero) shift  <- shift_new
    prob   <- prob_new
    if(count == 200) break
    
  }
  mod <- new("MixVARGaussian", prob=prob, shift=shift, vcov=Scale, arcoef=arcoef)
  list(model = mod, vallogf = cond_loglikV(mod, y))
}


