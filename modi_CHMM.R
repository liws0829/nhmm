library(CHMM)
coupledHMM <- function(X, Y, nb.states = 2, S = NULL, omega.list = c(0.3,0.5,0.7,.9),
                       var.equal = TRUE, exact = FALSE, meth.init = "mclust", viterbi = TRUE, itmax = 10,
                       threshold = 1e-4){
  # dim(X):t_*n_ -> t_*n_[q]
  # dim(Y): n_*t_*p_
  if (is.null(S))
    stop("Argument S must be specified.")
  nbI <- ncol(X)
  # inference ---------------------------
  RSS <- Inf
  if(exact == FALSE){
    resApprox <- apply(as.matrix(omega.list),c(1),FUN=function(y) CHMM_VEM(X, Y, nb.states, S , y, meth.init, var.equal, itmax, threshold))
    RSS <- unlist(lapply(resApprox,FUN = function(x) x$'RSS'))
    mod.sel <- resApprox[[which.min(RSS)]]
    omega.sel <- omega.list[which.min(RSS)]

    if (viterbi == TRUE){
      status <- apply(as.matrix(1:nbI),c(1),FUN=function(y) viterbi_algo(mod.sel$emisPrW[[y]], mod.sel$transPr, mod.sel$initPr))
    }else{
      status <-  apply(as.matrix(1:nbI),c(1),FUN=function(y) apply(mod.sel$postPr[[y]], 1, which.max))
    }
    colnames(status) <- colnames(X)
  }else{

    resExact <- apply(as.matrix(omega.list),c(1),FUN=function(y)
      CHMM_EM(X, nb.states, S , y, meth.init, var.equal, itmax, threshold))
    RSS <- unlist(lapply(resExact,FUN = function(x) x$'RSS'))
    mod.sel <- resExact[[which.min(RSS)]]
    omega.sel <- omega.list[which.min(RSS)]

    if (viterbi == TRUE){
      stsGb <- viterbi_algo(mod.sel$'emisGb', mod.sel$'transGb', mod.sel$'initGb')
      status <- mod.sel$ID.K[stsGb, ]
    }else{
      status <-  apply(as.matrix(1:nbI),c(1),FUN=function(y) apply(mod.sel$postPr[[y]], 1, which.max))
    }

    colnames(status) <- colnames(X)
  }#end if
  return(list(omega = omega.sel, model = mod.sel, status = status, RSS.omega = data.frame(RSS = RSS, omega = omega.list)))
}

CHMM_VEM <- function(X, Y, nb.states, S = NULL, omega = 0.7, meth.init = "mclust", var.equal = TRUE, itmax = 5e2,
                     threshold = 1e-4){
  if (is.null(S))
    stop("Argument S must be specified.")
  nbT <- nrow(X)
  nbI <- ncol(X)
  #dim(X) t_*(n_*q_)
  # initialisation   ------------------------------------------------------
  # if (is.null(init.esAvg)) {
  res.init <- init.VEM(X, Y, nb.states, meth.init, var.equal, nbI, nbT)
  esBeta <- res.init$esBeta
  esAvg <- res.init$esAvg
  esVar <- res.init$esVar
  transPr <- res.init$transPr
  initPr <- res.init$initPr
  postPr.last <- res.init$postPr

  # Variational E-M algorithm    ------------------------------------------------------
  Old.param = c(esBeta, esAvg, esVar, as.vector(transPr))
  for (iter in 1:itmax) {
    ## 1.VE-step    -----------------------------
    Ytemp = t(apply(Y,2,function(x) as.matrix(x)%*%esBeta))  #T*N
    postPr.list <- list()
    emisPr.list <- list()
    emisPrW.list <- list()
    trsTmp <- matrix(0, nb.states, nb.states)
    emisPr <- matrix(0, nbT, nb.states)
    for (ind in 1:nbI) {
      w <- matrix(0, nbT, nb.states)
      for (ij in c(1:nbI)[-ind])
        w <- w + S[ind, ij] * (1 - postPr.last[[ij]])
      emisPr <- Emis.Gauss(X[,ind]-Ytemp[,ind], esAvg, esVar)
      emisPr.list[[ind]] <- emisPr
      emisPrW <- omega^w * emisPr
      emisPrW.list[[ind]] <- omega^w * emisPr

      # Transform matrix to vectors     -----------------------------
      emisWVec <- as.vector(emisPrW)
      emisWVec <- pmax(emisWVec, .Machine$double.xmin)
      trsVec <- as.vector(transPr)
      initTmp <- as.vector(initPr * emisPrW[1,])
      initPr <- initTmp /sum(initTmp)
      # Forward-Backward recursion      -----------------------------
      resF <- ForwardR(emisWVec, initPr, trsVec)
      resB <- BackwardR(resF$Fpr, trsVec, nb.states)

      postPr.tmp <- matrix(resB$postPr, nbT, nb.states)
      postPr.tmp <- apply(postPr.tmp, 2, pmax, .Machine$double.xmin)
      postPr <- postPr.tmp / rowSums(postPr.tmp)
      postPr.list[[ind]] <- postPr
      Fpr <- matrix(c(resF$Fpr), nbT, nb.states)
      Gpr <- matrix(c(resB$Gpr), nbT, nb.states)

      trsTmp <- trsTmp + transPr * t(Fpr[-nbT,]) %*% (postPr[-1,] / Gpr[-1,])
    }
    postPr.last <- postPr.list
    #print(postPr.list[[1]])
    ## 2.M-step      ----------------------------
    # update transPr
    trsAvg <- trsTmp / nbI
    transPr <- trsAvg / rowSums(trsAvg)

    # update initPr       ----------------------------
    init.tmp <- rep(0, nb.states)
    for(i in 1:nbI) init.tmp <- init.tmp + postPr.list[[i]][1,]
    initPr <- init.tmp / sum(init.tmp)

    # update esAvg        ----------------------------
    esAvg <- NULL
    for (r in 1:nb.states) {
      nom <- 0
      den <- 0
      for (i in 1:nbI) {
        nom <- nom + sum(postPr.list[[i]][,r] * (X[,i]-Ytemp[,i]))
        den <- den + sum(postPr.list[[i]][,r])
      }
      esAvg <- c(esAvg, nom / den)
    }
    # cat("Avg",esAvg,"\n")

    # update esBeta         ----------------------------
    nom <- matrix(0,1,dim(Y)[3])
    den <- matrix(0,dim(Y)[3],dim(Y)[3])
    for(i in 1:nbI){
      temp = (postPr.list[[i]] * (X[,i] - rep(esAvg,each=nbT))) %*% matrix(1,nb.states,1)
      nom <- nom + apply(rep(temp,dim(Y)[3]) * Y[i,,],2,sum)
      den <- den + crossprod(Y[i,,])
    }

    esBeta <- matrix(t(nom%*%ginv(den)),dim(Y)[3],1)
    # print(esBeta)

    # update esVar         ----------------------------
    Ytemp = t(apply(Y,2,function(x) as.matrix(x)%*%esBeta))
    RSS <- 0
    esVar <- NULL
    if (var.equal == TRUE) {
      nom <- 0
      for(i in 1:nbI) nom <- nom + sum((X[,i]-Ytemp[,i])^2) - 2 * esAvg %*% t(postPr.list[[i]]) %*% (X[,i] - Ytemp[,i]) + sum(esAvg^2 %*% t(postPr.list[[i]]))
      esVar <- rep(nom/(nbI*nbT),nb.states)
      RSS <- nom
      #  } else if(var.model == "V") {
    } else {
      for (r in 1:nb.states) {
        nom <- 0
        den <- 0
        for (i in 1:nbI) {
          nom <- nom + sum(postPr.list[[i]][,r] * (X[,i]-Ytemp[,i]-esAvg[r])^2)
          den <- den + sum(postPr.list[[i]][,r])
        }
        esVar <- c(esVar, nom / den)
        RSS <- RSS + nom
      }
    }

    #  order the parameters          ----------------------------
    ordAvg <- order(esAvg)
    esAvg <- sort(esAvg)
    esVar <- esVar[ordAvg]
    transPr <- transPr[ordAvg,ordAvg]

    ## 3.Stop iteration          ----------------------------
    New.param <-  c(esBeta, esAvg, esVar, as.vector(transPr))
    crit <- New.param - Old.param
    esVar <- pmax(esVar, 1e-3)
    Old.param  <- New.param

    if (iter > 1 && max(abs(crit)) <= threshold) break()
  }
  return(list(postPr = postPr.list, initPr = initPr,  transPr = transPr,  esBeta = esBeta, esAvg = esAvg, esVar = esVar,
              emisPr = emisPr.list, emisPrW = emisPrW.list, RSS = as.numeric(RSS), iterstop = iter))
}

init.VEM <- function(X, Y, nb.states, meth.init, var.equal, nbI, nbT){

  esBeta <- matrix(1,dim(Y)[3],1)
  Ytemp = apply(Y,2,function(y) as.matrix(y)%*%esBeta)
  X = X - t(Ytemp)
  # mclust initialization ---------------------------------------------------------
  if(meth.init == "mclust") {
    if (var.equal == TRUE) {
      clust <- mclust::Mclust(data = as.vector(t(X)), G = nb.states, modelNames = "E",verbose = F)
      esVar <- rep(clust$parameters$variance$sigmasq, nb.states)
    } else {
      clust <- mclust::Mclust(data = as.vector(t(X)), G = nb.states, modelNames = "V",verbose = F)
      esVar <- clust$parameters$variance$sigmasq
    }
    esAvg <- clust$parameters$mean


    # kmeans initialization ---------------------------------------------------------
  } else if(meth.init == "kmeans") {
    clust <- kmeans(x = as.vector(t(X)), centers = nb.states)
    esAvg <- sort(as.vector(clust$centers))
    if (var.equal == TRUE) {
      esVar <- rep(clust$tot.withinss / (nbT * nbI - 1), nb.states)
    } else {
      esVar <- clust$withinss / (clust$size - 1)
    }
  }
  mat.tmp <- matrix(runif(nb.states^2), ncol = nb.states) + diag(rep(50, nb.states))
  transPr <- mat.tmp / rowSums(mat.tmp)

  # initial distribution -----------------------------------------
  eigenvalues <- round(eigen(t(transPr))$values, 3)
  pos  <- which(eigenvalues == 1.000)
  nuHMM <- eigen(t(transPr))$vectors[, pos]
  initPr <- pmax(as.numeric(nuHMM / sum(nuHMM)), 0)
  initPr <- initPr / sum(initPr)

  # postPr  ------------------------------------------------------
  postPr <- list()
  for(ind in 1:nbI){
    #tau.tmp <- data.frame(matrix(runif(nbT * nb.states, 0, 1), nbT, nb.states))
    tau.tmp <- matrix(runif(nbT * nb.states, 0, 1), nbT, nb.states)
    tau.tmp <- tau.tmp / rowSums(tau.tmp)
    postPr[[ind]] <- tau.tmp
  }
  return(list(esBeta = esBeta, esAvg = esAvg, esVar = esVar, transPr = transPr, initPr = initPr,
              postPr = postPr))

}

Emis.Gauss <- function(X, esAvg, esVar){
  nbS <- length(esAvg)
  apply(as.matrix(1:nbS), 1, function(r) dnorm(X, mean = esAvg[r], sd = sqrt(esVar[r])) )
}
