#' @title Frechet Ridge Selection Operator for Global Frechet Regression with SPD Matrix Output
#' 
#' @description This code implements the Frechet Ridge Selection Operator
#' on Symmetic Positive Definite Matrix type data equipped with the Cholesky decomposition distance 
#' 
#' @details The outline is as follows:
#' 
#' 1. Define the data generation function utilized for example 5.3.3.
#' 2. Define the global objective function for Individually penalized ridge regression 
#' 3. Define the function to iterate over which 
#'     a. generates validation/tuning data (seed is 2019),
#'     b. sets the grid for tau
#'     c. generates training data (seed is iteration #),
#'     d. fits an individually penalized global Frechet regression for each tau,
#'     e. calculates in sample error (on training data),
#'     f. calculates out of sample error (on validation/tuning data),
#'     g. collects the best tau from these results and the corresponding lambda hats,
#'     
#'     h. refits global Frechet regression on training data with only predictors whose lambda hat > 0.01,
#'     i. calculates in sample error
#'     j. calculates out of sample error  
#'     j. collects the best tau from these refitted model results and the lambda hats which produced it
#' 4. Iterate 100 times 

rm(list=ls())
###################################################
library(psych)
library(matlab)
library(ramify)
source("IndivRidgeGloCovReg.R") 
source("IndivRidgeGFRCovCholesky.R")
source("lambdaGlobalcoordesc.R")
source("GFRCovCholesky.R")
source("GloCovReg.R")
source("cholesky_error.R")
source("IndivRidgeGFRCovCholesky_v2.R")
Rcpp::sourceCpp("IndivRidgeGFRCovCholesky_vC.cpp")
####################################################

# Generates data from example 5.3.2
gendata=function(n,p,m,r,seed) {
  #' @param n: sample size 
  #' @param p: number of predictors
  #' @param m: dimension of SPD matrix outputs
  #' @param r: correlation between simulated predictors
  #' @param seed: seed set for the random generation
  
  # Recall the generation of Y=Min in example 5.3.2: Min = A^tA
  #                                                  A = (mu + sigma)I + (sigma)U                                                  
  
  x=2*pnorm(matrix(rnorm(n*p),n,p)%*%chol(r^(as.matrix(dist(1:p)))))
  
  M <- array(0,c(m,m,n))
  Vmu0=3
  Vsigma0=3 
  Vbeta=2 
  Vgamma=3
  
  Vnu1=1
  Vnu2=2  
  
  set.seed(seed)
  Vmu=Vbeta*(x[,1] + (x[,3]) )+rnorm(n, Vmu0,sqrt(Vnu1))
  GamP1=(Vsigma0+Vgamma*(x[,7] + x[,9]+ x[,5]))^2/Vnu2
  GamP2=Vnu2/(Vsigma0+Vgamma*(x[,7] + x[,9]+ x[,5]))
  I =diag(m)
  off_diag = matlab::ones(m)
  off_diag[lower.tri(off_diag, diag = TRUE)] <- 0
  
  for (j in 1:n) {
    set.seed(seed)
    Vsigma=rgamma(1,shape=GamP1[j], scale=GamP2[j]) 
    A = (Vmu[j]+Vsigma)*I + Vsigma*off_diag
    aux<-t(A)%*%A
    M[,,j]<-aux
  }
  
  return(list(x=x, Min=M))
}

GlobalObj4h <- function(x, Min, lambda) {
  Mres <- IndivRidgeGFRCovCholesky_vC(x = x, M = Min, xout = x, lambda = lambda, 
                                      optns = list(metric="cholesky"))$Mout
  error <- vector("double", nrow(x))
  for(i in 1:nrow(x)){
    error[i] <- cholesky_error_C(Mres[[i]], Min[,,i])^2
  }
  return(sum(error))
}

mymainfun <- function(irep){
  n=200
  p=10
  m=5
  r=0.5
  
  seed = 2019
  set.seed(2019)
  
  # Generate independent tuning/validation data
  wutest=gendata(n,p,m,r,seed)
  xtest=wutest$x
  Mintest=wutest$Min
  
  # Grid for tau = c(.5,1,1.5, 2,...,30)
  tauall=seq(1,80,1)*.5
  
  # Generate training data
  set.seed(irep)
  wu=gendata(n,p,m,r, irep)
  x=wu$x
  Min=wu$Min
  
  # Find optimal lambda given each tau
  hopt=lambdaGlobalcoordesc(x=x,Min=Min, tau=tauall[1],lambdainit=NULL)
  hopt=cbind(hopt, lambdaGlobalcoordesc(x,Min, tauall[2], lambdainit=hopt))
  for(k in 3:length(tauall)){
    hopt=cbind(hopt,lambdaGlobalcoordesc(x,Min, tauall[k], lambdainit=hopt[,k-1]) )
    print(hopt)
  }
  
  # Calculate error on training data 
  finalobj=NULL
  for (k in 1:length(tauall)) {
    Mres=IndivRidgeGloCovReg(x=x, M=Min, xout=x, lambda=hopt[,k], optns=list(metric="cholesky"))$Mout
    error = cholesky_error(Mres[[1]],Min[,,1])^2
    for(i in 2:nrow(x)){
      error = rbind(error,cholesky_error(Mres[[i]],Min[,,i])^2)
    }
    finalobj[k]=sum(error)
  }
  
  # Refit model on training data with only good predictors for each tau
  refit.finalobj=NULL
  for (k in 1:length(tauall)) {
    ind=which(hopt[,k]>1e-2)
    Mres=GloCovReg(x=as.matrix(x[,ind]), M = Min, xout=as.matrix(x[,ind]), optns=list(metric="cholesky"))$Mout
    error = cholesky_error(Mres[[1]],Min[,,1])^2
    for(i in 2:nrow(x)){
      error = rbind(error, cholesky_error(Mres[[i]],Min[,,i])^2)
    }
    refit.finalobj[k]=sum(error)
  }
  
  # Calculate error on independent validation data
  indep.finalobj=NULL
  for (k in 1:length(tauall)) {
    Mres=IndivRidgeGloCovReg(x=x, M=Min, xout=xtest, lambda=hopt[,k], optns=list(metric="cholesky"))$Mout
    error = cholesky_error(Mres[[1]],Mintest[,,1])^2
    for(i in 2:nrow(xtest)){
      error = rbind(error, cholesky_error(Mres[[i]],Mintest[,,i])^2)
    }
    indep.finalobj[k]=sum(error)
  }
  
  # Grab best model based on error on validation data
  indep.ind=which.min(indep.finalobj)
  indep.finallambda=hopt[,indep.ind]
  
  # Refit on model on independent validation data with only good predictors for each tau
  indep.refit.finalobj=NULL
  for (k in 1:length(tauall)) {
    ind=which(hopt[,k]>1e-2)
    Mres=GloCovReg(as.matrix(x[,ind]), M=Min, xout=as.matrix(xtest[,ind]), optns=list(metric="cholesky"))$Mout
    error = cholesky_error(Mres[[1]],Mintest[,,1])^2
    for(i in 2:nrow(xtest)){
      error = rbind(error, cholesky_error(Mres[[i]],Mintest[,,i])^2)
    }
    indep.refit.finalobj[k]=sum(error)
  }
  
  # Grab refitted model with best error in independent validation error
  indep.refit.ind=which.min(indep.refit.finalobj)
  indep.refit.finallambda=hopt[,indep.refit.ind]
  
  # Grab results
  result=list(hopt=hopt,                                      # all lambdas given tau
              finalobj=finalobj,                              # error on training data
              indep.ind=indep.ind,                            # best tau given error on validation data
              indep.finallambda=indep.finallambda,            # lambda given best tau from validation data
              refit.finalobj=refit.finalobj,                  # refitted error on training data
              indep.refit.finalobj=indep.refit.finalobj,      # refitted error on validation data
              indep.refit.ind=indep.refit.ind,                # refitted best tau on validation data
              indep.refit.finallambda=indep.refit.finallambda) # refitted lambda given best tau from validation data
  
  return(result)
  
} 

seed.start=0

result = list()

# Prepare parallel processing
library(foreach)
library(quadprog)
library(iterators)
library(parallel)
library(doParallel)
registerDoParallel(25)

# repeat main function 100 times
allresult=foreach(irep=1:100, .combine=c) %dopar% mymainfun(irep + seed.start)

save(file=sprintf("4Pcholeskycorr533simulationseed%d.Rdata",seed.start), list=ls())

stopImplicitCluster()
