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

#rm(list=ls())
###################################################
library(psych)
library(ramify)
file_folder1 <- "/gpfs/home/ly21a/projects/frechet_regression/shared_code_SPD_Matrices/"
source("IndivRidgeGloCovReg.R") 
source("IndivRidgeGFRCovCholesky.R")
source("lambdaGlobalcoordesc.R")
source("GFRCovCholesky.R")
source("GloCovReg.R")
source(paste0(file_folder1, "IndivRidgeGFRCovCholesky_v2.R"))
Rcpp::sourceCpp(paste0(file_folder1, "IndivRidgeGFRCovCholesky_vC.cpp"))
source("cholesky_error.R")
source(paste0(file_folder1, "ultility.R"))
file_folder2 <- "/gpfs/home/ly21a/projects/frechet_regression/ADNI/"
Rcpp::sourceCpp(paste0(file_folder2, "dcov_rcpp.cpp"))

####################################################


GlobalObj4h <- function(x, Min, lambda) {
  Mres <- IndivRidgeGFRCovCholesky_vC(x = x, M = Min, xout = x, lambda = lambda, 
                                      optns = list(metric="cholesky"))$Mout
  error <- vector("double", nrow(x))
  for(i in 1:nrow(x)){
    error[i] <- cholesky_error_C(Mres[[i]], Min[,,i])^2
  }
  return(sum(error))
}


foo <- function(irep){
  n=200
  p=20
  m <- 3
  r <- 0.5
  seed = 2022
  # Generate independent tuning/validation data
  wutest=gendata_SPD_cate2(n,p,m,r,seed) # same test set for different reps
  xtest=wutest$x
  Mintest=wutest$Min
  
  # Grid for tau = c(.5,1,1.5, 2,...,30)
  tauall=seq(1,60,1)*.5
  
  # Generate training data and compute DC
  train_data <- gendata_SPD_cate2(n, p, m, r, irep)
  x <- train_data$x
  Min <- train_data$Min
  dc_start <- Sys.time()
  dcov <- compute_dcov2(x, Min, "cate2", cholesky_error_C)
  dc_time <- difftime(Sys.time(), dc_start, units = "secs")
  
  # sorting dc time
  dc_sort_start <- Sys.time()
  dcov_ordered <- sort(dcov, decreasing = TRUE, index.return = TRUE)$ix
  dc_sort_time <- difftime(Sys.time(), dc_sort_start, units = "secs")

  # fitting time if we select 6 variables 
  ind <- dcov_ordered[1:6]
  ind <- ind[order(ind)]
  fit_start <- Sys.time()
  frechet_fit <- GloCovReg(x=as.matrix(x[, ind]), M = Min, xout = as.matrix(x[, ind]),
                           optns = list(metric="cholesky"))
  fit_time1 <- difftime(Sys.time(), fit_start, units = "secs")
  
  # fitting time if we select 8 variables 
  ind <- dcov_ordered[1:8]
  ind <- ind[order(ind)]
  fit_start <- Sys.time()
  frechet_fit <- GloCovReg(x=as.matrix(x[, ind]), M = Min, xout = as.matrix(x[, ind]),
                           optns = list(metric="cholesky"))
  fit_time2 <- difftime(Sys.time(), fit_start, units = "secs")
 
  
  friso_start <- Sys.time()
  # Find optimal lambda given each tau
  hopt=lambdaGlobalcoordesc(x=x,Min=Min, tau=tauall[1],lambdainit=NULL)
  hopt=cbind(hopt, lambdaGlobalcoordesc(x,Min, tauall[2], lambdainit=hopt))
  for(k in 3:length(tauall)){
    hopt=cbind(hopt,lambdaGlobalcoordesc(x,Min, tauall[k], lambdainit=hopt[,k-1]) )
    #print(hopt)
  }
  
  # Calculate error on training data 
  finalobj=NULL
  for (k in 1:length(tauall)) {
    Mres=IndivRidgeGloCovReg(x=x, M=Min, xout=x, lambda=hopt[,k], optns=list(metric="cholesky"))$Mout
    error = cholesky_error_C(Mres[[1]],Min[,,1])^2
    for(i in 2:nrow(x)){
      error = rbind(error,cholesky_error_C(Mres[[i]],Min[,,i])^2)
    }
    finalobj[k]=sum(error)
  }
  
  # Refit model on training data with only good predictors for each tau
  refit.finalobj=NULL
  for (k in 1:length(tauall)) {
    ind=which(hopt[,k]>1e-2)
    Mres=GloCovReg(x=as.matrix(x[,ind]), M = Min, xout=as.matrix(x[,ind]), optns=list(metric="cholesky"))$Mout
    error = cholesky_error_C(Mres[[1]],Min[,,1])^2
    for(i in 2:nrow(x)){
      error = rbind(error, cholesky_error_C(Mres[[i]],Min[,,i])^2)
    }
    refit.finalobj[k]=sum(error)
  }
  
  # Calculate error on independent validation data
  indep.finalobj=NULL
  for (k in 1:length(tauall)) {
    Mres=IndivRidgeGloCovReg(x=x, M=Min, xout=xtest, lambda=hopt[,k], optns=list(metric="cholesky"))$Mout
    error = cholesky_error_C(Mres[[1]],Mintest[,,1])^2
    for(i in 2:nrow(xtest)){
      error = rbind(error, cholesky_error_C(Mres[[i]],Mintest[,,i])^2)
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
    error = cholesky_error_C(Mres[[1]],Mintest[,,1])^2
    for(i in 2:nrow(xtest)){
      error = rbind(error, cholesky_error_C(Mres[[i]],Mintest[,,i])^2)
    }
    indep.refit.finalobj[k]=sum(error)
  }
  
  # Grab refitted model with best error in independent validation error
  indep.refit.ind=which.min(indep.refit.finalobj)
  indep.refit.finallambda=hopt[,indep.refit.ind]
  
  # record time of friso
  friso_time <- difftime(Sys.time(), friso_start, units = "secs")
  
  # Grab results
  result=list(hopt=hopt,                                      # all lambdas given tau
              finalobj=finalobj,                              # error on training data
              indep.ind=indep.ind,                            # best tau given error on validation data
              indep.finallambda=indep.finallambda,            # lambda given best tau from validation data
              refit.finalobj=refit.finalobj,                  # refitted error on training data
              indep.refit.finalobj=indep.refit.finalobj,      # refitted error on validation data
              indep.refit.ind=indep.refit.ind,                # refitted best tau on validation data
              indep.refit.finallambda=indep.refit.finallambda, # refitted lambda given best tau from validation data
              dcov = dcov, dc_time=dc_time, friso_time=friso_time,
              dc_sort_time = dc_sort_time, fit_time1 = fit_time1,
              fit_time2 = fit_time2) 
  
  return(result)
  
} 



# Prepare parallel processing
library(foreach)
library(quadprog)
library(iterators)
library(parallel)
library(doParallel)
registerDoParallel(25)
RNGkind("L'Ecuyer-CMRG")
set.seed(2022)
# repeat main function 100 times
friso_cate2=foreach(irep=1:100) %dopar% foo(irep)

save(file="friso_cate2.Rdata", friso_cate2)

stopImplicitCluster()
