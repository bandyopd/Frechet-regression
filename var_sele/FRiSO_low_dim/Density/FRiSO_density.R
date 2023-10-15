#' @title Frechet Ridge Selection Operator for Global Frechet Regression with Density Output
#' 
#' @description This code implements the Frechet Ridge Selection Operator
#' on Density type data equipped with the Wasserstein distance 
#' 
#' @details The outline is as follows:
#' 
#' 1. Define the data generation function utilized for example 5.2.2.
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

#####################################################################
#rm(list=ls())
#####################################################################
library(psych)
source("IndivRidgeGloWassReg.R")
source("lambdaGlobalcoordesc.R")
source("GloWassReg.R")
file_folder1 <- ""
file_folder2 <- ""
source(paste0(file_folder1, "utility.R"))
Rcpp::sourceCpp(paste0(file_folder2, "dcov_rcpp.cpp"))
####################################################################


GlobalObj4h=function(x,Qin,lambda) {
  Mres=IndivRidgeGloWassReg(x, Qin, x, lambda)
  return(sum((Mres-Qin)^2))
}

foo <- function(irep){
  n=200
  p=20
  m=20
  r=0.5

  seed = 2022
  
  # Generate independent tuning data
  wutest=gendata_prob(n,p,m,r,seed)
  xtest=wutest$x
  Qintest=wutest$Qin
  
  # Grid for tau = c(1,2,...,30)
  tauall=seq(0.5,30,0.5)*2
  
  # Generate training data and compute DC
  train_data <- gendata_prob(n, p, m, r, irep)
  x <- train_data$x
  Qin <- train_data$Qin
  dc_start <- Sys.time()
  dcov <- compute_dcov2(x, Qin, "density", d_vec)
  dc_time <- difftime(Sys.time(), dc_start, units = "secs")
  
  dcov_ordered <- sort(dcov, decreasing = TRUE, index.return = TRUE)$ix
  # fitting time if we select 6 variables 
  ind <- dcov_ordered[1:6]
  ind <- ind[order(ind)]
  fit_start <- Sys.time()
  frechet_fit <- GloWassReg(x[, ind], Qin, x[, ind])
  fit_time1 <- difftime(Sys.time(), fit_start, units = "secs")
  
  # fitting time if we select 12 variables 
  ind <- dcov_ordered[1:12]
  ind <- ind[order(ind)]
  fit_start <- Sys.time()
  frechet_fit <- GloWassReg(x[, ind], Qin, x[, ind])
  fit_time2 <- difftime(Sys.time(), fit_start, units = "secs")
  
  friso_start <- Sys.time()
  # Find optimal lambda given each tau
  hopt=lambdaGlobalcoordesc(x,Qin, tauall[1])
  hopt=cbind(hopt, lambdaGlobalcoordesc(x,Qin, tauall[2], lambdainit=hopt))
  for (k in 3:length(tauall)) {
    hopt=cbind(hopt, lambdaGlobalcoordesc(x,Qin, tauall[k], lambdainit=hopt[,k-1]))
    #print(hopt)}
  }
  # Calculate error on training data 
  finalobj=NULL
  for (k in 1:length(tauall)) {
    Mres=IndivRidgeGloWassReg(x, Qin, x, hopt[,k])
    finalobj[k]=base::sum((Qin-Mres)^2)
  }
  
  # Refit model on training data with only good predictors for each tau
  refit.finalobj=NULL
  for (k in 1:length(tauall)) {
    ind=which(hopt[,k]>1e-2)
    Mres=GloWassReg(x[,ind], Qin, x[,ind])
    refit.finalobj[k]=base::sum((Qin-Mres)^2)
  }
  
  # Calculate error on independent validation data
  indep.finalobj=NULL
  for (k in 1:length(tauall)) {
    Mres=IndivRidgeGloWassReg(x, Qin, xtest, hopt[,k])
    indep.finalobj[k]=base::sum((Qintest-Mres)^2)
  }
  
  # Grab best model based on error on validation data
  indep.ind=which.min(indep.finalobj)
  indep.finallambda=hopt[,indep.ind]
  

  # Refit on model on independent validation data with only good predictors for each tau
  indep.refit.finalobj=NULL
  for (k in 1:length(tauall)) {
     ind=which(hopt[,k]>1e-2)
     Mres=GloWassReg(x[,ind], Qin, xtest[,ind])
     indep.refit.finalobj[k]=base::sum((Qintest-Mres)^2)
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
              indep.refit.finallambda=indep.refit.finallambda,# refitted lambda given best tau from validation data
              dcov = dcov, dc_time = dc_time, friso_time = friso_time,
              fit_time1 = fit_time1, fit_time2 = fit_time2) 
  
  return(result)
} 


# prepare parallel processing
library(foreach)
library(quadprog)
library(iterators)
library(parallel)
library(doParallel)
registerDoParallel(25)
RNGkind("L'Ecuyer-CMRG")
set.seed(2022)
# repeat main function 100 times
FRiSO_density=foreach(irep=1:100) %dopar% foo(irep)

save(file="FRiSO_density.Rdata", FRiSO_density)

stopImplicitCluster()

