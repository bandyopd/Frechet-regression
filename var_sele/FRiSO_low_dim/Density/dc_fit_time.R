# this script compute fitting time after we select variables using dc

library(psych)
source("IndivRidgeGloWassReg.R")
source("lambdaGlobalcoordesc.R")
source("GloWassReg.R")
file_folder1 <- ""
file_folder2 <- ""
source(paste0(file_folder1, "utility.R"))
Rcpp::sourceCpp(paste0(file_folder2, "dcov_rcpp.cpp"))

GlobalObj4h=function(x,Qin,lambda) {
  Mres=IndivRidgeGloWassReg(x, Qin, x, lambda)
  return(sum((Mres-Qin)^2))
}


foo <- function(irep){
  n=200
  p=20
  m=20
  r=0.5
  
  # Generate training data and compute DC
  train_data <- gendata_prob(n, p, m, r, irep)
  x <- train_data$x
  Qin <- train_data$Qin
  dc_start <- Sys.time()
  dcov <- compute_dcov2(x, Qin, "density", d_vec)
  dc_time <- difftime(Sys.time(), dc_start, units = "secs")
  
  # record sort time
  dc_sort_start <- Sys.time()
  dcov_ordered <- sort(dcov, decreasing = TRUE, index.return = TRUE)$ix
  dc_sort_time <- difftime(Sys.time(), dc_sort_start, units = "secs")
  
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
  
  list(dc_time = dc_time, dc_sort_time = dc_sort_time, fit_time1 = fit_time1,
       fit_time2 = fit_time2)
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
density_fit_time=foreach(irep=1:100) %dopar% foo(irep)

save(file="density_fit_time.Rdata", density_fit_time)

stopImplicitCluster()
