file_folder1 <- "/gpfs/home/ly21a/projects/frechet_regression/shared_code_SPD_Matrices/"
source(paste0(file_folder1, "IndivRidgeGloCovReg.R"))
source(paste0(file_folder1, "IndivRidgeGFRCovCholesky.R"))
source(paste0(file_folder1, "lambdaGlobalcoordesc.R"))
source(paste0(file_folder1, "GFRCovCholesky.R"))
source(paste0(file_folder1, "GloCovReg.R"))
source(paste0(file_folder1, "cholesky_error.R"))
source(paste0(file_folder1, "IndivRidgeGFRCovCholesky_v2.R"))
source(paste0(file_folder1, "ultility.R"))
source(paste0(file_folder1, "ultility_smooth.R"))
file_folder2 <- "/gpfs/home/ly21a/projects/frechet_regression/ADNI/"
# Rcpp::sourceCpp("IndivRidgeGFRCovCholesky_vC.cpp")
Rcpp::sourceCpp(paste0(file_folder2, "dcov_rcpp.cpp"))
Rcpp::sourceCpp(paste0(file_folder2, "IndivRidgeGFRCovCholesky_vC.cpp"))
library(gglasso)
library(mgcv)
library(geoR)


gen_SWP <- function(x, correlation, grid_size, m, seed_) {
  # generate spatially correlated SPD matrices using spatial Wishart process
  # described in the paper "Geostatistical Modeling of Positive Definite 
  # Matrices and Its Applications to Diffusion Tensor Imaging".
  # x: a n*p design matrix, n is the number of subject, p is the number of variables
  # correlation: correlation parameter that control the strength of spatial
  #              correlation
  # grid_size: a grid_size*grid_size grid is generated 
  # m: degree of freedom of the Wishart distribution of U
  # seed_: seed for reproducing results
  set.seed(seed_)
  seqx <- seq_len(grid_size)
  seqy <- seq_len(grid_size)
  s <- expand.grid(x = seqx, y = seqy) # expand.grid in package 'base'
  N <- grid_size^2 # number of SPD matrices on the grid, number of voxels
  cov <- varcov.spatial(
    coords = s,
    cov.model = "exponential",
    cov.pars = c(1, correlation)
  )$`varcov`
  
  beta11 <- t(cbind(rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01)))
  
  beta22 <- t(cbind(rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01)))
  
  beta33 <- t(cbind(rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01)))
  
  beta21 <- t(cbind(rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01)))
  
  beta31 <- t(cbind(rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01)))
  
  beta32 <- t(cbind(rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01),
                    rmvn(n = 1, mu=rep(0, N),  cov*0.01)))
  
  n <- nrow(x) # number of row
  Y <- array(dim = c(3, 3, n, N)) # 3*3*n sample*n voxel
  for (sub in seq_len(n)) {
    # for each subject
    D1=rmvn(n = m, mu=rep(0, N), cov)
    D2=rmvn(n = m, mu=rep(0, N), cov)
    D3=rmvn(n = m, mu=rep(0, N), cov)
    
    for (v in seq_len(N)) {
      # for each voxel
      L <- matrix(0, 3, 3)
      L[1, 1] <- exp(x[sub, c(1, 5, 10)]%*%beta11[, v])
      L[2, 2] <- exp(x[sub, c(1, 5, 10)]%*%beta22[, v])
      L[3, 3] <- exp(x[sub, c(1, 5, 10)]%*%beta33[, v])
      L[2, 1] <- x[sub, c(15, 20, 25)]%*%beta33[, v]
      L[3, 1] <- x[sub, c(15, 20, 25)]%*%beta33[, v]
      L[3, 2] <- x[sub, c(15, 20, 25)]%*%beta33[, v]
      
      mm <- cbind(D1[, v], D2[, v], D3[, v])
      Y[, , sub, v] <- L%*%t(mm)%*%mm%*%t(L)/m
    }
  }
  Y
}

x <- gendata_SPD(200, 1000, 3, 0.5, 2022)$x
Y <- gen_SWP(x, 15, 20, 3, 2022)

foo <- function(i) {
  n <- 200
  p <- 1000
  type_data <- "SPD"
  var_group <- seq_len(p)
  # generate data x, y
  X_i <- x
  Y_i <- Y[,,, i]
  # ----- DC + Chol ------------
  d <- 40 # keep the first d variables
  pred_DC_chol <- matrix(0, n, 6)
  dcov <- compute_dcov2(X_i, Y_i, type_data, cholesky_error_C)
  y_chol_mat <- convert_3D2mat(Y_i) # n*6 matrix
  dc_selected <- sort(dcov, decreasing = TRUE, index.return = TRUE)$ix[1:d]
  dc_selected <- dc_selected[order(dc_selected)]
  for (i in seq_len(6)) {
    pred_DC_chol[, i] <- fitted(lm(y_chol_mat[, i] ~ X_i[, dc_selected]))
  }
  pred_DC_chol_mat <- convert_mat23D(pred_DC_chol) # 3*3*n array
  
  # ----- DC + Frechet ----------
  pred_DC_fr <- GloCovReg(as.matrix(X_i[, dc_selected]), M = Y_i,
                          xout = as.matrix(X_i[, dc_selected]),
                          optns = list(metric="cholesky"))$Mout
  pred_DC_fr <- convert_list23D(pred_DC_fr)
  # ------ LASSO + Chol ----------
  coef_mat_chol <- matrix(0, nrow = p+1, ncol = 6) # each column 
  for (i in seq_len(6)) {
    if (sum(y_chol_mat[, i] != 0) == 0) next
    chol_fit <- cv.gglasso(X_i, y_chol_mat[, i], group = var_group, 
                           loss = "ls", pred.loss = "L2",
                           nfolds = 5, eps = 1e-4, maxit = 1e+6)
    #print(length(coef(chol_fit, s = "lambda.min")))
    coef_mat_chol[, i] <- coef(chol_fit, s = "lambda.min")
  }
  pred_lasso <- cbind(rep(1, n), X_i) %*% coef_mat_chol
  pred_lasso_mat <- convert_mat23D(pred_lasso) # 3*3*n array
  chol_selected <- select_var(coef_mat_chol, var_group)
  
  # ------ compute prediction errors --------
  pred_error_lasso <- compute_diff(pred_lasso_mat, Y_i)
  pred_error_DC_chol <- compute_diff(pred_DC_chol_mat, Y_i)
  pred_error_DC_fr <- compute_diff(pred_DC_fr, Y_i)
 
  
  MSE_lasso <- mean(pred_error_lasso)
  MSE_DC_chol <- mean(pred_error_DC_chol)
  MSE_DC_fr <- mean(pred_error_DC_fr)
  
  list(MSE_lasso = MSE_lasso, MSE_DC_chol = MSE_DC_chol, 
       MSE_DC_fr = MSE_DC_fr, pred_DC_fr = pred_DC_fr, 
       dcov = dcov, chol_selected = chol_selected)
}




library(foreach)
library(quadprog)
library(iterators)
library(parallel)
library(doParallel)
registerDoParallel(25)

RNGkind("L'Ecuyer-CMRG")
set.seed(2022)
# repeat main function 500 times
n_vox <- dim(Y)[4]
pe_SPD_p1000 = foreach(irep=1:n_vox) %dopar% foo(irep)

save(pe_SPD_p1000, Y, file="pe_SPD_p1000.Rdata")

stopImplicitCluster()














