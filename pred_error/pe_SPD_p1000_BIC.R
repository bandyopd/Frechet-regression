file_folder1 <- ""
source(paste0(file_folder1, "IndivRidgeGloCovReg.R"))
source(paste0(file_folder1, "IndivRidgeGFRCovCholesky.R"))
source(paste0(file_folder1, "lambdaGlobalcoordesc.R"))
source(paste0(file_folder1, "GFRCovCholesky.R"))
source(paste0(file_folder1, "GloCovReg.R"))
source(paste0(file_folder1, "cholesky_error.R"))
source(paste0(file_folder1, "IndivRidgeGFRCovCholesky_v2.R"))
source(paste0(file_folder1, "ultility.R"))
source(paste0(file_folder1, "ultility_smooth.R"))
file_folder2 <- ""
# Rcpp::sourceCpp("IndivRidgeGFRCovCholesky_vC.cpp")
Rcpp::sourceCpp(paste0(file_folder2, "dcov_rcpp.cpp"))
Rcpp::sourceCpp(paste0(file_folder2, "IndivRidgeGFRCovCholesky_vC.cpp"))
library(gglasso)
library(mgcv)
library(geoR)

# get active predictor, the indices of zero rows of estimated beta
get_nonzero_row <- function(beta_hat) {
  # beta_hat: a p-by-r coefficient matrix
  which(apply(beta_hat, 1, function(x) any(x != 0)))
}


compute_BIC <- function(diff_square, d) {
  # Input 
  #   diff_square: the squared distance between Y and Y_hat
  #             d: model size, the number of variables used in Frechet regression
  m <- length(diff_square) # sample size
  log(mean(diff_square)) + log(m) * d / m
  #2*log(mean(diff_square))
  #log(m) * d / m
}


compute_diff2 <- function(Y, Y_hat) {
  # Y: a 3*3*m arrays
  # Y_hat: a list with length m
  m <- length(Y_hat)
  difference <- vector("double", m)
  for (i in seq_len(m)) {
    difference[i] <- cholesky_error_C(Y[, , i], Y_hat[[i]])^2
  }
  difference
}

# select model size d according to BIC
select_d_BIC <- function(X, Y, ranked_dcov) {
  # Input:
  # X: m*p matrix, each row is a sample
  # Y: 3*3*m array, where m is the number of samples
  # ranked_dcov: a p*1 vector
  # Output: 
  # d_BIC: the selected model size d based on BIC
  
  d_candidate <- seq(6, 37, 5) # all model size candidate
  BIC_val <- rep(0, length(d_candidate)) # initialize with 0
  for (i in seq_along(d_candidate)) {
    d <- d_candidate[i]
    selected_var <- ranked_dcov[1:d] # select d variables
    ind <- selected_var[order(selected_var)]
    Y_hat <- GloCovReg(as.matrix(X[, ind]), M = Y,
                       xout = as.matrix(X[, ind]), 
                       optns = list(metric="cholesky"))$Mout
    diff_square <- compute_diff2(Y, Y_hat)
    BIC_val[i] <- compute_BIC(diff_square, d)
  }
  list(d_selected = d_candidate[which.min(BIC_val)], # model size with the smallest BIC is selected
       BIC_val = BIC_val) 
}


#source("gen_SWP_modified.R")
#source("gen_SWP.R")
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
  
  mu_beta <- 0
  sigma_beta <- 0.01
  beta11 <- t(cbind(rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta),
                    rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta),
                    rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta)))
  
  beta22 <- t(cbind(rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta),
                    rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta),
                    rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta)))
  
  beta33 <- t(cbind(rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta),
                    rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta),
                    rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta)))
  
  beta21 <- t(cbind(rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta),
                    rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta),
                    rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta)))
  
  beta31 <- t(cbind(rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta),
                    rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta),
                    rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta)))
  
  beta32 <- t(cbind(rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta),
                    rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta),
                    rmvn(n = 1, mu=rep(mu_beta, N),  cov*sigma_beta)))
  
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
      L[2, 1] <- x[sub, c(15, 20, 25)]%*%beta21[, v]
      L[3, 1] <- x[sub, c(15, 20, 25)]%*%beta31[, v]
      L[3, 2] <- x[sub, c(15, 20, 25)]%*%beta32[, v]
      
      mm <- cbind(D1[, v], D2[, v], D3[, v])
      Y[, , sub, v] <- L%*%t(mm)%*%mm%*%t(L)/m
    }
  }
  Y
}


foo <- function(v) {
  type_data <- "SPD"
  var_group <- seq_len(p)
  # generate data x, y
  n_train <- nrow(x_train)
  n_test <- nrow(x_test)
  Y_v_train <- Y_train[,,, v]
  Y_v_test <- Y_test[,,, v]

  # ----- DC + Chol ------------
  pred_DC_chol <- matrix(0, n_test, 6)
  pred_DC_chol_train <- matrix(0, n_train, 6)
  dcov <- compute_dcov2(x_train, Y_v_train, type_data, cholesky_error_C)
  y_chol_mat <- convert_3D2mat(Y_v_train) # n*6 matrix

  # use d selected by BIC
  ranked_dc <- sort(dcov, decreasing = TRUE, index.return = TRUE)$ix
  d_BIC <- select_d_BIC(x_train, Y_v_train, ranked_dc)$d_selected
  dc_selected <- ranked_dc[1:d_BIC]
  dc_selected <- dc_selected[order(dc_selected)]
  #dc_selected <- c(1, 5, 10, 15, 20, 25)
  mean_abs_coef <- rep(0, 6) # record the mean of abs of coef
  for (i in seq_len(6)) {
    lm_fit <- lm(y_chol_mat[, i] ~ x_train[, dc_selected])
    coef_lm <- matrix(coef(lm_fit), ncol = 1)
    if (sum(is.na(coef_lm)) != 0) {
      coef_lm[is.na(coef_lm)] <- 0
    }
    pred_DC_chol[, i] <- cbind(rep(1, n_test), x_test[, dc_selected]) %*% coef_lm
    pred_DC_chol_train[, i] <- fitted(lm_fit)
    mean_abs_coef[i] <- mean(abs(coef_lm))
  }
  pred_DC_chol_mat <- convert_mat23D(pred_DC_chol) # 3*3*n array
  pred_DC_chol_train_mat <- convert_mat23D(pred_DC_chol_train)
  
  # ----- DC + Frechet ----------
  pred_DC_fr <- GloCovReg(as.matrix(x_train[, dc_selected]), M = Y_v_train,
                          xout = as.matrix(x_test[, dc_selected]),
                          optns = list(metric="cholesky"))$Mout
  pred_DC_fr <- convert_list23D(pred_DC_fr)

  pred_DC_fr_train <- GloCovReg(as.matrix(x_train[, dc_selected]), M = Y_v_train,
                                xout = as.matrix(x_train[, dc_selected]),
                                optns = list(metric="cholesky"))$Mout
  pred_DC_fr_train <- convert_list23D(pred_DC_fr_train)

  # ------ LASSO + Chol ----------
  coef_mat_chol <- matrix(0, nrow = p+1, ncol = 6) # each column 
  for (i in seq_len(6)) {
    if (sum(y_chol_mat[, i] != 0) == 0) next
    chol_fit <- cv.gglasso(x_train, y_chol_mat[, i], nlambda = 10, group = var_group, loss = "ls",
                           eps = 1e-3, maxit = 1e+3)
    #print(length(coef(chol_fit, s = "lambda.min")))
    coef_mat_chol[, i] <- coef(chol_fit, s = "lambda.min")
  }
  lasso_selected <- get_nonzero_row(coef_mat_chol[-1, ])

  # refit with selected variables
  n_element <- 6
  pred_lasso <- matrix(0, n_test, n_element)
  pred_lasso_train <- matrix(0, n_train, n_element)
  for (ii in seq_len(n_element)) {
    possible_error_lasso <- tryCatch({
      lm_fit <- lm(y_chol_mat[, ii] ~ x_train[, lasso_selected])
      coef_lm <- matrix(coef(lm_fit), ncol = 1)
      if (sum(is.na(coef_lm)) != 0) {
        coef_lm[is.na(coef_lm)] <- 0
      }
      pred_lasso[, ii] <- cbind(rep(1, n_test), x_test[, lasso_selected, drop = FALSE]) %*% coef_lm
      pred_lasso_train[, ii] <- fitted(lm_fit)
    }, error = function(e) e)
      
    if (inherits(possible_error_lasso, "error")) {
      pred_lasso[, ii] <- rep(mean(y_chol_mat[, ii]), n_test)
      pred_lasso_train[, ii] <- rep(mean(y_chol_mat[, ii]), n_train)
    } 
  }


  pred_lasso_mat <- convert_mat23D(pred_lasso) # 3*3*n array
  pred_lasso_mat_train <- convert_mat23D(pred_lasso_train) # 3*3*n array


  



  # ------ compute prediction errors --------
  pred_error_lasso <- compute_diff(pred_lasso_mat, Y_v_test)
  pred_error_DC_chol <- compute_diff(pred_DC_chol_mat, Y_v_test)
  pred_error_DC_fr <- compute_diff(pred_DC_fr, Y_v_test)

  # compute training error 
  pred_error_lasso_train <- compute_diff(pred_lasso_mat_train, Y_v_train)
  pred_error_DC_chol_train <- compute_diff(pred_DC_chol_train_mat, Y_v_train)
  pred_error_DC_fr_train <- compute_diff(pred_DC_fr_train, Y_v_train)


  MSE_lasso <- mean(pred_error_lasso)
  MSE_DC_chol <- mean(pred_error_DC_chol)
  MSE_DC_fr <- mean(pred_error_DC_fr)


  MSE_lasso_train <- mean(pred_error_lasso_train)
  MSE_DC_chol_train <- mean(pred_error_DC_chol_train)
  MSE_DC_fr_train <- mean(pred_error_DC_fr_train)

  # compute the abs of coef of lasso and dc
  mean_abs_coef_dc <- mean(mean_abs_coef)


  
  list(MSE_lasso = MSE_lasso, MSE_DC_chol = MSE_DC_chol, 
       MSE_DC_fr = MSE_DC_fr, pred_DC_fr = pred_DC_fr, 
       dc_selected = dc_selected, lasso_selected = lasso_selected,
       MSE_lasso_train = MSE_lasso_train,
       MSE_DC_chol_train = MSE_DC_chol_train, MSE_DC_fr_train = MSE_DC_fr_train,
       pred_DC_fr_train = pred_DC_fr_train,
       mean_abs_coef_dc = mean_abs_coef_dc, 
       pred_error_lasso = pred_error_lasso, pred_error_DC_chol = pred_error_DC_chol,
       pred_error_DC_fr = pred_error_DC_fr, 
       pred_error_lasso_train = pred_error_lasso_train, pred_error_DC_chol_train = pred_error_DC_chol_train,
       pred_error_DC_fr_train = pred_error_DC_fr_train, 
       d_BIC = d_BIC)
}




library(foreach)
library(quadprog)
library(iterators)
library(parallel)
library(doParallel)
library(glmnet)
registerDoParallel(25)

RNGkind("L'Ecuyer-CMRG")
set.seed(2022)

n <- 400 # sample size
p <- 1000   # number of predictor
r <- 0.5 # correlation between predictors
x <- 2*pnorm(matrix(rnorm(n*p),n,p)%*%chol(r^(as.matrix(dist(1:p)))))
#x <- matrix(rnorm(n*p),n,p)%*%chol(r^(as.matrix(dist(1:p))))
Y <- gen_SWP(x, 15, 20, 3, 2022)

# Create random indices
train_index <- sample(1:n, 200)
# Split the data
x_train <- x[train_index, ]
x_train <- scale(x_train)
x_test <- x[-train_index, ]
x_test <- scale(x_test)
Y_train <- Y[, , train_index, ]
Y_test <- Y[, , -train_index, ]

# repeat main function 500 times
n_vox <- dim(Y)[4]
pe_SPD_p1000 = foreach(irep=1:n_vox) %dopar% foo(irep)

save(pe_SPD_p1000, Y, train_index, file="pe_SPD_p1000_BIC.Rdata")

stopImplicitCluster()














