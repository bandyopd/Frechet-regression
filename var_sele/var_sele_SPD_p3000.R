file_folder1 <- "/gpfs/home/ly21a/projects/frechet_regression/shared_code_SPD_Matrices/"
source(paste0(file_folder1, "IndivRidgeGloCovReg.R"))
source(paste0(file_folder1, "IndivRidgeGFRCovCholesky.R"))
source(paste0(file_folder1, "lambdaGlobalcoordesc.R"))
source(paste0(file_folder1, "GFRCovCholesky.R"))
source(paste0(file_folder1, "GloCovReg.R"))
source(paste0(file_folder1, "cholesky_error.R"))
source(paste0(file_folder1, "IndivRidgeGFRCovCholesky_v2.R"))
source(paste0(file_folder1, "ultility.R"))
file_folder2 <- "/gpfs/home/ly21a/projects/frechet_regression/ADNI/"
# Rcpp::sourceCpp("IndivRidgeGFRCovCholesky_vC.cpp")
Rcpp::sourceCpp(paste0(file_folder2, "dcov_rcpp.cpp"))
library(gglasso)


foo <- function(i) {
  type_data <- "SPD"
  n <- 200
  p <- 3000
  d_y <- cholesky_error_C
  
  # first dcov, and obtain x, y
  dcov_result <- compute_dcov(n, p, type_data, i, d_y, r=0.5)
  x <- dcov_result$X
  y <- dcov_result$Y # y is a list
  var_group <- seq_len(p)
  # then FA + LASSO
  y_FA <- convert_3D2FA(convert_list23D(y)) # FA vector
  FA_fit <- cv.gglasso(x, y_FA, group = var_group, loss = "ls", pred.loss = "L2",
                       nfolds = 5, eps = 1e-4, maxit = 1e+6)
  coef_FA_fit <- coef(FA_fit, s = "lambda.min")
  FA_selected <- select_var_vec(coef_FA_fit, var_group)
  # coef includes the intercept, first element is intercept
  
  # last Cholesky + LASSO
  y_chol_mat <- convert_3D2mat(convert_list23D(y)) # n*6 matrix
  coef_mat_chol <- matrix(0, nrow = p+1, ncol = 6) # each column 
  for (i in seq_len(6)) {
    if (sum(y_chol_mat[, i] != 0) == 0) next
    chol_fit <- cv.gglasso(x, y_chol_mat[, i], group = var_group, 
                           loss = "ls", pred.loss = "L2",
                           nfolds = 5, eps = 1e-4, maxit = 1e+6)
    coef_mat_chol[, i] <- coef(chol_fit, s = "lambda.min")
  }
  chol_selected <- select_var(coef_mat_chol, var_group)
  list(dcov = dcov_result$dcov_SPD, FA_selected = FA_selected,
       chol_selected = chol_selected)
}


library(foreach)
library(quadprog)
library(iterators)
library(parallel)
library(doParallel)
registerDoParallel(25)

RNGkind("L'Ecuyer-CMRG")
set.seed(2022)
# repeat main function 100 times
var_sele_SPD_p3000 = foreach(irep=1:500) %dopar% foo(irep)

save(file="var_sele_SPD_p3000.Rdata", var_sele_SPD_p3000)

stopImplicitCluster()



#p <- 400
#n <- 100
#nrep <- n
#idx <- matrix(0, nrep, p)
#for (i in 1:nrep) {
#  idx[i,] <- 
#    sort(allresult_dcov[((i-1)*p+1):(i*p)], decreasing = TRUE, index.return = TRUE)$ix
#}
#selected_var <- idx[, 1:3]
#var_count <- rep(0, p)
#for (i in 1:nrep) {
#  var_count[selected_var[i, ]] = 1 + var_count[selected_var[i, ]]
#}
#var_count










