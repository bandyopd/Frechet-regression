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

library(foreach)
library(quadprog)
library(iterators)
library(parallel)
library(doParallel)


registerDoParallel(25)
RNGkind("L'Ecuyer-CMRG")
set.seed(2022)
# repeat main function 500 times
dcov_spherical_p1000=foreach(irep=1:500) %dopar% compute_dcov(200, 1000, "spherical", irep, d_geo, r=0.5)
dcov_spherical_p1000=lapply(dcov_spherical_p1000, "[[", 1)
save(file="dcov_spherical_p1000.Rdata", dcov_spherical_p1000)
stopImplicitCluster()

