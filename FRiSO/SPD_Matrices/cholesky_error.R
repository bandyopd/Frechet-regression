#' @title Cholesky distance function 
#' @description This is the function that computes Cholesky Errors (distance between two matrices)
#' @details 
#  For SPD matrices, P1 and P2, it decomposes P1=L1%*%t(L1), P2= L2%*%t(L2)
#' Where L1 and L2 are lower triangle matrices with positive diagonal components
#' It then takes the Frobenius norm:
#' sqrt(tr((L1-L2)%*%t(L1-L2)))
#' @param mat1: M by M symmetric positive definite matrix
#' @param mat2: M by M symmetric positive definite matrix

cholesky_error = function(mat1, mat2){
  chol1 = chol(mat1)
  chol2 = chol(mat2)
  error = sqrt(tr((chol1-chol2)%*%t(chol1-chol2)))
  return(error=error)
}