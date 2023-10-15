#include <RcppEigen.h>
#include <Rcpp.h>
// Eigen does not support implicit type casting between matrices using different types
// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
using Rcpp::List;

// [[Rcpp::export]]
MatrixXd chol_C(Map<MatrixXd> A) {
  //Eigen::MatrixXd P(3,3);
  //P << 6, 0, 0, 0, 4, 0, 0, 0, 7;
  MatrixXd L( A.llt().matrixL() );
  //std::cout << L.col(0) << std::endl;
  return(L.transpose());
}

// [[Rcpp::export]]
double cholesky_error_C(Map<MatrixXd> mat1, Map<MatrixXd> mat2){
  MatrixXd chol1 = chol_C(mat1);
  MatrixXd chol2 = chol_C(mat2);
  double error = ((chol1 - chol2) * (chol1 - chol2).transpose()).trace();
  error = sqrt(error);
  return(error);
}


// [[Rcpp::export]]
double vec_distance(Map<MatrixXd> X) {
  return((X.row(0)-X.row(1)).norm());
}

// [[Rcpp::export]]
double dcov_DTI_cpp(Map<MatrixXd> X, List Y) {
  // X is an n*p dimensional vector
  // Y is a list with length n, each element is a 3*3 matrix
  int n = X.rows();
  MatrixXd A(n, n);
  MatrixXd B(n, n);
  
  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      A(j, k) = (X.row(j)-X.row(k)).norm();
      B(j, k) = cholesky_error_C(Y[j], Y[k]);
    }
  }
  
  // row and col means
  VectorXd A_row_mean = A.rowwise().mean();
  VectorXd A_col_mean = A.colwise().mean();
  VectorXd B_row_mean = B.rowwise().mean();
  VectorXd B_col_mean = B.colwise().mean();
  double A_mean = A.sum() / (n * n);
  double B_mean = B.sum() / (n * n);
  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      A(j, k) = A(j, k) - A_row_mean(j) - A_col_mean(k) + A_mean;
      B(j, k) = B(j, k) - B_row_mean(j) - B_col_mean(k) + B_mean;
    }
  }
  
  return((A.array() * B.array()).sum() / (n * n));
}

