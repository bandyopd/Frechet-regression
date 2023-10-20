#include <RcppEigen.h>
#include <Rcpp.h>
// Eigen does not support implicit type casting between matrices using different types
// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
using Eigen::Map;
using Rcpp::List;


// [[Rcpp::export]]
MatrixXd mat_minus_vector(Map<MatrixXd> A, Map<VectorXd> v) {
  MatrixXd newmat = (A.transpose().colwise() - v).transpose();
  return(newmat);
}
// [[Rcpp::export]]
MatrixXd mat_minus_scala(Map<MatrixXd> A, double v) {
  MatrixXd newmat = A.array() - v;
  return(newmat);
}

// [[Rcpp::export]]
VectorXd mat_mul_vec(Map<MatrixXd> A, Map<MatrixXd> B, VectorXd v) {
  VectorXd newvec = A * (B.row(0) - v);
  return(newvec);
}

// [[Rcpp::export]]
MatrixXd Kron(MatrixXd A, int m) {
  // kronecker product between A and identity matrix with dim m
  // A is a column vector
  MatrixXd Im = MatrixXd::Identity(m, m);
  int n = A.rows();
  MatrixXd product(m * n, m * 1);
  product.setZero();
  for (int i = 0; i < n; i++) {
    product.block(i*m, 0, m, m) = Im.array() * A(i, 0);
    //std::cout << Im.array() * A(i, 0) << std::endl;
  }
  //MatrixXd temp = Im.array() * A(1, 1);
  return(product);
}


// [[Rcpp::export]]
List IRGFRCCCpp(Map<MatrixXd> chol_Y, Map<MatrixXd> x,
                                 Map<MatrixXd> xout, Map<MatrixXd> NewInvSigma,
                                 Map<VectorXd> mx, int m, int nout) {
  List Mout(nout);
  for(int j = 0; j < nout; j++) {
    VectorXd weight_Y_minus1 = (x.transpose().colwise() - mx).transpose() * NewInvSigma * (xout.row(j).transpose() - mx);
    //VectorXd weight_Y_minus1 = (x.transpose().colwise() - mx).transpose() * NewInvSigma * (xout.row(j-1) - mx);
    VectorXd weight_Y = weight_Y_minus1.array() + 1;
    MatrixXd weighted_mat = chol_Y * Kron(weight_Y, m);
    MatrixXd sqrt_mat = weighted_mat.array() / weight_Y.sum();
    Mout[j] = sqrt_mat.transpose() * sqrt_mat;
    //Mout[j] = weight_Y;
  }
  
  return(Mout);
}


// [[Rcpp::export]]
MatrixXd WeightedFrechetMeanCpp(Map<MatrixXd> chol_Y, VectorXd weight_Y, 
                            int m) {

  MatrixXd weighted_mat = chol_Y * Kron(weight_Y, m);
  MatrixXd sqrt_mat = weighted_mat.array() / weight_Y.sum();
  MatrixXd Mout = sqrt_mat.transpose() * sqrt_mat;

  
  return(Mout);
}

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
