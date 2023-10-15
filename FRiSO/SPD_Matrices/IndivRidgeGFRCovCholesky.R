#' @title Global Fréchet regression of covariance matrices with Log-Cholesky and Cholesky metric
#' @noRd
#' @description Global Fréchet regression of covariance matrices with Euclidean predictors.
#' 
#' @param x an n by p matrix of predictors.
#' @param M an q by q by n array (resp. a list of q by q matrices) where \code{M[,,i]} (resp. \code{M[[i]]}) contains the i-th covariance matrix of dimension q by q.
#' @param xout an m by p matrix of output predictor levels.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are 
#' \describe{
#' \item{corrOut}{Boolean indicating if output is shown as correlation or covariance matrix. Default is \code{FALSE} and corresponds to a covariance matrix.}
#' \item{metric}{Metric type choice, "log_cholesky", "cholesky" - default: \code{log_cholesky} for log Cholesky metric}
#' }
#' 
#' @return A list containing the following fields:
#' \item{xout}{An m by p matrix of output predictor levels.}
#' \item{Mout}{A list of estimated conditional covariance matrices at \code{xout}.}
#' \item{optns}{A list containing the \code{optns} parameters utilized.}
#' 

IndivRidgeGFRCovCholesky <- function(x, M, xout, lambda, optns = list()){
  
  if(!is.matrix(x)&!is.vector(x)){
    stop('x must be a matrix or vector')
  }
  if(!is.matrix(xout)&!is.vector(xout)){
    stop('xout must be a matrix or vector')
  }
  if(is.vector(x)){x<- matrix(x,length(x)) }
  
  if(is.vector(xout)){xout<- matrix(xout,length(xout)) }
  
  if(ncol(x) != ncol(xout)){
    stop('x and xout must have same number of columns')
  }
  
  if(is.null(M)){
    stop("M must be provided")
  }
  if(class(M) == 'list'){
    M=array(as.numeric(unlist(M)), dim=c(dim(M[[1]])[1],dim(M[[1]])[1],length(M)))
  }else{
    if(!class(M)=="array"){
      stop('M must be an array or a list')
    }
  }
  
  if(nrow(x)!=dim(M)[3]){
    stop("the number of rows of x must be the same as the number of covariance matrices in M")
  }
  
  if(is.null(optns$corrOut)){
    corrOut = FALSE
  } else {
    corrOut = optns$corrOut
  }
  
  if(is.null(optns$metric)){
    metric = 'log_cholesky'
  } else {
    metric =  optns$metric
  }
  
  n = nrow(x)
  p = ncol(x)
  nout = nrow(xout)
  invVa = solve(var(x))
  mx = apply(x,2,mean)
  
  # Compute the covariance matrix weighted by lambda parameter
  diagsqrtlambda=diag(sqrt(lambda))
  NewInvSigma=diagsqrtlambda%*%solve(diagsqrtlambda%*%var(x)%*%diagsqrtlambda+diag(length(lambda)))%*%diagsqrtlambda
  
  # Convert matrix data into list type
  MM = list()
  if(class(M) == 'array'){
    for (i in 1:n) {
      MM[[i]] = M[,,i]
    }
  } else {MM = M}
  
  # Ensure symmetric matrix
  M = lapply(MM, function(X) (X+t(X))/2)
  m = nrow(M[[1]])
  Mout = list()
  if(metric == 'log_cholesky'){
    LL = lapply(M, chol)
    L = lapply(LL, function(X) X - diag(diag(X)))
    D = lapply(LL, function(X) diag(X))
    
    for (j in 1:nout){
      ss = 0
      U = 0
      E = 0
      s = array(0,n)
      for (i in 1:n) {
        s[i] = 1+(x[i,]-mx)%*%NewInvSigma%*%(xout[j,]-mx)
        ss = ss + s[i]
        U = U + s[i]*L[[i]]
        E = E + s[i]*log(D[[i]])
      }
      SS = U/ss + diag(exp(E/ss))
      Mout[[j]] = t(SS)%*%SS
    }
  } else{
    # Fit individually weighted ridge Frechet regression
    # nout x n matrix of numerator weights
    smat = t(sapply(1:nout, function(j){
      # n x 1 vector
      sj = as.vector(1 + t(t(x) - mx) %*% NewInvSigma %*% (xout[j,] - mx))
      sout=ramify::flatten(matrix(rep(1,n*m),nrow=n)*sj)
      return(as.vector(sout))
    }))
    
    # example, row 1 from 1:10 corresponds to s[1]*L[[1]] for the first xout observation
    # row 1 from 11:20 corresponds to s[2]*L[[2]] for first xout observation
    smat2 = t(sapply(1:nout, function(j){
      # n x 1 vector
      sj = 1 + t(t(x) - mx) %*% NewInvSigma %*% (xout[j,] - mx)
      return(as.vector(sj))
    }))
    
    # nout x 1 vector of denominator weights
    ss=rowSums(smat2)
    # create large matrix of all Choleskys
    L=do.call(rbind, lapply(M,chol))
    
    # For each xout, weights all n L's and then sums them up
    Mout = lapply(1:nout, function(j){
      temp =L*smat[j,]
      temp2 = split.data.frame(temp, cut(seq_len(n*m), n))
      numeratorj = Reduce("+", temp2)
      return((t(numeratorj)/ss[j]) %*% (numeratorj/ss[j]))
    })
  }
  if(corrOut){
    for(j in 1:nrow(xout)){
      D=diag(1/sqrt(diag(Mout[[j]])))
      Mout[[j]]=D%*%Mout[[j]]%*%D
      Mout[[j]]=as.matrix(Matrix::forceSymmetric(Mout[[j]]))
    }
  }
    
  out = list(xout=xout, Mout=Mout, optns=list(corrOut=corrOut,metric=metric))
  return(out)
}