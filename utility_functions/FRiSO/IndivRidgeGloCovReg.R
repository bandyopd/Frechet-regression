#'@title Global Fréchet regression of covariance matrices
#'@description Global Fréchet regression of covariance matrices with Euclidean predictors.
#'@param x An n by p matrix of predictors.
#'@param y An n by l matrix, each row corresponds to an observation, l is the length of time points where the responses are observed. See 'metric' option in 'Details' for more details.
#'@param M A q by q by n array (resp. a list of q by q matrices) where \code{M[,,i]} (resp. \code{M[[i]]}) contains the i-th covariance matrix of dimension q by q.  See 'metric' option in 'Details' for more details.
#'@param xout An m by p matrix of output predictor levels.
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{corrOut}{Boolean indicating if output is shown as correlation or covariance matrix. Default is \code{FALSE} and corresponds to a covariance matrix.}
#' \item{metric}{Metric type choice, \code{"frobenius"}, \code{"power"}, \code{"log_cholesky"}, \code{"cholesky"} - default: \code{"frobenius"} which corresponds to the power metric with \code{alpha} equal to 1.
#' For power (and Frobenius) metrics, either \code{y} or \code{M} must be input; \code{y} would override \code{M}. For Cholesky and log-Cholesky metrics, \code{M} must be input and \code{y} does not apply.}
#' \item{alpha}{The power parameter for the power metric. Default is 1 which corresponds to Frobenius metric.}
#' }
#' @return A \code{covReg} object --- a list containing the following fields:
#' \item{xout}{An m by p matrix of output predictor levels.}
#' \item{Mout}{A list of estimated conditional covariance or correlation matrices at \code{xout}.}
#' \item{optns}{A list containing the \code{optns} parameters utilized.}
#' @examples

IndivRidgeGloCovReg= function(x,y=NULL,M=NULL,xout,lambda,optns = list()){
  if (is.null(optns$metric)){
    metric="frobenius"
  } else {
    metric=optns$metric
  }
  if(!metric%in%c("frobenius","power","cholesky","log_cholesky")){
    stop("metric choice not supported.")
  }
  if(metric=="frobenius"){
    res <- IndivRidgeGFRCov(x=x, y=NULL,M=M,xout=xout,lambda=lambda,optns = optns)
  } else if(metric=="power"){
    res <- IndivRidgeGFRCovPower(x=x, y=y,M=M,xout=xout,lambda=lambda,optns = optns)
  } else {
    if (is.null(M))
      stop("M must be input for Cholesky and log-Cholesky metrics; y does not apply.")
    res <- IndivRidgeGFRCovCholesky(x=x, M=M, xout=xout, lambda,optns = optns)
  }
  class(res) <- "covReg"
  return(res)
}
