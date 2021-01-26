#' @title Numeric Delta Method approximation for the variance-covariance matrix
#'
#' @description
#'  Computes delta method variance-covariance matrix of results of any generic 
#'  function fct that computes a vector of estimates as a function of a set of 
#'  estimated parameters par.
#'  
#' @details 
#' This function is copied vebatim from package \link{mrds}.
#' 
#' The delta method (aka propagation of errors is based on Taylor series 
#' approximation). It uses 
#' the first derivative of fct with respect to par which is computed in this 
#' function numerically using the central-difference formula. It also uses the 
#' variance-covariance matrix of the estimated parameters which is derived in 
#' estimating the parameters and is an input argument.
#' 
#' The first argument of fct should be par which is a vector of parameter 
#' estimates. It should return a single value (or vector) of estimate(s). The 
#' remaining arguments of fct if any can be passed to fct by including them at 
#' the end of the call to DeltaMethod as name=value pairs.
#' 
#' @param par vector of parameter values at which estimates should be constructed.
#' @param fct function that constructs estimates from parameters \code{par}.
#' @param vcov variance-covariance matrix of the parameters.
#' @param  delta proportional change in parameters used to numerically estimate 
#' first derivative with central-difference formula.
#' @param ... any additional arguments needed by \code{fct}.
#' 
#' @return Returns a list with values
#' \itemize{
#'  \item{variance}{ Estimated variance-covariance matrix of estimates derived
#'  by \code{fct}.}
#'  \item{partial}{ Matrix (or vector) of partial derivatives of\code{fct} with
#'  respect to the parameters  \code{par}.}
#' }
#' @examples 
#' 
#' par = c(-4.77, 2.64)
#' vcv = matrix(c(0.067,-0.054, -0.054, 0.068),nrow=2,byrow=TRUE)
#' fn = function(par) invlogit(sum(par))
#' DeltaMethod(par=beta.g0.d, fct=invlogit, vcov=vcv)
#' 
#' @export 
#' 
DeltaMethod = function (par, fct, vcov, delta=0.01, ...) {
  theta <- function(par) fct(par, ...)
  savepar <- par
  value1 <- theta(par)
  partial <- matrix(0, nrow = length(par), ncol = length(value1))
  for (i in seq_along(par)) {
    if (savepar[i] != 0) {
      deltap <- delta * savepar[i]
    }
    else {
      deltap <- delta
    }
    par <- savepar
    par[i] <- savepar[i] + deltap
    value1 <- theta(par)
    par <- savepar
    par[i] <- savepar[i] - deltap
    value2 <- theta(par)
    partial[i, ] <- (value1 - value2)/(2 * deltap)
  }
  variance <- t(partial) %*% vcov %*% partial
  return(list(variance = variance, partial = partial))
}
