#' GPvecchia: fast, scalable Gaussian process approximations
#'
#' The package can be used for parameter inference and prediction for Gaussian and non-Gaussian spatial data using many popular GP approximation methods.
#'
#'
#' @docType package
#' @name GPvecchia
#' @useDynLib GPvecchia
#' @importFrom Rcpp sourceCpp
#'
#' @example
#' z=rnorm(10); locs=matrix(1:10,ncol=1); nuggets = rep(.1, 10)
#' vecchia.approx=vecchia_specify(locs,m=5)
#' preds=vecchia_prediction(z,vecchia.approx,covparms=c(1,2,.5), nuggets)
NULL



