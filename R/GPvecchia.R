#' GPvecchia: fast, scalable Gaussian process approximations
#'
#' The foo package provides three categories of important functions:
#' foo, bar and baz.
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



