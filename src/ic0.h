#define ARMA_WARN_LEVEL 0

#include <iostream>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "dist.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;



double dot_prod(int l1, int u1, int l2, int u2, NumericVector row_inds, NumericVector cells);
  
NumericVector ic0(NumericVector ptrs, NumericVector inds, NumericVector vals);

NumericVector createUcppM(NumericVector ptrs, NumericVector inds, NumericVector cov_vals);

NumericVector createUcpp(NumericVector ptrs, NumericVector inds, mat locsord, vec covparams);
