#define ARMA_DONT_PRINT_ERRORS

#include <iostream>
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;
using namespace std;



// cov fun: exponential + square exponential
// [[Rcpp::export]]
arma::mat EsqeFun( arma::mat distmat, arma::vec covparms ){ //covparms=c(sig2_1,r1,sig2_2,r2)

  int d1 = distmat.n_rows;
  int d2 = distmat.n_cols;
  int j1;
  int j2;
  arma::mat covmat(d1,d2);
  double scaledist;
  double scaledist2;

  for (j1 = 0; j1 < d1; j1++){
    for (j2 = 0; j2 < d2; j2++){
      if ( distmat(j1,j2) == 0 ){
        covmat(j1,j2) = covparms(0) + covparms(2);
      } else {
        scaledist = distmat(j1,j2)/covparms(1);
        scaledist2 = pow(distmat(j1,j2)/covparms(3),2);
        covmat(j1,j2) = covparms(0)*exp(-scaledist)+covparms(2)*exp(-scaledist2);//exponential + square exponential
      }
    }
  }
  return covmat;
}
