#define ARMA_WARN_LEVEL 0

#include <iostream>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace std;

double dist(rowvec l1, rowvec l2){
  double ssq = 0.0;
  for(arma::uword k=0; k<l1.size(); ++k){
    ssq += (l1[k] - l2[k])*(l1[k] - l2[k]);
  }
  return sqrt(ssq);
}



arma::mat calcPWD( arma::mat x) {
  arma::uword outrows = x.n_rows ;
  arma::uword outcols = x.n_rows ;
  arma::mat out(outrows, outcols) ;
  for (arma::uword arow = 0 ; arow < outrows ; arow++) {
    for (arma::uword acol = 0 ; acol < outcols ; acol++) {
      out(arow, acol) = dist(x.row(arow), x.row(acol) );
    }
  }
  return (out) ;
}
