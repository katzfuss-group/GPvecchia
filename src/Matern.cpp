#define ARMA_DONT_PRINT_ERRORS

#include <iostream>
#include <RcppArmadillo.h>
#include <omp.h>
#include <Rcpp.h>
#include "Matern.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
mat MaternFun( mat distmat, vec covparms ){ //covparms=c(sig2,range,smooth)

int d1 = distmat.n_rows;
int d2 = distmat.n_cols;
int j1;
int j2;
mat covmat(d1,d2);
double scaledist;
if (covparms(2)==0.5) { // Exponential cov function
  for (j1 = 0; j1 < d1; j1++){
    for (j2 = 0; j2 < d2; j2++){
      if ( distmat(j1,j2) == 0 ){
        covmat(j1,j2) = covparms(0);
      } else {
        scaledist = distmat(j1,j2)/covparms(1);
        covmat(j1,j2) = covparms(0)*exp(-scaledist);
      }
    }
  }
} else if(covparms(2)==1.5) {
  for (j1 = 0; j1 < d1; j1++){
    for (j2 = 0; j2 < d2; j2++){
      if ( distmat(j1,j2) == 0 ){
        covmat(j1,j2) = covparms(0);
      } else {
        scaledist = distmat(j1,j2)/covparms(1);
        covmat(j1,j2) = covparms(0)*(1+sqrt(3)*scaledist)*exp(-sqrt(3)*scaledist);
      }
    }
  }
} else if(covparms(2)==2.5){
  for (j1 = 0; j1 < d1; j1++){
    for (j2 = 0; j2 < d2; j2++){
      if ( distmat(j1,j2) == 0 ){
        covmat(j1,j2) = covparms(0);
      } else {
        scaledist = distmat(j1,j2)/covparms(1);
        covmat(j1,j2) = covparms(0)*(1+sqrt(5)*scaledist+(5.0/3.0)*scaledist*scaledist)*exp(-sqrt(5)*scaledist);
      }
    }
  }
} else {// Matern cov with bessel function
  double normcon = covparms(0)/(pow(2.0,covparms(2)-1)*boost::math::tgamma(covparms(2))); //Rf_gammafn(covparms(2)));//
  for (j1 = 0; j1 < d1; j1++){
    for (j2 = 0; j2 < d2; j2++){
      if ( distmat(j1,j2) == 0 ){
        covmat(j1,j2) = covparms(0);
      } else {
        scaledist = distmat(j1,j2)/covparms(1);
        covmat(j1,j2) = normcon*pow( scaledist, covparms(2) )*boost::math::cyl_bessel_k(covparms(2),scaledist); //Rf_bessel_k(scaledist,covparms(2),1.0);
      }
    }
  }
}
return covmat;
}
