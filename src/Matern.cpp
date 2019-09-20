#define BOOST_DISABLE_ASSERTS
#define ARMA_DONT_PRINT_ERRORS
#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
#include <omp.h>
#include <Rcpp.h>
#include "Matern.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace Rcpp;
using namespace arma;
using namespace std;

//' Calculate Matern covariance function
//'
//' @param distmat A matrix with distances between points
//' @param covparms A vector with parameters (marg. variance, range, smoothness)
//' @return A matrix with covariance values corrsponding to the distance matrix
//' @export
// [[Rcpp::export]]
arma::mat MaternFun( arma::mat distmat, arma::vec  covparms ){ //covparms=c(sig2,range,smooth)

int d1 = distmat.n_rows;
int d2 = distmat.n_cols;
int j1;
int j2;
arma::mat covmat(d1,d2);
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
  double alpha = covparms(1);
  double normcon = (pow(alpha, -3)/sqrt(M_PI))*(boost::math::tgamma(2)/boost::math::tgamma(1.5));
  for (j1 = 0; j1 < d1; j1++){
    for (j2 = 0; j2 < d2; j2++){
      if ( distmat(j1,j2) == 0 ){
        covmat(j1,j2) = covparms(0);
      } else {
        scaledist = distmat(j1,j2)/covparms(1);
        //covmat(j1,j2) = covparms(0)*(1+sqrt(3)*scaledist)*exp(-sqrt(3)*scaledist);
        double alpha = covparms(1);
        covmat(j1,j2) = normcon*covparms(0)*0.5*M_PI*pow(alpha,3)*exp(-scaledist)*(1+scaledist);
      }
    }
  }
} else if(covparms(2)==2.5){
  double alpha = covparms(1);
  double normcon = (pow(alpha, -5)/sqrt(M_PI))*(boost::math::tgamma(3)/boost::math::tgamma(2.5));
  for (j1 = 0; j1 < d1; j1++){
    for (j2 = 0; j2 < d2; j2++){
      if ( distmat(j1,j2) == 0 ){
        covmat(j1,j2) = covparms(0);
      } else {
        scaledist = distmat(j1,j2)/covparms(1);
        covmat(j1,j2) = normcon*covparms(0)*0.125*M_PI*pow(alpha,5)*exp(-scaledist)*(3+3*scaledist+scaledist*scaledist);
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
