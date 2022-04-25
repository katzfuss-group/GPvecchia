#define ARMA_DONT_PRINT_ERRORS

#include <iostream>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <boost/math/special_functions/bessel.hpp>
#include "Sphere.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

mat Rx(double theta) {
  mat res(3, 3, fill::eye);
  res(1, 1) = cos(theta);
  res(2, 2) = cos(theta);
  res(1, 2) = -sin(theta);
  res(2, 1) = sin(theta);
  return res;
}

mat Ry(double theta) {
  mat res(3, 3, fill::eye);
  res(0, 0) = cos(theta);
  res(2, 2) = cos(theta);
  res(0, 2) = sin(theta);
  res(2, 0) = -sin(theta);
  return res;
}

mat Rz(double theta) {
  mat res(3, 3, fill::eye);
  res(0, 0) = cos(theta);
  res(1, 1) = cos(theta);
  res(0, 1) = -sin(theta);
  res(1, 0) = sin(theta);
  return res;
}

double Matern(double d, double mu, double range, double sigma2)
{
	if (abs(d) < 1.0E-10)
	{
		return sigma2;
	}
	return sigma2 * pow(2.0, 1 - mu) / tgamma(mu) * pow(d / range, mu) * boost::math::cyl_bessel_k(mu, d / range);
}

double truncate(double x) {
  return x < -1.0 ? -1.0 : (x > 1.0 ? 1.0 : x);
}

// [[Rcpp::export]]
arma::mat SphereFun( arma::mat distmat, arma::vec covparms, arma::mat locs ){ 

  double alpha1 = covparms(0);
  double alpha2 = covparms(1);
  double alpha3 = covparms(2);
  double beta1 = covparms(3);
  double beta2 = covparms(4);
  double beta3 = covparms(5);
  double kappa = covparms(6);
  double mu = covparms(7);
  double range = covparms(8);
  
  int d1 = distmat.n_rows;
  int d2 = distmat.n_cols;
  int j1;
  int j2;
  arma::mat covmat(d1,d2);

  // compute lon/lat
  // x = cos(lat) * cos(lon)
  // y = cos(lat) * sin(lon)
  // z = sin(lat)
  vector<double> lon(locs.n_rows, 0.0), lat(locs.n_rows, 0.0);
  for (unsigned int i = 0; i < locs.n_rows; i++) {
    lat[i] = asin(truncate(locs(i, 2)));
    lon[i] = acos(truncate(locs(i, 0) / cos(lat[i])));
    if (locs(i, 1) < 0) {
      lon[i] = 2 * 3.14159265359 - lon[i];
    }
  }

  vector<mat> Sigma(locs.n_rows);
  vec gamma(3);
  gamma(0) = 1;
  for (unsigned int i = 0; i < locs.n_rows; i++) {
    gamma(1) = exp(alpha1 + alpha2 * sin(lon[i]) + alpha3 * lat[i]);
    gamma(2) = exp(beta1 + beta2 * sin(lon[i]) + beta3 * lat[i]);
    mat rzry = Rz(lon[i]) * Ry(lat[i]);
    Sigma[i] =  rzry * (Rx(kappa) * diagmat(gamma) * Rx(kappa).t()) * rzry.t();
  }

  for (j1 = 0; j1 < d1; j1++) {
    for (j2 = 0; j2 <= j1; j2++) {
      rowvec diff = (locs.row(j1) - locs.row(j2));
      mat dis_m = 2.0 * diff * (Sigma[j1] + Sigma[j2]).i() * diff.t();
      double dis = pow(dis_m(0,0), 0.5);
      double c = pow(det(Sigma[j1]), 0.25) * pow(det(Sigma[j2]), 0.25) * pow(det((Sigma[j1] + Sigma[j2]) / 2.0), -0.5);
      covmat(j1, j2) = c * Matern(dis, mu, range, 1.0);
      covmat(j2, j1) = covmat(j1, j2);
      if(j1 == j2)
          covmat(j2, j1) += 0.0025;
    }
  }
  return covmat;
}
