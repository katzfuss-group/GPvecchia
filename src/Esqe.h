#define ARMA_WARN_LEVEL 0

#include <iostream>
#include <RcppArmadillo.h>
#include <Rcpp.h>



using namespace Rcpp;
using namespace arma;
using namespace std;

mat EsqeFun( mat distmat, vec covparms );
