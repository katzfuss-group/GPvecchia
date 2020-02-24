#define ARMA_DONT_PRINT_ERRORS

#include <iostream>
#include <RcppArmadillo.h>
#include <Rcpp.h>



using namespace Rcpp;
using namespace arma;
using namespace std;

mat EsqeFun( mat distmat, vec covparms );
