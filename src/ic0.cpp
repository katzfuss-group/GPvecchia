#define ARMA_DONT_PRINT_ERRORS

#include <iostream>
#include <RcppArmadillo.h>
#include <omp.h>
#include <Rcpp.h>
#include "Matern.h"
#include "dist.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;



double dot_prod(int l1, int u1, int l2, int u2, NumericVector row_inds, NumericVector cells){

  double result = 0.0;
  while(l1<=u1 && l2<=u2){
    if(row_inds[l1]==row_inds[l2]) {
      result += cells[l1]*cells[l2];
      l1++; l2++;
    }
    else if(l1<l2)
      l1++;
    else
      l2++;
  }
  return result;
}



//' Incomplete Cholesky decomposition of a sparse matrix passed in
//' the compressed sparse row format
//'
//' @param ptrs pointers to the beginning of the row
//' @param inds indices of nonzero elements in a row
//' @param vals nonzero values
//' @return vector of the values of the incomplete Cholesky factor
//' @export
// [[Rcpp::export]]
NumericVector ic0(NumericVector ptrs, NumericVector inds, NumericVector vals){

  const int N = ptrs.size()-1;

  for( int i = 0; i<N; ++i ){
    for( int j = ptrs[i]; j<ptrs[i+1]; ++j ){
      int u1 = ptrs[i];
      int u2 = ptrs[inds[j]];
      double dp = dot_prod( u1, ptrs[i+1]-2, u2, ptrs[inds[j] + 1]-2, inds, vals );

      if( inds[j] < i ){
        vals[j] = (vals[j] - dp) / vals[ ptrs[inds[j] + 1] - 1 ];
      }
      else if( inds[j]==i ){
        vals[j] = sqrt( vals[j] - dp );
      }
      else
        Rcout << "ERROR" << endl;
    }
  }
  return vals;
}


// [[Rcpp::export]]
NumericVector createUcppM(NumericVector ptrs, NumericVector inds, NumericVector cov_vals){
  ic0(ptrs, inds, cov_vals);
  return cov_vals;
}




// [[Rcpp::export]]
NumericVector createUcpp(NumericVector ptrs, NumericVector inds, arma::mat locsord, arma::vec covparams){

  const int nvals = inds.size();
  const int N = ptrs.size()-1;
  NumericVector vals(nvals);

  for(int i=0; i<N; ++i){
    for(int j=ptrs[i]; j<ptrs[i+1]; ++j){
      vec D = vec{dist(locsord.row(i), locsord.row(inds[j]))};
      vals[j] = MaternFun(D, covparams)(0,0);
    }
  }

  ic0(ptrs, inds, vals);
  return vals;
}
