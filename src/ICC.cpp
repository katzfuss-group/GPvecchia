#include <iostream>
#include <chrono>
#include <RcppArmadillo.h>
#include <unordered_map>
#include <omp.h>
#include <Rcpp.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace Rcpp;
using namespace arma;
using namespace std;

#define ARMA_DONT_PRINT_ERRORS

#include <iostream>
#include <math.h>       /* sqrt */

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


double dist(rowvec l1, rowvec l2){
  double ssq = 0.0;
  for(int k=0; k<l1.size(); ++k){
    ssq += (l1[k] - l2[k])*(l1[k] - l2[k]);
  }
  return sqrt(ssq);
}




void ic0(NumericVector ptrs, NumericVector inds, NumericVector vals){

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
        cout << "ERROR" << endl;
    }
  }

}


// [[Rcpp::export]]
NumericVector createUcppM(NumericVector ptrs, NumericVector inds, NumericVector cov_vals){

  //cout << "passing values" << endl;
  //cout << cov_vals << endl;

  auto start = chrono::high_resolution_clock::now();
  ic0(ptrs, inds, cov_vals);
  auto finish = chrono::high_resolution_clock::now();

  chrono::duration<double> elapsed = finish - start;
  cout << "IC0 Elapsed time: " << elapsed.count() << " s\n";
  return cov_vals;
}




// [[Rcpp::export]]
NumericVector createUcpp(NumericVector ptrs, NumericVector inds, mat locsord){

  const int nvals = inds.size();
  const int N = ptrs.size();
  NumericVector vals(nvals);

  //auto start = chrono::high_resolution_clock::now();
  for(int i=0; i<N; ++i){
    for(int j=ptrs[i]; j<ptrs[i+1]; ++j){
      double d = dist(locsord.row(i), locsord.row(inds[j]));
      double v = exp(-d);
      vals[j] = v;
    }
  }

  ic0(ptrs, inds, vals);
  //auto finish = chrono::high_resolution_clock::now();

  //chrono::duration<double> elapsed = finish - start;
  //cout << "IC0 Elapsed time: " << elapsed.count() << " s\n";

  return vals;
}





// [[Rcpp::export]]
arma::mat slowIC0( arma::mat A, arma::mat S){
  int N = A.n_rows;

  for( int i = 0; i < N; i++ ){
    for( int j = 0; j < i; j++ ){
      if(S(i,j)==0) {
        A(i,j)=0;
      } else if(i==0 || j==0) {
        A(i,j) = A(i,j)/A(j,j);
      } else {
        A(i,j) = as_scalar((A(i,j) - sum(A(i,span(0,j-1)) % A(j,span(0,j-1))))/A(j,j));
      }
    }

    if( i>0 ){
      A(i,i) = as_scalar(sqrt(A(i,i) - sum(A(i,span(0,i-1)) % A(i,span(0,i-1)))));
    }
  }

  // zero the upper triangle
  for( int i = 0; i < N; i++ ){
    for( int j = i+1; j < N; j++ ){
      A(i,j) = 0;
    }
  }
  return(A);
}

