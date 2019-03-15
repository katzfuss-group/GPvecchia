#include <iostream>
#include <RcppArmadillo.h>
#include <unordered_map>
#include <omp.h>
#include <Rcpp.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


#define ARMA_DONT_PRINT_ERRORS

//#include <armadillo>
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
        cout << "ERROR" << endl;
    }
  }

  return vals;
}





// [[Rcpp::export]]
List createUcpp(NumericVector ptrs, NumericVector inds, mat locsord){



  const int nvals = inds.size();
  const int N = ptrs.size();
  NumericVector vals(nvals);
  NumericVector vals2(nvals);

  for(int i=0; i<N; ++i){
    for(int j=ptrs[i]; j<ptrs[i+1]; ++j){
      //cout << "(row,col)=("  << i << "," << inds[j] << ")" << endl;
      //cout << locsord.row(i) << endl;
      //cout << locsord.row(inds[j]) << endl;
      double d = dist(locsord.row(i), locsord.row(inds[j]));
      double v = exp(-d);
      vals[j] = v;
      //cout << v << endl;
    }
  }

  vals2 = ic0(ptrs, inds, vals);

  //cout << "vals" << endl;
  //cout << vals2 << endl;

  List Udata;
  Udata["ptrs"]=ptrs;
  Udata["inds"]=inds;
  Udata["vals"]=vals;
  return Udata;
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
  //cout << A << endl;
  return(A);
}

