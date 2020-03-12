#define ARMA_DONT_PRINT_ERRORS

#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;


//' Calculate the covariance values required by HV for
//' matrix factors passed as sparse matrices
//'
//' @param F factor of a matrix in a sparse format
//' @param revNNarray array with the neighbourhood structure
//' @return matrix with covariance values
//' @export
// [[Rcpp::export]]
arma::mat getMatCovFromFactorCpp(arma::sp_mat F, arma::umat revNNarray){
  
  arma::mat sigSel = arma::zeros<arma::mat>(revNNarray.n_rows, revNNarray.n_cols);

  for(int i=0; i < revNNarray.n_rows; i++) {
    
    arma::uvec r = revNNarray.row( i ).t();
    arma::uvec inds = find( r );
    arma::uvec cols = r.elem( inds ) - 1;
    arma::sp_mat thisCol = F.col( i ).t();
     
    for(int colnum=0; colnum<cols.n_rows; colnum++){
      arma::sp_mat cl = F.col( cols(colnum) ); 
      arma::sp_mat val = thisCol * cl;
      sigSel( i, inds(colnum) ) = val( 0, 0 );
    }
  }
  return sigSel;
}