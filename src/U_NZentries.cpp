#define ARMA_WARN_LEVEL 0

#include <iostream>
#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <Rcpp.h>
#include "Matern.h"
#include "Esqe.h"
#include "dist.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;
using namespace std;





// [[Rcpp::export]]
List U_NZentries (const int Ncores, const arma::uword n, const arma::mat& locs, const arma::umat& revNNarray, const arma::mat& revCondOnLatent, const arma::vec& nuggets, const arma::vec& nuggets_obsord, const std::string covType, const arma::vec covparms){

  if ((covType!="matern") && (covType!="esqe")){
    Rcerr << "Error message: " << covType << " covariance is not implemented"<< endl;
  }
  
  const uword m = revNNarray.n_cols - 1;
  const uword Nlocs = locs.n_rows;
  arma::mat Lentries = zeros( Nlocs, m + 1 );

  #ifdef _OPENMP

    #pragma omp parallel for num_threads(Ncores) shared(Lentries) schedule(static)
  
    for (uword k = 0; k < Nlocs; k++) {
    
      arma::uvec inds = revNNarray.row( k ).t();
      arma::vec revCon_row = revCondOnLatent.row(k).t();
    
      arma::uvec inds00 = inds.elem( find( inds ) ) - 1;
      uword n0 = inds00.n_elem;
    
      arma::vec nug = nuggets.elem(inds00) % (ones(n0)-revCon_row(arma::span(m+1-n0,m)));
      arma::mat dist = calcPWD(locs.rows( inds00 ));
      arma::mat covmat;
    
      if( covType=="matern"){
        covmat = MaternFun(dist,covparms) + diagmat(nug);
      } else if(covType=="esqe") {
        covmat = EsqeFun(dist,covparms) + diagmat(nug);
      }
    
      arma::vec onevec = zeros(n0);
      onevec[n0-1] = 1;

      try {
        arma::mat R = chol(covmat, "upper");
        arma::vec M = solve(R, onevec);
        Lentries(k,arma::span(0,n0-1)) = M.t();
      } catch (const runtime_error& error) {
        Rcerr << "Error message: Cholesky decomposition failed:" << error.what() << endl;
      }
      //arma::vec M = solve(chol(covmat,"upper"),onevec);
      //Lentries(k,arma::span(0,n0-1)) = M.t();
    }

  #else

    for (uword k = 0; k < Nlocs; k++) {
      
      arma::uvec inds = revNNarray.row( k ).t();
      arma::vec revCon_row = revCondOnLatent.row(k).t();
      
      arma::uvec inds00 = inds.elem( find( inds ) ) - 1;
      uword n0 = inds00.n_elem;
      
      arma::vec nug = nuggets.elem(inds00) % (ones(n0)-revCon_row(arma::span(m+1-n0,m)));
      arma::mat dist = calcPWD(locs.rows( inds00 ));
      arma::mat covmat;
      
      if( covType=="matern"){
        covmat = MaternFun(dist,covparms) + diagmat(nug);
      } else if(covType=="esqe") {
        covmat = EsqeFun(dist,covparms) + diagmat(nug);
      }
      
      arma::vec onevec = zeros(n0);
      onevec[n0-1] = 1;


      try {
        arma::mat R = chol(covmat, "upper");
        arma::vec M = solve(R, onevec);
        Lentries(k,arma::span(0,n0-1)) = M.t();
      } catch (const runtime_error& error) {
        Rcerr << "Error message: Cholesky decomposition failed:" << error.what() << endl;
      }
      //arma::vec M = solve(chol(covmat,"upper"),onevec);
      //Lentries(k,arma::span(0,n0-1)) = M.t();
    }

  #endif
  
  


  arma::mat Zentries=zeros(2*n);
  for (uword i = 0; i < n; i++){
    Zentries[2*i] = (-1)/sqrt(nuggets_obsord[i]);
    Zentries[2*i+1] = 1/sqrt(nuggets_obsord[i]);
  }

  return List::create( _["Lentries"] = Lentries, _["Zentries"] = Zentries );
}






// [[Rcpp::export]]
List U_NZentries_mat (int Ncores, const arma::uword n, const arma::mat& locs, const arma::umat& revNNarray, const arma::mat& revCondOnLatent, const arma::vec& nuggets, const arma::vec& nuggets_obsord, arma::mat& covVals, const arma::vec covparms){
  
  const uword m = revNNarray.n_cols - 1;
  const uword Nlocs = locs.n_rows;
  arma::mat Lentries = zeros( Nlocs, m + 1 );

  #ifdef _OPENMP
  
    #pragma omp parallel for num_threads(Ncores) shared(Lentries) schedule(static)
  
    for (uword k = 0; k < Nlocs; k++) {
       
      arma::uvec inds = revNNarray.row( k ).t();
      arma::vec revCon_row = revCondOnLatent.row(k).t();

      arma::uvec inds00 = inds.elem( find( inds ) ) - 1;
      uword n0 = inds00.n_elem;

      arma::mat covmat = covVals.submat(inds00, inds00);

      arma::vec onevec = zeros(n0);
      onevec[n0-1] = 1;


      try {
        arma::mat R = chol(covmat, "upper");
        arma::vec M = solve(R, onevec);
        Lentries(k,arma::span(0,n0-1)) = M.t();
      } catch (const runtime_error& error) {
        Rcerr << "Error message: Cholesky decomposition failed:" << error.what() << endl;
      }
      //arma::vec M = solve(chol(covmat,"upper"),onevec);
      //Lentries(k,arma::span(0,n0-1)) = M.t();
    }
    
  #else

    for (uword k = 0; k < Nlocs; k++) {
      
      arma::uvec inds = revNNarray.row( k ).t();
      arma::vec revCon_row = revCondOnLatent.row(k).t();
      
      arma::uvec inds00 = inds.elem( find( inds ) ) - 1;
      uword n0 = inds00.n_elem;
      
      arma::mat covmat = covVals.submat(inds00, inds00);
      
      arma::vec onevec = zeros(n0);
      onevec[n0-1] = 1;

      try {
        arma::mat R = chol(covmat, "upper");
        arma::vec M = solve(R, onevec);
        Lentries(k,arma::span(0,n0-1)) = M.t();
      } catch (const runtime_error& error) {
        Rcerr << "Error message: Cholesky decomposition failed:" << error.what() << endl;
      }
      //arma::vec M = solve(chol(covmat,"upper"),onevec);
      //Lentries(k,arma::span(0,n0-1)) = M.t();
    }

  #endif
    
  arma::mat Zentries=zeros(2*n);
  for (uword i = 0; i < n; i++){
    Zentries[2*i] = (-1)/sqrt(nuggets_obsord[i]);
    Zentries[2*i+1] = 1/sqrt(nuggets_obsord[i]);
  }

  return List::create( _["Lentries"] = Lentries, _["Zentries"] = Zentries );

}
