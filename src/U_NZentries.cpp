#define BOOST_DISABLE_ASSERTS
#define ARMA_DONT_PRINT_ERRORS

#include <iostream>
#include <RcppArmadillo.h>
#include <omp.h>
#include <Rcpp.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "Matern.h"
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;
using namespace std;




double dist(rowvec l1, rowvec l2){
  double ssq = 0.0;
  for(int k=0; k<l1.size(); ++k){
    ssq += (l1[k] - l2[k])*(l1[k] - l2[k]);
  }
  return sqrt(ssq);
}



arma::mat calcPWD( arma::mat x) {
  int outrows = x.n_rows ;
  int outcols = x.n_rows ;
  arma::mat out(outrows, outcols) ;
  for (int arow = 0 ; arow < outrows ; arow++) {
    for (int acol = 0 ; acol < outcols ; acol++) {
      out(arow, acol) = dist(x.row(arow), x.row(acol) );
    }
  }
  return (out) ;
}




double besselK_boost(double v, double x) {
  return boost::math::cyl_bessel_k(v,x);//*exp(x);
}

double gamma_fn_boost(double x) {
  return boost::math::tgamma(x);
}

// cov fun: exponential + square exponential
// [[Rcpp::export]]
arma::mat EsqeFun( arma::mat distmat, arma::vec covparms ){ //covparms=c(sig2_1,r1,sig2_2,r2)

  int d1 = distmat.n_rows;
  int d2 = distmat.n_cols;
  int j1;
  int j2;
  arma::mat covmat(d1,d2);
  double scaledist;
  double scaledist2;

  for (j1 = 0; j1 < d1; j1++){
    for (j2 = 0; j2 < d2; j2++){
      if ( distmat(j1,j2) == 0 ){
        covmat(j1,j2) = covparms(0) + covparms(2);
      } else {
        scaledist = distmat(j1,j2)/covparms(1);
        scaledist2 = pow(distmat(j1,j2)/covparms(3),2);
        covmat(j1,j2) = covparms(0)*exp(-scaledist)+covparms(2)*exp(-scaledist2);//exponential + square exponential
      }
    }
  }
  return covmat;
}

int get_nonzero_count(int k, int m){
  int n0;
  if (k < m){
    n0=k+1;
  } else {
    n0=m+1;
  }
  return n0;
}


int get_nonzero_count_general(const arma::uvec inds){
  //Rcpp defaults NA to 0, so look for values !=0
  int nonzero_counter = 0;
  for(uword idx = 0; idx <inds.n_elem; idx++ ){
    if(inds.at(idx)!=0){nonzero_counter++;}
  }
  return nonzero_counter;
}


arma::uvec get_idx_vals_general(int n0, const arma::uvec inds){
  //Rcpp defaults NA to 0, so look for values !=0
  arma::uvec inds00(n0);
  int nonzero_counter = 0;
  for(uword idx = 0; idx <inds.n_elem; idx++ ){
    if(inds.at(idx)!=0){
      inds00.at(nonzero_counter) = inds.at(idx)-1;// shift the indices by -1
      nonzero_counter++;
      }
  }
  return inds00;
}



arma::uvec get_idx_vals(int n0, int m, const arma::uvec inds){
  arma::uvec inds00;//
  inds00=inds(span(m+1-n0,m))-ones<uvec>(n0);// shift the indices by -1
  return inds00;
}



// [[Rcpp::export]]
List U_NZentries_mat (int Ncores,int n, const arma::mat& locs, const arma::umat& revNNarray,const arma::mat& revCondOnLatent,const arma::vec& nuggets,const arma::vec& nuggets_obsord, arma::mat& COV, const arma::vec covparms){
  // initialize the output matrix
  int m= revNNarray.n_cols-1;
  int nnp=locs.n_rows;
  arma::mat Lentries=zeros(nnp,m+1);
  int n0; //number of !is_na elements
  arma::uvec inds;//
  arma::vec revCon_row;//
  arma::uvec inds00;//
  arma::vec nug;//
  arma::mat covmat;//
  arma::vec onevec;//
  arma::vec M;//
  arma::mat dist;//
  int k;//
  mat Zentries=zeros(2*n);
  int attempt;
  bool succ;


  omp_set_num_threads(Ncores);// selects the number of cores to use.
  // initialized all elements outside of omp part, and claim them as private
#pragma omp parallel for shared(locs,revNNarray,revCondOnLatent,nuggets, nnp,m,Lentries,COV) private(k,M,dist,onevec,covmat,nug,n0,inds,revCon_row,inds00,succ,attempt) default(none) schedule(static)
  for (k = 0; k < nnp; k++) {
    // extract a row to work with

    inds=revNNarray.row(k).t();

        n0 = get_nonzero_count_general(inds); // for general case
    inds00 = get_idx_vals_general(n0, inds);
    // covmat = zeros(n0, n0);
    // for(int l=n0-1; l>=0; --l){
    //   vec auxrow = COV.row(inds[l]);
    //   covmat( l, span(0, auxrow.size()-1) ) = auxrow;
    // }

    covmat = COV.submat(inds00, inds00);
    // covmat = symmatl(covmat);
    // get Cholesky decomposition : upper triagular
    // cholmat = chol(covmat,"upper");
    // get last row of inverse Cholesky
    onevec.resize(n0);
    onevec = zeros(n0);
    onevec[n0-1] = 1;
    M=solve(chol(covmat,"upper"),onevec);

    // save the entries to matrix
    Lentries(k,span(0,n0-1)) = M.t();
  }

//
//   uvec inf = find_nonfinite(nuggets_obsord);
//   uvec fin = find_finite(nuggets_obsord);
//
//   for(const auto& i : inf){
//     Zentries[2*i] = 0;
//     Zentries[2*i+1] = 0;
//   }
//
//   for(const auto& i : fin){
//     Zentries[2*i] = (-1)/sqrt(nuggets_obsord[i]);
//     Zentries[2*i+1] = 1/sqrt(nuggets_obsord[i]);
//   }


  for (int i = 0; i < n; i++){
    Zentries[2*i] = (-1)/sqrt(nuggets_obsord[i]);
    Zentries[2*i+1] = 1/sqrt(nuggets_obsord[i]);
  }

  List LZentries;
  LZentries["Lentries"]=Lentries;
  LZentries["Zentries"]=Zentries;
  return LZentries;
}

// [[Rcpp::export]]
List U_NZentries (int Ncores,int n, const arma::mat& locs, const arma::umat& revNNarray,const arma::mat& revCondOnLatent,const arma::vec& nuggets,const arma::vec& nuggets_obsord, std::string COV, const arma::vec covparms){

  // initialize the output matrix
  int m= revNNarray.n_cols-1;
  int nnp=locs.n_rows;
  arma::mat Lentries=zeros(nnp,m+1);
  int n0; //number of !is_na elements
  arma::uvec inds;
  arma::vec revCon_row;
  arma::uvec inds00;
  arma::vec nug;
  arma::mat covmat;
  arma::vec onevec;
  arma::vec M;
  arma::mat dist;
  int k;
  mat Zentries=zeros(2*n);
  int attempt;
  bool succ;

  if ((COV!="matern")&(COV!="esqe")){
    Rcerr << "Error message: That covariance is not implemented"<< endl;
  }

  omp_set_num_threads(Ncores);// selects the number of cores to use.
  // initialized all elements outside of omp part, and claim them as private

  #pragma omp parallel for shared(locs,revNNarray,revCondOnLatent,nuggets, nnp,m,Lentries,COV) private(k,M,dist,onevec,covmat,nug,n0,inds,revCon_row,inds00,succ,attempt) default(none) schedule(static)

  for (k = 0; k < nnp; k++) {

    inds=revNNarray.row(k).t();
    revCon_row=revCondOnLatent.row(k).t();
    n0 = get_nonzero_count_general(inds); // for general case
    inds00 = get_idx_vals_general(n0, inds);

    nug=nuggets.elem(inds00) % (ones(n0)-revCon_row(span(m+1-n0,m))); // vec is vec, cannot convert to mat
    dist = calcPWD(locs.rows(inds00));

    // Rcout << "k: " << k << endl;
    // Rcout << nuggets.elem(inds00) << endl;
    // Rcout << (ones(n0)-revCon_row(span(m+1-n0,m))) << endl;
    // Rcout << inf * 0 << endl;

    #pragma omp critical
    {
    // add nugget if cond on observed i.e., not in CondOnLatent

    //uvec nnf = find_nonfinite(nug);
    //if(nnf.size()>0){
    //  nug(nnf) = 1e3 * ones<vec>(nnf.size());
    //}


    if( COV=="matern"){
      covmat= MaternFun(dist,covparms) + diagmat(nug) ; // summation from arma
    } else if(COV=="esqe") {
      covmat= EsqeFun(dist,covparms) + diagmat(nug) ; // summation from arma
    }

    }

    onevec.resize(n0);
    onevec = zeros(n0);
    onevec[n0-1] = 1;

    M=solve(chol(covmat,"upper"),onevec);

    // save the entries to matrix
    Lentries(k,span(0,n0-1)) = M.t();
  }

    // uvec inf = find_nonfinite(nuggets_obsord);
    // uvec fin = find_finite(nuggets_obsord);

    // for(const auto& i : inf){
    //   Zentries[2*i] = 0;
    //   Zentries[2*i+1] = 0;
    // }

    // for(const auto& i : fin){
    //   Zentries[2*i] = (-1)/sqrt(nuggets_obsord[i]);
    //   Zentries[2*i+1] = 1/sqrt(nuggets_obsord[i]);
    // }

   //nuggets_obsord << endl;

  for (int i = 0; i < n; i++){
    Zentries[2*i] = (-1)/sqrt(nuggets_obsord[i]);
    Zentries[2*i+1] = 1/sqrt(nuggets_obsord[i]);
  }

  List LZentries;
  LZentries["Lentries"]=Lentries;
  LZentries["Zentries"]=Zentries;
  return LZentries;
}


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
  const int N = ptrs.size();
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
