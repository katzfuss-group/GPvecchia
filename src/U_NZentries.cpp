#define ARMA_DONT_PRINT_ERRORS

#include <iostream>
#include <RcppArmadillo.h>
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

// distance function for 1 pair of locs
// [[Rcpp::export]]
double dist2(arma::vec l1,arma::vec l2) {
  double dist = norm(l1-l2,2) ;
  return (dist) ;
}
// [[Rcpp::export]]
double dist1(double x, double y){
  double d = sqrt( pow(x - y, 2) );
  return d;
}
//calculate distance matrix for multiple pairs of locs: have to use two separate functions because 1D is vec locs 2D is mat locs
//for 2D
// [[Rcpp::export]]
arma::mat calcPWD2( arma::mat x) {//Rcpp::NumericMatrix
  int outrows = x.n_rows ;
  int outcols = x.n_rows ;
  arma::mat out(outrows, outcols) ;
  for (int arow = 0 ; arow < outrows ; arow++) {
    for (int acol = 0 ; acol < outcols ; acol++) {
      out(arow, acol) = dist2(x.row(arow).t(),x.row(acol).t()) ; //extract row from mat, have to transpose to colvec
    }
  }
  return (out) ;
}
//for 1D
// [[Rcpp::export]]
arma::mat calcPWD1( arma::vec x) {//Rcpp::NumericMatrix
  int outrows = x.size() ;
  int outcols = x.size() ;
  arma::mat out(outrows, outcols) ;
  for (int arow = 0 ; arow < outrows ; arow++) {
    for (int acol = 0 ; acol < outcols ; acol++) {
      out(arow, acol) = dist1(x[arow],x[acol]) ; //extract element from vec
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


// rewrite MaternFun with arma object input/output
// [[Rcpp::export]]
arma::mat MaternFun( arma::mat distmat, arma::vec covparms ){ //covparms=c(sig2,range,smooth)

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
  }else{ // Matern cov with bessel function
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


// [[Rcpp::export]]
int get_nonzero_count(int k, int m){
  int n0;
  if (k < m){
    n0=k+1;
  } else {
    n0=m+1;
  }
  return n0;
}
// [[Rcpp::export]]
int get_nonzero_count_general(const arma::uvec inds){
  //Rcpp defaults NA to 0, so look for values !=0
  int nonzero_counter = 0;
  for(uword idx = 0; idx <inds.n_elem; idx++ ){
    if(inds.at(idx)!=0){nonzero_counter++;}
  }
  return nonzero_counter;
}
// [[Rcpp::export]]
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
// [[Rcpp::export]]
arma::uvec get_idx_vals(int n0, int m, const arma::uvec inds){
  arma::uvec inds00;//
  inds00=inds(span(m+1-n0,m))-ones<uvec>(n0);// shift the indices by -1
  return inds00;
}






// [[Rcpp::export]]
List U_NZentries (int Ncores,int n, const arma::mat& locs, const arma::umat& revNNarray,const arma::mat& revCondOnLatent,const arma::vec& nuggets,const arma::vec& nuggets_obsord, std::string COV, const arma::vec covparms){
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
  //vec revCond;//
  //mat cholmat;//
  //mat locs0;//
  if ((COV!="matern")&(COV!="esqe")){
    cerr << "Error message: That covariance is not implemented"<< endl;
  }

  omp_set_num_threads(Ncores);// selects the number of cores to use.
  // initialized all elements outside of omp part, and claim them as private
#pragma omp parallel for shared(locs,revNNarray,revCondOnLatent,nuggets, nnp,m,Lentries,COV) private(k,M,dist,onevec,covmat,nug,n0,inds,revCon_row,inds00,succ,attempt) default(none) schedule(static)
   for (k = 0; k < nnp; k++) {
// extract a row to work with
     inds=revNNarray.row(k).t();
     revCon_row=revCondOnLatent.row(k).t();


     n0 = get_nonzero_count_general(inds); // for general case
     inds00 = get_idx_vals_general(n0, inds);



// extract locations
     //locs0=locs.rows(inds00); // to extract multiple rows from matrix
     //revCond = revCon_row(span(m+1-n0,m));
     // "%" indicates element-wise multiplication
     nug=nuggets.elem(inds00) % (ones(n0)-revCon_row(span(m+1-n0,m))); // vec is vec, cannot convert to mat
    if (locs.n_cols==1){
      dist=calcPWD1(locs.rows(inds00));
    } else {
      dist=calcPWD2(locs.rows(inds00));
    }


#pragma omp critical
{
// add nugget if cond on observed i.e., not in CondOnLatent
    if( COV=="matern"){
      covmat= MaternFun(dist,covparms) + diagmat(nug) ; // summation from arma
    }else if(COV=="esqe") {
      covmat= EsqeFun(dist,covparms) + diagmat(nug) ; // summation from arma
    }
}

// get Cholesky decomposition : upper triagular
    //cholmat = chol(covmat,"upper");
// get last row of inverse Cholesky
    onevec.resize(n0);
    onevec = zeros(n0);
    onevec[n0-1] = 1;
    try {
      M=solve(chol(covmat,"upper"),onevec);
    } catch(...) {
        M=zeros(n0);
        succ=false;
        attempt=1;
        while(succ==false && attempt<n0-1) {
            // remove the farthest conditioning locs, i.e.the first element corresponding to revNNarray (last elem in NNarray)
            // to avoid resize nug, covmat and dist, use directly in chol() witout pre-defining
            onevec.resize( onevec.size()-1 );
            onevec = zeros(n0-attempt);
            onevec.tail(1) = 1;
            try {
              if (COV=="matern"){
                if (locs.n_cols==1){
                  M.tail(n0-attempt)=solve(chol(MaternFun(calcPWD1(locs.rows(inds00(span(attempt,n0-1)))),covparms)+ diagmat(nuggets.elem(inds00(span(attempt,n0-1)))),"upper"),onevec);
                } else {
                  M.tail(n0-attempt)=solve(chol(MaternFun(calcPWD2(locs.rows(inds00(span(attempt,n0-1)))),covparms)+ diagmat(nuggets.elem(inds00(span(attempt,n0-1)))),"upper"),onevec);
                }
              }else if(COV=="esqe"){
                if (locs.n_cols==1){
                  M.tail(n0-attempt)=solve(chol(EsqeFun(calcPWD1(locs.rows(inds00(span(attempt,n0-1)))),covparms)+ diagmat(nuggets.elem(inds00(span(attempt,n0-1)))),"upper"),onevec);
                } else {
                  M.tail(n0-attempt)=solve(chol(EsqeFun(calcPWD2(locs.rows(inds00(span(attempt,n0-1)))),covparms)+ diagmat(nuggets.elem(inds00(span(attempt,n0-1)))),"upper"),onevec);
                }
              }
              succ=true;
            } catch(...){
              attempt=attempt+1;
            }

        }

    }
// save the entries to matrix
    Lentries(k,span(0,n0-1)) = M.t();
  }

   for (int i = 0; i < n; i++){
     Zentries[2*i] = (-1)/sqrt(nuggets_obsord[i]);
     Zentries[2*i+1] = 1/sqrt(nuggets_obsord[i]);
   }

  List LZentries;
   LZentries["Lentries"]=Lentries;
   LZentries["Zentries"]=Zentries;
  return LZentries;
}

