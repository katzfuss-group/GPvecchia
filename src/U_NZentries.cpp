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
double dist2(vec l1,vec l2) {
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
mat calcPWD2( mat x) {//Rcpp::NumericMatrix
  int outrows = x.n_rows ;
  int outcols = x.n_rows ;
  mat out(outrows, outcols) ;
  for (int arow = 0 ; arow < outrows ; arow++) {
    for (int acol = 0 ; acol < outcols ; acol++) {
      out(arow, acol) = dist2(x.row(arow).t(),x.row(acol).t()) ; //extract row from mat, have to transpose to colvec
    }
  }
  return (out) ;
}
//for 1D
// [[Rcpp::export]]
mat calcPWD1( vec x) {//Rcpp::NumericMatrix
  int outrows = x.size() ;
  int outcols = x.size() ;
  mat out(outrows, outcols) ;
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


// rewrite MaternFun with arma object input/output
// [[Rcpp::export]]
mat MaternFun( mat distmat, vec covparms ){ //covparms=c(sig2,range,smooth)

  int d1 = distmat.n_rows;
  int d2 = distmat.n_cols;
  int j1;
  int j2;
  mat covmat(d1,d2);
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
List U_NZentries (int Ncores,int n, const mat& locs, const umat& revNNarray,const mat& revCondOnLatent,const vec& nuggets, const vec covparms){
  // initialize the output matrix
  int m= revNNarray.n_cols-1;
  int nnp=locs.n_rows;
  mat Lentries=zeros(nnp,m+1);
  int n0; //number of !is_na elements
  uvec inds;//
  vec revCon_row;//
  uvec inds00;//
  vec nug;//
  mat covmat;//
  vec onevec;//
  vec M;//
  mat dist;//
  int k;//
  mat Zentries=zeros(2*n);
  int attempt;
  bool succ;
  //vec revCond;//
  //mat cholmat;//
  //mat locs0;//

  omp_set_num_threads(Ncores);// selects the number of cores to use.
  // initialized all elements outside of omp part, and claim them as private
  #pragma omp parallel for shared(locs,revNNarray,revCondOnLatent,nuggets,nnp,m,Lentries) private(k,M,dist,onevec,covmat,nug,n0,inds,revCon_row,inds00,succ,attempt) default(none) schedule(static)
   for (k = 0; k < nnp; k++) {
// extract a row to work with
     inds=revNNarray.row(k).t();
     revCon_row=revCondOnLatent.row(k).t();

    if (k < m){
      n0=k+1;
    } else {
      n0=m+1;
    }
     inds00=inds(span(m+1-n0,m))-ones<uvec>(n0);// shift the indices by -1

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
    covmat= MaternFun(dist,covparms) + diagmat(nug) ; // summation from arma
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
              if (locs.n_cols==1){
                M.tail(n0-attempt)=solve(chol(MaternFun(calcPWD1(locs.rows(inds00(span(attempt,n0-1)))),covparms)+ diagmat(nuggets.elem(inds00(span(attempt,n0-1)))),"upper"),onevec);
              } else {
                M.tail(n0-attempt)=solve(chol(MaternFun(calcPWD2(locs.rows(inds00(span(attempt,n0-1)))),covparms)+ diagmat(nuggets.elem(inds00(span(attempt,n0-1)))),"upper"),onevec);
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
     Zentries[2*i] = (-1)/sqrt(nuggets[i]);
     Zentries[2*i+1] = 1/sqrt(nuggets[i]);
   }

  List LZentries;
   LZentries["Lentries"]=Lentries;
   LZentries["Zentries"]=Zentries;
  return LZentries;
}

