#' Vecchia prediction 
#' 
#' @param vecchia.approx @@ MK
#' @param covparms: covariance parameters as a vector
#' @param nuggets: nugget
#' @param var.exact  @@ MK
#' @param covmodel: covariance model, 'matern' by default.
#' 
#' @return Posterior mean, obverved mean, posterior variance, observed variance and V matrix.
#' @examples
#' vecchia_prediction=function(vecchia.approx,covparms,nuggets,var.exact,covmodel='matern')
#' @export



vecchia_prediction=function(vecchia.approx,covparms,nuggets,var.exact,covmodel='matern') {

  # create the U matrix
  U=createU(vecchia.approx,covparms,nuggets,covmodel)

  # compute V for posterior inference
  V.ord=U2V(U,vecchia.approx)

  # compute the posterior mean
  vecchia.mean=vecchia_mean(vecchia.approx,U,V.ord)

  # compute posterior variances
  if(missing(var.exact)) var.exact = (sum(!vecchia.approx$obs)<2*1e4)
  vars.vecchia=vecchia_var(vecchia.approx,V.ord,exact=var.exact)

  # return everything
  return(list(mu.pred=vecchia.mean$mu.pred,mu.obs=vecchia.mean$mu.obs,
               var.pred=vars.vecchia$vars.pred,var.obs=vars.vecchia$vars.obs,
               V.ord=V.ord))

}





######  compute V for posterior inference   #######

U2V=function(U,vecchia.approx){

  U.y=U[vecchia.approx$U.prep$y.ind,]

  if(vecchia.approx$ord.pred!='obspred'){

    W=Matrix::tcrossprod(U.y)
    W.rev=rev.mat(W)
    V.ord=t(chol(W.rev))

  } else {  # for obspred ordering

    n=sum(vecchia.approx$obs)
    n.p=sum(!vecchia.approx$obs)

    # pred columns are unchanged
    V.pr=rev.mat(U.y[,2*n+(1:n.p)]) # in reverse order

    # have to compute cholesky for obs block
    U.oo=U.y[1:n,1:(2*n)]
    A=Matrix::tcrossprod(U.oo)
    A.rev=rev.mat(A)
    V.oor=t(chol(A.rev))

    # combine the blocks into one matrix
    zeromat.sparse=sparseMatrix(c(),c(),dims=c(n.p,n))
    V.or=rbind(zeromat.sparse,V.oor)
    V.ord=as(cbind(V.pr,V.or),'dtCMatrix')

  }

  return(V.ord)
}




######  posterior mean (predictions)   #######

vecchia_mean=function(vecchia.approx,U,V.ord){

  # compute entire posterior mean vector
  #z1=t(U[-vecchia.approx$U.prep$y.ind,])%*%vecchia.approx$zord      ## t(U[...,]) cause error when put in package
  z1=Matrix::crossprod(U[-vecchia.approx$U.prep$y.ind,],vecchia.approx$zord) ## crossprod(x,y) = t(x)%*%y
  z2=as.numeric(U[vecchia.approx$U.prep$y.ind,]%*%z1)
  temp=solve(V.ord,rev(z2))
  mu.rev=-solve(Matrix::t(V.ord),temp) ## base::t() cause error when put in package: Error in t.default(V.ord) : argument is not a matrix
  mu.ord=rev(mu.rev)

  # extract obs and pred parts; return to original ordering
  orig.order=order(vecchia.approx$ord)
  mu=mu.ord[orig.order]
  mu.obs=mu[1:sum(vecchia.approx$obs)]
  if(sum(vecchia.approx$obs)<length(mu)) {
    mu.pred=mu[(sum(vecchia.approx$obs)+1):length(mu)]
  } else mu.pred=c()

  return(list(mu.obs=mu.obs,mu.pred=mu.pred))
}




#' linear combination of predictions
#' 
#' @param H: @@MK
#' @param vecchia.approx @@ MK
#' @param V.ord ordered V matrix
#' @param cov.mat logical TRUE or FALSE
#' 
#' @return Variance of linear combination of predictions.
#' @examples
#' vecchia_lincomb=function(H,vecchia.approx,V.ord,cov.mat=FALSE)
#' @export



######  linear combination   #######

vecchia_lincomb=function(H,vecchia.approx,V.ord,cov.mat=FALSE) {
  H.tt=Matrix::t(H[,rev(vecchia.approx$ord),drop=FALSE])
  temp=Matrix::solve(V.ord,H.tt)
  if(cov.mat){
    lincomb.cov=as.matrix(Matrix::t(temp)%*%temp)
    return(lincomb.cov)
  } else {
    lincomb.vars=as.numeric(Matrix::t(temp*temp)%*%rep(1,ncol(H)))
    return(lincomb.vars)
  }
}


######  selected inverse of a sparse matrix   #######

SelInv=function(cholmat){
  n.all=nrow(cholmat)
  Takahashi_Davis(Q=sparseMatrix(c(),c(),dims=c(n.all,1)),
                  cholQp=cholmat,P=sparseMatrix(i=1:n.all,j=1:n.all,x=1))
}



######  posterior variances   #######

vecchia_var=function(vecchia.approx,V.ord,exact=FALSE){

  # joe added these two lines (4/8/2018)
  n <- length(vecchia.approx$zord)
  n.p <- length(vecchia.approx$ord) - n

  # compute selected inverse and extract variances
  inv.sparse=SelInv(V.ord)
  vars.ord=rev(Matrix::diag(inv.sparse))

  # extract obs and pred parts; return to original ordering
  orig.order=order(vecchia.approx$ord)
  vars=vars.ord[orig.order]
  vars.obs=vars[1:sum(vecchia.approx$obs)]

  if(sum(vecchia.approx$obs)<length(vars)) {
    # if exact prediction variances desired, compute using lincomb
    if(exact & vecchia.approx$ord.pred=='obspred'){
      H=sparseMatrix(i=1:(n+n.p),j=1:(n+n.p),x=1)[(n+1):(n+n.p),]
      vars.pred=vecchia_lincomb(H,vecchia.approx,V.ord)
    } else vars.pred=vars[(sum(vecchia.approx$obs)+1):length(vars)]
  } else vars.pred=c()


  return(list(vars.obs=vars.obs,vars.pred=vars.pred))
}




#' compute covariance matrix from V.ord
#' 
#' @param preds: Results from vecchia_prediction @@MK
#' @param vecchia.approx Results from vecchia_specify
#' 
#' @return Covariance matrix of prediction. Do not do this for large n or n.p!!!
#' @examples
#' V2covmat=function(preds,vecchia.approx)
#' @export

######  compute covariance matrix from V.ord   #######
# do not do this for large n or n.p!!!

V2covmat=function(preds,vecchia.approx){

  orig.order=order(vecchia.approx$ord)
  W=as.matrix(rev.mat(preds$V.ord%*%Matrix::t(preds$V.ord))[orig.order,orig.order])
  Sigma=solve(W)

  n=sum(vecchia.approx$obs)
  n.p=sum(!vecchia.approx$obs)
  Sigma.obs=Sigma[1:n,1:n]
  Sigma.pred=Sigma[n+(1:n.p),n+(1:n.p)]

  return(list(Sigma.obs=Sigma.obs,Sigma.pred=Sigma.pred))
}
