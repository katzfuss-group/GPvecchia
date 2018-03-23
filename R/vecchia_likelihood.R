##  evaluation of the likelihood

vecchia_likelihood=function(vecchia.approx,covparms,nuggets) {
  U=createU_matern(vecchia.approx,covparms,nuggets)
  vecchia_likelihood_U(vecchia.approx,U)
}



## evaluate vecchia likelihood based on U

vecchia_likelihood_U=function(vecchia.approx,U) {
  ### output: loglikelihood (for z)  
  
  y.ind=vecchia.approx$U.prep$y.ind
  
  # constants
  n.y=length(y.ind)
  n.z=nrow(U)-n.y
  const=n.z*log(2*pi)
  
  # numerator
  #z1=t(U[-y.ind,])%*%vecchia.approx$zord      ## t(U[-y.ind,]) does not work when put in package
  z1=Matrix::crossprod(U[-y.ind,],vecchia.approx$zord) ## crossprod(x,y) = t(x)%*%y
  quadform.num=sum(z1^2)
  logdet.num=-2*sum(log(spam::diag(U))) ## base::diag(U) cause error when put in package: Error in diag(U) : no method for coercing this S4 class to a vector
  
  # denominator
  U.y=U[y.ind,]
  z2=as.numeric(U.y%*%z1)
  W=Matrix::tcrossprod(U.y)
  V.ord=t(chol(rev.mat(W)))
  z3=solve(V.ord,rev(z2),system='L')
  quadform.denom=sum(z3^2)
  logdet.denom=-2*sum(log(diag(V.ord)))
  
  # putting everything together
  neg2loglik=logdet.num-logdet.denom+quadform.num-quadform.denom+const
  loglik=-neg2loglik/2
  return(loglik)
  
}


## function to reverse-order a matrix
rev.mat=function(mat) mat[nrow(mat):1,ncol(mat):1]
