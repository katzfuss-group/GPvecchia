#to compare the original part with rcpp version
U_NZentries_R=function(n,locs,revNNarray, revCondOnLatent,nuggets,covparms){
  m= ncol(revNNarray)-1
  nnp=nrow(locs);
  nentries <- sum( !is.na(revNNarray) )
  Lentries <- matrix(NA,nrow(revNNarray),m+1)
  sig2=covparms[1]
  range=covparms[2]
  smooth=covparms[3]
  # loop over all locations
  for( k in 1:nnp){ ## need input: nnp, revNNarray, revCondOnLatent, nuggets,covfun
    ##using pre-reverse-ordered NNarray and CondOnLatent: revNNarray,revCondOnLatent
    inds    <- revNNarray[k,]
    inds0   <- inds[!is.na(inds)]
    n0 <- length(inds0)
    revCond <- revCondOnLatent[k,!is.na(inds)]

    # subset locations
    locs0   <- locs[inds0,]
    # compute latent covariance matrix
    ##i.e.,C(y_i,y_qy(i)), no nugget added for any, why? because no one condition on observed for now. It occurs in 2D.
    nugmat <- diag(n0)
    diag(nugmat) <- nuggets[inds0]*(!revCond) # has nugget if cond on observed i.e., not in CondOnLatent
    covmat  <- covfun( locs0 ) + nugmat
    # get Cholesky decomposition
    cholmat <- t(chol(covmat))#tryCatch(t(chol(covmat)) , error = function(a) numeric(0) )
    # get last row of inverse Cholesky ### was it just regular Cholesky of PA(row-reversed A) here?
    onevec     <- rep(0,n0)
    onevec[n0] <- 1
    M   <- backsolve( t(cholmat) , onevec )
    # save the entries
    Lentries[k,1:n0] <- t(M)
  }

  Zentries <- rep(NA, 2*n )
  for(i in 1:n){
    Zentries[(2*i-1):(2*i)] <- 1/sqrt(nuggets[i])*c(-1,1)
  }
  return(list(Lentries=Lentries,Zentries=Zentries))
}

