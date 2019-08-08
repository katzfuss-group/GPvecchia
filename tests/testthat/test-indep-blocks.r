set.seed(1988)
spatial.dim=2
n=100; m=20

# simulate locations
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs = cbind(runif(n),runif(n))
}

# define covariance function
sig2=1; range=1; smooth=0.5
me.var = 1e-8
covparms =c(sig2,range,smooth)
covfun = function(locs1, locs2=NULL) {
  if(is.null(locs2)){
    sig2*MaternFun(fields::rdist(locs1),covparms)
  } else {
    c(sig2*MaternFun(matrix(fields::rdist.vec(locs1, locs2), ncol=1),covparms))
  }
}
nuggets=rep(me.var,n)

# simulate observations
Sigma = covfun(locs) + diag(nuggets)
Sigma.c = chol(Sigma)
z=as.numeric(t(Sigma.c)%*%rnorm(n))

mra.options = list(plots=FALSE, r=c(0,m)) # independent blocks of size m
V = vecchia_specify(locs, m, conditioning = 'mra', verbose=FALSE, mra.options = mra.options)
vecchia_loglik = vecchia_likelihood(z,V,covparms,nuggets,covmodel='matern')

test_that("hierarchical Vecchia allows for independent blocks", {
  expect_equal(vecchia_loglik, -37.1934687)
})
