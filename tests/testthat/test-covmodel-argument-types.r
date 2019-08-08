library(GPvecchia)

set.seed(1988)
spatial.dim=2# number of spatial dimensions
n=20**2  # number of observed locs
m=20

# simulate locations
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs = cbind(runif(n),runif(n))
}

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
if(n <= 1e4) {
  Sigma = covfun(locs) + diag(nuggets)
  Sigma.c = chol(Sigma)
  z=as.numeric(t(Sigma.c)%*%rnorm(n))
} else z=rnorm(n)

V = vecchia_specify(locs, m, conditioning = 'mra')

##### likelihood evaluation #####
Sig.sel = getMatCov(V, covfun(locs))

vecchia_loglik1 = vecchia_likelihood(z,V,covparms,nuggets,covmodel=Sig.sel)
vecchia_loglik2 = vecchia_likelihood(z,V,covparms,nuggets,covmodel=covfun)
vecchia_loglik3 = vecchia_likelihood(z,V,covparms,nuggets,covmodel='matern')

test_that("cov. model can be passed as matern, matrix with only the required entries or R function", {
  expect_equal(vecchia_loglik1-vecchia_loglik2, 0)
  expect_equal(vecchia_loglik1-vecchia_loglik2, 0)
  expect_equal(vecchia_loglik1-vecchia_loglik2, 0)
})
