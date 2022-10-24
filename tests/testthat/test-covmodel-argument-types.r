set.seed(1988)
spatial.dim=2# number of spatial dimensions
n=9  # number of observed locs
m=4

# simulate locations
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs = cbind(runif(n),runif(n))
}

sig2=1; range=1; smooth=0.5
me.var = 1e-4
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


### MRA conditioning; in this case covmodel can be an n x (m + 1) matrix
V = vecchia_specify(locs, m, conditioning = 'mra')

##### likelihood evaluation #####
Sig.sel = getMatCov(V, covfun(locs))

vecchia_loglik1 = vecchia_likelihood(z,V,covparms,nuggets,covmodel=Sig.sel)
vecchia_loglik2 = vecchia_likelihood(z,V,covparms,nuggets,covmodel=covfun)
vecchia_loglik3 = vecchia_likelihood(z,V,covparms,nuggets,covmodel='matern')
vecchia_loglik4 = vecchia_likelihood(z,V,covparms,nuggets,covmodel=covfun(locs))

test_that("likelihood is the same for all covariance argument types", {
  expect_equal(vecchia_loglik1-vecchia_loglik2, 0)
  expect_equal(vecchia_loglik3-vecchia_loglik2, 0)
  expect_equal(vecchia_loglik4-vecchia_loglik2, 0)
})


#### prediction ####
vecchia_pred1 = vecchia_prediction(z, V, covparms, nuggets, covmodel = Sig.sel, return.values = 'mean')$mu.obs
vecchia_pred2 = vecchia_prediction(z, V, covparms, nuggets, covmodel = covfun, return.values = 'mean')$mu.obs
vecchia_pred3 = vecchia_prediction(z, V, covparms, nuggets, covmodel = 'matern', return.values = 'mean')$mu.obs
vecchia_pred3 = vecchia_prediction(z, V, covparms, nuggets, covmodel = covfun(locs), return.values = 'mean')$mu.obs


test_that("prediction is the same for all covariance argument types", {
  expect_equal(sum(abs(vecchia_pred1-vecchia_pred2)), 0)
  expect_equal(sum(abs(vecchia_pred3-vecchia_pred2)), 0)
  expect_equal(sum(abs(vecchia_pred4-vecchia_pred2)), 0)
})
