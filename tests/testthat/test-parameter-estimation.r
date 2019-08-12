set.seed(1988)
spatial.dim=2
n=50
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs <- cbind(runif(n),runif(n))
}

beta=2
sig2=1; range=.1; smooth=1.5
covparms =c(sig2,range,smooth)
covfun <- function(locs) sig2*MaternFun(fields::rdist(locs),covparms)
nuggets=rep(.1,n)

Om0 <- covfun(locs)+diag(nuggets)
z=as.numeric(t(chol(Om0))%*%rnorm(n))
data=z+beta

vecchia.est=vecchia_estimate(data,locs)


#test_that("VL Posterior mean for fixed Poisson data returns expected values for first 3", {
#  expect_equal(posterior$mean[1], 0.445597969)
#  expect_equal(posterior$mean[2], 0.814020221)
#  expect_equal(posterior$mean[3], -1.001667008)
#})
