spatial.dim=2
n=50

beta=2
sig2=1; range=.1; smooth=1.5
covparms =c(sig2,range,smooth)
covfun <- function(locs) sig2*MaternFun(fields::rdist(locs),covparms)
nuggets=rep(.1,n)

## data that gives too high smoothness
set.seed(1989)
locs_bad = cbind(runif(n),runif(n))

Om0 <- covfun(locs_bad)+diag(nuggets)
z=as.numeric(t(chol(Om0))%*%rnorm(n))
data_bad=z+beta


## data for which everything works
set.seed(1988)
locs_good = cbind(runif(n),runif(n))

Om0 <- covfun(locs_good)+diag(nuggets)
z=as.numeric(t(chol(Om0))%*%rnorm(n))
data_good=z+beta


test_that("VL Posterior mean for fixed Poisson data returns expected values for first 3", {
  expect_error(vecchia_estimate(data_bad,locs_bad), "The smoothness parameter did not converge. Try a custom optimization routine.")
  expect_output(vecchia_estimate(data_good,locs_good), "Exiting from Nelder Mead minimizer")
})
