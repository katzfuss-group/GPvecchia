### This is to test three special MRA settings. The first one si the independent blocks case.
### The second one corresponds to the full scale approximation. The third one is the modified predictive process
rm(list=ls())
set.seed(1988)
spatial.dim=1 # number of spatial dimensions
n=200  # number of observed locs

# simulate locations
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs = cbind(runif(n),runif(n))
}

sig2=1; range=.1; smooth=0.5
me.var = 1e-8
covparms =c(sig2,range,smooth)
covfun <- function(locs) sig2*MaternFun(fields::rdist(locs),covparms)
nuggets=rep(me.var,n)

# simulate observations
if(n < 1e4) {
  Sigma = covfun(locs) + diag(nuggets)
  Sigma.c = chol(Sigma)
  z=as.numeric(t(Sigma.c)%*%rnorm(n))
} else z=rnorm(n)




m = 20
mra.options.fulls=list(M=1)# Block-full-scale approximation: Pick the first r_0 = 20/2 = 10 locations in maxmin ordering as knots at resolution 0, then split the remaining locations into groups of maximum size 10+1=11.
mra.options.mpproc=list(r=c(m,1))# Modified predictive process: Pick the first 20 locations in maxmin ordering as knots at resolution
mra.options.indep=list(r=c(0,m))# Independent blocks: The locations would be split into groups of maximum size 20+1=21, and each group is then assumed to be independent.


test_that("normal MRA", {
  expect_warning(vecchia_specify(locs, m, conditioning='mra'),"ordering for the selected conditioning scheme changed to required 'maxmin'")
})
test_that("hierarchical Vecchia allows block full-scale", {
   expect_output(vecchia_specify(locs, m, 'maxmin', conditioning='mra', mra.options=mra.options.fulls, verbose=TRUE), paste("r=",m/2,sep=""))
})
test_that("hierarchical Vecchia allows MPP", {
   expect_output(vecchia_specify(locs, m, 'maxmin', conditioning='mra', mra.options=mra.options.mpproc, verbose=TRUE), paste("r=",m,sep=""))
})
test_that("hierarchical Vecchia allows for independent blocks", {
   expect_output(vecchia_specify(locs, m, 'maxmin', conditioning='mra', mra.options=mra.options.indep, verbose=TRUE), "r=0,")
})

