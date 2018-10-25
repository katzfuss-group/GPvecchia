rm(list=ls())
library(GpGp)
library(parallel)
source("R/vecchia_specify.R")
source("R/RcppExports.R")
source("R/ordering_functions.R")
source("MRA/mraNN.r")
source("R/whichCondOnLatent.R")
source("R/U_sparsity.R")

spatial.dim=2 # number of spatial dimensions
n=200  # number of observed locs
m=4


# simulate locations
set.seed(10)
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs = cbind(runif(n),runif(n))
}


sig2=1; range=.1; smooth=1.5
covparms =c(sig2,range,smooth)
covfun <- function(locs) sig2*MaternFun(fields::rdist(locs),covparms)
nuggets=rep(.1,n)

# simulate observations
if(n < 1e4) {
  Om0 <- covfun(locs)+diag(nuggets)
  z=as.numeric(t(chol(Om0))%*%rnorm(n))
} else z=rnorm(n)


V = vecchia_specify(z, locs, m, conditioning='mra')
