rm(list=ls())
library(GpGp)
library(parallel)
library(Matrix)

setwd("/home/marcin/GPvecchia")
source("R/vecchia_specify.R")
source("R/createU.R")
source("R/vecchia_likelihood.R")
source("R/vecchia_prediction.R")
source("R/RcppExports.R")
source("R/ordering_functions.R")
source("MRA/mraNN.r")
source("R/whichCondOnLatent.R")
source("R/U_sparsity.R")
source("R/NN_kdtree.R")
source("MRA/mra-tree.r")

Rcpp::sourceCpp('src/U_NZentries.cpp')



spatial.dim=2 # number of spatial dimensions
n=15  # number of observed locs
m=2

# simulate locations
#set.seed(1988)
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
  #ord = order_coordinate((locs))
  #locs = matrix(locs[ord])
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

mra.options = list(r=c(10), M=1)
V = vecchia_specify(locs, 20, conditioning='mra', mra.options=mra.options)


##### likelihood evaluation #####
covparms=c(sig2,range,smooth)
vecchia_loglik = vecchia_likelihood(z,V,covparms,nuggets,covmodel=Sigma)
#vecchia_loglik = vecchia_likelihood(z,V,covparms,nuggets,covmodel='matern')

# exact likelihood
const = dim(locs)[1]*log(2*pi)
logdet = as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
quad.form2 = as.numeric(t(z) %*% solve(Sigma) %*% z)
neg2loglike = const + logdet + quad.form2
loglik = -neg2loglike/2


print(loglik)# - vecchia_loglik)
print(vecchia_loglik)
