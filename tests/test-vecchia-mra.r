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

Rcpp::sourceCpp('src/U_NZentries.cpp')

spatial.dim=2 # number of spatial dimensions
n=8  # number of observed locs
m=4

# simulate locations
set.seed(10)
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs = cbind(runif(n),runif(n))
}

sig2=1; range=.1; smooth=1.5
me.var = 0.1

covparms =c(sig2,range,smooth)
covfun <- function(locs) sig2*MaternFun(fields::rdist(locs),covparms)

nuggets=rep(me.var,n)

# simulate observations
if(n < 1e4) {
  Sigma = covfun(locs)+ diag(nuggets)
  Sigma.c = chol(Sigma)
  z=as.numeric(t(Sigma.c)%*%rnorm(n))
} else z=rnorm(n)


V = vecchia_specify(z, locs, m, conditioning='mra')

##### likelihood evaluation #####
covparms=c(sig2,range,smooth)
vecchia_loglik = vecchia_likelihood(V,covparms,nuggets)

# exact likelihood
const = log(2*pi)
logdet = as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
quad.form = as.numeric(t(z) %*% t(Sigma.c) %*% Sigma.c %*% z)
neg2loglik = const + logdet + quad.form
loglik = -neg2loglik/2


print(loglik)
print(vecchia_loglik)
