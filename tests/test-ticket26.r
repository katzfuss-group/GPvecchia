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
source("R/MRA/mraNN.r")
source("R/whichCondOnLatent.R")
source("R/U_sparsity.R")
source("R/NN_kdtree.R")

Rcpp::sourceCpp('src/U_NZentries.cpp')


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



##### specify tests #####
test.scenarios = vector("list", 6)
test.scenarios[[1]] = list(m=10, r=4, J=2, M=2, plot=TRUE) # normal MRA
test.scenarios[[2]] = list(r=c(5,1), plot=TRUE) # low-rank
test.scenarios[[3]] = list(m = 10, M=1, plot=TRUE) # FSA
test.scenarios[[4]] = list(m = 10, plot=TRUE) # the simplest the user can do
test.scenarios[[5]] = list(m = 4, r=c(0, 4), plot=TRUE) # indep. blocks
test.scenarios[[6]] = list(m=3, r=2, J=2, M=2, plot=TRUE) # an MRA case that didn't work


# values obtained from verified runs of the tests
# answers = c(
#   -16.4636835006451,
#   -22.4847597626423,
#   -16.5221542862629,
#   -19.1256294397091,
#   -16.4087289718365,
#   -16.4654724744534,
#   -16.4087289718365
# )

problems = c()
##### run tests #####
for(scen.no in 1:length(test.scenarios)){
  print(paste("===== scenario ", scen.no," =====", sep=""))
  scen = test.scenarios[[scen.no]]
  m = scen$m; mra.options = list(J=scen$J, r=scen$r, M=scen$M)
  V = vecchia_specify(locs, m, conditioning='mra', mra.options=mra.options)

  ##### likelihood evaluation #####
  covparms=c(sig2,range,smooth)
  #vecchia_loglik = vecchia_likelihood(z,V,covparms,nuggets,covmodel=Sigma)
  vecchia_loglik = vecchia_likelihood(z,V,covparms,nuggets,covmodel='matern')
  # if(abs(vecchia_loglik-answers[scen.no])>1e-10) {
  #   #print(paste("Test ", scen.no, " passed!", sep=""))
  #   problems = c(problems, scen.no)
  # }

  # exact likelihood
  const = dim(locs)[1]*log(2*pi)
  logdet = as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  quad.form2 = as.numeric(t(z) %*% solve(Sigma) %*% z)
  neg2loglike = const + logdet + quad.form2
  loglik = -neg2loglike/2

  print("Likelihood")
  print(paste("True: ",loglik,sep=""))
  print(paste("Vecchia: ", vecchia_loglik, sep=""))
}
# if(length(problems)==0) {
#   print("all tests passed!")
# } else print(paste("problems in tests:", problems, sep=" "))
