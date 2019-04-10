rm(list=ls())
setwd("/home/marcin")
library(GpGp); library(Matrix); library(RcppParallel)
library(parallel); library(sparseinv); library(fields)
for (nm in list.files('GPvecchia/R',pattern = "\\.[RrSsQq]$")) {
  #cat(nm,":");
  source(file.path('GPvecchia/R',nm))#; cat("\n")
}
Rcpp::sourceCpp('GPvecchia/src/U_NZentries.cpp')
Rcpp::sourceCpp('GPvecchia/src/MaxMin.cpp')
Rcpp::sourceCpp('~/GPvecchia/src/ICC.cpp')

set.seed(1988)
spatial.dim=2 # number of spatial dimensions
n=100000  # number of observed locs
m=


# simulate locations
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs = cbind(runif(n),runif(n))
}


# set covariance params
sig2=1; range=1; smooth=0.5
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


#mra.options = list(plots=FALSE, r=c(0,16))
#mra.options = list(plots=FALSE, M=4)
#mra.options = list(plots=FALSE, r=c(0,11))
#mra.options = list(plots=FALSE, r=c(0,5), J=10, M=1)
#mra.options = list(plots=TRUE)

#V = vecchia_specify(locs, m=20, conditioning='mra',cond.yz = 'SGV')
#V = vecchia_specify(locs, m=3, conditioning='mra')#, mra.options=mra.options)
#V = vecchia_specify(locs, m=20, conditioning = 'NN', cond.yz="SGV")
V = vecchia_specify(locs, m, conditioning = 'mra')



#Sig.sel = t(apply(V$U.prep$revNNarray, 1, function(r) Sigma[V$ord,V$ord][r[m+1],r]))

##### likelihood evaluation #####
covparms=c(sig2,range,smooth)
#vecchia_loglik = vecchia_likelihood(z,V,covparms,nuggets,covmodel=Sig.sel)
#vecchia_loglik = vecchia_likelihood(z,V,covparms,nuggets,covmodel=Sigma)
lltime = proc.time()
vecchia_loglik.obj= vecchia_likelihood(z,V,covparms,nuggets,covmodel='matern')
print("ll time:")
print(proc.time() - lltime)
print(vecchia_loglik.obj[[1]])
