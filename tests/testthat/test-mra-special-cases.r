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
mra.options.mpproc=list(r=c(20,1))# Modified predictive process: Pick the first 20 locations in maxmin ordering as knots at resolution
mra.options.indep=list(r=c(0,m))# Independent blocks: The locations would be split into groups of maximum size 20+1=21, and each group is then assumed to be independent.


test_that("hierarchical Vecchia allows block full-scale", {
  expect_message(vecchia_specify(locs, m, 'maxmin', conditioning='mra', mra.options=mra.options.fulls, verbose=TRUE), "r=0,")
})
test_that("hierarchical Vecchia allows MPP", {
  expect_message(vecchia_specify(locs, m, 'maxmin', conditioning='mra', mra.options=mra.options.mpproc, verbose=TRUE), "r=0,")
})
test_that("hierarchical Vecchia allows for independent blocks", {
  expect_message(vecchia_specify(locs, m, 'maxmin', conditioning='mra', mra.options=mra.options.indep, verbose=TRUE), "r=0,")
})

# V = vecchia_specify(locs, m, 'maxmin', conditioning='mra', mra.options=mra.options.fulls, verbose=TRUE)
# V = vecchia_specify(locs, m, 'maxmin', conditioning='mra', mra.options=mra.options.mpproc, verbose=TRUE)
# V = vecchia_specify(locs, m, 'maxmin', conditioning='mra', mra.options=mra.options.indep, verbose=TRUE)
#
#
#
#
# ##### specify tests #####
# test.scenarios = vector("list", 6)
# test.scenarios[[1]] = list(m=10, r=4, J=2, M=2, plot=TRUE) # normal MRA
# test.scenarios[[2]] = list(r=c(5,1), plot=TRUE) # low-rank
# test.scenarios[[3]] = list(m = 10, M=1, plot=TRUE) # FSA
# test.scenarios[[4]] = list(m = 10, plot=TRUE) # the simplest the user can do
# test.scenarios[[5]] = list(m = 4, r=c(0, 4), plot=TRUE) # indep. blocks
# test.scenarios[[6]] = list(m=3, r=2, J=2, M=2, plot=TRUE) # an MRA case that didn't work
#
#
# problems = c()
# ##### run tests #####
# for(scen.no in 1:length(test.scenarios)){
#   print(paste("===== scenario ", scen.no," =====", sep=""))
#   scen = test.scenarios[[scen.no]]
#   m = scen$m; mra.options = list(J=scen$J, r=scen$r, M=scen$M)
#   V = vecchia_specify(locs, m, conditioning='mra', mra.options=mra.options)
#
#   ##### likelihood evaluation #####
#   covparms=c(sig2,range,smooth)
#   #vecchia_loglik = vecchia_likelihood(z,V,covparms,nuggets,covmodel=Sigma)
#   vecchia_loglik = vecchia_likelihood(z,V,covparms,nuggets,covmodel='matern')
#   # if(abs(vecchia_loglik-answers[scen.no])>1e-10) {
#   #   #print(paste("Test ", scen.no, " passed!", sep=""))
#   #   problems = c(problems, scen.no)
#   # }
#
#   # exact likelihood
#   const = dim(locs)[1]*log(2*pi)
#   logdet = as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
#   quad.form2 = as.numeric(t(z) %*% solve(Sigma) %*% z)
#   neg2loglike = const + logdet + quad.form2
#   loglik = -neg2loglike/2
#
#   print("Likelihood")
#   print(paste("True: ",loglik,sep=""))
#   print(paste("Vecchia: ", vecchia_loglik, sep=""))
# }
# # if(length(problems)==0) {
# #   print("all tests passed!")
# # } else print(paste("problems in tests:", problems, sep=" "))
