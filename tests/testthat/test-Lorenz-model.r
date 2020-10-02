#setwd("~/HVLF")
rm(list = ls())
source("~/HVLF/aux-functions.r")
source("~/HVLF/simulations-lorenz/aux-functions-Lorenz.r")
source("~/HVLF/scores.r")


library(Matrix)
library(GPvecchia)


AllParamsAsString = list()



######### set parameters #########
#set.seed(1988)
n = 960 
m = 50
frac.obs = 1.0
Tmax = 1
b = 0.2


## likelihood settings
me.var = 0.01;
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  if (!(args[1] %in% c("gauss", "poisson", "logistic", "gamma"))) {
    stop("One of the models has to be passed as argument")
  } else {
    data.model = args[1]
  }
} else {
  #data.model = "poisson"
  data.model = "gauss"
}
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)


## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(n))
locs = matrix(grid.oneside, ncol = 1)


## set initial state
cat("Loading the moments of the long-run Lorenz\n")
moments = getLRMuCovariance(n, 10, 0.005, 32)
Sig0 = (b**2)*moments[["Sigma"]] + diag(1e-12, n)
mu = b*moments[["mu"]]
x0 = b*getX0(n, 10, 32, 0.005)
#x0 = t(chol(Sig0)) %*% matrix(rnorm(n), ncol=1) + mu

## define Vecchia approximation
cat("Calculating the approximations\n")
mra = GPvecchia::vecchia_specify(locs, m, conditioning = 'mra', ordering = 'maxmin', verbose=TRUE)
exact = GPvecchia::vecchia_specify(locs, nrow(locs)-1, ordering = 'maxmin', conditioning = 'firstm')
low.rank = GPvecchia::vecchia_specify(locs, ncol(mra$U.prep$revNNarray) - 1, ordering = 'maxmin', conditioning = 'firstm', verbose=TRUE)
approximations = list(mra = mra, low.rank = low.rank, exact = exact)


XY = simulate.xy(x0, evolFun, Sigt, frac.obs, lik.params, Tmax)
obs = as.numeric(XY$y[[1]])


predsMRA = list()
covmodel = GPvecchia::getMatCov(mra, as.matrix(Sig0))
preds.aux = GPvecchia::calculate_posterior_VL( obs, mra, prior_mean = mu,
                                       likelihood_model = data.model, covmodel = covmodel,
                                       covparms = NULL, likparms = lik.params, return_all = TRUE)
  
mu.tt = matrix(preds.aux$mean, ncol = 1)
predsMRA[[1]] = list(state = mu.tt, W = preds.aux$W)


predsLR = list()
covmodel = GPvecchia::getMatCov(low.rank, as.matrix(Sig0))
preds.aux = GPvecchia::calculate_posterior_VL( obs, low.rank, prior_mean = mu,
                                      likelihood_model = data.model, covmodel = covmodel,
                                      covparms = NULL, likparms = lik.params, return_all = TRUE)
  
mu.tt = matrix(preds.aux$mean, ncol = 1)
predsLR[[1]] = list(state = mu.tt, W = preds.aux$W)


predsE = list()
covmodel = GPvecchia::getMatCov(exact, as.matrix(Sig0))
preds.aux = GPvecchia::calculate_posterior_VL( obs, exact, prior_mean = mu,
                                      likelihood_model = data.model, covmodel = covmodel,
                                      covparms = NULL, likparms = lik.params, return_all = TRUE)
  
mu.tt = matrix(preds.aux$mean, ncol = 1)
predsE[[1]] = list(state = mu.tt, W = preds.aux$W)


LogSc = calculateLSs(predsMRA, predsLR, predsE, XY$x)
print(LogSc)

