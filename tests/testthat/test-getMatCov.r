setwd("~/HVLF")
rm(list=ls())
source("aux-functions.r")


######### set parameters #########
set.seed(1996)
n = 34**2
m = 50
diffusion = 0.00004
advection = 0.01
frac.obs = 0.1
Tmax = 1


## covariance parameters
sig2=0.5; range=.15; smooth=0.5; 
covparms =c(sig2,range,smooth)
covfun <- function(locs) MaternFun(fields::rdist(locs),covparms)

## likelihood settings
me.var=0.1;
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 1) {
  if (args[1] %in% c("gauss", "poisson", "logistic", "gamma")){
    stop("One of the models has to be passed as argument")
  } else {
    data.model = args[1]
  }
} else {
  data.model = "gauss"  
}

lik.params = list(data.model = data.model, sigma = sqrt(me.var))

## generate grid of pred.locs
grid=seq(0,1,length=sqrt(n))
locs=as.matrix(expand.grid(grid,grid)) 

## set initial state
Q = covfun(locs)
x0 = t(chol(Q)) %*% matrix(rnorm(n), ncol=1);
evolFun = function(x) evol(x, diff=diffusion, adv=advection)
XY = simulate.xy(x0, evolFun, Q, frac.obs, lik.params, Tmax)

vecchia.approx=vecchia_specify(locs, m, conditioning = 'mra')


########## filtering ##########
L.tt = getL00(vecchia.approx, covfun, locs)
mu.tt = matrix(rep(0, n), ncol=1)

obs.aux = as.numeric(XY$y[[1]])

E = evol(diag(n))
Fmat = E %*% L.tt
covmodel = getMatCov(vecchia.approx, Fmat %*% Matrix::t(Fmat) + Q, factor=TRUE)

preds.aux.vl = calculate_posterior_VL( obs.aux, vecchia.approx, prior_mean = mu.tt, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)

L.tt = getLtt(vecchia.approx, preds.aux.vl)
mu.tt = matrix(preds.aux.vl$mean, ncol=1)

preds = list(state=mu.tt, L=L.tt)