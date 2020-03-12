rm(list = ls())
source("~/HVLF/aux-functions.r")
library(RandomFields)
library(GPvecchia)

######### set parameters #########
set.seed(1996)
n = 150**2
m = 200
diffusion = 0.0000001
advection = 0.001
#diffusion = 0.00004
#advection = 0.01
frac.obs = 0.1
Tmax = 2


## covariance parameters
sig2 = 0.5; range = .15; smooth = 0.5; 
covparms = c(sig2,range,smooth)
covfun.d = function(D) MaternFun(D, covparms)
covfun <- function(locs) covfun.d(fields::rdist(locs))


## likelihood settings
me.var = 1e-4;
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  if (args[1] %in% c("gauss", "poisson", "logistic", "gamma")) {
    stop("One of the models has to be passed as argument")
  } else {
    data.model = args[1]
  }
} else {
  data.model = "gauss"  
}

lik.params = list(data.model = data.model, sigma = sqrt(me.var))



## generate grid of pred.locs
grid = seq(0,1,length = sqrt(n))
locs = as.matrix(expand.grid(grid,grid)) 

## set initial state

x0 = matrix(sig2*RandomFields::RFsimulate(model = RMmatern(nu = smooth, scale = range),
                                          x = locs[,1], y = locs[,2], spConform=FALSE), ncol = 1)
evolFun = function(x) evol(x, diff = diffusion, adv = advection)
XY = simulate.xy(x0, evolFun, NULL, frac.obs, lik.params, Tmax, sig2 = sig2, smooth = smooth, range = range, locs = locs)


## define Vecchia approximation
vecchia.approx = vecchia_specify(locs, m, conditioning = 'mra')


########## filtering ##########

predsVL = list()

covmodel = getMatCov(vecchia.approx, covfun.d)
mu.tt1 = rep(0, n)
obs.aux = as.numeric(XY$y[[1]])
preds.aux.vl = calculate_posterior_VL( obs.aux, vecchia.approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)

L.tt = getLtt(vecchia.approx, preds.aux.vl)
mu.tt = matrix(preds.aux.vl$mean, ncol = 1)

predsVL[[1]] = list(state = mu.tt, L = L.tt)


#for (t in 2:Tmax) {
t = 2

cat(paste("filtering: t=", t, "\n", sep = ""))
obs.aux = as.numeric(XY$y[[t]])
 
E = evolFun(Matrix::Diagonal(n))
Fmat = E %*% L.tt
 
cat("Calculate covariance elements from factor\n")
M1 = getMatCov(vecchia.approx, Matrix::t(Fmat), factor = TRUE)
cat("... from function\n")
M2 = getMatCov(vecchia.approx, covfun.d)
covmodel = M1 + M2
 
mu.tt1 = E %*% mu.tt

cat("calculate posterior\n")
preds.aux.vl = calculate_posterior_VL( obs.aux, vecchia.approx, prior_mean = mu.tt1, likelihood_model = data.model, covmodel = covmodel, covparms = covparms, likparms = lik.params, return_all = TRUE)
cat("done")
 
L.tt = getLtt(vecchia.approx, preds.aux.vl)
mu.tt = matrix(preds.aux.vl$mean, ncol = 1)

predsVL[[t]] = list(state = mu.tt, L = L.tt)

#}


########## plot results ########## 
## m = M = 0
## for (t in 1:Tmax) {
##   #zrange = range(c(predsVL[[t]][["state"]], unlist(lapply(XY$x, function(t) range(t, na.rm = TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm = TRUE)))))
##   zrange = range(c(unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm=TRUE)))))
##   m = min(m, zrange[1])
##   M = max(M, zrange[2])
## }
## for (t in 1:Tmax) {
##   defpar = par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
##   nna.obs = which(!is.na(XY$y[[t]]))
##   fields::quilt.plot( locs[nna.obs,], XY$y[[t]][nna.obs], zlim = c(m, M), nx = sqrt(n), ny = sqrt(n), main = "obs" )
##   fields::quilt.plot( locs, as.numeric(XY$x[[t]]), zlim = c(m, M), nx = sqrt(n), ny = sqrt(n), main = "truth" )
##   fields::quilt.plot( locs, predsVL[[t]]$state, zlim = zrange, nx = sqrt(n), ny = sqrt(n), main = "prediction" )
##   par(defpar)
## }

