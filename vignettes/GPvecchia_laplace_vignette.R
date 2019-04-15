##### example for calculating likelihood and predictions for general vecchia

rm(list = ls())

###  load GPvecchia package
# library(GPvecchia)
library(devtools)
install_github("katzfuss-group/GPvecchia")

source("../GPVecchia/server/importer.R")

#####################   simulate data    #######################

data.distr = 'gamma' #options: "gaussian","logistic", "poisson", "gamma"
spatial.dim = 1 # number of spatial dimensions
n=100 # number of observed locs

default_lh_params = list("alpha"=2, "sigma"=sqrt(.1))

# simulate locations
set.seed(12)
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs <- cbind(runif(n),runif(n))
  #locs = create_locs(n, dimen = spatial.dim, dom = 1)
}

# covariance parameters
sig2=1; range=.4; smooth = 1.5
covparms=c(sig2, range, smooth)
covfun <- function(locs) sig2*Matern(fields::rdist(locs),range=range,smoothness=smooth)

# simulate latent process
if(n < 1e4) {
  Om0 <- covfun(locs)
  y=as.numeric(t(chol(Om0))%*%rnorm(n))
} else y=rnorm(n)

# simulate data
if(data.distr=='gaussian'){
  nuggets=rep(default_lh_params$sigma,n)
  z = rnorm(n,y,sqrt(nuggets))
} else if(data.distr=='poisson'){
  z = rpois(n, exp(y))
} else if(data.distr=='logistic'){
  z = rbinom(n,1,prob = exp(y)/(1+exp(y)))
} else if(data.distr=='gamma_alt'){
  z = rgamma(n, shape = default_lh_params$alpha, rate = exp(y))
} else if(data.distr=='gamma'){
  z = rgamma(n, shape = default_lh_params$alpha, rate = default_lh_params$alpha*exp(-y))
} else{
  print('Error: Distribution not implemented yet.')
}


# plot simulated data
par(mfrow=c(1,2))
if(spatial.dim==1) {
  plot(locs,y, main = "latent")
  plot(locs,z, main = "observed")
} else {
  quilt.plot(locs,y, main = "Latent")
  quilt.plot(locs,z, main = "Observed")
}


#####################   specify Vecchia approx    #######################
# (this only has to be run once)
m=10
if(spatial.dim==1){
  vecchia.approx=vecchia_specify(locs,m)
} else {
  vecchia.approx=vecchia_specify(locs,m,cond.yz='zy')
  vecchia.approx.lr=vecchia_specify(locs,m,conditioning = "firstm")
}

#####################   prediction at observed locations    ######################
# Perform inference on latent mean with Vecchia Laplace approximation

posterior = calculate_posterior_VL(z,vecchia.approx,likelihood_model=data.distr,
                                   covparms,likparms = default_lh_params, prior_mean = 40)


# Laplace approximation for comparison
nugget = diag(rep(1e-5, n))  #stability
laplace_cov = covfun(locs)+nugget
if(m==0)laplace_cov = diag(diag(laplace_cov))
post_lap = calculate_posterior_laplace(z, data.distr, C =laplace_cov,
                                       return_all = TRUE)


if (spatial.dim==1){
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  ord = order(locs) # order so that lines appear correctly
  y_limits = c(min(y, posterior$mean[ord]), max(y, posterior$mean[ord]))
  plot(locs[ord], y[ord], type = "l", ylim = y_limits )
  lines(locs[ord], posterior$mean[ord], type = "l", col=3, lwd=3)
  lines(locs[ord], post_lap$mean[ord], type = "l", col=2)
  legend("topright", legend = c("Latent", "VL", "Laplace"), col= c(1,3,2), lwd=c(1,3,1))
} else if (spatial.dim==2){
  par(mfrow=c(1,3), mar=c(2,3,2,3))
  # ordering unnecessary; we are using a scatter plot rather than lines
  quilt.plot(locs, y, main= "Truth")
  quilt.plot(locs, posterior$mean,  main= "VL")
  quilt.plot(locs, array(post_lap$mean),  main= "Laplace")
}


#####################   evaluation of integrated likelihood   #######################

######  Calculate likelihood   #######
vecchia.approx=vecchia_specify(locs,m=mv,cond.yz='zy')
vecchia_laplace_likelihood(z,vecchia.approx,likelihood_model=data.distr,covparms,
                             likparms = default_lh_params, prior_mean =0)


#####################   prediction at unobserved locations    #######################

######  specify prediction locations   #######
n.p=20^2
if(spatial.dim==1){  #  1-D case
  locs.pred=matrix(seq(0,1,length=n.p),ncol=1)
} else {   # 2-D case
  grid.oneside=seq(0,1,length=round(sqrt(n.p)))
  locs.pred=as.matrix(expand.grid(grid.oneside,grid.oneside)) # grid of pred.locs
}
n.p=nrow(locs.pred)


# get pseudodata and nuggets from the latent y discovered by VL
z_VLpseudo = posterior$t
nuggets_VLpseudo = posterior$D

######  specify Vecchia approximation   #######
vecchia.approx_pred = vecchia_specify(locs, m, locs.pred=locs.pred)

######  carry out prediction   #######
preds=vecchia_prediction(z_VLpseudo, vecchia.approx_pred,covparms,nuggets_VLpseudo)
# returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord

if (spatial.dim==1){
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  ord = order(locs) # order so that lines appear correctly
  plot(locs[ord], y[ord], type = "l")
  lines(locs.pred, preds$mu.pred, type = "l", col=3, lwd=3, lty=1)
  lines(locs.pred,preds$mu.pred+sqrt(preds$var.pred), type = "l", lty = 3, col=3)
  lines(locs.pred,preds$mu.pred-sqrt(preds$var.pred), type = "l", lty = 3, col=3)
  legend("topright", legend = c("Latent", "VL: Pred", "VL: 1 stdev"), col= c(1,3,3), lwd=c(1,2,1), lty = c(1,1,3))
} else if (spatial.dim==2){
  par(mfrow=c(1,3), mar=c(2,3,2,3))
  # ordering unnecessary; we are using a scatter plot rather than lines
  quilt.plot(locs, y, main= "True Latent")
  quilt.plot(locs.pred,preds$mu.pred,  main= "VL Prediction",nx = 20, ny=20)
  quilt.plot(locs.pred,sqrt(preds$var.pred)*1.96,  main= "VL Uncertainty, 1 stdev")
}




#####################   parameter estimation via Nelder-Mead   #######################

#setup for 1D, not adjusted for ZY ordering for 2D problems
vecchia.approx=vecchia_specify(matrix(locs, ncol=1),m)
vl_likelihood = function(x0){
  theta = exp(x0)
  covparms=c(theta[1], theta[2], theta[3]) # sigma range smoothness
  default_lh_params = list("alpha"=2, "sigma"=sqrt(.1))
  prior_mean = 0#log(theta[4])
  # Perform inference on latent mean with Vecchia Laplace approximation
  vll = vecchia_laplace_likelihood(z,vecchia.approx, likelihood_model=data.distr,
                                   covparms, return_all = TRUE, likparms = default_lh_params, prior_mean = prior_mean)
  return(-vll)

}
x0 = log(c(.07,1.88, 1.9))
vl_likelihood(x0)
res = optim(x0, vl_likelihood, method = "Nelder-Mead", control = list("trace" = 4))
exp(res$par[1:3])
vl_likelihood(x0)
