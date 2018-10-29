##### example for calculating likelihood and predictions for general vecchia

rm(list = ls())

###  load GPvecchia package
# library(GPvecchia)
library(devtools)
install_github("katzfuss-group/GPvecchia")

source("server/importer.R")

#####################   simulate data    #######################

data.distr = 'logistic' #options: "gaussian","logistic", "poisson", "gamma"
spatial.dim =1 # number of spatial dimensions
n=225 # number of observed locs

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
sig2=1; range=.2; smooth = .5
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
} else if(data.distr=='gamma'){
  z = rgamma(n, shape = default_lh_params$alpha, rate = exp(y))
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
m=5
if(spatial.dim==1){
  vecchia.approx=vecchia_specify(locs,m)
} else {
  vecchia.approx=vecchia_specify(locs,m,cond.yz='zy')
}

#####################   prediction at observed locations    ######################
# Perform inference on latent mean with Vecchia Laplace approximation

posterior = calculate_posterior_VL(z,vecchia.approx,likelihood_model=data.distr,
                                   covparms,likparms = default_lh_params)

# Laplace approximation for comparison
post_lap = calculate_posterior_laplace(z, data.distr, C =covfun(locs),
                                       return_all = TRUE)


if (spatial.dim==1){
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  ord = order(locs) # order so that lines appear correctly
  plot(locs[ord], y[ord], type = "l")
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


# get pseudodata and nuggets from the latent y discovered by VL
z_VLpseudo = posterior$t
nuggets_VLpseudo = posterior$D

######  specify Vecchia approximation   #######
vecchia.approx = vecchia_specify(locs, m) # use SGV for likelihood, even if 2D

######  Calculate likelihood   #######
vecchia_likelihood(z_VLpseudo, vecchia.approx,covparms,nuggets_VLpseudo)





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
# m specified earlier
vecchia.approx_pred =vecchia_specify(locs, m, locs.pred=locs.pred)



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

vl_likelihood = function(theta){
  covparms=c(1, theta[1], theta[2]) # sigma range smoothness
  # Perform inference on latent mean with Vecchia Laplace approximation
  posterior.sgv = calculate_posterior_VL(z,vecchia.approx, likelihood_model=data.distr,
                                     covparms, return_all = TRUE, likparms = default_lh_params)
  return(posterior.sgv$vec_lh)
}
x0 = c(.1,.5)
optim(x0,vl_likelihood, method = "Nelder-Mead", control = list("trace" = 4, maxit = 100))

