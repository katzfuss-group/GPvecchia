##### example for calculating likelihood and predictions for general vecchia

rm(list = ls())

###  load GPvecchia package
# library(GPvecchia)
library(devtools)
install_github("katzfuss-group/GPvecchia")


#####################   simulate data    #######################

data.distr = 'logistic' # options: "gaussian","logistic", "poisson", "gamma"
spatial.dim = 1 # number of spatial dimensions
n=10^2 # number of observed locs

default_lh_params = list("alpha"=2, "sigma"=.1)

# simulate locations
set.seed(12)
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs <- cbind(runif(n),runif(n))
}

# covariance parameters
sig2=1; range=.2; smooth=1.5
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
  print('Error: Not implemented yet.')
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
m=2
vecchia.approx=vecchia_specify(z, locs, m)

#####################   prediction at observed locations    #######################

covparms=c(sig2,range,smooth)
# Perform inference on latent mean with Vecchia Laplace approximation
posterior = calculate_posterior_VL(vecchia.approx, likelihood_model=data.distr, covparms,likparms = default_lh_params)
# Laplace approximation for comparison
post_lap = calculate_posterior_laplace(z, data.distr, C = covfun(locs))

if (spatial.dim==1){
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  ord = order(locs) # order so that lines appear correctly
  plot(locs[ord], y[ord], type = "l")
  lines(locs[ord], posterior$mean[ord], type = "l", col=3, lwd=3, lty=1)
  lines(locs[ord], post_lap$mean[ord], type = "l", col=2)
  legend("topright", legend = c("Latent", "VL", "Laplace"), col= c(1,3,2), lwd=c(1,3,1))
} else if (spatial.dim==2){
  par(mfrow=c(1,3), mar=c(2,3,2,3))
  # ordering unnecessary; we are using a scatter plot rather than lines
  quilt.plot(locs, y, main= "Truth")
  quilt.plot(locs,posterior$mean,  main= "VL")
  quilt.plot(locs,array(post_lap$mean),  main= "Laplace")
}




#####################   prediction at unobserved locations    #######################

# not implemented yet



#####################   evaluation of integrated likelihood   #######################

# not implemented yet
