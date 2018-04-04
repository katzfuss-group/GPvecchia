##### example for calculating likelihood and predictions for general vecchia

rm(list = ls())
setwd("/Users/danielzilber/Desktop/Katzfuss/GPVecchia/")


###  load GPvecchia package
# library(GPvecchia)
for (nm in list.files('R',pattern = "\\.[RrSsQq]$")) {
  cat(nm,":"); source(file.path('R',nm)); cat("\n")
}
Rcpp::sourceCpp('src/U_NZentries.cpp')



#####################   simulate data    #######################

spatial.dim = 1 # number of spatial dimensions
n=10^2  # number of observed locs

require(fields)

# simulate locations
set.seed(10)
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs <- cbind(runif(n),runif(n))
}

# covariance parameters (only matern implemented so far)
sig2=1; range=.05; smooth=3
covfun <- function(locs) sig2*Matern(rdist(locs),range=range,smoothness=smooth)
nuggets=rep(.1,n) #.001+(locs[,1]<.5)

# simulate observations
if(n < 1e4) {
  Om0 <- covfun(locs) #+diag(nuggets) # include nugget for Gaussian obs
  y=as.numeric(t(chol(Om0))%*%rnorm(n))
} else y=rnorm(n)

# z = y #Gaussian obs
#z = rpois(n, exp(y)) #link for poisson obs
z = rbinom(n,1,prob = exp(y)/(1+exp(y))) # logisitic link for binary obs

# plot simulated data
par(mfrow=c(1,2))
if(spatial.dim==1) {
  plot(locs,y, main = "latent")
  plot(locs,z, main = "observed")

} else {
  quilt.plot(locs,y, main = "Latent")
  quilt.plot(locs,z, main = "Observed")
}

lk_m = define_likelihood_model(model_type ="logistic",locs = locs, obs = z)

#####################   specify Vecchia approx    #######################
# (this only has to be run once)
m=5
vecchia.approx=vecchia_specify(lk_m$z,lk_m$locs, m)#, cond.yz = "z"  )


#####################   likelihood evaluation    #######################

covparms=c(sig2,range,smooth)
posterior = calculate_posterior_VL(vecchia.approx, likelihood_model=lk_m, covparms)
post_lap = calculate_posterior_laplace(lk_m, C = covfun(locs))

if (spatial.dim==1){
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  ord = vecchia.approx$ord
  plot(locs[ord], y[ord], type = "l")
  points(locs[ord], posterior$mean[ord], type = "l", col=3, lwd=3)
  points(locs[ord], post_lap$mean[ord], type = "l", col=2)
  legend("topright", legend = c("Latent", "VL", "Laplace"), col= c(1,3,2), lwd=c(1,3,1))
}else if (spatial.dim==2){
  par(mfrow=c(1,3), mar=c(2,3,2,3))
  quilt.plot(locs, y, main= "Truth")
  quilt.plot(locs,posterior$mean,  main= "VL")
  quilt.plot(locs,array(post_lap$mean),  main= "Laplace")
}
dev.off()

vecchia_likelihood(vecchia.approx,covparms,nuggets)
# currently, only isotropic matern is implemented


# ## compare to exact likelihood
# library(mvtnorm)
# dmvnorm(z,mean=rep(0,n),sigma=Om0,log=TRUE)




#####################   prediction:  Missing for now   #######################
