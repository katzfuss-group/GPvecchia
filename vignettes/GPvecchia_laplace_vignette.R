##### example for calculating likelihood and predictions for general vecchia

rm(list = ls())

###  load GPvecchia package
library(devtools)
install_github("katzfuss-group/GPvecchia")
library(GPvecchia)
require(fields)


#####################   simulate data    #######################

data.distr = 'gamma' #options: "gaussian","logistic", "poisson", "gamma"
spatial.dim = 2 # number of spatial dimensions
n = 30^2 # number of observed locs

default_lh_params = list("alpha"=2, "sigma"=sqrt(.1), "beta"=4, "phi"=1)

# simulate locations
# set.seed(14)
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs <- cbind(runif(n),runif(n))
}

# covariance parameters
sig2=1; range=.2; smooth = 1.4
covparms=c(sig2, range, smooth)
covfun <- function(locs) sig2*Matern(fields::rdist(locs),range=range,smoothness=smooth)

# simulate latent process
if(n < 1e4) {
  Om0 <- covfun(locs)
  y=as.numeric(t(chol(Om0))%*%rnorm(n))
} else y=rnorm(n)

# simulate data
if(data.distr=='gaussian'){
  nuggets=rep(default_lh_params$sigma^2,n)
  z = rnorm(n,y,sqrt(nuggets))
} else if(data.distr=='poisson'){
  z = rpois(n, exp(y))
} else if(data.distr=='logistic'){
  z = rbinom(n,1,prob = exp(y)/(1+exp(y)))
} else if(data.distr=='gamma'){
  z = rgamma(n, shape = default_lh_params$alpha, rate = default_lh_params$alpha*exp(-y))
} else if(data.distr=='beta'){
  z = rbeta(n, shape1 =default_lh_params$beta*exp(y), shape2 = default_lh_params$beta)
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
  vecchia.approx=vecchia_specify(locs,m) #IW ordering
  vecchia.approx.lr=vecchia_specify(locs,m,conditioning='firstm') #low rank/mra
} else {
  vecchia.approx=vecchia_specify(locs,m,cond.yz='zy') #RF ordering
  vecchia.approx.lr=vecchia_specify(locs,m,conditioning='firstm') #low rank/mra
}

#####################   prediction at observed locations    ######################
# Perform inference on latent mean with Vecchia Laplace approximation

posterior = calculate_posterior_VL(z, vecchia.approx, likelihood_model=data.distr,
                                   covparms, likparms = default_lh_params, prior_mean = 0)
posterior.lr = calculate_posterior_VL(z,vecchia.approx.lr,likelihood_model=data.distr,
                                   covparms,likparms = default_lh_params, prior_mean = 0)

# Laplace approximation, for comparison (slow)
laplace_cov = covfun(locs)
if(m==0)laplace_cov = diag(diag(laplace_cov))
post_lap = calculate_posterior_laplace(z, data.distr, C = laplace_cov,
                                       likparms = default_lh_params,
                                       return_all = TRUE)

if (spatial.dim==1){
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  ord = order(locs) # order so that lines appear correctly
  y_limits = c(min(y, posterior$mean[ord]), max(y, posterior$mean[ord]))
  plot(locs[ord], y[ord], type = "l", ylim = y_limits )
  lines(locs[ord], posterior$mean[ord], type = "l", col=3, lwd=3)
  lines(locs[ord], post_lap$mean[ord], type = "l", col=2)
  legend("bottomright", legend = c("Latent", "VL", "Laplace"), col= c(1,3,2), lwd=c(1,3,1))
} else if (spatial.dim==2){
  par(mfrow=c(1,3), mar=c(2,2,2,2))
  # ordering unnecessary; we are using a scatter plot rather than lines
  quilt.plot(locs, posterior$mean,  main= "VL m=10")
  quilt.plot(locs, y, main= "Truth")
  quilt.plot(locs, array(post_lap$mean),  main= "Laplace")#, nx=50, ny=50)
  #quilt.plot(locs, posterior.lr$mean,  main= "LR m=10")#, nx=50, ny=50)
}

#####################   evaluation of integrated likelihood   #######################

######  Calculate likelihood   #######
apprx_log_likelihood = vecchia_laplace_likelihood(z,vecchia.approx,likelihood_model=data.distr,covparms,
                             likparms = default_lh_params, prior_mean = locs[,1])
apprx_log_likelihood   # integrated or marginal likelihood of the data, p(z|theta)

### RMSE
VL_MSE = mean((y - posterior$mean)^2)
Lap_MSE = mean((y - post_lap$mean)^2)

#####################   prediction at unobserved locations    #######################

######  specify prediction locations   #######
n.p=30^2
if(spatial.dim==1){  #  1-D case
  locs.pred=matrix(seq(0,1,length=n.p),ncol=1)
} else {   # 2-D case
  grid.oneside=seq(0,1,length=round(sqrt(n.p)))
  locs.pred=as.matrix(expand.grid(grid.oneside,grid.oneside)) # grid of pred.locs
}
n.p=nrow(locs.pred)

######  specify Vecchia approximation   #######
vecchia.approx.pred = vecchia_specify(locs, m=20, locs.pred=locs.pred)
###  carry out prediction
preds = vecchia_laplace_prediction(posterior, vecchia.approx.pred, covparms, pred.mean = 0)
### Repeat for low rank for comparison, if desired
vecchia.approx.pred.lr = vecchia_specify(locs, m=20, locs.pred=locs.pred, conditioning = "firstm")
preds.lr = vecchia_laplace_prediction(posterior.lr, vecchia.approx.pred.lr, covparms, pred.mean = 0)
# returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord

# plotting predicitions
if (spatial.dim==1){
  par(mfrow=c(1,1),mar=c(2,2,2,2))
  ord = order(locs) # order so that lines appear correctly
  plot(locs[ord], y[ord], type = "l", xlim=c(0,1.2), ylim = c(-1,3))
  lines(locs, posterior$mean, type = "p", col=4, lwd=3, lty=1)
  lines(locs.pred, preds$mu.pred, type = "l", col=3, lwd=3, lty=1)
  lines(locs.pred,preds$mu.pred+sqrt(preds$var.pred), type = "l", lty = 3, col=3)
  lines(locs.pred,preds$mu.pred-sqrt(preds$var.pred), type = "l", lty = 3, col=3)
  legend("topleft", legend = c("Latent", "VL: Pred", "VL: 1 stdev"), col= c(1,3,3), lwd=c(1,2,1), lty = c(1,1,3))
} else if (spatial.dim==2){
  par(mfrow=c(1,3), mar=c(2,3,2,3))
  # ordering unnecessary; we are using a scatter plot rather than lines
  quilt.plot(locs, array(post_lap$mean), main= "True Latent", xlim = c(0,1), ylim = c(0,1), nx=64, ny=64)
  quilt.plot(locs.pred, preds$mu.pred,  main= "VL Prediction",nx = 30, ny=30)
  quilt.plot(locs.pred, preds.lr$mu.pred,  main= "LR Prediction",nx = 64, ny=64)
}


#####################   parameter estimation via Nelder-Mead   #######################

# currently set up for gamma distributed data
vecchia.approx=vecchia_specify(locs, m=10, cond.yz = "zy") # for posterior
vecchia.approx.IW = vecchia_specify(locs, m=10) # for integrated likelihood
vl_likelihood = function(x0){
  theta = exp(x0)
  covparms=c(theta[1], theta[2], theta[3]) # sigma range smoothness
  default_lh_params = list("alpha"=theta[4], "sigma"=sqrt(.1))
  prior_mean = 0#log(theta[4])
  # Perform inference on latent mean with Vecchia Laplace approximation
  vll = vecchia_laplace_likelihood(z,vecchia.approx, likelihood_model=data.distr,
                                   covparms, return_all = FALSE,
                                   likparms = default_lh_params, prior_mean = prior_mean,
                                   vecchia.approx.IW = vecchia.approx.IW)
  return(-vll)

}
x0 = log(c(.07,1.88, 1.9, .5))
vl_likelihood(x0)
res = optim(x0, vl_likelihood, method = "Nelder-Mead", control = list("trace" = 4))
exp(res$par[1:4])
vl_likelihood(x0)

