rm(list=ls())

set.seed(1988)
nx = 50; ny = 50
n=nx*ny
m=42

# covariance function
sig2=1; range=.15; smooth=1.5
covparms =c(sig2,range,smooth)
covfun <- function(locs) sig2*GPvecchia::MaternFun(fields::rdist(locs),covparms)
nugget_size = 0
models = c("poisson", "logistic", "gamma", "gauss")
frac.obs = 0.1



# generate grid
if(ny==1){
  locs=matrix(runif(n),ncol=1)
} else {
  grid.x=seq(0,1,length=nx)
  grid.y=seq(0,1,length=ny)
  locs=as.matrix(expand.grid(grid.x,grid.y)) # grid of pred.locs
}


# simulate latent process
#nuggets = rep(nugget_size, n)
Om0 <- covfun(locs)#+diag(nuggets)
y=as.numeric(t(chol(Om0))%*%rnorm(n))


for( data.model in models ){
  # simulate data
  if(data.model=='poisson'){
    z = rpois(n, exp(y))
  } else if(data.model=='logistic'){
    z = rbinom(n,1,prob = exp(y)/(1+exp(y)))
  } else if(data.model=='gamma'){
    default_lh_params = list("alpha"=2, "sigma"=sqrt(.1), "beta"=.9, "phi"=1.5)
    z = rgamma(n, shape = default_lh_params$alpha, rate = default_lh_params$alpha*exp(-y))
  } else if(data.model=='gauss'){
    me.std=1.0
    z = rnorm(n, mean=y, sd=me.std)
  } else {
    print('Error: Distribution not implemented yet.')
  }
  obs.inds = sample(1:n, round(n*frac.obs), replace = FALSE)
  z[-obs.inds] = NA


  ######  specify Vecchia approximation   #######
  vecchia.approx = GPvecchia::vecchia_specify(locs, 'maxmin', conditioning='mra', verbose=TRUE, cond.yz='y', mra.options = list(r=c(25, 16, 16, 9, 9)))
  posterior = GPvecchia::calculate_posterior_VL(z, vecchia.approx, likelihood_model=data.model, covparms = covparms)
}

# plot simulated data, 1 or 2D
defpar = par(mfrow=c(1,2))
if(ny==1) {
  ord = order(locs) # order so that lines appear correctly
  plot(locs[ord], y[ord], main = "latent", type="l")
  plot(locs[ord], z[ord], main = "observed")
} else {
  fields::quilt.plot(locs,y, main = "Latent", nx=nx, ny=ny)
  fields::quilt.plot(locs,z, main = "Observed", nx=nx, ny=ny)
}
par(defpar)



# plotting predicitions
if (ny==1){
  defpar = par(mfrow=c(1,1))
  ord = order(locs) # order so that lines appear correctly
  ylim = c(min(z, y, posterior$mean, na.rm=TRUE), max(z, y, posterior$mean, na.rm=TRUE))
  plot(locs[ord], y[ord], xlim=c(0,1), type="l", ylim=ylim)
  points(locs[ord], z[ord], col="red", pch=16)
  lines(locs[ord], posterior$mean[ord], col="blue", lwd=2)
  par(defpar)
} else {
  defpar =  par(mfrow=c(1,2))
  # ordering unnecessary; we are using a scatter plot rather than lines
  fields::quilt.plot(locs, y, main= "True Latent", nx=nx, ny=ny)
  fields::quilt.plot(locs, posterior$mean,  main= "VL Prediction",nx=nx, ny=ny)
  par(defpar)
}