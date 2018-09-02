
# import VL code somehow
#source("server/importer.R")

data_distr = "poisson"
spatial.dim =2 # number of spatial dimensions
n=100^2 # high res locs
n_c=25^2 # coarse locs

dom = 1 # [0, domn]

# covariance parameters
sig2=1; range=.1; smooth = 1
covfun <- function(locs) sig2*Matern(fields::rdist(locs),range=range,smoothness=smooth)

# simulate locations
set.seed(16)
locs = create_locs(n, dimen = spatial.dim, dom = dom)

# simulate latent process
if(n <= 1e4) {
  Om0 <- covfun(locs)
  y=as.numeric(t(chol(Om0))%*%rnorm(n))
} else y=rnorm(n)





###  Simulate poisson directly
pt_density = n^(1/spatial.dim)
c = dom / pt_density
z = rpois(n, c*exp(y))


# aggregate and plot simulated data
pdf("LGCP_2D_sample.pdf", height = 4, width = 10)
par(mfrow=c(1,3))
if(spatial.dim==1) {

  plot(locs,y, main = "latent")
  plot(locs,z, main = "observed")



} else {
  # aggregate simulated events
  down_sampled = as.image(z, locs, nx = pt_density/2, ny = pt_density/2, FUN = sum)
  coarse_z = as.vector(t(down_sampled$z))
  coarse_idx = unique(down_sampled$ind)
  coarse_locs = cbind(down_sampled$x[coarse_idx[,1]], down_sampled$y[coarse_idx[,2]])

  quilt.plot(locs,y, main = "Latent", nx =pt_density, ny = pt_density )
  quilt.plot(locs,z, main = "Observed", nx=pt_density, ny = pt_density)
  quilt.plot(coarse_locs, coarse_z, main = "Down-Sampled", nx=pt_density/2, ny = pt_density/2)
}
dev.off()

#####################   specify Vecchia approx    #######################
# (this only has to be run once)
m=10
vecchia.approx = vecchia_specify(coarse_z, coarse_locs, m, cond.yz = "zy")
vecchia.approx.sgv = vecchia_specify(coarse_z, coarse_locs, m)
#vecchia.approx.ll=vecchia_specify(z, locs, m, conditioning = "firstm")

#####################   prediction at observed locations    ######################
covparms=c(1, .1, 1.5)
# Perform inference on latent mean with Vecchia Laplace approximation


## Check:  W  = inv(Cov)
posterior = calculate_posterior_VL(vecchia.approx, likelihood_model=data_distr,
                                   covparms, likparms = default_lh_params,
                                   max.iter = 50, return_all = TRUE, prior_mean = 0)
posterior.sgv = calculate_posterior_VL(vecchia.approx.sgv, likelihood_model=data_distr,
                                       covparms, likparms = default_lh_params,
                                       max.iter = 50, return_all = TRUE, prior_mean = 0)

par(mfrow=c(1,2))
quilt.plot(locs,y, main = "Latent", nx =pt_density, ny = pt_density )
quilt.plot(coarse_locs,posterior$mean, main = "Posterior", nx =pt_density/2, ny = pt_density/2 )
