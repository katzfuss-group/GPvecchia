data.distr = 'poisson' #options: "gaussian","logistic", "poisson", "gamma"
spatial.dim = 2 # number of spatial dimensions
n = 20^2 # number of observed locs

default_lh_params = list("alpha"=2, "sigma"=sqrt(.1), "beta"=.9, "phi"=1.5)

# simulate locations
set.seed(256)
locs <- cbind(runif(n),runif(n))

# covariance parameters
sig2=1; range=.2; smooth = 1.4
covparms=c(sig2, range, smooth)
covfun <- function(locs) sig2*MaternFun(fields::rdist(locs),c(sig2,range,smooth))

# simulate latent process and data
Om0 <- covfun(locs)
y=as.numeric(t(chol(Om0))%*%rnorm(n))
z = rpois(n, exp(y))


#####################   specify Vecchia approx    #######################
# (this only has to be run once)
m=10
vecchia.approx=vecchia_specify(locs,m,cond.yz='zy') #RF ordering

#####################   prediction at observed locations    ######################
# Perform inference on latent mean with Vecchia Laplace approximation

posterior = calculate_posterior_VL(z, vecchia.approx, likelihood_model=data.distr,
                                   covparms, likparms = default_lh_params, prior_mean = 0)

test_that("VL Posterior mean for fixed Poisson data returns expected values for first 3", {
  expect_equal(posterior$mean[1], 0.445597969)
  expect_equal(posterior$mean[2], 0.814020221)
  expect_equal(posterior$mean[3], -1.001667008)
})
