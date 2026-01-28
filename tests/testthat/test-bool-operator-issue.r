n = 10^2
locs = cbind(runif(n), runif(n))
covparms = c(1, 0.1, 0.5)
nuggets = rep(0.1, n)
Sigma = exp(-fields::rdist(locs) / covparms[2]) + diag(nuggets)
z = as.numeric(t(chol(Sigma)) %*% rnorm(n));
data = z + 1
vecchia.est = vecchia_estimate(data, locs, theta.ini = c(covparms, nuggets[1]))
