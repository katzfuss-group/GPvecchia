suppressPackageStartupMessages(library(RandomFields))
######### set parameters #########
set.seed(1988)
spatial.dim = 2
n = 10**2
frac.obs = 0.3
m = 10


## covariance parameters
sig2 = 0.5; range = .15; smooth = 0.5; 
covparms = c(sig2,range,smooth)
#covfun = function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
#covfun.d = function(D) GPvecchia::MaternFun(D, covparms)


## likelihood settings
me.var = 0.1;
data.model = "gauss"  
lik.params = list(data.model = data.model, sigma = sqrt(me.var), alpha=2)


## generate grid of pred.locs
grid.oneside = seq(0,1,length = round(sqrt(n)))
locs = as.matrix(expand.grid(grid.oneside,grid.oneside)) 


## set initial state
x = matrix(RFsimulate(model = RMmatern(nu = smooth, scale = range, var = sig2),
                      x = locs[,1], y = locs[,2], spConform = FALSE), ncol = 1)
n.obs = round(frac.obs*n)
obs.inds = sample(n, n.obs)
eps = me.var*rnorm(length(obs.inds))
y = rep(NA, n)
y[obs.inds] = x[obs.inds] + eps

nuggets = rep(me.var, n.obs)



######  specify Vecchia approximation   #######
vecchia.approx = vecchia_specify(locs, m, 'maxmin', conditioning='mra', verbose=FALSE)
posterior = vecchia_likelihood(y, vecchia.approx, covparms = covparms, nuggets)

