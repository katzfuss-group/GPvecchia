n=10

## covariance parameters
sig2=0.5; range=.15; smooth=1.5; 
covparms =c(sig2,range,smooth)
covfun <- function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)


## generate grid of pred.locs
grid.oneside=seq(0,1,length=round(sqrt(n)))
locs=as.matrix(expand.grid(grid.oneside,grid.oneside)) 
Sigma = covfun(locs)
covmodel = GPvecchia::getMatCov(vecchia.approx, covfun)