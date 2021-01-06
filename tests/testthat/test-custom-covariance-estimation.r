set.seed(1988)
### set parameters
sig2=1
nugget=1
smooth=.5
range=.5
covparms = c(sig2,range,smooth)

### simulate some data
n.all=20^2
grid.oneside=seq(0,1,length=round(sqrt(n.all)))
locs.all=as.matrix(expand.grid(grid.oneside,grid.oneside))
y.all=sqrt(sig2)*RandomFields::RFsimulate(
  model=RandomFields::RMmatern(nu=smooth,scale=range),
  x=grid.oneside,y=grid.oneside)[[1]]

### obs vs pred data and locs
n=round(n.all*.8)
n.p=n.all-n
obs.ind=sample(1:n.all,n)
z=y.all[obs.ind]+sqrt(nugget)*rnorm(n)
locs=locs.all[obs.ind,]
locs.pred=locs.all[-obs.ind,]


### specify latent Vecchia + IC0 for likelihood approx
m.lik=50
v.approx=vecchia_specify(locs,m.lik,cond.yz='y', conditioning='mra')

covfun = function(locs1, locs2=NULL) {
  if(is.null(locs2)){
    sig2*MaternFun(fields::rdist(locs1),covparms)
  } else {
    c(sig2*MaternFun(matrix(fields::rdist.vec(locs1, locs2), ncol=1),covparms))
  }
}


### evaluate likelihood (need to optimize this for MLE)


negloglik.vecchia=function(logparms) {
    -vecchia_likelihood(z, v.approx, NULL ,exp(logparms)[n.par], covmodel=covmodel)
}
