rm(list=ls())
#library(GPvecchia); library(Matrix); library(fields)

######### set parameters #########
set.seed(1988)
spatial.dim=2
n=20**2
m=25
frac.obs = 0.3
Tmax = 10

## covariance parameters
sig2=0.5; range=.15; smooth=1.5; me.var=1e-4;
covparms =c(sig2,range,smooth)
covfun <- function(locs) GPvecchia::MaternFun(fields::rdist(locs),covparms)
nuggets=rep(me.var,n)




######### define the temporal evolution function #########
evol = function(state){
  n = dim(as.matrix(state))[1]
  diags = list(rep(0, n), rep(1, n))
  #E = Matrix::bandSparse(n, k=0, diag=list(rep(1,n)))
  E = Matrix::bandSparse(n, k=c(0, 2), diag=diags)
  if(dim(state)[2]>1) return( E %*% state )
  else as.numeric(E %*% as.matrix(state))
}


######### simulate and plot the data #########
simulate.xy = function(x0, E, Q, me, Tmax, seed=NULL){

  if(!is.null(seed)) set.seed(seed)
  n=nrow(x0); nobs=round(frac.obs*n)
  x = list(); y = list()
  
  ind.obs = sample(1:n, nobs)
  x[[1]] = x0; y[[1]] = rep(NA, n)
  y[[1]][ind.obs] = x0[ind.obs] + me*rnorm(nobs)
  
  Qc = chol(Q)
  for(t in 2:Tmax){
    x[[t]] = E(x[[t-1]]) + t(Qc) %*% matrix(rnorm(n), ncol=1)
    ind.obs = sample(1:n, nobs)
    y[[t]] = rep(NA, n)
    y[[t]][ind.obs] = matrix(x[[t]][ind.obs] + me*rnorm(nobs), ncol=1)
  }
  
  return(list(x=x, y=y))
}


## grid of pred.locs
grid.oneside=seq(0,1,length=round(sqrt(n)))
locs=as.matrix(expand.grid(grid.oneside,grid.oneside)) 
## initial state
x0 = matrix(rep(0, n), ncol=1); Sig0 = covfun(locs)
XY = simulate.xy(x0, evol, covfun(locs), me.var, Tmax)


########## filtering ##########


### This function is used mainly for testing.
### It takes the entire covariance matrix and creates
### a matrix of covariances


getL00 = function(vecchia.approx, covfun, locs){
  Sig.sel = getMatCov(vecchia.approx, covfun(locs))
  inds = Filter(function(i) !is.na(i), as.vector(t(vecchia.approx$U.prep$revNNarray - 1)))
  ptrs = c(0, cumsum(apply(vecchia.approx$U.prep$revNNarray, 1, function(r) sum(!is.na(r)))))
  cov.vals = Filter(function(i) !is.na(i), c(t(Sig.sel)))
  vals = GPvecchia::ic0(ptrs, inds, cov.vals)
  Laux = Matrix::sparseMatrix(j=inds, p=ptrs, x=vals, index1=FALSE)
  ro = order(vecchia.approx$ord)
  #return(Laux[ro,ro])
  return(Laux)
}


getLtt = function(vecchia.approx, preds){
  orig.order=order(vecchia.approx$ord)
  V = preds$V
  L.tt = (Matrix::solve(Matrix::t(V), sparse=TRUE)[seq(n, 1), ])[orig.order,]
  return(L.tt)
}


vecchia.approx=GPvecchia::vecchia_specify(locs, m, conditioning='mra', verbose=TRUE)#, mra.options=list(r=c(32)))



preds = list()
L.00 = getL00(vecchia.approx, covfun, locs)
preds[[1]] = list(state=x0, L=L.00)


forecast=x0; Fmat=L.00
for(t in 1:Tmax){
  cat(paste("filtering: t=", t, "\n", sep=""))
  obs.aux = as.numeric(XY$y[[t]])# - forecast)
  covmodel = getMatCov(vecchia.approx, Fmat %*% Matrix::t(Fmat) + covfun(locs))
  preds.aux = calculate_posterior_VL( obs.aux, vecchia.approx, likelihood_model = 'gaussian', covmodel=covmodel, covparms = covparms, likparms = list("sigma"=sqrt(me.var)), return_all = TRUE)
  L.tt = getLtt(vecchia.approx, preds.aux)
  mu.tt = matrix(preds.aux$mean, ncol=1)
  preds[[t]] = list(state=mu.tt, L=L.tt)
  if(t<Tmax){
    forecast = evol(mu.tt)
    Fmat = evol(L.tt)
  }
}


## plot results ##
for(t in 2:Tmax){
  zrange = range(c(preds[[t]][["state"]], unlist(lapply(XY$x, function(t) range(t, na.rm=TRUE))), unlist(lapply(XY$y, function(t) range(t, na.rm=TRUE)))))
  defpar = par(mfrow=c(1, 3), oma=c(0, 0, 2, 0))
  nna.obs = which(!is.na(XY$y[[t]]))
  fields::quilt.plot( locs[nna.obs,], XY$y[[t]][nna.obs], main="data",zlim=zrange, nx=sqrt(n), ny=sqrt(n) )
  fields::quilt.plot( locs, as.numeric(XY$x[[t]]), main="latent", zlim=zrange, nx=sqrt(n), ny=sqrt(n) )
  fields::quilt.plot( locs, preds[[t]]$state, main="Vecchia pred.", zlim=zrange, nx=sqrt(n), ny=sqrt(n) )
  #fields::quilt.plot(locs, preds[[t]]$state, main="prediction", zlim=zrange, nx=sqrt(n), ny=sqrt(n))
  mtext(paste("t=", t, sep=""), outer = TRUE, cex = 1.5)
  par(defpar)
}