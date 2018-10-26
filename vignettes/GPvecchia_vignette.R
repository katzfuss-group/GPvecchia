##### example for calculating likelihood and predictions for general vecchia


###  load GPvecchia package
# library(GPvecchia)
library(devtools)
install_github("katzfuss-group/GPvecchia")


#####################   simulate data    #######################

spatial.dim=2 # number of spatial dimensions
n=20^2  # number of observed locs

# simulate locations
set.seed(10)
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs <- cbind(runif(n),runif(n))
}

# covariance parameters (only matern implemented so far)
sig2=1; range=.1; smooth=1.5
covparms =c(sig2,range,smooth)
covfun <- function(locs) sig2*MaternFun(fields::rdist(locs),covparms)
nuggets=rep(.1,n)

# simulate observations
if(n < 1e4) {
  Om0 <- covfun(locs)+diag(nuggets)
  z=as.numeric(t(chol(Om0))%*%rnorm(n))
} else z=rnorm(n)

# plot simulated data
if(spatial.dim==1) {
  plot(locs,z)
} else {
  quilt.plot(locs,z)
}


#####################   specify Vecchia approx    #######################
# (this only has to be run once)
m=20
vecchia.approx=vecchia_specify(z,locs,m)



#####################   likelihood evaluation    #######################

covparms=c(sig2,range,smooth)
vecchia_likelihood(vecchia.approx,covparms,nuggets)

# ## compare to exact likelihood
# library(mvtnorm)
# dmvnorm(z,mean=rep(0,n),sigma=Om0,log=TRUE)




#####################   prediction   #######################

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
m=20
vecchia.approx=vecchia_specify(z,locs,m,locs.pred=locs.pred)



######  carry out prediction   #######
preds=vecchia_prediction(vecchia.approx,covparms,nuggets)
# returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord



#####################   plot and compare to exact pred   #######################

##  exact prediction
mu.exact=as.numeric(covfun(rbind(locs,locs.pred))[,1:n]%*%solve(Om0,z))
cov.exact=covfun(rbind(locs,locs.pred))-
  covfun(rbind(locs,locs.pred))[,1:n]%*%solve(Om0,t(covfun(rbind(locs,locs.pred))[,1:n]))
var.exact=diag(cov.exact)
cov.exact.pred=cov.exact[n+(1:n.p),n+(1:n.p)]


### plot Vecchia and exact predictions
if(spatial.dim==1) {
  plot(locs,z)
  lines(locs.pred,preds$mu.pred,col='blue')
  lines(locs.pred,preds$mu.pred-1.96*sqrt(preds$var.pred),col='blue',lty=2)
  lines(locs.pred,preds$mu.pred+1.96*sqrt(preds$var.pred),col='blue',lty=2)
  lines(locs.pred,mu.exact[n+(1:n.p)],col='red')
  lines(locs.pred,mu.exact[n+(1:n.p)]-1.96*sqrt(var.exact[n+(1:n.p)]),col='red',lty=2)
  lines(locs.pred,mu.exact[n+(1:n.p)]+1.96*sqrt(var.exact[n+(1:n.p)]),col='red',lty=2)
} else {
  sdrange=range(sqrt(c(preds$var.pred,var.exact[n+(1:n.p)])))
  par(mfrow=c(2,3))
  quilt.plot(locs,z)
  quilt.plot(locs.pred,preds$mu.pred)
  quilt.plot(locs.pred,sqrt(preds$var.pred),zlim=sdrange)
  quilt.plot(locs,z)
  quilt.plot(locs.pred,mu.exact[n+(1:n.p)])
  quilt.plot(locs.pred,sqrt(var.exact[n+(1:n.p)]),zlim=sdrange)
  par(mfrow=c(1,1))
}


### plot entire predictive covariance matrix
Sigma=V2covmat(preds,vecchia.approx)$Sigma.pred
cov.range=quantile(rbind(Sigma,cov.exact.pred),c(.01,.99))
par(mfrow=c(1,2))
image.plot(cov.exact.pred,zlim=cov.range)
image.plot(Sigma,zlim=cov.range)
par(mfrow=c(1,1))



#####################   linear combinations   #######################

### example: subset pred locs

H=sparseMatrix(i=1:(n+n.p),j=1:(n+n.p),x=1)[(n+1):(n+n.p),]

# compute variances of Hy
lincomb.vars=vecchia_lincomb(H,vecchia.approx,preds$V.ord)


### example: overall mean over space

mean(preds$mu.pred)
H=sparseMatrix(i=rep(1,n.p),j=n+(1:n.p),x=1/n.p)

# compute entire covariance matrix of Hy (here, 1x1)
lincomb.cov=vecchia_lincomb(H,vecchia.approx,preds$V.ord,cov.mat=TRUE)


