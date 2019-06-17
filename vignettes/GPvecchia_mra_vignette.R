library(GPvecchia)

#####################   simulate data    #######################


set.seed(1988)
spatial.dim=2# number of spatial dimensions
n=20  # number of observed locs

# simulate locations
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs = cbind(runif(n),runif(n))
}

# define covariance function
sig2=1; range=1; smooth=0.5
me.var = 1e-8
covparms =c(sig2,range,smooth)
covfun = function(locs1, locs2=NULL) {
  if(is.null(locs2)){
    sig2*MaternFun(fields::rdist(locs1),covparms)
  } else {
    c(sig2*MaternFun(matrix(fields::rdist.vec(locs1, locs2), ncol=1),covparms))
  }
}
nuggets=rep(me.var,n)

# simulate observations
if(n <= 1e4) {
  Sigma = covfun(locs) + diag(nuggets)
  Sigma.c = chol(Sigma)
  z=as.numeric(t(Sigma.c)%*%rnorm(n))
} else z=rnorm(n)



#####################   specify Vecchia approx    #######################
m=5


#mra.options = list(plots=FALSE, M=1, J=2) # only the number of levels and splits specified
mra.options = list(plots=FALSE, r=c(0,10)) # independent blocks of size m
#mra.options = list(plots=FALSE, r=c(0,3), J=10, M=1) # independent blocks fo size m
#mra.options = list(plots=TRUE)
browser()
V = vecchia_specify(locs, m, conditioning = 'mra', verbose=TRUE, mra.options = mra.options)



#####################   likelihood evaluation    #######################
# exact likelihood
const = dim(locs)[1]*log(2*pi)
logdet = as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
quad.form2 = as.numeric(t(z) %*% solve(Sigma) %*% z)
neg2loglike = const + logdet + quad.form2
loglik = -neg2loglike/2
print("Likelihood")
print(paste("True: ",loglik,sep=""))

# covariance passed as a matrix with only the required elements
Sig.sel = getMatCov(V, covfun(locs))
vecchia_loglik1 = vecchia_likelihood(z,V,covparms,nuggets,covmodel=Sig.sel)
print(paste("Vecchia: ", vecchia_loglik1, sep=""))

# covariance passed as a function
vecchia_loglik2 = vecchia_likelihood(z,V,covparms,nuggets,covmodel=covfun)
print(paste("Vecchia: ", vecchia_loglik2, sep=""))

# using a predermined covariance
vecchia_loglik3 = vecchia_likelihood(z,V,covparms,nuggets,covmodel='matern')
print(paste("Vecchia: ", vecchia_loglik3, sep=""))



#####################   prediction   #######################

######  specify prediction locations   #######
# n.p=30^2
# if(spatial.dim==1){  #  1-D case
#   locs.pred=matrix(seq(0,1,length=n.p),ncol=1)
# } else {   # 2-D case
#   grid.oneside=seq(0,1,length=round(sqrt(n.p)))
#   locs.pred=as.matrix(expand.grid(grid.oneside,grid.oneside)) # grid of pred.locs
# }
# n.p=nrow(locs.pred)
#
#
# ######  specify Vecchia approximation   #######
# m=20
# vecchia.approx=vecchia_specify(locs,m,locs.pred=locs.pred, conditioning = 'mra')
#
#
#
# ######  carry out prediction   #######
# # returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord
# preds=vecchia_prediction(z,vecchia.approx,covparms,nuggets)
#
#
# #####################   plot and compare to exact pred   #######################
#
# ##  exact prediction
# mu.exact=as.numeric(covfun(rbind(locs,locs.pred))[,1:n]%*%solve(Sigma,z))
# cov.exact=covfun(rbind(locs,locs.pred))-
#   covfun(rbind(locs,locs.pred))[,1:n]%*%solve(Sigma,t(covfun(rbind(locs,locs.pred))[,1:n]))
# var.exact=diag(cov.exact)
# cov.exact.pred=cov.exact[n+(1:n.p),n+(1:n.p)]
#
#
# ### plot Vecchia and exact predictions
# if(spatial.dim==1) {
#   plot(locs,z)
#   lines(locs.pred,preds$mu.pred,col='blue')
#   lines(locs.pred,preds$mu.pred-1.96*sqrt(preds$var.pred),col='blue',lty=2)
#   lines(locs.pred,preds$mu.pred+1.96*sqrt(preds$var.pred),col='blue',lty=2)
#   lines(locs.pred,mu.exact[n+(1:n.p)],col='red')
#   lines(locs.pred,mu.exact[n+(1:n.p)]-1.96*sqrt(var.exact[n+(1:n.p)]),col='red',lty=2)
#   lines(locs.pred,mu.exact[n+(1:n.p)]+1.96*sqrt(var.exact[n+(1:n.p)]),col='red',lty=2)
# } else {
#   sdrange=range(sqrt(c(preds$var.pred,var.exact[n+(1:n.p)])))
#   par(mfrow=c(2,3))
#   quilt.plot(locs,z)
#   quilt.plot(locs.pred,preds$mu.pred)
#   quilt.plot(locs.pred,sqrt(preds$var.pred),zlim=sdrange)
#   quilt.plot(locs,z)
#   quilt.plot(locs.pred,mu.exact[n+(1:n.p)])
#   quilt.plot(locs.pred,sqrt(var.exact[n+(1:n.p)]),zlim=sdrange)
#   par(mfrow=c(1,1))
# }
#
#
