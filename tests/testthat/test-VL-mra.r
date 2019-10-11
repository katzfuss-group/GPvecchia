rm(list=ls())

spatial.dim=1

set.seed(1988)

n=10
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs <- cbind(runif(n),runif(n))
}

beta=2
sig2=1; range=.1; smooth=1.5
covparms =c(sig2,range,smooth)
covfun <- function(locs) sig2*MaternFun(fields::rdist(locs),covparms)
nuggets=rep(.1,n)

Om0 <- covfun(locs)+diag(nuggets)
z=as.numeric(t(chol(Om0))%*%rnorm(n))
data=z+beta

# # plot simulated data
# if(spatial.dim==1) {
#   plot(locs,data)
# } else {
#   fields::quilt.plot(locs,data, nx=n, ny=n)
# }


# simulate latent process
y=as.numeric(t(chol(Om0))%*%rnorm(n))

data.model = "logistic"

z = rbinom(n,1,prob = exp(y)/(1+exp(y)))

# # plot simulated data, 1 or 2D
# defpar = par(mfrow=c(1,2))
# if(spatial.dim==1) {
#   plot(locs,y, main = "latent")
#   plot(locs,z, main = "observed")
# } else {
#   fields::quilt.plot(locs,y, main = "Latent")
#   fields::quilt.plot(locs,z, main = "Observed")
# }
# par(defpar)


m=10


######  specify prediction locations   #######
n.p=4^2
locs.pred=matrix(seq(0,1,length=n.p),ncol=1)
n.p=nrow(locs.pred)

######  specify Vecchia approximation   #######
vecchia.approx.pred = vecchia_specify(locs, m, 'maxmin', locs.pred=locs.pred, conditioning='mra', verbose=TRUE)
posterior = calculate_posterior_VL(z,vecchia.approx.pred,likelihood_model=data.model,
                                   covparms = covparms)
###  carry out prediction
#browser()
preds = vecchia_laplace_prediction(posterior, vecchia.approx.pred, covparms)

# plotting predicitions
# if (spatial.dim==1){
#   defpar = par(mfrow=c(1,1))
#   ord = order(locs) # order so that lines appear correctly
#   plot(locs[ord], y[ord], type = "l", xlim=c(0,1.2), ylim = c(-1,3))
#   lines(locs, posterior$mean, type = "p", col=4, lwd=3, lty=1)
#   lines(locs.pred, preds$mu.pred, type = "l", col=3, lwd=3, lty=1)
#   lines(locs.pred,preds$mu.pred+sqrt(preds$var.pred), type = "l", lty = 3, col=3)
#   lines(locs.pred,preds$mu.pred-sqrt(preds$var.pred), type = "l", lty = 3, col=3)
#   legend("topleft", legend = c("Latent", "VL: Pred", "VL: 1 stdev"), 
#          col= c(1,3,3), lwd=c(1,2,1), lty = c(1,1,3))
#   par(defpar)
# } else if (spatial.dim==2){
#   defpar =  par(mfrow=c(1,2))
#   # ordering unnecessary; we are using a scatter plot rather than lines
#   fields::quilt.plot(locs, y, main= "True Latent", 
#                      xlim = c(0,1), ylim = c(0,1), nx=64, ny=64)
#   fields::quilt.plot(locs.pred, preds$mu.pred,  main= "VL Prediction",nx = 30, ny=30)
#   par(defpar)
# }