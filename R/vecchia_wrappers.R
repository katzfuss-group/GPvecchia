#' estimate mean and covariance parameters using Vecchia
#'
#' @param data data vector of length n
#' @param locs n x d matrix of spatial locations
#' @param X n x p matrix of trend covariates. default is vector of ones (constant trend).
#'      set to NULL if data are already detrended
#' @param m number of neighbors for vecchia approximation. default is 20
#' @param covmodel covariance model. default is Matern.
#'    see \code{\link{vecchia_likelihood}} for details.
#' @param theta.ini initial values of covariance parameters. nugget variance must be last.
#' @param ... additional input parameters for \code{\link{vecchia_specify}}
#'
#' @return object containing detrended data z, trend coefficients beta.hat,
#'    covariance parameters theta.hat, and other quantities necessary for prediction
#' @examples
#' n=10^2; locs=cbind(runif(n),runif(n))
#' covparms=c(1,.1,.5); nuggets=rep(.1,n)
#' Sigma=exp(-fields::rdist(locs)/covparms[2])+diag(nuggets)
#' z=as.numeric(t(chol(Sigma))%*%rnorm(n)); data=z+1
#' vecchia.est=vecchia_estimate(data,locs,theta.ini=c(covparms,nuggets[1]))
#' @export
vecchia_estimate=function(data,locs,X,m=20,covmodel='matern',theta.ini,...) {

  ## default trend is constant over space (intercept)
  if(missing(X)){

    beta.hat=mean(data)
    z=data-beta.hat
    trend='constant'

  } else if(is.null(X)){
  ## if X=NULL, do not estimate any trend

    beta.hat=c()
    z=data
    trend='none'

  } else {
  ## otherwise, estimate and de-trend

    beta.hat=solve(crossprod(X),crossprod(X,data))
    z=data-X%*%beta.hat
    trend='userspecified'

  }

  ## specify vecchia approximation
  vecchia.approx=vecchia_specify(locs,m,...)

  ## initial covariance parameter values
  if(missing(theta.ini)){
    if(covmodel=='matern'){
      var.res=var(z)
      n=length(z)
      dists.sample=fields::rdist(locs[sample(1:n,min(n,300)),])
      theta.ini=c(.9*var.res,mean(dists.sample)/4,.8,.1*var.res) # var,range,smooth,nugget
    } else {
      stop("Initial cov. parameter values must be specified if cov.model!='matern'")
    }
  }

  ## specify vecchia loglikelihood
  n.par=length(theta.ini)
  negloglik.vecchia=function(logparms)
    -vecchia_likelihood(z,vecchia.approx,exp(logparms)[-n.par],exp(logparms)[n.par],
                        covmodel=covmodel)

  ## find MLE of theta (given beta.hat)
  opt.result=optim(par=log(theta.ini),fn=negloglik.vecchia,
                   control=list(trace=1,maxit=300)) # trace=1 outputs iteration counts
  theta.hat=exp(opt.result$par)

  ## return estimated parameters
  cat('estimated trend coefficients: beta.hat=',beta.hat,"\n",sep='')
  cat('estimated covariance parameters: theta.hat=',theta.hat,"\n",sep=',')
  return(list(z=z,beta.hat=beta.hat,theta.hat=theta.hat,
              trend=trend,locs=locs,covmodel=covmodel))

}



#' make spatial predictions using Vecchia based on estimated parameters
#'
#' @param vecchia.est object returned by \code{\link{vecchia_estimate}}
#' @param locs.pred n.p x d matrix of prediction locations
#' @param X.pred n.p x p matrix of trend covariates at prediction locations.
#'      does not need to be specified if constant or no trend was used in
#'      \code{\link{vecchia_estimate}}
#' @param m number of neighbors for vecchia approximation. default is 30.
#' @param ... additional input parameters for \code{\link{vecchia_specify}}
#'
#' @return object containing prediction means mean.pred and variances var.pred
#' @examples
#' n.p=30^2; grid.oneside=seq(0,1,length=round(sqrt(n.p)))
#' locs.pred=as.matrix(expand.grid(grid.oneside,grid.oneside))
#' vecchia.pred=vecchia_pred(vecchia.est,locs.pred)
#' @export
vecchia_pred=function(vecchia.est,locs.pred,X.pred,m=30,...) {

  ## specify vecchia approximation
  vecchia.approx=vecchia_specify(vecchia.est$locs,m,locs.pred=locs.pred,...)

  ## compute predictions
  theta.hat=vecchia.est$theta.hat
  n.par=length(theta.hat)
  preds=vecchia_prediction(vecchia.est$z,vecchia.approx,
                           theta.hat[-n.par],theta.hat[n.par])

  ## add back the trend if possible
  if(!missing(X.pred)){
    mu.pred=preds$mu.pred+X.pred%*%vecchia.est$beta.hat
  } else if(vecchia.est$trend=='none'){
    mu.pred=preds$mu.pred
  } else if(vecchia.est$trend=='constant'){
    mu.pred=preds$mu.pred+vecchia.est$beta.hat
  } else {
    mu.pred=preds$mu.pred
    warning(paste0('X.pred was not specified, so no trend was ',
                   'added back to the predictions'))
  }

  ## return mean and variance at prediction locations
  return(list(mean.pred=mu.pred,var.pred=preds$var.pred))

}

