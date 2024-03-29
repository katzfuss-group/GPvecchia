### Vecchia-Laplace approximation using efficient CPP packages ###


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
######################    Vecchia + Laplace     #####################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###




#' Vecchia Laplace extension of GPVecchia for non-Gaussian data
#'
#' @param z  an array of real numbers representing observations
#' @param vecchia.approx a vecchia object as generated by vecchia_specify()
#' @param likelihood_model text describing likelihood model to be used for observations.  Can be "gaussian","logistic", "poisson", "gamma", or "beta"
#' @param covparms covariance parameters as a vector
#' @param covmodel type of the model covariance or selected elements of the covariance matrix
#' @param likparms likelihood parameters for the likelihood_model, as a list.  Default values are sqrt(.1) for Gaussian noise and 2 for the alpha parameter for Gamma data.
#' @param max.iter maximum iterations to perform
#' @param convg  convergence criteria.  End iterations if the Newton step is this small
#' @param return_all  Return additional posterior covariance terms, TRUE or FALSE
#' @param y_init Specify initial guess for posterior mode
#' @param prior_mean  specify the prior latent mean
#' @param verbose if TRUE messages about the posterior estimation will be displayed
#'
#' @return multivariate normal posterior parameters calculated by the Vecchia-Laplace approximation
#' @examples
#' z=rnorm(10); locs=matrix(1:10,ncol=1); vecchia.approx=vecchia_specify(locs,m=5)
#' calculate_posterior_VL(z,vecchia.approx,"gaussian",covparms=c(1,2,.5))
#' @export
calculate_posterior_VL = function(z,vecchia.approx,
                                  likelihood_model=c("gaussian","logistic", "poisson", "gamma", "beta", "gamma_alt"),
                                  covparms, covmodel='matern', likparms = list("alpha"=2, "sigma"=sqrt(.1)),
                                  max.iter=50, convg = 1e-6, return_all = FALSE, y_init = NA,
                                  prior_mean = rep(0,length(z)), verbose=FALSE){

  likelihood_model <- match.arg(likelihood_model)

    if(is.character(covmodel) && covmodel=='matern' && length(covparms)!=3) {
        stop(sprintf("Matern kernel requires 3 parameters but %d were passed", length(covparms)))
    }
  # Avoid crashes due to bad data
  obs.inds = which(!is.na(z))
  z.obs = z[obs.inds]
  invalid_data_support = switch(likelihood_model,
                      "gaussian" = FALSE,
                      "logistic" = .logistic_data_req(z.obs),
                      "poisson" = .poisson_data_req(z.obs),
                      "gamma" = .gamma_data_reqs(z.obs),
                      "gamma_alt" = .gamma_data_reqs(z.obs),
                      "beta" = .beta_data_reqs(z.obs))
  if(invalid_data_support){
    stop("Data invalid for likelihood type. Make sure that your data lies in the support of the likelihood function.")
  }

  # pull out score and second derivative for readability
  model_funs = switch(likelihood_model,
                      "gaussian" = .gauss_model(likparms),
                      "logistic" = .logistic_model(),
                      "poisson" = .poisson_model(),
                      "gamma" = .gamma_model(likparms),
                      "gamma_alt" = .gamma_model_alt(likparms),
                      "beta" = .beta_model(likparms))


  ell_dbl_prime = model_funs$hess
  ell_prime = model_funs$score
  link_fun = model_funs$link

  # for logging purposes, output scenario

  if(verbose){
      log_comment = cat(paste("Running VL-NR for",likelihood_model, "with m=",
                              ncol(vecchia.approx$U.prep$revNNarray)-1," and sample size",length(z.obs)))
  }

  # record duration of NR
  t_start = Sys.time()

  # init latent variable
  y_o = y_init
  if(any(is.na(y_o))) y_o = prior_mean

  if(length(y_o)>1) y_o = y_o[obs.inds]
  
  convgd = FALSE
  tot_iters = 0
  
  for( i in 1:max.iter){

    y_prev = y_o    # save y_prev for convergence test

    # update pseudo-data and -variances
    D_inv = ell_dbl_prime(y_o, z.obs)
    
    if(any(D_inv<0)) {
      stop("Negative variances occurred, check parameters")
      D_inv = abs(D_inv)
    }

    D = 1/D_inv
      u = ell_prime(y_o,z.obs)
      if (any(!is.finite(u))) stop("Derivative of the loglikehood is infinite. Try different parameter values")
      pseudo.data = rep(NA, length(z))

    pseudo.data[obs.inds] = D * u + y_o - prior_mean[obs.inds]

    nuggets = rep(Inf, length(z))
    nuggets[obs.inds] = D
    # make prediction

    preds=vecchia_prediction(pseudo.data,vecchia.approx,covparms,covmodel=covmodel,
                             nuggets,return.values='meanmat')

    y_o = preds$mu.obs[obs.inds] + prior_mean[obs.inds]
    
    if(is.na(max(abs(y_o-y_prev)))){
      # convergence failed due to machine precision?
      fail_comment = message(paste("VL-NR hit NA on iteration ",tot_iters,", convergence failed."))

      y_o = y_prev
      break
    }
    if (max(abs(y_o-y_prev))<convg){
      convgd = TRUE
      tot_iters = i
      break
    }
    tot_iters = tot_iters +1
  }
  t_end = Sys.time()
  LV_time = as.double(difftime(t_end, t_start, units = "secs"))

  optional_data = list()
  if(return_all){
    # return additional information if needed, can be slow
    orig.order=order(vecchia.approx$ord)
    V.ord=preds$V.ord
    # if ZY_liklihood works, dont need W or V?
    W = methods::as(revMat(V.ord%*%Matrix::t(V.ord))[orig.order,orig.order], 'dgCMatrix')
    if (vecchia.approx$cond.yz=="zy"){
      n = length(y_o)
      V.ord = V.ord[1:n, 1:n]
      #W = W[(n+1):(2*n), (n+1):(2*n)]
    }
    optional_data = list( "V"=V.ord, "W" = W)
  }

  posterior_data = list("mean" = preds$mu.obs + prior_mean, "cnvgd" = convgd, "runtime" = LV_time,
                        "iter" = tot_iters, "t"=pseudo.data+prior_mean, "D" = D,
                        "prediction" = preds, "data_link"=link_fun, "model_llh"=  model_funs$llh,
                        "prior_mean"=prior_mean)

  return ( c(posterior_data, optional_data))
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
######################    Laplace + True Covar ######################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


.calculate_posterior_laplace = function(z, likelihood_model, C,  likparms = list("alpha"=2, "sigma"=sqrt(.1)),
                                       convg = 1e-6, return_all = FALSE, prior_mean = 0,verbose=FALSE){
  # pull out score and second derivative for readability
  model_funs = switch(likelihood_model,
                      "gaussian" = .gauss_model(likparms),
                      "logistic" = .logistic_model(),
                      "poisson" = .poisson_model(),
                      "gamma" = .gamma_model(likparms),
                      "gamma_alt" = .gamma_model_alt(likparms),
                      "beta" = .beta_model(likparms))
  ell_dbl_prime = model_funs$hess
  ell_prime = model_funs$score

  log_comment = paste("Running Laplace for",likelihood_model, "with sample size", length(z) )
  if(verbose){
      cat(log_comment)
  }

  t_start = Sys.time()
  y_o = rep(1, length(z))
  tot_iters=0
  # begin NR iteration
  for( i in 1:50){
    D_inv = ell_dbl_prime(y_o, z)
    D = Matrix::sparseMatrix(i=1:length(y_o), j = 1:length(y_o), x= 1/D_inv)
    u =  ell_prime(y_o,z)
    t = D%*%u+y_o - prior_mean
    y_prev = y_o
    y_o = t - D%*%Matrix::solve(D+C,t) + prior_mean # woodbury morrison of previous line
    #W = D_inv +  C_inv
    tot_iters = i
    if (max(abs(y_o-y_prev))<convg) break
  } # end iterate
  t_end = Sys.time()
  Lap_time = as.double(difftime(t_end, t_start, units = "secs"))
  if(return_all){
    # Caclulating sd is expensive
    D_inv_mat = Matrix::sparseMatrix(i=1:length(y_o), j = 1:length(y_o), x= D_inv)
    W = D_inv_mat +  Matrix::solve(C)
    sd_posterior = sqrt(Matrix::diag(Matrix::solve(W)))
    return (list("mean" = y_o, "W"=W,"sd" = sd_posterior, "iter"=tot_iters,
                 "C"=C, "t" = t, "D" = D, "runtime"=Lap_time))
  }
  return (list("mean" = y_o, "iter"=tot_iters))
}



#################  Logistic   #########################
.logistic_model = function(){
  logistic_llh = function(y_o, z) sum(z*y_o-log(1+exp(y_o)))
  logistic_hess = function(y_o, z) exp(y_o)/(1+exp(y_o))^2
  logistic_score = function(y_o, z) z - exp(y_o)/(1+exp(y_o))
  logistic_link = function(y) exp(y)/(1+exp(y))
  # return object with all components of model
  return(list("hess" = logistic_hess, "score"=logistic_score,"llh" = logistic_llh, "link" = logistic_link))
}

.logistic_data_req = function(z){
  return(!all(z %in% c(0,1)))
}

#################  Poisson  #########################
.poisson_model = function(){
  pois_llh = function(y_o, z) sum(z*y_o -exp(y_o)-lfactorial(z))
  pois_hess =function(y_o, z) exp(y_o)
  pois_score = function(y_o, z) z-exp(y_o)
  pois_link = function(y) exp(y)
  return(list("hess" = pois_hess, "score"=pois_score,"llh" = pois_llh, "link" = pois_link))
}

.poisson_data_req = function(z){
  negative = any(z<0)
  non_integral = any(z%%1>0)
  return(negative | non_integral)
}

#################  Gaussian  #########################
.gauss_model = function(likparms){
  # default nugget sd = .3
  sigma = ifelse("sigma" %in% names(likparms),likparms$sigma, sqrt(.1))
  gauss_llh = function(y_o, z) sum(-.5*(z-y_o)^2/sigma^2) -length(y_o)*(log(sigma)+log(2*pi)/2)
  gauss_hess = function(y_o, z)  rep(1/sigma^2, length(y_o))
  gauss_score = function(y_o, z) (z-y_o)/sigma^2
  gauss_link = function(y) y
  return(list("hess" = gauss_hess, "score"=gauss_score, "llh"=gauss_llh, "link" = gauss_link))
}

#################  Gamma  #########################
.gamma_model_alt = function(likparms){
  alpha = ifelse("alpha" %in% names(likparms),likparms$alpha, 2)
  gamma_hess = function(y_o, z)  z*exp(y_o)
  gamma_score = function(y_o, z) -z*exp(y_o)+ alpha
  #gamma_llh = function(y_o, z) sum(-y_o*z + (alpha-1)*log(z) +alpha*log(y_o)-n*log(gamma(alpha))) # canonical link
  gamma_llh = function(y_o, z) sum(-exp(y_o)*z + (alpha-1)*log(z) +alpha*y_o - lgamma(alpha)) # log link
  gamma_link = function(y) alpha/exp(y)
  return(list("hess" = gamma_hess, "score"=gamma_score, "llh" = gamma_llh, "link" = gamma_link))
}


# This alternate parameterization encodes the latent y directly as the mean.
# Assuming a fixed alpha implies a beta
.gamma_model = function(likparms){
  alpha = ifelse("alpha" %in% names(likparms),likparms$alpha, 2)
  gamma_hess = function(y_o, z)  alpha*z*exp(-y_o)
  gamma_score = function(y_o, z) alpha*(z*exp(-y_o)-1)
  gamma_llh = function(y_o, z) sum(-alpha*z*exp(-y_o) + (alpha-1)*log(z) -alpha*y_o +alpha*log(alpha) - lgamma(alpha)) # log link

  #gamma_score_alpha = function(a, y_o, z) sum(-exp(-y_o)*z+log(z)-y_o+log(a)+1-digamma(a))
  #gamma_hess_alpha = function(a, y_o, z) length(z)*(1/a-trigamma(a))

  gamma_link = function(y) exp(y)
  return(list("hess" = gamma_hess, "score"=gamma_score, "llh" = gamma_llh, "link" = gamma_link))
}

.gamma_data_reqs = function(z){
  return(any(z<=0))
}

################# Beta #########################

.beta_model = function(likparms){
  beta = ifelse("beta" %in% names(likparms),likparms$beta, .5)
  # store score and hessian functions, update
  beta_hess = function(y_o, z)  -exp(y_o)*beta*(log(z)-digamma(exp(y_o)*beta) + digamma(beta*(1+exp(y_o))))+
                                -(exp(y_o)*beta)^2*(-trigamma(exp(y_o)*beta)  + trigamma(beta*(1+exp(y_o))))

  beta_score = function(y_o, z) exp(y_o)*beta*(log(z)-digamma(exp(y_o)*beta) + digamma(beta*(1+exp(y_o))))
  beta_llh = function(y_o, z) sum((exp(y_o)*beta-1)*log(z) + (beta-1)*log(1-z)-
                                    log(beta(beta*exp(y_o), beta)))

  beta_link = function(y) 1/(1+exp(-y)) #logit link
  # return object with all components of model
  return(list("hess" = beta_hess, "score"=beta_score, "llh" = beta_llh, "link" = beta_link))

}

.beta_model_alt = function(likparms){
  phi = ifelse("phi" %in% names(likparms),likparms$phi, .5)
  # store score and hessian functions, update
  beta_hess = function(y_o, z)
    -exp(y_o)*phi*(log(z/(1-z))-digamma(exp(y_o)*phi) + digamma(beta*(1-exp(y_o))))+
    -exp(2*y_o)*phi^2*(-trigamma(exp(y_o)*phi) - trigamma(beta*(1-exp(y_o))))

  beta_score = function(y_o, z) exp(y_o)*phi*(log(z/(1-z))-digamma(exp(y_o)*phi) + digamma(beta*(1-exp(y_o))))

  beta_llh = function(y_o, z) sum((exp(y_o)*phi-1)*log(z) + ((1-exp(y_o))*phi-1)*log(1-z)-
                                    log(beta(phi*exp(y_o), phi*(1-exp(y_o)))))

  beta_link = function(y) 1/(1+exp(-y)) #logit link
  # return object with all components of model
  return(list("hess" = beta_hess, "score"=beta_score, "llh" = beta_llh, "link" = beta_link))
}

.beta_data_reqs = function(z){
  negative = any(z<0)
  large = any(z>1)
  return(negative | large)
}


################# Negative Binomial (NOT FINISHED) #########################

.negbin_model = function(likparms){

  # store score and hessian functions, update
  negbin_hess = 0#function(y_o, z) diag(array(exp(y_o)/(1+exp(y_o))^2))
  negbin_score = 0 # function(y_o, z) z - exp(y_o)/(1+exp(y_o))

  # return object with all components of model
  return(list("hess" = negbin_hess, "score"=negbin_score))

}



#' Wrapper for VL version of vecchia_likelihood
#'
#' @param z  an array of real numbers representing observations
#' @param vecchia.approx a vecchia object as generated by vecchia_specify()
#' @param likelihood_model text describing likelihood model to be used for observations
#' @param covparms covariance parameters as a vector
#' @param likparms likelihood parameters for the likelihood_model, as a list
#' @param covmodel describes the covariance model, "matern" by default
#' @param max.iter maximum iterations to perform
#' @param convg  convergence criteria.  End iterations if the Newton step is this small
#' @param return_all  Return additional posterior covariance terms
#' @param y_init Specify initial guess for posterior mode
#' @param prior_mean  specify the prior latent mean
#' @param vecchia.approx.IW an optional vecchia approximation object, can reduce computation if method is called repeatedly
#'
#'
#' @return (multivariate normal) loglikelihood implied by the Vecchia approximation
#' @examples
#' z=rnorm(10); locs=matrix(1:10,ncol=1); vecchia.approx=vecchia_specify(locs,m=5)
#' vecchia_laplace_likelihood(z,vecchia.approx,"gaussian",covparms=c(1,2,.5))
#' @export
vecchia_laplace_likelihood <- function(z,vecchia.approx,likelihood_model, covparms,
                                      likparms = list("alpha"=2, "sigma"=sqrt(.1)),
                                      covmodel='matern', max.iter=50, convg = 1e-5,
                                      return_all = FALSE, y_init = NA,
                                      prior_mean = rep(0,length(z)),
                                      vecchia.approx.IW = NA) {


  # vecchia.approx.IW can be passed in for parameter estimation to reduce cpu time,
  posterior = calculate_posterior_VL(z,vecchia.approx, likelihood_model, covparms, covmodel, likparms,
                                     max.iter, convg, return_all, y_init, prior_mean)

  if(!posterior$cnvgd){warning("Convergence Failed, returning -Inf"); return(-Inf)}


    
  m = ncol(vecchia.approx$U.prep$revNNarray)-1
  locs = as.matrix(vecchia.approx$locsord[order(vecchia.approx$ord.z),])

  # get pseudodata and nuggets from the latent y discovered by VL
  z_pseudo = posterior$t - prior_mean
  if( length(posterior$D) <  length(z_pseudo) & any(is.na(z_pseudo)) ){
      nuggets_pseudo = rep(NA, length(z_pseudo))
      nuggets_pseudo[ !is.na(z_pseudo) ] = posterior$D
  } else {
      nuggets_pseudo = posterior$D
  }

  # create an approximation to llh using interweaved ordering.
  if(all(is.na(vecchia.approx.IW))){
    vecchia.approx.IW = vecchia.approx
    if(vecchia.approx$cond.yz == "zy"){
      vecchia.approx.IW = vecchia_specify(locs, m)
    }
  }
  pseudo_marginal_loglik_vecchia = vecchia_likelihood(z_pseudo, vecchia.approx.IW,
                                                      covparms,nuggets_pseudo, covmodel)


  # get true model log likelihood
  ind.obs = which(!is.na(z))
  true_llh = posterior$model_llh(posterior$mean[ind.obs], z[ind.obs])

  # get gaussian (pseudo-data) approximate log likelihood
  pseudo_cond_loglik = sum(stats::dnorm(z_pseudo, mean = posterior$mean - prior_mean, sd =sqrt(nuggets_pseudo), log = TRUE), na.rm=TRUE)

    # combine three log likelihood terms
    loglik_vecchia = pseudo_marginal_loglik_vecchia - pseudo_cond_loglik
    loglik_vecchia = loglik_vecchia + true_llh

    if(any(is.na(y_init))){
        return(loglik_vecchia)
    } else {
        return(list("llv"=loglik_vecchia, "mean" = posterior$mean))
    }
}








#' Wrapper for VL version of vecchia_likelihood
#'
#' @param z  an array of real numbers representing observations
#' @param posterior posterior distribution obtained from calculate_posterior_VL()
#' @param vecchia.approx a vecchia object as generated by vecchia_specify()
#' @param likelihood_model text describing likelihood model to be used for observations
#' @param covparms covariance parameters as a vector
#' @param likparms likelihood parameters for the likelihood_model, as a list
#' @param covmodel describes the covariance model, "matern" by default
#' @param max.iter maximum iterations to perform
#' @param convg  convergence criteria.  End iterations if the Newton step is this small
#' @param return_all  Return additional posterior covariance terms
#' @param y_init Specify initial guess for posterior mode
#' @param prior_mean  specify the prior latent mean
#' @param vecchia.approx.IW an optional vecchia approximation object, can reduce computation if method is called repeatedly
#'
#'
#' @return (multivariate normal) loglikelihood implied by the Vecchia approximation
#' @export
vecchia_laplace_likelihood_from_posterior <- function(z, posterior, vecchia.approx,likelihood_model, covparms,
                                                      likparms = list("alpha"=2, "sigma"=sqrt(.1)),
                                                      covmodel='matern', max.iter=50, convg = 1e-5,
                                                      return_all = FALSE, y_init = NA,
                                                      prior_mean = rep(0,length(z)),
                                                      vecchia.approx.IW = NA) {

    
  m = ncol(vecchia.approx$U.prep$revNNarray)-1
  locs = as.matrix(vecchia.approx$locsord[order(vecchia.approx$ord.z),])

  # get pseudodata and nuggets from the latent y discovered by VL
  z_pseudo = posterior$t - prior_mean
  if( length(posterior$D) <  length(z_pseudo) & any(is.na(z_pseudo)) ){
      nuggets_pseudo = rep(NA, length(z_pseudo))
      nuggets_pseudo[ !is.na(z_pseudo) ] = posterior$D
  } else {
      nuggets_pseudo = posterior$D
  }

  # create an approximation to llh using interweaved ordering.
  if(all(is.na(vecchia.approx.IW))){
    vecchia.approx.IW = vecchia.approx
    if(vecchia.approx$cond.yz == "zy"){
      vecchia.approx.IW = vecchia_specify(locs, m)
    }
  }
  pseudo_marginal_loglik_vecchia = vecchia_likelihood(z_pseudo, vecchia.approx.IW,
                                                      covparms,nuggets_pseudo, covmodel)


  # get true model log likelihood
  ind.obs = which(!is.na(z))
  true_llh = posterior$model_llh(posterior$mean[ind.obs], z[ind.obs])

  # get gaussian (pseudo-data) approximate log likelihood
  pseudo_cond_loglik = sum(stats::dnorm(z_pseudo, mean = posterior$mean - prior_mean, sd =sqrt(nuggets_pseudo), log = TRUE), na.rm=TRUE)

    # combine three log likelihood terms
    loglik_vecchia = pseudo_marginal_loglik_vecchia - pseudo_cond_loglik
    loglik_vecchia = loglik_vecchia + true_llh

    if(any(is.na(y_init))){
        return(loglik_vecchia)
    } else {
        return(list("llv"=loglik_vecchia, "mean" = posterior$mean))
    }
}










######  wrapper for VL version w/pseudo-data   #######

#' Wrapper for VL version of vecchia_prediction
#'
#' @param vl_posterior  a posterior estimate object produced by calculate_posterior_VL
#' @param vecchia.approx a vecchia object as generated by vecchia_specify()
#' @param covparms covariance parameters as a vector
#' @param pred.mean  provides the prior latent mean for the prediction locations
#' @param var.exact should prediction variances be computed exactly, or is a (faster) approximation acceptable
#' @param covmodel covariance model, 'matern' by default.
#' @param return.values either 'mean' only, 'meanvar', 'meanmat', or 'all'
#'
#'
#' @return (multivariate normal) loglikelihood implied by the Vecchia approximation
#' @examples
#' z=rnorm(10); locs=matrix(1:10,ncol=1); vecchia.approx=vecchia_specify(locs,m=5)
#' vl_posterior = calculate_posterior_VL(z,vecchia.approx,"gaussian",covparms=c(1,2,.5))
#' locs.pred=matrix(1:10+.5,ncol=1)
#' vecchia.approx.pred = vecchia_specify(locs, m=5, locs.pred=locs.pred )
#' vecchia_laplace_prediction(vl_posterior,vecchia.approx.pred,covparms=c(1,2,.5))
#' @export
vecchia_laplace_prediction=function(vl_posterior, vecchia.approx, covparms, pred.mean = 0, var.exact = FALSE,
                                    covmodel='matern',return.values='all') {
  # perform vecchia_prediction with pseudo-data
  z_pseudo = vl_posterior$t -vl_posterior$prior_mean
  nuggets_pseudo = vl_posterior$D


  preds=vecchia_prediction(z_pseudo, vecchia.approx, covparms, nuggets_pseudo,
                           var.exact, covmodel, return.values)
  preds$mu.pred = preds$mu.pred + pred.mean
  preds$mu.obs = preds$mu.obs+vl_posterior$prior_mean
  data_preds = list()

  # Convert predicted mean (median) to data scale
  data_preds$data.pred <- vl_posterior$data_link(preds$mu.pred)
  data_preds$data.obs <- vl_posterior$data_link(preds$mu.obs)

  # Convert predicted variance (via quantiles) to data scale
  data_preds$data_pred_upper_quantile = vl_posterior$data_link(stats::qnorm(p=.95, mean = preds$mu.pred,
                                                                     sd = sqrt(preds$var.pred)))
  data_preds$data_pred_lower_quantile = vl_posterior$data_link(stats::qnorm(p=.05, mean = preds$mu.pred,
                                                                     sd = sqrt(preds$var.pred)))
  data_preds$data_obs_upper_quantiles = vl_posterior$data_link(stats::qnorm(p=.95, mean = preds$mu.obs,
                                                                     sd = sqrt(preds$var.obs)))
  data_preds$data_obs_lower_quantiles = vl_posterior$data_link(stats::qnorm(p=.05, mean = preds$mu.obs,
                                                                     sd = sqrt(preds$var.obs)))

  return(c(preds, data_preds))
}



#### May be useful:  backtracking version of NR
### Does not help with cycles.  Ex usage:
##    tau = backtrack(model_funs, y_o, y_prev, z)
##    y_o = y_prev - tau*(y_prev-y_o)

.backtrack = function(model_funs, y_o, y_prev,z, beta = .9, alpha = .4){
  score = model_funs$score
  llh = model_funs$llh
  v = y_o - y_prev
  tau = 1
  for(i in 1:50){
    #cat(llh(y_prev + tau*v,z), llh(y_prev,z) + alpha*tau*t(score(y_prev,z))%*%v)
    if (llh(y_prev + tau*v,z) > llh(y_prev,z) + alpha*tau*t(score(y_prev,z))%*%v){
      return(tau)
    }else tau = tau*beta
  }
  return(tau)
}

