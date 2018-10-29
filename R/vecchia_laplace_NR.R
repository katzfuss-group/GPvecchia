### Vecchia-Laplace approximation using efficient CPP packages ###


#####################################################################
######################    Vecchia + Laplace     #####################
#####################################################################




# algorithm to find latent GP for non-gaussian likelihood
calculate_posterior_VL = function(z,vecchia.approx,
          likelihood_model=c("gaussian","logistic", "poisson", "gamma"),
          covparms, likparms = list("alpha"=2, "sigma"=sqrt(.1)),
          max.iter=50, convg = 1e-5, return_all = FALSE, y_init = NA,
          prior_mean = rep(0,length(z))){

  zy.conditioning = (vecchia.approx$cond.yz=="zy")
  likelihood_model <- match.arg(likelihood_model)

  # pull out score and second derivative for readability
  model_funs = switch(likelihood_model,
                      "gaussian" = .gauss_model(likparms),
                      "logistic" = .logistic_model(),
                      "poisson" = .poisson_model(),
                      "gamma" = .gamma_model(likparms))
  ell_dbl_prime = model_funs$hess
  ell_prime = model_funs$score

  # for logging purposes, output scenario

  log_comment = print(paste("Running VL-NR for",likelihood_model, "with m=",
                ncol(vecchia.approx$U.prep$revNNarray)-1," and sample size",length(z)))

  # record duration of NR
  t_start = Sys.time()

  # init latent variable
  y_o = y_init
  if(is.na(y_o)) y_o = prior_mean

  #points(locs[order(locs)], y_o[order(locs)], type = "l", col = alpha("black", .4))
  convgd = FALSE
  tot_iters = 1
  for( i in 1:max.iter){

    y_prev = y_o    # save y_prev for convergence test

    # update pseudo-data and -variances
    D_inv = ell_dbl_prime(y_o, z)
    D = 1/diag(D_inv)
    u = ell_prime(y_o,z)
    pseudo.data = D * u + y_o - prior_mean
    nuggets = D

    # make prediction
    preds=vecchia_prediction(pseudo.data,vecchia.approx,covparms,
                                    nuggets,return.values='meanmat')
    if( zy.conditioning){
      y_o = preds$mu.pred + prior_mean # y treated as prediction
    }else{
      y_o = preds$mu.obs + prior_mean # SGV, firstm, etc
    }
    #print(sum(abs(V.ord%*%(y_o-y_prev)))) #newton increment-> reveals cycling
    #print(y_o)
    #points(locs[order(locs)], y_o[order(locs)], type = "l", col = alpha("black", .2))
    if(is.na(max(abs(y_o-y_prev)))){
      # convergence failed due to machine precision?
      fail_comment = print(paste("VL-NR hit NA on iteration ",tot_iters,", convergence failed."))

      y_o = y_prev
      tot_iters = -1
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

  if(return_all){
    # return additional information if needed, can be slow
    orig.order=order(vecchia.approx$ord)
    V.ord=preds$V.ord
    W = as(rev.mat(V.ord%*%t(V.ord))[orig.order,orig.order], 'dgCMatrix')
    if (vecchia.approx$cond.yz=="zy"){
      n = length(y_o)
      V.ord = V.ord[1:n, 1:n]
      W = W[(n+1):(2*n), (n+1):(2*n)]
      vec_likelihood = NA # needs to be adjusted: run with SGV
    } else {
      vec_likelihood = vecchia_likelihood(pseudo.data,vecchia.approx,covparms,nuggets)
    }

    return (list("mean" = y_o, "sd" =sqrt(diag(solve(W))), "iter"=tot_iters,
                 "cnvgd" = convgd, "D" = D, "t"=pseudo.data, "V"=V.ord,
                 "W" = W, "vec_lh"=vec_likelihood, "runtime" = LV_time))
  }
  return (list("mean" = y_o, "cnvgd" = convgd, "runtime" = LV_time,
               "iter" = tot_iters, "t"=pseudo.data, "D" = D))
}



#####################################################################
######################    Laplace + True Covar ######################
#####################################################################


calculate_posterior_laplace = function(z, likelihood_model, C,  likparms = list("alpha"=2, "sigma"=sqrt(.1)),
                                       convg = 1e-6, return_all = FALSE, prior_mean = 0){
  # pull out score and second derivative for readability
  model_funs = switch(likelihood_model,
                      "gaussian" = .gauss_model(likparms),
                      "logistic" = .logistic_model(),
                      "poisson" = .poisson_model(),
                      "gamma" = .gamma_model(likparms))
  ell_dbl_prime = model_funs$hess
  ell_prime = model_funs$score

  log_comment = paste("Running Laplace for",likelihood_model, "with sample size", length(z) )
  print(log_comment)

  t_start = Sys.time()
  y_o = rep(1, length(z))
  tot_iters=0
  # begin NR iteration
  for( i in 1:50){
    #b= exp(y_o) #scale, g''(y)
    D_inv = ell_dbl_prime(y_o, z)
    D = sparseMatrix(i=1:length(y_o), j = 1:length(y_o), x= 1/diag(D_inv))
    u =  ell_prime(y_o,z)
    t = D%*%u+y_o - prior_mean
    y_prev = y_o
    #y_o = solve(W , D_inv) %*% t # prior mean 0 update
    y_o = t - D%*%solve(D+C,t) + prior_mean # woodbury morrison of previous line
    #W = D_inv +  C_inv
    #print((t(y_o-y_prev)%*%W%*%(y_o-y_prev))[1,1]) #newton increment, llh
    tot_iters = i
    if (max(abs(y_o-y_prev))<convg) break
  } # end iterate
  t_end = Sys.time()
  Lap_time = as.double(difftime(t_end, t_start, units = "secs"))
  if(return_all){
    # Caclulating sd is expensive
    W = D_inv +  solve(C)
    sd_posterior = sqrt(diag(solve(W)))
    return (list("mean" = y_o, "W"=W,"sd" = sd_posterior, "iter"=tot_iters,
                 "C"=C, "t" = t, "D" = D, "runtime"=Lap_time))
  }
  return (list("mean" = y_o, "iter"=tot_iters))
}



#################  Logistic   #########################
.logistic_model = function(){
  logistic_llh = function(y_o, z) sum(z*y_o-log(1+exp(y_o)))
  logistic_hess = function(y_o, z) sparseMatrix(i=1:length(y_o), j = 1:length(y_o), x=exp(y_o)/(1+exp(y_o))^2)
  logistic_score = function(y_o, z) z - exp(y_o)/(1+exp(y_o))
  # return object with all components of model
  return(list("hess" = logistic_hess, "score"=logistic_score,"llh" = logistic_llh))
}

#################  Poisson  #########################
.poisson_model = function(){
  pois_llh = function(y_o, z) sum(z*y_o -exp(y_o)-log(factorial(z)))
  pois_hess =function(y_o, z)  sparseMatrix(i=1:length(y_o), j = 1:length(y_o), x=exp(y_o))
  pois_score = function(y_o, z) z-exp(y_o)
  return(list("hess" = pois_hess, "score"=pois_score,"llh" = pois_llh))
}

#################  Gaussian  #########################
.gauss_model = function(likparams){
  # default nugget sd = .3
  sigma = ifelse("sigma" %in% names(likparams),likparams$sigma, sqrt(.1))
  gauss_llh = function(y_o, z) sum(-.5*(z-y_o)^2/sigma^2) -n*(log(sigma)+log(2*pi)/2)
  gauss_hess = function(y_o, z)  sparseMatrix(i=1:length(y_o), j = 1:length(y_o), x=rep(1/sigma^2, length(y_o)))
  gauss_score = function(y_o, z) (z-y_o)/sigma^2
  return(list("hess" = gauss_hess, "score"=gauss_score, "llh"=gauss_llh))
}

#################  Gamma  #########################
.gamma_model = function(likparams){
  alpha = ifelse("alpha" %in% names(likparams),likparams$alpha, 2)
  gamma_hess = function(y_o, z)  sparseMatrix(i=1:length(y_o), j = 1:length(y_o), x=z*exp(y_o))
  gamma_score = function(y_o, z) -z*exp(y_o)+ alpha
  #gamma_llh = function(y_o, z) sum(-y_o*z + (alpha-1)*log(z) +alpha*log(y_o)-n*log(gamma(alpha))) # canonical link
  gamma_llh = function(y_o, z) sum(-exp(y_o)*z + (alpha-1)*log(z) +alpha*y_o-n*log(gamma(alpha))) # log link
  return(list("hess" = gamma_hess, "score"=gamma_score, "llh" = gamma_llh))
}




################# Negative Binomial (NOT FINISHED) #########################

.negbin_model = function(likparams){

  # store score and hessian functions, update
  negbin_hess = 0#function(y_o, z) diag(array(exp(y_o)/(1+exp(y_o))^2))
  negbin_score = 0 # function(y_o, z) z - exp(y_o)/(1+exp(y_o))

  # return object with all components of model
  return(list("hess" = logistic_hess, "score"=logistic_score))

}


#### May be useful:  backtracking version of NR
### Does not help with cycles.  Ex usage:
##    tau = backtrack(model_funs, y_o, y_prev, z)
##    y_o = y_prev - tau*(y_prev-y_o)

backtrack = function(model_funs, y_o, y_prev,z, beta = .9, alpha = .4){
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



# Bringing loglikhd calculation into this class requires recalculating vecchia_specify in the class

VL_loglik = function(data_gen, pred_lv, ms = 1, ps = 0){
  pseudo_marginal_loglik_vecchia=pred_lv$vec_lh
  # for testing approximation, allow perturbation of mean
  n=length(pred_lv$mean)
  perturbed = rnorm(n, sd = ps)
  post_mean = pred_lv$mean*ms + perturbed
  pseudo_cond_loglik = dmvnorm(x=t(pred_lv$t), mean = post_mean, sigma = diag(pred_lv$D), log = TRUE)
  true_loglike = data_gen$llh(post_mean, data_gen$z)
  loglik_vecchia = pseudo_marginal_loglik_vecchia + true_loglike - pseudo_cond_loglik
  return(loglik_vecchia)
}
