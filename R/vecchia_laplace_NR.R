### Vecchia-Laplace approximation using efficient CPP packages ###


#####################################################################
######################    Vecchia + Laplace     #####################
#####################################################################




# algorithm to find latent GP for non-gaussian likelihood
calculate_posterior_VL = function(vecchia.approx, likelihood_model=c("gaussian","logistic", "poisson", "gamma"),
                                  covparms, likparms = list("alpha"=2, "sigma"=.1),
                                  max.iter=50, convg = 1e-5, return_all = FALSE){

  zy.conditioning = (vecchia.approx$cond.yz=="zy")
  likelihood_model <- match.arg(likelihood_model)

  # undo ordering - Vecchia code reorders
  orig_ord = order(vecchia.approx$ord)

  z = vecchia.approx$zord[orig_ord]
  z = z[!is.na(z)] # for 'zy' ordering, z and loc duplicated
  #locs = vecchia.approx$locsord[orig_ord,][1:length(z)]


  # pull out score and second derivative for readability
  model_funs = switch(likelihood_model,
                      "gaussian" = .gauss_model(likparms),
                      "logistic" = .logistic_model(),
                      "poisson" = .poisson_model(),
                      "gamma" = .gamma_model(likparms))
  ell_dbl_prime = model_funs$hess
  ell_prime = model_funs$score

  # for logging purposes, output scenario

  log_comment = print(paste("Running VL-NR for",likelihood_model, "with",
                            ncol(vecchia.approx$U.prep$revNNarray)-1,"nbrs and sample size",length(z)))

  # record duration of NR
  t_start = Sys.time()

  # init latent variable
  y_o = rep(1, length(vecchia.approx$zord))
  #points(locs[order(locs)], y_o[order(locs)], type = "l", col = alpha("black", .4))
  convgd = FALSE
  tot_iters = max.iter
  for( i in 1:max.iter){
    y_prev = y_o    # save y_prev for convergence test
    D_inv =ell_dbl_prime(y_o, z)
    D = 1/diag(D_inv)
    u = ell_prime(y_o,z)
    pseudo.data =matrix(D * u + y_o, ncol=1)
    nuggets = D
    # Update the pseudo data stored in the approximation
    vecchia.approx$zord=pseudo.data[vecchia.approx$ord]
    vecchia.approx$zord = vecchia.approx$zord[!is.na(vecchia.approx$zord)]
    # Update U matrix with new nuggets, make the prediction
    U=createU(vecchia.approx,covparms,nuggets)
    V.ord=U2V(U,vecchia.approx)
    vecchia.mean=vecchia_mean(vecchia.approx,U,V.ord)
    if( zy.conditioning){
      y_o = vecchia.mean$mu.pred # y treated as prediction
    }else{
      y_o = vecchia.mean$mu.obs  # SGV, firstm, etc
    }
    #print(sum(abs(V.ord%*%(y_o-y_prev)))) #newton increment-> reveals cycling
    #print(y_o)
    #points(locs[order(locs)], y_o[order(locs)], type = "l", col = alpha("black", .2))
    if (max(abs(y_o-y_prev))<convg){
      convgd = TRUE
      tot_iters = i
      break
    }
  }
  t_end = Sys.time()
  LV_time = as.double(difftime(t_end, t_start, units = "secs"))

  if(return_all){
    # return additional information if needed
    orig.order=order(vecchia.approx$ord)
    vec_likelihood = vecchia_likelihood(vecchia.approx,covparms,nuggets)
    W = as.matrix(rev.mat(V.ord%*%t(V.ord))[orig.order,orig.order])
    if (vecchia.approx$cond.yz=="zy"){
      n = length(y_o)
      W = W[(n+1):(2*n), (n+1):(2*n)]
      vec_likelihood = NA # needs to be corrected?
    }

    return (list("mean" = y_o, "sd" =sqrt(diag(solve(W))), "iter"=tot_iters,
                 "cnvgd" = convgd, "D" = D, "t"=pseudo.data, "V"=V.ord,
                 "W" = W, "vec_lh"=vec_likelihood, "runtime" = LV_time, "U" = U))
  }
  return (list("mean" = y_o, "cnvgd" = convgd, "runtime" = LV_time,
               "iter" = tot_iters, "t"=pseudo.data, "D" = D))
}



#####################################################################
######################    Laplace + True Covar ######################
#####################################################################


calculate_posterior_laplace = function(z, likelihood_model, C,  likparms = list("alpha"=2, "sigma"=.1),
                                       convg = 1e-6, return_all = FALSE){
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
    t = D%*%u+y_o
    y_prev = y_o
    #y_o = solve(W , D_inv) %*% t # prior mean 0 update
    y_o = t - D%*%solve(D+C,t) # woodbury morrison of previous line
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
  sigma = ifelse("sigma" %in% names(likparams),likparams$sigma, .1)
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
  gamma_llh = function(y_o, z) sum(-y_o*z + (alpha-1)*log(z) +alpha*log(y_o)-n*log(gamma(alpha)))
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

