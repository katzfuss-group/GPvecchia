### Collection of methods to wrap data with appropriate
### likelihood model for using Vecchia-Laplace algorithm



#################  Logistic   #########################
logistic_model = function(){
  logistic_llh = function(y_o, z) sum(z*y_o-log(1+exp(y_o)))
  logistic_hess = function(y_o, z) diag(array(exp(y_o)/(1+exp(y_o))^2))
  logistic_score = function(y_o, z) z - exp(y_o)/(1+exp(y_o))
  # return object with all components of model
  return(list("hess" = logistic_hess, "score"=logistic_score,"llh" = logistic_llh))
}

#################  Poisson  #########################
pois_model = function(){
  pois_llh = function(y_o, z) sum(z*y_o -exp(y_o)-log(factorial(z)))
  pois_hess =function(y_o, z) diag(array(exp(y_o)))
  pois_score = function(y_o, z) z-exp(y_o)
  return(list("hess" = pois_hess, "score"=pois_score,"llh" = pois_llh))
}

#################  Gaussian  #########################
gauss_model = function(sigma = .3){

  gauss_llh = function(y_o, z) sum(-.5*(z-y_o)^2/sigma^2) -n*(log(sigma)+log(2*pi)/2)
  gauss_hess = function(y_o, z) diag(array(rep(1/sigma^2, length(y_o))))
  gauss_score = function(y_o, z) (z-y_o)/sigma^2
  return(list("hess" = gauss_hess, "score"=gauss_score, "llh"=gauss_llh))
}

#################  Gamma  #########################
gamma_sample = function(alpha = 2 ){
  gamma_hess = function(y_o, z) diag(array(z*exp(y_o)))
  gamma_score = function(y_o, z) -z*exp(y_o)+ alpha
  gamma_llh = function(y_o, z) sum(-y_o*z + (alpha-1)*log(z) +alpha*log(y_o)-n*log(Gamma(alpha)))
  return(list("hess" = gamma_hess, "score"=gamma_score, "llh" = gamma_llh))
}




################# Negative Binomial (NOT FINISHED) #########################

negbin_sample = function(n, covfun, seed = 125, dom = 1){

  # store score and hessian functions, update
  negbin_hess = 0#function(y_o, z) diag(array(exp(y_o)/(1+exp(y_o))^2))
  negbin_score = 0 # function(y_o, z) z - exp(y_o)/(1+exp(y_o))

  # return object with all components of model
  return(list("hess" = logistic_hess, "score"=logistic_score))

}




#### General wrapper

define_likelihood_model = function(model_type = c("gaussian","logistic", "poisson", "gamma"),
                                   locs, obs){

  model_type <- match.arg(model_type)

  if (!is.matrix(locs)) stop("Locations must be a matrix")

  model_funs = switch(model_type,
                      "gaussian" = gauss_model(),
                      "logistic" = logistic_model(),
                      "poisson" = poisson_model(),
                      "gamma" = gamma_model())

  return(list("type"=model_type, "locs" = locs, "z"=obs,
              "hess" = model_funs$hess,
              "score"=model_funs$score,
              "llh" = model_funs$llh))
}


