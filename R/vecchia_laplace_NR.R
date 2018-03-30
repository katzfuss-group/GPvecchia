### Vecchia-Laplace approximation using efficient CPP packages ###

#  stopping criteria for NR interations
stopping_criteria_met = function(y_o, y_prev, cnvg){
  if (max(abs(y_o-y_prev))<cnvg) return(TRUE)
  return(FALSE)
}


#####################################################################
######################    Vecchia + Laplace     #####################
#####################################################################

# algorithm to find latent GP for non-gaussian likelihood
calculate_posterior_VL = function(likelihood_model, covparms, m = 3, use_low_rank = FALSE, convg = 1e-6, return_all = FALSE){
  # pull out constants for readability
  z = likelihood_model$z
  locs = likelihood_model$locs
  # pull out score and second derivative for readability
  ell_dbl_prime = likelihood_model$hess
  ell_prime = likelihood_model$score

  # for logging purposes, output scenario
  output_algo =  ifelse(use_low_rank, "Low Rank",  "VL")
  l_type = likelihood_model$type
  log_comment = paste("Running",output_algo, "(cpp) for",l_type, "with",m,"nbrs and sample size", length(z) )
  message(log_comment)


  # make initial approximation and Cond, NNArray
  cond_type = ifelse(use_low_rank, 'firstm', 'NN')
  vecchia.approx=vecchia_specify(z, locs, m, conditioning = cond_type)#, cond.yz = "z" )

  # init latent variable
  y_o = rep(1, length(z))
  convgd = FALSE
  for( i in 1:50){
    D_inv = ell_dbl_prime(y_o, z)
    D = solve(D_inv)
    u = ell_prime(y_o,z)
    #Calculate pseudo obs
    pseudo.data = D %*% u + y_o
    pseudo.vars = diag(D) # nuggets
    # Update the pseudo data stored in the approximation
    vecchia.approx$zord=pseudo.data[vecchia.approx$ord]
    # update the approx and make the prediction
    updated=vecchia_prediction(vecchia.approx,covparms,pseudo.vars)
    # save y_prev for convergence test
    y_prev = y_o
    y_o = updated$mu.obs
    if (stopping_criteria_met(y_o,y_prev,convg)){
      convgd = TRUE
      tot_iters = i
      break
    }
  }

  if(return_all){
    # return additional information if needed
    orig.order=order(vecchia.approx$ord)
    vec_likelihood = vecchia_likelihood(vecchia.approx,covparms,pseudo.vars)
    W = as.matrix(rev.mat(updated$V.ord%*%t(updated$V.ord))[orig.order,orig.order])
    return (list("mean" = y_o, "sd" =sqrt(updated$var.obs), "iter"=tot_iters,
                 "cnvgd" = convgd, "D" = D, "t"=pseudo.data, "V"=updated$V.ord, "W" = W, "vec_lh"=vec_likelihood))
  }

  return (list("mean" = y_o, "sd" =sqrt(updated$var.obs), "cnvgd" = convgd, "iter" = tot_iters))

}


#####################################################################
######################    Laplace + True Covar ######################
#####################################################################


calculate_posterior_laplace = function(likelihood_model, C, convg = 1e-6, return_all = FALSE){

  l_type = likelihood_model$type
  locsord = likelihood_model$locs
  z = likelihood_model$z
  ell_dbl_prime = likelihood_model$hess
  ell_prime = likelihood_model$score

  log_comment = paste("Running Laplace for",l_type, "with sample size", length(z) )
  message(log_comment)

  C_inv = solve(C)

  y_o = rep(1,length(z))
  tot_iters=0
  # begin NR iteration
  for( i in 1:50){
    #b= exp(y_o) #scale, g''(y)
    D_inv = ell_dbl_prime(y_o, z)
    D = solve(D_inv)  # d is diagonal (hessian), cheap to invert
    u =  ell_prime(y_o,z)
    t = D%*%u+y_o
    W=D_inv+C_inv
    y_prev=y_o
    y_o = solve(W , D_inv) %*% t
    if (stopping_criteria_met(y_o,y_prev,convg)){
      tot_iters = i
      break
    }
  } # end iterate

  if(return_all){
    # Caclulating sd is expensive, avoid for fair comparison
    sd_posterior = sqrt(diag(solve(W)))
    return (list("mean" = y_o, "W"=W,"sd" = sd_posterior, "iter"=tot_iters, "C"=C, "t" = t, "D" = D))
  }
  return (list("mean" = y_o, "W"=W, "iter"=tot_iters))
}


