EllipticalSliceBeta = function(y, X, xi, beta_init, u, sigma_sq_eps, sigma_sq_beta = 10000.0){

  # Perform Compatibility Checks

  # Call C++ update_beta function to implement the algorithm
  newbeta = update_beta(y, X, xi, beta_init, u, sigma_sq_eps, sigma_sq_beta)

  return(newbeta$beta)
}
