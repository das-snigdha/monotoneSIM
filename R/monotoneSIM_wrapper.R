#' Bayesian Estimation of Monotone Single Index Models
#'
#' @param y
#' @param X
#' @param beta.init
#' @param xi.init
#' @param Sigma.xi
#' @param knots
#' @param monotone
#' @param iter.HMC
#' @param sigma.sq.beta
#' @param sigma.sq.eps
#' @param a.eps
#' @param b.eps
#' @param Burn.in
#' @param M
#' @param grid.x
#' @param size.grid.x
#'
#' @return
#' @export
#'
#' @examples
monotoneSIM = function(y, X, beta.init, xi.init, Sigma.xi, knots = NULL, monotone = TRUE, iter.HMC = 10,
                       sigma.sq.beta = 1, sigma.sq.eps = 1, a.eps = 1.0, b.eps = 1.0,
                       Burn.in = 100, M = 1000){
  #Compatibility Checks

  # Standardize x and beta.init
  # scale beta.init so that starting value of beta has unit euclidean norm
  beta.std = beta.init/norm(beta.init, "2")

  # Scale the covariate matrix X so that every row of X has euclidean norm <= 1.
  X.scaled = scale(X, center = FALSE, scale = TRUE)
  weights = max(sqrt(rowSums(X.scaled*X.scaled)))
  X.std = X.scaled/weights

  # Weights to back scale beta corresponding to the original covariates.
  weights.beta = attr(X.scaled, "scaled:scale")*weights

  # If knots are not provided, take equispaced knots between -1 and 1 of length as that of xi.init.
  if(is.null(knots)){
    knots = seq(-1, 1, length.out = length(xi.init))
  }
  else{
    knots = as.vector(knots)
    #If provided, throw an error if length does not match the length of xi.init.
    if(length(knots) != length(xi)){
      stop("Number of knots should be same as that of basis coefficients.")
    }

    #sort them in ascending order and throw an error if min(knots) != -1 and max(knots) != 1.
    knots = sort(knots)
    if(min(knots) != -1 || max(knots) != 1)
      stop("Knots should be lie in [-1, 1].")
  }

  # Call C++ function monotoneSIM_c to perform an MCMC algorithm to get samples from the posterior distribution of xi, beta and sigma_sq
  out = monotoneSIM_c(y, X.std, beta.std, xi.init, Sigma.xi, knots, monotone, iter.HMC, sigma.sq.beta,
                      sigma.sq.eps, a.eps, b.eps, Burn.in, M)

  # Back scale beta corresponding to the original covariates
  beta.backscaled = t( t(out$beta) / weights.beta )

  #Adjust the back-scaled betas to have unit euclidean norm
  norms.beta = sqrt(rowSums(beta.backscaled * beta.backscaled))
  beta.backscaled = beta.backscaled / norms.beta

  # return xi, beta.backscaled, beta, sigma.sq.eps, X.std and knots
  return(list("xi" = out$xi, "beta" = beta.backscaled, "beta.Xscaled" = out$beta, "sigma.sq.eps" = out$sigma_sq_eps, "X.scaled" = X.std, "knots" = knots))
}
