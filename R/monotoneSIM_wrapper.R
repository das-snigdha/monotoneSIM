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
monotoneSIM = function(y, X, beta.init, xi.init, Sigma.xi, knots, monotone = TRUE, iter.HMC = 10,
                       sigma.sq.beta = 10, sigma.sq.eps = 1, a.eps = 1.0, b.eps = 1.0,
                       Burn.in = 100, M = 1000, grid.x = NULL, size.grid.x = 100){
  #Compatibility Checks

  # Standardize x and beta.init

  # Call C++ function monotoneSIM_c to perform an MCMC algorithm to get samples from the posterior distribution of xi, beta and sigma_sq
  out = monotoneSIM_c(y, X, beta.init, xi.init, Sigma.xi, knots, monotone, iter.HMC, sigma.sq.beta,
                      sigma.sq.eps, a.eps, b.eps, Burn.in, M)

  # If grid values of x are not provided, take a grid of size = size.grid.x with equispaced values from -1 and 1.
  if (is.null(grid.x)){
    grid.x = seq(-1, 1, length.out = size.grid.x)
  }

  # Call c++ function g_mtx to calculate g(x) for the grid.x using xi's obtained from the MCMC.
  g.sample = g_mtx(out$xi, grid.x, knots)

  # return xi, beta, sigma.sq.eps, estimated g(x) and grid.x
  return(list("xi" = out$xi, "beta" = out$beta, "sigma.sq.eps" = out$sigma_sq_eps, "g.estimated" = g.sample, "grid.x" = grid.x))
}
