#' Extract Fitted values after Bayesian estimation of a Monotone Single Index Model.
#'
#' @param mono.sim A returned object of \code{monotoneSIM} function.
#' @param size.grid.x length of vector of values for which the monotone function is estimated. Takes the value \eqn{100} by default.
#' @param grid.x (Optional) Vector of user supplied values for which the monotone function is estimated. Takes \code{NULL} by default. If \code{NULL}, vector (of length size.grid.x) of equispaced values between -1 and 1 are chosen.
#'
#' @return A list with the following elements.
#' \item{Y.fitted}{ \code{M} \eqn{x n} Matrix of Fitted values of the response. \code{M} is the size of the MCMC sample obtained in \code{monotoneSIM}. }
#' \item{g.estimated}{ \code{M} \eqn{x} \code{size.grid.x} Matrix OR (if \code{grid.x} is supplied) \code{M} \eqn{x} \code{length(grid.x)} Matrix of estimated values of the link function. Each row represents estimated value of the link function evaluated on grid.x, corresponding to each posterior sample of basis coefficients.}
#' \item{grid.x}{ Grid of x values at which the link function is estimated. }
#' @export
#'
#' @examples
monotoneFIT = function(mono.sim, size.grid.x = 100, grid.x = NULL){

  # If grid values of x are not provided, take a grid of size = size.grid.x with equispaced values from -1 and 1.
  if (is.null(grid.x)){
    grid.x = seq(-1, 1, length.out = size.grid.x)
  }
  else{
    # If provided, convert grid.x into a vector
    grid.x = as.vector(grid.x)
  }
  xi.mtx = as.matrix(mono.sim$xi)

  # Call c++ function g_mtx to calculate g(x) for the grid.x using xi's obtained from the MCMC.
  g.sample = g_mtx(xi_mtx = xi.mtx, grid_x = grid.x, u = mono.sim$knots)

  # Estimate beta.Xscaled by the posterior sample mean.
  est.beta.Xscaled = colMeans(mono.sim$beta.Xscaled)

  #Use estimated beta and X.scaled to calculate Xbeta
  Xbeta = mono.sim$X.scaled %*% est.beta.Xscaled

  # Call c++ function g_mtx to calculate g(Xbeta) i.e. estimate Y using xi's obtained from the MCMC.
  y.fitted = g_mtx(xi_mtx = xi.mtx, grid_x = Xbeta, u = mono.sim$knots)

  # return y.fitted, g.sample and grid.x
  return(list("Y.fitted" = y.fitted, "g.estimated" = g.sample, "grid.x" = grid.x))

}
