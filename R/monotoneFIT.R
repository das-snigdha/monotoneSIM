#' Extract Fitted values after Bayesian estimation of a Monotone Single Index Model.
#'
#' @description
#' The function \code{monotoneFIT} calculates fitted values of the response using posterior mean of the single index parameter after a Markov Chain Monte Carlo sample has been drawn using the function \code{\link{monotoneSIM}}. It also returns the estimated link function of the Single Index Model based on basis coefficients for a grid of supplied values.
#'
#' @param mono.sim A returned object of \code{monotoneSIM} function.
#' @param size.grid.x length of vector of values for which the link function is estimated. Takes the value \eqn{100} by default.
#' @param grid.x Vector of user supplied values for which the link function is estimated. Takes \code{NULL} by default. If \code{NULL}, vector (of length \code{size.grid.x}) of equispaced values between -1 and 1 are chosen.
#'
#' @return A list with the following elements.
#' \item{Y.fitted}{ \code{M} \eqn{x n} Matrix of Fitted values of the response. \code{M} is the size of the MCMC sample obtained in \code{monotoneSIM}. Each row represents estimated value of the response corresponding to each posterior sample of basis coefficients.}
#' \item{link.estimated}{ \code{M} \eqn{x} \code{length(grid.x)} Matrix of estimated values of the link function. Each row represents estimated value of the link function evaluated on grid.x, corresponding to each posterior sample of basis coefficients.}
#' \item{grid.x}{ Grid of x values at which the link function is estimated. }
#' @export
#'
#' @examples
#' n = 100; p = 3; L = 20
#' # We take 2 continuous variables and 1 dichotomous attribute as predictors.
#' X = matrix(rnorm(n*(p-1)), nrow = n, ncol = (p-1))
#' X = cbind(X, rbinom(n, 1, 0.5))
#'
#' # True Value of the parameter (having unit euclidean norm).
#' true.beta = rnorm(p); true.beta = true.beta/ norm(true.beta, "2")
#'
#' beta.start = rnorm(p)  #Starting value of beta
#' xi = abs(rnorm((L+1), 0, 5))   #Starting value of xi
#' S_xi = 5*diag(L+1) #Prior Variance of xi
#' sigma.sq.eps.start = 0.01 #Starting value of sigma.sq.eps
#'
#' # True monotone link function
#' true.g = function(x){
#'   y = (x+1)/2
#'   5*(pnorm(y, mean=0.5, sd=0.1) - pnorm(0, mean = 0.5, sd = 0.1))
#'   }
#'
#' # Generate the response
#' y.true = true.g(X%*%true.beta) + rnorm(n, 0, sqrt(sigma.sq.eps.start))
#'
#' MCMC.sample = monotoneSIM(y = y.true, X = X, beta.init = beta.start,
#'   xi.init = xi, Sigma.xi =  S_xi, monotone = TRUE,
#'   sigma.sq.eps.init = sigma.sq.eps.start , Burn.in = 100, M = 500)
#'
#' fit = monotoneFIT(MCMC.sample)
#'
#' # Obtain the fitted response
#' y.fit = colMeans(fit$Y.fitted)
#'
#' # Plot the fitted and true values of the response
#' Y = data.frame(y.true, y.fit); Y = Y[order(y.fit), ]
#' plot(Y$y.true, pch = 16, ylab = "Y (Response)", xlab = "",
#'   main = "Plot of true responses and the fitted model.")
#' lines(Y$y.fit, col = "red", lwd = 2)
#' legend("topleft", c("True response", "Fitted model"), col = c("black", "red"),
#'   lwd = 2, lty = c(0,1), pch = c(16, NA))
#'
#' # Obtain the estimated value of the link function g(x)
#' est.func = colMeans(fit$link.estimated)
#'
#' # Calculate the true value of g(x) for each x in grid.x
#' true.func = rep(0,length(est.func))
#' for(i in 1:length(est.func))
#'   true.func[i] = true.g(fit$grid.x[i])
#'
#' # Plot the estimated and the true link function vs grid.x
#' plot(fit$grid.x, est.func, type="l", xlab="grid.x", ylab="Link function, g(x)",
#'   main = "Plot of the link function, g(x) vs x")
#' lines(fit$grid.x, true.func, lwd=2, col="red")
#' legend("topleft", c("Estimated link function", "True link function"),
#'   col = c("black", "red"), lwd = c(1,3))

monotoneFIT = function(mono.sim, size.grid.x = 100, grid.x = NULL){

  #Check if size.grid.x is a positive integer
  if(size.grid.x <= 0 || size.grid.x != as.integer(size.grid.x)){
    stop("Number of values at which the link function is estimated should be a positive integer.")
  }

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
  return(list("Y.fitted" = y.fitted, "link.estimated" = g.sample, "grid.x" = grid.x))

}
