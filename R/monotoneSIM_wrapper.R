#' Bayesian Estimation of Monotone Single Index Models
#'
#' @param y \eqn{n x 1} Vector of Response Variable.
#'
#' @param X \eqn{n x p} Matrix of Predictors. Each row represents a \eqn{p x 1} predictor vector.
#'
#' @param beta.init \eqn{p x 1} Vector; Starting value of beta for the algorithm.
#'
#' @param xi.init \eqn{(L+1) x 1} Vector; Starting value of Basis Coefficients for the algorithm.
#'
#' @param Sigma.xi \eqn{(L+1) X (L+1)} Matrix; Hyperparameter specifying prior Variance of xi.
#'
#' @param knots \eqn{(L+1) X 1} Vector of user supplied knots. Takes \code{NULL} by default. If \code{NULL}, equispaced knots in \eqn{[-1, 1]} are chosen for the algorithm.
#'
#' @param monotone  Logical; Takes \code{TRUE} by default. If \code{TRUE}, the link function is taken to be monotone and the posterior of xi is drawn from a Truncated Normal distribution by an exact Hamiltonian Monte Carlo algorithm. If \code{FALSE},  the posterior of xi is drawn from a Normal distribution.
#'
#' @param iter.HMC  Number of iterations of the Hamiltonian Monte Carlo algorithm in Truncated Normal sampling if \code{monotone = TRUE}. Takes the value \eqn{10} by default.
#'
#' @param sigma.sq.beta Scalar; Hyperparameter specifying prior Variance of beta. Takes the value \eqn{1} by default.
#'
#' @param sigma.sq.eps  Scalar; Starting of Error Variance. Takes the value \eqn{1} by default.
#'
#' @param a.eps Scalar; Hyperparameter specifying prior distirbution of sigma.sq.eps. Takes the value \eqn{1} by default.
#'
#' @param b.eps Scalar; Hyperparameter specifying prior distirbution of sigma.sq.eps. Takes the value \eqn{1} by default.
#'
#' @param Burn.in Non-negative Integer; Burn-in period of the Markov Chain Monte Carlo algorithm. Takes the value \eqn{100} by default.
#'
#' @param M Positive Integer; required size of the  Markov Chain Monte Carlo sample. Takes the value \eqn{1000} by default.
#'
#' @return A list with the following elements.
#' \item{xi}{ \code{M} \eqn{x (L+1)} Matrix of Basis Coefficients. Each row represent one sample from the Conditional posterior of basis coefficients.}
#' \item{beta.Xscaled}{ \code{M} \eqn{x p} Matrix of parameters corresponding to scaled covariates. Each row represent one sample from the Conditional posterior of single index parameter.}
#' \item{beta}{ \code{M} \eqn{x p} Matrix of parameters corresponding to original covariates.}
#' \item{sigma.sq.eps}{ \code{M} \eqn{x 1} Vector; Sample of size \eqn{M} drawn from the Conditional posterior of the Error Variance.}
#' \item{X.scaled}{ \eqn{n x p} Matrix of scaled Covariates. Each row has euclidean norm less than or equal to \eqn{1}.}
#' \item{knots}{  \eqn{(L+1) x 1} Vector of knots used for estimation of the link function.}
#' @export
#'
#' @examples
#' n=100; p=3; L = 20
#'
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
#' y = (x+1)/2
#' 5*(pnorm(y, mean=0.5, sd=0.1) - pnorm(0, mean = 0.5, sd = 0.1))
#' }
#'
#' # Generate the response
#' y.true = true.g(X%*%true.beta) + rnorm(n, 0, sqrt(sigma.sq.eps.start))
#'
#' MCMC.sample = monotoneSIM(y = y.true, X = X, beta.init = beta.start,
#'   xi.init = xi, Sigma.xi =  S_xi, monotone = TRUE,
#'   sigma.sq.eps = sigma.sq.eps.start, Burn.in = 100, M = 500)
#'
#' #Posterior mean of beta
#' beta.estimated = colMeans(MCMC.sample$beta); beta.estimated
#' true.beta  #Compare with true beta
#'
#' #Posterior Standard Deviation of beta
#' sd.beta = apply(MCMC.sample$beta, 2, sd); sd.beta

monotoneSIM = function(y, X, beta.init, xi.init, Sigma.xi, knots = NULL, monotone = TRUE, iter.HMC = 10,
                       sigma.sq.beta = 1, sigma.sq.eps = 1, a.eps = 1.0, b.eps = 1.0,
                       Burn.in = 100, M = 1000){

  # Convert y, beta.init and xi.init into vectors
  y = as.vector(y);
  beta.init = as.vector(beta.init);
  xi.init = as.vector(xi.init)

  # Convert X and Sigma.xi into matrices
  X = as.matrix(X);
  Sigma.xi = as.matrix(Sigma.xi)

  n = nrow(X); p = ncol(X); L = length(xi.init)

  # Compatibility Checks
  # Check if dimensions of y and X match
  if(length(y) != n){
    stop("Number of responses and predictors do not match.")
  }
  # Check if dimensions of beta.init and X match
  if(length(beta.init) != p){
    stop("Dimension of predictors and number of parameters do not match.")
  }
  # Check if dimensions of Sigma_xi and xi.init match
  if(nrow(Sigma.xi) != L || ncol(Sigma.xi) != L){
    stop("Dimensions of basis coefficient vector and its Variance matrix do not match.")
  }
  else if(isSymmetric(Sigma.xi) == FALSE){
    # Check if Sigma_xi is symmetric
    stop("Variance matrix of basis coefficients is not Symmetric.")
  }

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
    if(length(knots) != length(xi.init)){
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
  return(list("xi" = out$xi, "beta.Xscaled" = out$beta, "beta" = beta.backscaled,  "sigma.sq.eps" = out$sigma_sq_eps, "X.scaled" = X.std, "knots" = knots))
}
