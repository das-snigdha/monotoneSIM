#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "CppBasics.h"

using namespace Rcpp;

// Triangular basis function h(x) = (1-|x|)_+
// x : value at which h(.) is evaluated
double h(double x){
  if(x > -1 && x < 1)
    return(1.0 - abs(x));   // return (1-|x|) if -1 <= x <= 1
  else
    return(0.0);    // return 0 otherwise
}

// h_l(.) = B-spline basis function of order 2
// x : value at which h_l(.) is evaluated
// l : index ranging from 0 to (#knots - 1)
// u : vector of specified knots from -1 to 1
double h_l(double x, int l, const NumericVector& u){
  if(x < u[l]){
    // if x < u_l  then h_l((x - u_l)/(u_l - u_{l-1}))

    double u1 ;
    if(l==0)
      u1 = 2.0*u[0] - u[1] ;
    else
      u1 = u[l-1] ;
    return( h((x - u[l])/(u[l]-u1)) ) ;
  }
  else{
    // if x >= u_l  then h_l((x - u_l)/(u_{l+1} - u_l))

    double u3 ;
    int L=u.size()-1 ;
    if(l==L)
      u3 = 2.0*u[L] - u[L-1] ;
    else
      u3 = u[l+1] ;
    return( h((x - u[l])/(u3-u[l])) ) ;
  }
}

// \psi(x) = \int_{-1}^x h(t) dt
// x : value at which \psi(.) is evaluated
double psi(double x){
  if( x < -1.0 )
    return(0.0) ;   // return 0 for x < -1
  else if( x < 0.0 )
    return( 0.5*(x+1.0)*(x+1.0) ) ;   // return 1/2*(x+1)^2 if -1<x<0
  else if(x < 1.0)
    return( 0.5 + 0.5*(2.0-x)*x ) ;   // return 1/2 * 1/2*(2-x)x if 0<=x<1
  else
    return(1.0) ;   // return 1 if x >= 1
}

// psi_l(x) = \int_{-1}^x h_l(t) dt
// l : index ranging from 0 to (#knots - 1)
// u : specified knots from -1 to 1
double psi_l(double x, int l, const NumericVector& u){
  if(l==0){
    // specific algerba for psi_0
    if(x <= -1)
      return (0.0);
    else if(x < u[1]){
      return( ( u[1]*(1+x) + 0.5 - 0.5*x*x )/(u[1]+1) ) ;
    }
    else
      return(0.5*(u[1]+1)) ;
  }
  else{
    double u1, tmp ;
    u1 = u[l-1] ;

    if(x < u[l]){
      // if x < u_l, return (u_l - u_{l-1})* \psi( (x - u_l)/(u_l - u_{l-1}) )

      tmp = u[l] - u1 ;
      return( tmp*psi((x-u[l])/tmp) ) ;
    }
    else{
      // if x >= u_l, return 0.5(u_l - u_{l-1})* (u_{l-1} - u_l){\psi( (x - u_l)/(u_{l-1} - u_l) ) - 0.5}

      double u3 ;
      int L=u.size()-1 ;
      if(l==L)
        u3 = 2.0*u[L] - u[L-1] ;
      else
        u3 = u[l+1] ;
      tmp = u3 - u[l] ;
      return( 0.5*(u[l]-u1) + tmp*(psi((x-u[l])/tmp)-0.5) ) ;
    }
  }
}

// Monotone single index function as a B-Spline basis expansion.
// g(x) = \sum_{l=0}^L \xi_l * \psi_l(x)
// x : value at which g(.) is evaluated
// xi : vector of basis coefficients
// u : vector of specified knots from -1 to 1
double g(double x, const NumericVector& xi, const NumericVector& u){
  if(x > 1.0)
    return(g(1.0, xi, u)) ;   // If x>1, g(x) = g(1)
  else if(x < -1.0)
    return(0.0) ;     // If x<-1, g(x) = g(0)
  else{
    double res=0.0 ;
    for(int l=0; l<xi.size();l++)
      res += xi[l]*psi_l(x,l,u) ;     // For -1<=x<=1, g(x) = \sum_{l=0}^L \xi_l * \psi_l(x)
    return(res) ;
  }
}

// log likelihood function of beta
// log_L =  -(1/2*\sigma_sq_eps) \sum_{i=1}^n (y_i - g(X_i^{\top}\beta))^2
// y : nx1 response variable
// X : nxp matrix of covariates
// Xbeta : nx1 vector containing X*beta, X: nxp matrix of covariates, beta: px1 predictors
// xi : (L+1)x1 vector of coefficients to contruct the single index function
// u : vector of knots ranging from -1 to 1
// sigma_sq_eps : error variance of the model
double log_L(const NumericVector& y, const NumericVector& Xbeta,
             const NumericVector& xi, double sigma_sq_eps, const NumericVector& u){
  double res=0.0, tmp ;
  for(int i=0; i<y.size(); i++){
    // res = \sum_{i=1}^n (y_i - g(X_i^{\top}\beta))^2

    tmp = y[i] - g(Xbeta[i], xi, u) ;
    res += tmp*tmp ;
  }
  res = -0.5*res/sigma_sq_eps ;     // multiplying res by -(1/2*\sigma_sq_eps) to get log_L as defined
  return(res) ;
}

// Function to Draw a sample from the Posterior Distribution of beta using Elliptical Slice Sampler
// y : nx1 response variable
// X : nxp matrix of covariates
// xi : (L+1)x1 vector of coefficients to contruct the single index function
// beta_init : px1 starting value of beta, having unit Euclidean norm
// u : vector of knots ranging from -1 to 1
// sigma_sq_eps : error variance of the model
// sigma_sq_beta : Prior variance of beta
List update_beta(const NumericVector& y, const NumericMatrix& X, const NumericVector& xi,
                 const NumericVector&  beta_init, const NumericVector& u, double sigma_sq_eps,
                 double sigma_sq_beta){

  // Initialize the required variables
  int n=X.nrow(), p=X.ncol() ;
  double uu, log_y, theta, theta_min, theta_max ;
  NumericVector beta(p), beta_tilde(p), beta_new(p), beta_tilde_new(p), Xbeta(n), Xbeta_new(n), nu(p) ;

  // Set the starting values of beta, beta_tilde and Xbeta = X * beta
  vec_copy(beta_init, beta) ;
  beta_tilde = sqrt((double)p*sigma_sq_beta)*beta ;
  product(X, beta, Xbeta) ;

  // Perform Elliptical Slice Sampling Algorithm

  // Draw an observation from Uniform(0,1)
  uu = runif(1)[0] ;

  // Draw an observation from p-variate Normal(0, sigma_sq_beta*I_p)
  nu = rnorm(p, 0.0, sqrt(sigma_sq_beta)) ;

  // Calculate log_y = Log_L(beta) + log(u) that defines a slice where log likelihood is atleast y
  log_y = log_L(y, Xbeta, xi, sigma_sq_eps, u) + log(uu) ;

  // Draw theta from Uniform(0, 2PI)
  theta = 2.0*PI*runif(1)[0] ;

  // Shrink the bracket of drawing theta
  theta_min = theta - 2.0*PI ;
  theta_max = theta ;

  // Draw a new beta_tilde from the ellipse passing through current beta_tilde and nu
  beta_tilde_new = cos(theta)*beta_tilde + sin(theta)*nu ;
  beta_new = beta_tilde_new / norm(beta_tilde_new) ;    // beta_new with unit norm
  product(X, beta_new, Xbeta_new) ;

  // If Log_L(beta_new) > log y, accept the value
  // Else draw a new observation from the slice
  while( log_L(y, Xbeta_new, xi, sigma_sq_eps, u) <= log_y ){
    if(theta < 0)
      theta_min = theta ;
    else
      theta_max = theta ;
    theta = (theta_max-theta_min)*runif(1)[0] + theta_min ;   // shrink the bracket
    beta_tilde_new = cos(theta)*beta_tilde + sin(theta)*nu ;  // try a new slice
    beta_new = beta_tilde_new / norm(beta_tilde_new) ;
    product(X, beta_new, Xbeta_new) ;
  }

  // Set beta, beta_tilde and Xbeta as the final values obtained from the algorithm
  vec_copy(beta_new, beta) ;
  vec_copy(beta_tilde_new, beta_tilde) ;
  vec_copy(Xbeta_new, Xbeta) ;

  return(List::create(Named("beta") = beta, Named("beta_tilde") = beta_tilde, Named("Xbeta") = Xbeta));
}

// Function to draw a sample of size n from truncated Normal distribution using an exact Hamiltonian Monte Carlo Algorithm
// x_init ~ d-variate Normal(mu, Sigma)
// constraints: <ff(j,_), x_init> + gg[j] \geq 0, j=1,...,m
// n: required sample size
// mu, Sigma : Parameters of the normal distribution
// x_init : initial value of the observation that follows the d-variate Normal distribution
// ff, gg : functions that define constraints
// n_burn : Burn in period of the Monte Carlo Algorithm
NumericMatrix rtmvnormHMC(int n, const NumericVector& mu, const NumericMatrix& Sigma,
                          const NumericVector& x_init, const NumericMatrix& ff,
                          const NumericVector& gg, int n_burn){

  // Initialize required varaibles
  int d=mu.size(), m=ff.nrow(), h ;
  double u, phi, tmp, tmp2, T_end, T_h, alpha_h ;
  NumericVector x(d), s(d), a(d), b(d), T(m), x_dot(d), g(m) ;
  NumericMatrix f(m,d), res(n,d), Sigma_chol_L(d,d), Sigma_chol_U(d,d) ;

  // Sigma = Sigma_col_L * Sigma_col_U by Cholesky Decomposition
  chol(Sigma, Sigma_chol_U) ;
  transpose(Sigma_chol_U, Sigma_chol_L) ;

  // x = (Sigma_chol_L)^-1 (x_init-mu) ~ d-variate Normal(0, I_d)
  product(solve(Sigma_chol_L), x_init-mu, x) ;

  // Adjust the constraints on x_init to apply them on x
  for(int j=0; j<m; j++){
    g[j] = innerProduct(ff(j,_), mu) + gg[j] ;
    for(int i=0; i<d; i++)
      f(j,i) = innerProduct(Sigma_chol_U(i,_), ff(j,_)) ;
  }

  // Perform the HMC algorithm
  for(int nn=0; nn<n+n_burn; nn++){

    s = rnorm(d) ;    // Draw s ~ Normal(0, I_d)
    vec_copy(s, a) ;   // a = s
    vec_copy(x, b) ;   // b = s

    T_end = PI/2.0 ;    // T_end = PI/2, end time point

    while(1){
      for(int j=0; j<m; j++){

        // u = sqrt((<f_j,a>)^2 + (<f_j,b>)^2)
        tmp = innerProduct(f(j,_), a) ;
        tmp2 = innerProduct(f(j,_), b) ;
        u = sqrt(tmp*tmp + tmp2*tmp2) ;

        // phi = tan^-1 (-tmp/tmp2)
        phi = atan2(-1.0*tmp, tmp2) ;

        if((u < g[j]) || (u < -1.0*g[j]))     // condition to check if the jth particle has hit the wall
          T[j] = T_end ;    // Set the time point of the jth particle as the end point
        else{
          // if the jth particle is still along trajectory, calculate the time to hit the wall using
          // u * cos(T_j + phi) + g_j = 0
          T[j] = acos(-1.0*g[j]/u) - phi ;
          if(T[j] < 0.0){
            T[j] += 2.0*PI ;
            // Rprintf("1\n") ;
          }
        }
      }

      // T_h = min {T_1, T_2, ..., T_m}
      h = 0 ;
      T_h = T[0] ;
      for(int j=1; j<m; j++){
        if(T[j] < T_h){
          T_h = T[j] ;
          h = j ;
        }
      }

      if(T_h < T_end){
        // if atleast one of the particles has not hit the wall, Set:
        // x_i = a_i sin(T_h) + b_i cos(T_h)
        // x_dot_i = a_i sin(T_h) - b_i cos(T_h)
        for(int i=0; i<d; i++){
          tmp = sin(T_h - 1e-10) ;
          tmp2 = cos(T_h - 1e-10) ;
          x[i] = a[i]*tmp + b[i]*tmp2 ;
          x_dot[i] = a[i]*tmp2 - b[i]*tmp ;
        }

        // set a_i as the reflected velocity at time T_h
        alpha_h = innerProduct(f(h,_), x_dot) / norm(f(h,_),2,true) ;
        for(int i=0; i<d; i++)
          a[i] = x_dot[i] - 2.0*alpha_h*f(h,i) ;

        // Set b_i = x_i
        vec_copy(x, b) ;

        // Decrease the end time point by T_h
        T_end -= T_h ;
      }
      else{

        // if all particles have hit the wall, Set x_i = a_i sin(T_h) + b_i cos(T_h)
        // and stop the algorithm
        for(int i=0; i<d; i++){
          tmp = sin(T_end) ;
          tmp2 = cos(T_end) ;
          x[i] = a[i]*tmp + b[i]*tmp2 ;
        }
        break ;
      }
    }
    // Store samples after burn-in period
    if(nn >= n_burn)
      for(int i=0; i<d; i++)
        res(nn-n_burn,i) = innerProduct(Sigma_chol_L(i,_), x) + mu[i] ;
    // Obtained x is a sample from truncated Normal(0, I_d)
    // Transform it to get the required truncated Normal(mu, Sigma)
  }
  return(res) ;
}

// Function to Draw a sample from the Posterior Distribution of xi
// y : nx1 response variable
// Xbeta : nx1 vector containing X*beta, X: nxp matrix of covariates, beta: px1 predictors
// xi_init: starting value of xi
// Sigma_xi : Prior variance of xi
// u : vector of knots ranging from -1 to 1
// sigma_sq_eps : error variance of the model
// monotone : TRUE OR FALSE denoting whether the single index function is monotone or not
// n_HMC : number of the HMC iterations in truncated normal sampling
NumericVector update_xi(const NumericVector& y, const NumericVector& Xbeta, const NumericVector& xi_init,
               const NumericMatrix& Sigma_xi, const NumericVector& u, double sigma_sq_eps,
               bool monotone, int n_HMC){

  // Initialize required variables
  int n = y.size(), L = xi_init.size() - 1 ;
  NumericVector nu_xi(L+1), tmp_L(L+1), xi(L+1);
  NumericMatrix Psi(n,L+1), PsiPsi(L+1,L+1), Omega_xi(L+1, L+1), B_xi(L+1,L+1), B_xi_inv(L+1,L+1), ff(L+1,L+1) ;

  //Psi_il = psi_l(Xbeta[i])
  for(int i=0; i<n; i++){
    for(int l=0; l<L+1; l++){
      Psi(i,l) = psi_l(Xbeta[i], l, u) ;
    }
  }
  AtA(Psi, PsiPsi) ;    //PsiPsi = Psi^{\top}*Psi
  solve(Sigma_xi, Omega_xi) ;   //Omega_xi = Sigma_xi^(-1)

  // Defining the constraints of xi. Here xi follows a truncated Normal distribution N_+(0, Sigma_xi).
  // i.e. xi_l >= 0 for all l i.e. <F_l,xi_l> + g_l = 0
  // F = (F_0,...,F_L) = I_{L+1} and g = (0,...,0)
  for(int l=0; l<L+1; l++)
    ff(l,l) = 1.0 ;

  // Posterior variance of xi is B_xi = [Omega_xi + (Psi^\top Psi)/sigma_sq_eps]^(-1)
  product(PsiPsi, 1.0/sigma_sq_eps, B_xi) ; // B_xi: temporary mtx
  sum(B_xi, Omega_xi, B_xi_inv) ;
  solve(B_xi_inv, B_xi) ;

  // Posterior mean of xi is nu_xi = (B_xi Psi^\top y)/ sigma_sq_eps
  Atx(Psi, y, tmp_L) ;
  product(B_xi, tmp_L, nu_xi) ;
  for(int l=0; l<L+1; l++)
    nu_xi[l] /= sigma_sq_eps ;

  // if g is not monotone, then xi's need not be non negative.
  // then draw xi from Normal distribution with mean = nu_xi and variance = B_xi
  if(monotone==false){
    rmvnorm(nu_xi, B_xi, xi) ;
  }
  else{
    // if g is monotone, then xi_l >= 0 for all l
    // draw xi from truncated Normal distribution N_+(nu_xi, B_xi) using HMC algorithm
    xi = rtmvnormHMC(1+n_HMC, nu_xi, B_xi, xi_init, ff, rep(0.0, L+1), 0)(n_HMC,_) ;
  }
  return(xi);
}

// Function to Draw a sample from the Posterior Distribution of xi
// y : nx1 response variable
// Xbeta : nx1 vector containing X*beta, X: nxp matrix of covariates, beta: px1 predictors
// xi : (L+1)x1 vector of coefficients to construct the single index function
// u : vector of knots ranging from -1 to 1
// a_eps, b_eps : Hyperparameters specifying Prior distribution of sigma_sq_eps which is Inv Gamma(a_eps, b_eps)
double update_sigma_sq_eps(const NumericVector& y, const NumericVector& Xbeta, const NumericVector& xi,
                           const NumericVector& u, double a_eps, double b_eps){

  // Initialize required variables
  int n = y.size();
  NumericVector eps(n);

  // eps = y - g(X^\top beta)
  for(int i=0; i<n; i++){
    eps[i] = y[i] - g(Xbeta[i], xi, u) ;
  }
  // sigma_sq_eps ~ Inv Gamma (a_eps + n/2, 1/(b_eps + (eps^\top eps)/2))
  double sigma_sq_eps = 1.0/rgamma(1, a_eps+0.5*n, 1.0/(b_eps + 0.5*innerProduct(eps,eps)))[0] ;
  return(sigma_sq_eps);
}

// MCMC algorithm to draw samples from the conditional posterior distributions of xi, beta and sigma_sq_eps
// y : nx1 response variable
// X: nxp matrix of covariates
// beta_init: px1 vector containing starting value of predictors
// xi_init: (L+1)x1 vector containing starting value of coefficients that construct the single index function
// Sigma_xi : Prior variance of xi
// u : vector of knots ranging from -1 to 1
// monotone : TRUE OR FALSE denoting whether the single index function is monotone or not
// n_HMC : number of the HMC iterations in truncated normal sampling
// sigma_sq_beta : Prior variance of beta
// sigma_sq_eps_init : Starting value of error variance of the model
// a_eps, b_eps : Hyperparameters specifying Prior distribution of sigma_sq_eps which is Inv Gamma(a_eps, b_eps)
// M_burn : Burn-in period of the MCMC algorithm
// M : required size of the  MCMC sample
// [[Rcpp::export]]
List monotoneSIM_c(const NumericVector& y, const NumericMatrix& X, const NumericVector& beta_init,
                   const NumericVector& xi_init, const NumericMatrix& Sigma_xi, const NumericVector& u,
                   bool monotone, int n_HMC, double sigma_sq_beta, double sigma_sq_eps_init,
                   double a_eps, double b_eps, int M_burn, int M){

  // Initialize required variables
  int n=y.size(), p=X.ncol(), L=xi_init.size()-1 ;
  double sigma_sq_eps ;
  NumericVector beta(p), xi(L+1), Xbeta(n) ;
  NumericVector sigma_sq_eps_vec(M) ;
  NumericMatrix xi_mtx(M,L+1), beta_mtx(M,p) ;

  // Copy the initial provided values of xi, beta and sigma_sq_eps and calculate Xbeta
  vec_copy(xi_init, xi) ;
  vec_copy(beta_init, beta) ;
  product(X, beta, Xbeta) ;
  sigma_sq_eps = sigma_sq_eps_init ;

  // MCMC algorithm
  for(int m=0; m<M+M_burn; m++){
    // update xi
    xi = update_xi(y, Xbeta, xi, Sigma_xi, u, sigma_sq_eps, monotone, n_HMC) ;
    // update beta
    // beta_ES stores beta updated using Elliptical Slice Sampler
    List beta_ES = update_beta(y, X, xi, beta, u, sigma_sq_eps, sigma_sq_beta) ;
    beta = beta_ES["beta"] ;
    Xbeta = beta_ES["Xbeta"] ;
    // update sigma_sq_eps
    sigma_sq_eps = update_sigma_sq_eps(y, Xbeta, xi, u, a_eps, b_eps) ;

    if(m >= M_burn){
      // store the samples after Burn-in
      xi_mtx(m-M_burn,_) = xi ;
      beta_mtx(m-M_burn,_) = beta ;
      sigma_sq_eps_vec[m-M_burn] = sigma_sq_eps ;
    }
  }
  // return the MCMC sample
  return(List::create(Named("xi") = xi_mtx, Named("beta") = beta_mtx, Named("sigma_sq_eps") = sigma_sq_eps_vec));
}

// Function to calculate the value of g(x) for a grid of x values and xi vectors.
// xi_mtx : Matrix whose rows contain (L+1)x1 vector xi_j, j = 1(1)M
// grid_x : Nx1 vector contaning the x-values
// u : (L+1)x1 vector of knots from -1 to 1
// [[Rcpp::export]]
NumericMatrix g_mtx(const NumericMatrix& xi_mtx,const NumericVector& grid_x, const NumericVector& u){
  // Initialize variables
  int M = xi_mtx.nrow(), N = grid_x.size() ;
  NumericMatrix res(M, N) ;

  // Calculate g(x) for every x and every vector xi
  for(int i = 0; i<M; i++){
    for(int j = 0; j<N; j++){
      res(i, j) = g(grid_x(j), xi_mtx(i,_), u) ;
    }
  }
  return(res) ;
}
