#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Triangular basis function h(x) = (1-|x|)_+
double h(double x){
  if(x > -1 && x < 1)
    return(1.0 - abs(x));   // return (1-|x|) if -1 <= x <= 1
  else
    return(0.0);    // return 0 otherwise
}

// h_l(.) = B-spline basis function of order 2
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
// Here, u is a vector of specified knots from -1 to 1
// [[Rcpp::export]]
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
