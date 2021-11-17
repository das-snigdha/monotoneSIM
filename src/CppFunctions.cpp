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
