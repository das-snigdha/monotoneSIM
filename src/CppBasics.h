using namespace Rcpp;

// Copy elements of one vector containing real numbers to another vector
void vec_copy(NumericVector from, NumericVector to){
  for(int i=0; i<to.size(); i++)
    to[i] = from[i] ;
}

// Calculate the p-norm of a vector as ||x||_p = (\sum |x_i|^p)^(1/p)
double norm(NumericVector x, int p=2, bool power=false){
  double res=0.0 ;

  // res = \sum |x_i|^p
  for(int i=0; i<x.size(); i++)
    res += pow(abs(x[i]), p) ;
  if(power == true)
    return(res) ;   //if power = true, return \sum |x_i|^p
  else
    return(pow(res, 1.0/(double)p)) ;   // otherwise return (\sum |x_i|^p)^(1/p)
}

// Calculate the product of 2 numeric matrices and store it in "res" matrix passed as an argument
void product(NumericMatrix A, NumericMatrix B, NumericMatrix res){
  for(int i=0; i<A.nrow(); i++)	{
    for(int j=0; j<B.ncol(); j++)	{
      res(i,j) = 0.0 ;
      for(int k=0; k<A.ncol(); k++)
        res(i,j) += A(i,k)*B(k,j) ;   // perform matrix multiplication
    }
  }
}

// Multiply a matrix with a vector and store it in "res" vector passed as an argument
void product(NumericMatrix A, NumericVector x, NumericVector res){
for(int i=0; i<A.nrow(); i++)	{
  res[i] = 0.0 ;
  for(int j=0; j<A.ncol(); j++)
    res[i] += A(i,j)*x[j] ; // res = A*x
  }
}

// Multiply a matrix with a scalar and store it in "res" matrix passed as an argument
void product(NumericMatrix A, double x, NumericMatrix res){
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      res(i,j) = x * A(i,j) ;   // Multply every element of A with x
}

// Calculate the transpose of a matrix and store in a matrix passed as an argument
void transpose(NumericMatrix A, NumericMatrix t_A){
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      t_A(j,i) = A(i,j) ;
}

// Calculate Cholesky Decomposition of a matrix and store in a matrix passed as an argument
void chol(NumericMatrix A, NumericMatrix res){
  int n = A.nrow() ;
  arma::mat tmp = arma::chol(as<arma::mat>(A)) ;
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      res(i,j) = tmp(i,j) ;
}

// Calculate the inverse of a matrix and store in a matrix passed as an argument
void solve(NumericMatrix A, NumericMatrix A_inv){
  arma::mat tmp = as<arma::mat>(A).i() ;

  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      A_inv(i,j) = tmp(i,j) ;
}

// Calculate the inverse of a matrix and return it
NumericMatrix solve(NumericMatrix A){
  arma::mat res = as<arma::mat>(A).i() ;
  return(wrap(res)) ;
}

// Calculate A^{\top} * x, A: matrix and X: vector
void Atx(NumericMatrix A, NumericVector x, NumericVector res){
  for(int i=0; i<A.ncol(); i++)	{
    res[i] = 0.0 ;
    for(int j=0; j<A.nrow(); j++)
      res[i] += A(j,i)*x[j] ;
  }
}

// Calculate A^{\top} * A, where A is a matrix
void AtA(NumericMatrix A, NumericMatrix res){
  // res = t(A) %*% A
  for(int i=0; i<A.ncol(); i++){
    for(int j=0; j<A.ncol(); j++){
      res(i,j) = 0.0 ;
      for(int k=0; k<A.nrow(); k++)
        res(i,j) += A(k,i)*A(k,j) ;
    }
  }
}

// Calculate the innerproduct of 2 vectors x and y
double innerProduct(NumericVector x, NumericVector y){
  double res=0.0 ;
  for(int i=0; i<x.size(); i++)
    res += x[i]*y[i] ;
  return(res) ;
}

// Draw a sample from multivariate Normal (mu, Sigma) and return it
void rmvnorm(NumericVector mu, NumericMatrix Sigma, NumericVector res){
  arma::mat Y = arma::randn(1,mu.size());
  arma::mat Z = arma::repmat(as<arma::vec>(mu),1,1).t() + Y * arma::chol(as<arma::mat>(Sigma)) ;
  for(int i=0; i<res.size(); i++)
    res[i] = Z(0,i) ;
}

// Add 2 matrices and store it in a matrix passed as an argument
void sum(NumericMatrix A, NumericMatrix B, NumericMatrix res){
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      res(i,j) = A(i,j) + B(i,j) ;
}
