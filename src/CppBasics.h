using namespace Rcpp;

// Copy elements of one vector containing real numbers to another vector
void vec_mem_cpy(NumericVector from, NumericVector to){
  for(int i=0; i<to.size(); i++)
    to[i] = from[i] ;
}

// Copy elements of one vector to another upto the n-th element
void vec_mem_cpy(NumericVector from, NumericVector to, int n){
  for(int i=0; i<n; i++)
    to[i] = from[i] ;
}

// Copy elements of one vector containing integers to another vector
void vec_mem_cpy(IntegerVector from, IntegerVector to){
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

// Calculate the product of 2 numeric matrices and return it as a matrix
NumericMatrix product(NumericMatrix A, NumericMatrix B){
  NumericMatrix res(A.nrow(), B.ncol()) ;
  for(int i=0; i<A.nrow(); i++)	{
    for(int j=0; j<B.ncol(); j++)	{
      res(i,j) = 0.0 ;
      for(int k=0; k<A.ncol(); k++)
        res(i,j) += A(i,k)*B(k,j) ;    // perform matrix multiplication
    }
  }
  return(res) ;   // return the resulting matrix
}

