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
