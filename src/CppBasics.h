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
