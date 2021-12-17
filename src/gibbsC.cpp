#include <Rcpp.h>
using namespace Rcpp;
//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param thin the number of between-sample random numbers
//' @param n the parameter of Binomial distribution
//' @param a the parameter of Beta distribution
//' @param b the parameter of Beta distribution
//' @return a random sample 
//' @examples
//' \dontrun{
//' rnC <- gibbsC(100,10,20,4,1) 
//' }
//' @export

// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int thin, int n, int a, int b) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rbinom(1, n, y)[0];
      y = rbeta(1, (x+a), (n-x+b))[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}
