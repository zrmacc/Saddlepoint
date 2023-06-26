// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Score ECGF
//' 
//' Evaluate \eqn{K_{n}^{S}(t)}.
//' 
//' @param g Genotype.
//' @param mu Mean of Y under the null.
//' @param t Evaluation point.
//' @param y Phenotype.
//' @return Numeric evaluation.
// [[Rcpp::export]]
SEXP ECGFScpp(
  const arma::colvec g,
  const arma::colvec mu,
  const double t,
  const arma::colvec y
){
  const int n = g.size();
  const arma::colvec r = (y - mu);
  double outer = 0; 
  for(int i=0; i<n; i++){
    double inner = 0;
    for(int j=0; j<n; j++) {
      inner += std::exp(t * g(i) * r(j));
    }
    outer += std::log(inner / n);
  }
  return Rcpp::wrap(outer);
}


//' First Derivative of Score ECGF
//' 
//' Evaluate \eqn{\dot{K}_{n}^{S}(t)}.
//' 
//' @param g Genotype.
//' @param mu Mean of Y under the null.
//' @param t Evaluation point.
//' @param y Phenotype.
//' @return Numeric evaluation.
// [[Rcpp::export]]
SEXP dECGFScpp(
   const arma::colvec g,
   const arma::colvec mu,
   const double t,
   const arma::colvec y
){
 const int n = g.size();
 const arma::colvec r = (y - mu);
 double outer = 0;
 for(int i=0; i<n; i++){
   double sum_0 = 0;
   double sum_1 = 0;
   for(int j=0; j<n; j++) {
     double wij = std::exp(t * g(i) * r(j));
     sum_0 += wij;
     sum_1 += wij * g(i) * r(j);
   }
   outer += (sum_1 / sum_0);
 }
 return Rcpp::wrap(outer);
}


//' Second Derivative of Score ECGF
//' 
//' Evaluate \eqn{\ddot{K}_{n}^{S}(t)}.
//' 
//' @param g Genotype.
//' @param mu Mean of Y under the null.
//' @param t Evaluation point.
//' @param y Phenotype.
//' @return Numeric evaluation.
// [[Rcpp::export]]
SEXP ddECGFScpp(
   const arma::colvec g,
   const arma::colvec mu,
   const double t,
   const arma::colvec y
){
 const int n = g.size();
 const arma::colvec r = (y - mu);
 double outer = 0;
 for(int i=0; i<n; i++){
   double sum_0 = 0;
   double sum_1 = 0;
   double sum_2 = 0;
   for(int j=0; j<n; j++) {
     double wij = std::exp(t * g(i) * r(j));
     sum_0 += wij;
     sum_1 += wij * g(i) * r(j);
     sum_2 += wij * g(i) * g(i) * r(j) * r(j);
   }
   double num = sum_2 * sum_0 - sum_1 * sum_1;
   double denom = sum_0 * sum_0;
   outer += (num / denom);
 }
 return Rcpp::wrap(outer);
}
