// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Genotypic Residual
//' 
//' Residualize genotype with respect to covariates.
//' 
//' @param g Genotype.
//' @param x Covariates.
//' @return Numeric evaluation.
//' @export
// [[Rcpp::export]]
SEXP GenoResid(
  const arma::colvec g,
  const arma::mat x
){
  arma::colvec out = g - x * arma::solve(x.t() * x, x.t() * g, arma::solve_opts::likely_sympd);
  return Rcpp::wrap(out);
}

