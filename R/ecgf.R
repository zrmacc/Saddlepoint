#' Empirical CGF of R
#' 
#' Calculate the empirical cumulant generating function of the residual.
#' 
#' @param mu Mean of Y under the null.
#' @param y Phenotype.
#' @return Function, ECGF(t).
ECGFR <- function(mu, y) {
  ecgf <- function(t) {
    resid <- (y - mu)
    n <- length(resid)
    mi <- mean(exp(t * resid))
    ki <- log(mi)
    return(ki)
  }
  return(ecgf)
}


#' Empirical CGF of S
#' 
#' Calculate the empirical cumulant generating function of the score.
#' 
#' @param g Genotype.
#' @param mu Mean of Y under the null.
#' @param y Phenotype.
#' @return Function, ECGF(t).
ECGFS <- function(g, mu, y) {
  out <- function(t) {
    k_gt <- sapply(g * t, ECGFR(mu, y))
    sum_k_gt <- sum(k_gt)
    return(sum_k_gt)
  }
  return(out)
}
