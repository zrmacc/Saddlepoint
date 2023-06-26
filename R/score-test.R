#' Calculate P-value
#' 
#' Saddlepoint p-value calculation.
#'
#' @param sp Saddlepoint object.
#' @return Numeric p-value
CalcPval <- function(sp) {
  s <- sp@score
  w <- sp@w
  v <- sp@v
  arg <- w + (1 / w) * log(v / w)
  if (s < 0) {
    p <- 2 * stats::pnorm(arg, lower.tail = TRUE)
  } else {
    p <- 2 * stats::pnorm(arg, lower.tail = FALSE)
  }
  return(p)
}


#' Score Test
#' 
#' Perform a score test of the null hypothesis that genotype is unrelated
#' to the phenotype, where the p-value is calculated by empirical 
#' saddlepoint approximation. 
#' 
#' @param g Genotypic residuals.
#' @param mu Mean of Y under the null.
#' @param sigma2 Residual variance under the null.
#' @param y Phenotype.
#' @param alpha Type I error.
#' @param radius Radius about the origin to search for saddlepoint. 
#' @param tol Tolerance.
#' @param x Null model covariates, if not already regressed from genotype.
#' @return Data.frame containing the score and p-value.
#' @export 
ScoreTest <- function(
    g, 
    mu, 
    sigma2,
    y,
    alpha = 0.05,
    radius = 10,
    tol = 1e-6,
    x = NULL
) {
  
  # Residualize genotype.
  if (!is.null(x)) {
    g <- GenoResid(g = g, x = cbind(1, x))
  }
  
  # Score statistic.
  score <- sum(g * (y - mu))
  score_se <- sqrt(sigma2 * sum(g^2))
  score_stat <- abs(score) / score_se
  
  # Choose between saddlepoint and normal approximation. 
  z <- stats::qnorm(1 - alpha / 2)
  if (score_stat > z) {
    
    # Case of saddlepoint approximation.
    sp <- FindSaddlepoint(
      g = g,
      mu = mu,
      y = y,
      radius = radius,
      tol = tol
    )
    pval <- CalcPval(sp)
    used_sp <- TRUE
    
  } else {
    
    # Case of normal approximation.
    pval <- 2 * stats::pnorm(q = score_stat, lower.tail = FALSE)
    used_sp <- FALSE
    
  }

  # Output.
  out <- data.frame(
    score_stat = sign(score) * score_stat,
    p = pval,
    used_sp = used_sp
  )
  
  return(out)
}

