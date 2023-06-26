#' Data Generating Process
#' 
#' Simulates example data. 
#'
#' @param n Sample size.
#' @param maf Genotype minor allele frequency.
#' @param pve_g Proportion of variation explained by genotype.
#' @param pve_x Proportion of variation explained by covariates.
#' @return Data.frame.
#' @export
DGP <- function(
    n,
    maf = 0.25,
    pve_g = 0.0, 
    pve_x = 0.20
  ) {
  
  pve_e <- 1 - (pve_g + pve_x)
  if (pve_e < 0) {
    msg <- glue::glue("Environmental PVE {pve_e} is negative.")
    stop(msg)
  }
  
  # Genotype.
  g <- stats::rbinom(n = n, size = 2, prob = maf)
  bg <- sqrt(pve_g / (2 * maf * (1 - maf)))
  
  # Covariate.
  x <- stats::rnorm(n = n)
  bx <- sqrt(pve_x)
  
  # Residual.
  e <- stats::rnorm(n = n)
  be <- sqrt(pve_e)
  
  # Phenotype.
  y <- g * bg + x * bx + e * be
  
  out <- data.frame(
    g = g,
    x = x,
    y = y
  )
  return(out)
}

