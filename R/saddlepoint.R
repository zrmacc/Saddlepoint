# Purpose: Find the saddlepoint.
# Updated: 2023-06-23

#' Approximation Arguments
#' 
#' Calculate the quantities \eqn{w} and \eqn{v} utilized to calculate the
#' p-value.
#'
#' @param ddksp Second derivative of ECGF at the saddlepoint.
#' @param ksp ECGF at the saddlepoint.
#' @param score Score statistic.
#' @param sp Saddlepoint.
#' @return List containing the saddlepoint arguments.
#' @noRd
ApproxArgs <- function(ddksp, ksp, score, sp) {
  w <- sign(sp) * sqrt(2 * (sp * score - ksp))
  v <- sp * sqrt(ddksp)
  out <- list(w = w, v = v)
  return(out)
}


#' Find Saddlepoint
#' 
#' Calculate the saddlepoint: that \eqn{t} such that the first derivative of
#' the empirical cumulant generating function at \eqn{t} equals the observed
#' value of the score statistic.
#' 
#' @param g Genotype.
#' @param mu Mean of Y under the null.
#' @param y Phenotype.
#' @param radius Radius about the origin to search for saddlepoint. 
#' @param tol Tolerance.
#' @return Saddlepoint.
#' @export 
FindSaddlepoint <- function(
    g, 
    mu, 
    y,
    radius = 10,
    tol = 1e-6
) {
  
  # Score.
  r <- (y - mu)
  score <- sum(g * r)
  
  # Find saddlepoint.
  aux <- function(t){
    delta <- dECGFScpp(g = g, mu = mu, t = t, y = y) - score
    return(delta)
  }
  zero <- stats::uniroot(
    f = aux,
    lower = -radius,
    upper = +radius,
    extendInt = "upX",
    tol = tol
  )
  sp <- zero$root
  
  # Evaluations.
  ksp <- ECGFScpp(g = g, mu = mu, t = sp, y = y)
  dksp <- dECGFScpp(g = g, mu = mu, t = sp, y = y)
  ddksp <- ddECGFScpp(g = g, mu = mu, t = sp, y = y)
  
  # Saddlepoint arguments.
  args <- ApproxArgs(
    ddksp = ddksp, 
    ksp = ksp, 
    score = score, 
    sp = sp
  )
  
  # Evaluations.
  out <- methods::new(
    Class = "saddlepoint",
    ddksp = ddksp,
    dksp = dksp,
    ksp = ksp,
    score = score,
    sp = sp,
    w = args$w,
    v = args$v
  )
  return(out)
}

