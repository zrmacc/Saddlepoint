# Purpose: Define class returned by the parametric fitting and contrast functions.
# Updated: 2023-06-22


# -----------------------------------------------------------------------------
# Saddlepoint class.
# -----------------------------------------------------------------------------


#' Saddlepoint
#'
#' @slot ddksp Second derivative of ECGF at the saddlepoint.
#' @slot dksp First derivative of ECGF at the saddlepoint.
#' @slot ksp ECGF at the saddlepoint.
#' @slot score Score statistic.
#' @slot sp Saddlepoint.
#' @slot w Argument for saddlepoint approximation.
#' @slot v Argument for saddlepoint approximation.
#' @name saddlepoint-class
#' @rdname saddlepoint-class
#' @exportClass saddlepoint
setClass(
  Class = "saddlepoint", 
  representation = representation(
    ddksp = "numeric",
    dksp = "numeric",
    ksp = "numeric",
    score = "numeric",
    sp = "numeric",
    w = "numeric",
    v = "numeric"
    )
  )


# -----------------------------------------------------------------------------

#' Print Method for Saddlepoint
#'
#' Print method for objects of class \code{saddlepoint}.
#'
#' @param x An object of class \code{saddlepoint}.
#' @param digits Display digits.
#' @param ... Unused.
#' @export

print.saddlepoint <- function(x, digits = 3, ...) {

  # Display.
  cat("Saddlepoint:\n")
  print(signif(x@sp, digits = digits))
  cat("\n")
  
  cat("Score:\n")
  print(signif(x@score, digits = digits))
  cat("\n")

}


#' Show Method for Saddlepoint Class
#'
#' @param object An object of class \code{saddlepoint}.
#' @importFrom methods show
#' @rdname saddlepoint-method
setMethod(
  f = "show",
  signature = c(object = "saddlepoint"),
  definition = function(object) {
    print.saddlepoint(x = object)
  }
)


# -----------------------------------------------------------------------------
