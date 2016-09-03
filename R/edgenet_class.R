#' edgenetfit
#' 
#' Class that represents an edgenet fit.
#' 
#' @noRd
#' 
#' @slot coefficients  the estimated coefficients
#' @slot intercept  the estimated intercept
#' @slot lambda the shrinkage parameter for the LASSO
#' @slot psigy  the penalization for G.Y
#' @slot psigx  the penalization for G.Y
#' @slot family  family of the response variable
#' @slot call  match.call() 
setClass(
  "edgenetfit",
   representation = representation(
     coefficients="numeric",
     intercept="numeric",
     lambda="numeric",
     psigy="numeric",
     psigx="numeric",
     family="character",
     call="call"),
  validity = function(object)
  {
    length(coefficients) == length(intercept)
  }
)
