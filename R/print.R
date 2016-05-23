#' Print method for edgenet objects
#' 
#' @export
#' 
#' @author Simon Dirmeier, \email{netreg@@simon-dirmeier.net}
#' @description Print method for objects of class \emph{edgenet} for easy vizualization of fitted edge-regularized models
#' 
#' @param x  an object of class \emph{netreg.graph}
#' @param ...  further arguments
print.edgenet <- 
function
(
 x,
 ...
)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nIntercept:\n")
  print(x$intercept)
}