#' @noRd
#' @export
#' @author Simon Dirmeier, \email{netreg@@simon-dirmeier.net}
print.gaussian.edgenet <- 
function
(
 x,
 ...
)
{
  cat("\nCall: ")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nIntercept:\n")
  print(x$intercept)
  cat("Lambda:\n")
  print(x$lambda)
  cat("Psi_gx:\n")
  print(x$psigx)
  cat("Psi_gy:\n")
  print(x$psigy)
}

#' @noRd
#' @export
#' @author Simon Dirmeier, \email{netreg@@simon-dirmeier.net}
print.gaussian.cv.edgenet <-
function
(
  x, 
  ...
)
{
  cat("\nCall: ")
  print(x$call)
  cat("\nParameters: ")
  cat(paste("lambda=", x$lambda, ", psi_gx=", x$psigx, ", psi_gy=", x$psigy, "\n", sep=""))
  cat("\nFamily: ")
  cat(x$family, "\n")
}
