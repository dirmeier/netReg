#' @noRd
#' @export
#' @author Simon Dirmeier, \email{netreg@@simon-dirmeier.net}
print.edgenet <- 
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
  cat("\nParameters: ")
  cat(paste("lambda=", x$lambda, ", psi_gx=", x$psigx, ", psi_gy=", x$psigy, "\n", sep=""))
  cat("\nFamily: ")
  cat(x$family, "\n")
}

#' @noRd
#' @export
#' @author Simon Dirmeier, \email{netreg@@simon-dirmeier.net}
print.cv.edgenet <-
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
