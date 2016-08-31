#' @noRd
new.edgenet.obj <- function(res, lambda, psigx, psigy, X, Y)
{
  coefficients <- res$coefficients
  intr         <- res$intercept
  rownames(coefficients) <- colnames(X)
  colnames(coefficients) <- colnames(Y)
  ret <- list(coefficients=coefficients, 
              intercept=intr,
              lambda=lambda,
              psigx=psigx,
              psigy=psigy)
  ret
}