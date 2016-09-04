#' @noRd
.fit <-function(X, Y, G.X, G.Y, 
                lambda, psigx, psigy, 
                thresh, maxit, family)
{
  n <- dim(X)[1]                              
  p <- dim(X)[2]     
  q <- dim(Y)[2]
  res <- .Call("edgenet_", X, Y, G.X, G.Y, 
               as.integer(n), as.integer(p), as.integer(q),
               as.double(lambda), as.double(psigx),  as.double(psigy),
               as.integer(maxit), as.double(thresh), 
               as.character(family),
               PACKAGE="netReg")
  
  coefficients <- res$coefficients
  intr         <- res$intercept
  rownames(coefficients) <- colnames(X)
  colnames(coefficients) <- colnames(Y)
  ret <- list(coefficients=coefficients, 
              intercept=intr,
              lambda=lambda,
              psigx=psigx,
              psigy=psigy)
  ret$family <- family
  class(ret) <- paste(family, ".edgenet", sep="")
  ret
}