#' @noRd
binom.edgenet <- function(X, Y, G.X, G.Y, 
                          n, p, q, 
                          lambda, psigx, psigy, 
                          maxit, thresh)
{
  res  <-  .Call("binom_edgenet", 
                 X, Y,
                 G.X, G.Y, 
                 as.integer(n), as.integer(p), as.integer(q),
                 as.double(lambda), 
                 as.double(psigx),  as.double(psigy),
                 as.integer(maxit), as.double(thresh), 
                 PACKAGE="netReg")
  ret        <- new.edgenet.obj(res, lambda, psigx, psigy, X, Y)
  class(ret) <- "binomial.edgenet"
  ret
}