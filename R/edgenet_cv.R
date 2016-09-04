#' @noRd
.cv <- function(X, Y, G.X, G.Y, 
                psigx, psigy, thresh, maxit, family,
                nfolds, foldid)
{
  n <- dim(X)[1]                              
  p <- dim(X)[2]     
  q <- dim(Y)[2]
  cv <- .Call("cv_edgenet_", X, Y,G.X, G.Y, 
               as.integer(n), as.integer(p), as.integer(q),
               as.double(psigx),  as.double(psigy),
               as.integer(maxit), as.double(thresh),
               as.integer(nfolds), as.integer(foldid),
               as.integer(length(foldid)),
               as.character(family),
               PACKAGE="netReg")
  ret <- list(lambda=cv$shrinkage_parameters[1],
              psigx=cv$shrinkage_parameters[2], 
              psigy=cv$shrinkage_parameters[3],
              foldids=cv$fold_ids)
  ret$family <- family
  ret
}

