#' @noRd
.gaussian.cv.edgenet <-
function
(
  X, Y, G.X, G.Y,
  n, p, q, 
  psigx, psigy, 
  maxit, thresh,
  nfolds, foldid
)
{
  print(psigx)
  print(psigy)
  cv <- .Call("gauss_cv_edgenet", 
                X, Y,
                G.X, G.Y, 
                as.integer(n), as.integer(p), as.integer(q),
                as.double(psigx),  as.double(psigy),
                as.integer(maxit), as.double(thresh),
                as.integer(nfolds), as.integer(foldid),
                as.integer(length(foldid)),
                PACKAGE="netReg")
  ret <- list(lambda=cv$shrinkage_parameters[1],
              psigx=cv$shrinkage_parameters[2], 
              psigy=cv$shrinkage_parameters[3],
              foldids=cv$fold_ids)
  class(ret) <- "gaussian.cv.edgenet"
  ret
}

  