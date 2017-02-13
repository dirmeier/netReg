# netReg: graph-regularized linear regression models.
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of netReg.
#
# netReg is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# netReg is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with netReg. If not, see <http://www.gnu.org/licenses/>.


#' 
#' Demo script of an R implementation of edgenet for benchmarking
#' @noRd
#' 
edgenet.pure.R <- function(X, Y, G.X, G.Y, lambda, psigx, psigy, thresh, maxit)
{
  n     <- nrow(X)
  p     <- ncol(X)
  q     <- ncol(Y)
  B     <- matrix(rnorm(p*q), p, q)
  B.old <- matrix(0, p, q)
  
  LX   <- .laplacian(G.X)
  LY   <- .laplacian(G.Y)
  TXX <- t(X) %*% X
  TXY <- t(X) %*% Y
  
  niter <- 0
  repeat 
  {
    niter <- niter + 1
     for (i in seq(q))
     {
        .uccd(p, q,
              thresh, maxit,
              lambda,
              psigx, psigy,
              TXX, TXY,LX, LY,
              B, B.old,
              i);
     }
    
    if (sum(abs(B-B.old)) <= thresh | niter < maxit) break
  } 
}

.uucd <- function(p, q, thresh, maxit, lambda, psigx, psigy,
                  TXX, TXY,LX, LY, B, B.old, qi)
{
  niter <- 0
  repeat
  {
    niter <- niter + 1
    for (pi in seq(p))
    {
      Bold[pi, qi] <- B[pi, qi]
      s <-  (TXY[pi, qi] + (TXX[pi, pi] * B[pi, qi])) - sum(TXX[pi, ] %*% cfs[,qi])
      norm <- TXX[pi, pi]
      if (psigx != 0) 
      {
        
      }
    }
    if (sum(abs(B[,qi] - B.old[,qi])) <= thresh | niter < maxit) break  
  }
    
}

    set_params
    (s, norm, TXX, TXY, coef,
      LX, LY, P, Q, pi, qi, psigx, psigy, false);
    //                // soft-thresholded version of estimate
    coef(pi, qi) = softnorm(s, lalph, enorm * norm)



.laplacian <- function(X)
{
  m          <- nrow(X)
  degrees    <- apply(X, 1, sum)
  laplacian  <- matrix(0, m, m)
  for (i in 1:m)
  {
    for (j in 1:m)
    {
      if (i == j & degrees[i] != 0) {
        laplacian[i, j] = 1 - (X[i, j] / degrees[i])
      } else if (i != j & X[i, j] != 0) {
        laplacian[i, j] = - X[i, j] / sqrt(degrees[i] * degrees[j]);
      } else {
        laplacian[i,  j] = 0.0
      }
    }
  }
  laplacian
}
