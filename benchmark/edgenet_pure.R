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


 
#' Demo script of an R implementation of edgenet for benchmarking
#' 
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
        l <- .uccd(p, q,
              thresh, maxit,
              lambda,
              psigx, psigy,
              TXX, TXY,LX, LY,
              B, B.old,
              i);
        B     <- l$B
        B.old <- l$B.old
     }
    
    if (sum(abs(B - B.old)) <= thresh | niter > maxit) break
  } 
  
  list(B=B, mu=netReg:::intercept(Y, X, B, n))
}

.uccd <- function(p, q, thresh, maxit, lambda, psigx, psigy,
                  TXX, TXY,LX, LY, B, B.old, qi)
{
  niter <- 0
  repeat
  {
    niter <- niter + 1
    for (pi in seq(p))
    {
      B.old[pi, qi] <- B[pi, qi]
      s <-  (TXY[pi, qi] + (TXX[pi, pi] * B[pi, qi])) - sum(TXX[pi, ] %*% B[,qi])
      norm <- TXX[pi, pi]
      if (psigx > 0.00001)
      {
        x.penalty <- .x.graph.penalize(LX, B, pi, qi)
        s         <- s - 2 * psigx * x.penalty;
        norm      <- norm + 2 * psigx * LX[pi, pi];
      }
      if (psigy > 0.00001)
      {
        y.penalty <- .y.graph.penalize(LY, B, pi, qi)
        s         <- s - 2 * psigy * y.penalty
        norm      <- norm + 2 * psigy * LY[qi, qi]
      }
      B[pi, qi] <- .softnorm(s, lambda, norm)
    }
    if (sum(abs(B[,qi] - B.old[,qi])) <= thresh | niter < maxit) break  
  }
  list(B=B, B.old=B.old)
}

.softnorm <- function(s, lalph, norm)
{
  sabs  <- abs(s)
  r <- 0.0
  if (lalph < sabs)
  {
    if (s > 0)
      r <- (s - lalph) / norm
    else
      r <- (s + lalph) / norm
  }
  r
}

.x.graph.penalize <- function(G, cfs, pi, qi)
{
  penalty <- 0
  if (pi <= nrow(G) && pi <= ncol(G))
  {
    penalty <- -G[pi, pi] * cfs[pi, qi] + sum(G[pi, ] * cfs[ ,qi])
  }
  penalty
}

.y.graph.penalize <- function(G, cfs, pi, qi)
{
  penalty <- 0
  if (qi <= nrow(G) && qi <= ncol(G))
  {
    penalty <- - cfs[pi, qi] * G[qi, qi] + sum(cfs[pi, ] %*% G[ ,qi])
  }
  penalty
}

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
