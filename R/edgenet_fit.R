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


#' @noRd
.fit <-function(X, Y, G.X, G.Y, 
                lambda, psigx, psigy, 
                thresh, maxit, family)
{
  n <- dim(X)[1]                              
  p <- dim(X)[2]     
  q <- dim(Y)[2]
  res <- .Call("edgenet", X, Y, G.X, G.Y, 
               as.integer(n), as.integer(p), as.integer(q),
               as.double(lambda), as.double(psigx),  as.double(psigy),
               as.integer(maxit), as.double(thresh), 
               as.character(family),
               PACKAGE="netReg")
  # finalize output
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
  class(ret) <- paste0(family, ".edgenet")
  ret
}