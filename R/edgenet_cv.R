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
.cv <- function(X, Y, G.X, G.Y, 
                psigx, psigy, thresh, maxit, family,
                nfolds, foldid)
{
  n <- dim(X)[1]                              
  p <- dim(X)[2]     
  q <- dim(Y)[2]
  cv <- .Call("cv_edgenet", X, Y,G.X, G.Y, 
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
  class(ret) <- paste0(family, ".cv.edgenet")
  ret
}

