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
