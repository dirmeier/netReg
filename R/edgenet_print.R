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


#' @export
#' @method print edgenet
print.edgenet <- function(x,...)
{
    cat("\nCall: ")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients)
    cat("\nIntercept:\n")
    print(x$intercept)
    cat("\nParameters:\n")
    print(c(lambda=x$lambda, psi_gx=x$psigx, psi_gy=x$psigy))
    cat("\nFamily:\n")
    print(x$family)
}

#' @export
#' @method print cv.edgenet
print.cv.edgenet <- function(x, ...)
{
    cat("\nCall: ")
    print(x$call)
    cat("\nParameters:\n")
    print(c(x$lambda, x$psigx, x$psigy))
    cat("Family:\n")
    print(x$family)
}
