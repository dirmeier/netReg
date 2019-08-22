# netReg: graph-regularized linear regression models.
#
# Copyright (C) 2015 - 2019 Simon Dirmeier
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
summary.edgenet <- function(x,...)
{
    cat(sprintf("'%s' object", class(x)[1]))
    cat("\n\ncall:\n")
    print(x$call)
    cat("\nintercept:\n")
    print(x$alpha)
    cat("\ncoefficients:\n")
    print(x$beta)

}


#' @export
#' @method print cv.edgenet
summary.cv.edgenet <- function(x, ...)
{
    cat(sprintf("'%s' object", class(x)[1]))
    cat("\nOptimal parameters:\n")
    print(c(x$lambda, x$psigx, x$psigy))
}
