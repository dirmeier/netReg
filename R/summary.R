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



.summary.edgenet <- function(x, obj.name)
{
    cat(sprintf("'%s' object", class(x)[1]))
    cat("\n\ncall:\n")
    print(x$call)
    prmstr <- "parameters"
    if (class(x)[-1] == "cv.edgenet")
        prmstr <- paste("optimal", prmstr)
    cat(sprintf("\n%s:\n", prmstr))
    print(x$parameters)
    cat("\nfamily: ", x$family$family)
    cat("\nlink: ", x$family$link, "\n")
    cat(sprintf("\n-> call coef(%s) for coefficients", obj.name))
}


#' @export
#' @method summary edgenet
summary.edgenet <- function(object, ...)
    .summary.edgenet(object, deparse(substitute(object)))



#' @export
#' @method summary cv.edgenet
summary.cv.edgenet <- function(object, ...)
    .summary.edgenet(object, deparse(substitute(object)))


#' @export
#' @method summary netReg.family
summary.netReg.family <- function(object, ...) print(object)

