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
print.edgenet <- function(x, ...) {
  cat(sprintf("'%s' object", class(x)[1]))
}

#' @export
#' @method print group.lasso
print.group.lasso <- function(x, ...) {
  print.edgenet(x, ...)
}


#' @export
#' @method print cv.edgenet
print.cv.edgenet <- function(x, ...) {
  cat(sprintf("'%s' object", class(x)[1]))
  cat("\n\noptimal parameters:\n")
  print(x$parameters)
}


#' @export
#' @method print cv.group.lasso
print.cv.group.lasso <- function(x, ...) {
  print.cv.edgenet(x, ...)
}


#' @export
#' @method print netReg.family
print.netReg.family <- function(x, ...) {
  cat(sprintf("family: %s\n", x$family))
  cat(sprintf("link: %s\n", x$link))
}
