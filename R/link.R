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


#' @noRd
identity <- function(x) x


#' @noRd
#' @importFrom tensorflow tf
exp <- function(x) tf$maximum(tf$exp(x), tf$float32$min)


#' @noRd
#' @importFrom tensorflow tf
inverse <- function(x) 1 / x


#' @noRd
#' @importFrom tensorflow tf
logistic <- function(x) 1 / (1 + tf$exp(-x))


#' @noRd
gcdf <- function(x)
{
    std <- tfp$distributions$Normal(0, 1)
    std$cdf(x)
}

