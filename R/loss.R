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
#' @importFrom tensorflow tf
gaussian.loss <- function(y, eta, ...)
{
    obj <- tf$reduce_sum(tf$square(y - eta))
    obj
}


#' @noRd
#' @importFrom tensorflow tf
binomial.loss <- function(y, eta, ...)
{
    obj <- 0
    for (j in seq(ncol(y))) {
        prob <- tfp$distributions$Bernoulli(logits = eta[ ,j])
        obj <- obj + tf$reduce_sum(prob$log_prob(y[,j]))
    }

    -obj
}


#' @noRd
#' @importFrom tensorflow tf
poisson.loss <- function(y, eta, ...)
{
    obj <- 0
    for (j in seq(ncol(y))) {
        prob <- tfp$distributions$Poisson(log_rate = eta[ ,j])
        obj <- obj + tf$reduce_sum(prob$log_prob(y[,j]))
    }

    -obj
}


#' @noRd
#' @importFrom tensorflow tf
inverse.gaussian.loss <- function(y, eta, ...)
{
    obj <- 0

    for (j in seq(ncol(y))) {
        # loc := mu
        # concentration := lambda (shape)
        prob <- tfp$distributions$InverseGaussian(validate_args=TRUE,
            loc = 1 / sqrt(eta[ ,j]), concentration = 1)
        obj <- obj + tf$reduce_sum(prob$log_prob(y[,j]))
    }

    -obj
}


#' @noRd
#' @importFrom tensorflow tf
lasso <- function(lambda, beta)
{
    lambda * tf$reduce_sum(tf$abs(beta))
}


#' @noRd
#' @importFrom tensorflow tf
ridge<- function(lambda, beta)
{
    lambda * tf$reduce_sum(tf$square(beta))
}


#' @noRd
elastic <- function(alpha, lambda, beta)
{
    lambda * (ridge((1 - alpha) / 2, beta) +
              lasso(alpha, beta))
}
