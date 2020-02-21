# netReg: graph-regularized linear regression models.
#
# Copyright (C) 2015 - 2020 Simon Dirmeier
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
binomial.loss <- function(y, eta, invlink, ...)
{
    obj <- 0
    for (j in seq(ncol(y))) {
        prob <- tfp$distributions$Bernoulli(probs = invlink(eta[ ,j]))
        obj <- obj + tf$reduce_sum(prob$log_prob(y[,j]))
    }

    -obj
}


#' @noRd
#' @importFrom tensorflow tf
poisson.loss <- function(y, eta, invlink, ...)
{
    obj <- 0
    for (j in seq(ncol(y))) {
        prob <- tfp$distributions$Poisson(rate = invlink(eta[ ,j]))
        obj <- obj + tf$reduce_sum(prob$log_prob(y[,j]))
    }

    -obj
}


#' @noRd
#' @importFrom tensorflow tf
gamma.loss <- function(y, eta, invlink, ...)
{
    obj <- 0
    for (j in seq(ncol(y))) {
        prob <- tfp$distributions$Gamma(
            rate = invlink(tf$exp(eta[ ,j])), concentration = 1)
        obj <- obj + tf$reduce_sum(prob$log_prob(y[,j]))
    }

    -obj
}


#' @noRd
#' @importFrom tensorflow tf
beta.loss <- function(y, eta, invlink, ...)
{
    obj <- 0
    eps <- .Machine$double.eps * 1e9
    for (j in seq(ncol(y))) {
        mu <- invlink(eta[ ,j])
        phi <- 1 # TODO: replace this with tf$variable

        # reparametrize
        # concentration1 := alpha = mu * phi
        p <- mu * phi
        # concentration0 := beta = (1. - mu) * phi
        q <- (1 - mu) * phi

        # deal with numerical instabilities
        p.trans <- tf$math$maximum(p, eps)
        q.trans <- tf$math$maximum(q, eps)

        prob <- tfp$distributions$Beta(
            concentration1 = p.trans, concentration0 = q.trans)
        obj <- obj + tf$reduce_sum(prob$log_prob(y[,j]))
    }

    -obj
}


#' @noRd
#' @importFrom tensorflow tf
inverse.gaussian.loss <- function(y, eta, invlink, ...)
{
    obj <- 0
    for (j in seq(ncol(y))) {
        # loc := mu
        # concentration := lambda (shape)
        prob <- tfp$distributions$InverseGaussian(
            loc = invlink(tf$exp(eta[ ,j])), concentration = 1)
        obj <- obj + tf$reduce_sum(prob$log_prob(y[,j]))
    }

    -obj
}
