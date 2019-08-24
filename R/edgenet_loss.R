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
#' @import tensorflow
edgenet.loss <- function(gx, gy, family, q)
{
    family <- family$family
    loss.function <- switch(
        family,
        "gaussian" = .edgenet.gaussian.loss,
        "binomial" = .edgenet.binomial.loss,
        not.supported.yet(family))

    loss <- function(alpha, beta, lambda, psigx, psigy, x, y)
    {
        eta <- tf$matmul(x, beta) + tf$ones(shape(q, 1), tf$float32) * tf$transpose(alpha)
        obj <- loss.function(y, eta, ncol=q) + .lasso(lambda, beta)

        if (!is.null(gx)) {
            x.penalty <- tf$trace(tf$matmul(tf$transpose(beta), tf$matmul(gx, beta)))
            obj <- obj + psigx * x.penalty
        }
        if (!is.null(gy)) {
            y.penalty <- tf$trace(tf$matmul(beta, tf$matmul(gy, tf$transpose(beta))))
            obj <- obj + psigy * y.penalty
        }

        obj
    }

    loss
}


#' @noRd
#' @import tensorflow
.edgenet.gaussian.loss <- function(y, mean, ...)
{
    obj <- tf$reduce_sum(tf$square(y - mean))
    obj
}


#' @noRd
#' @import tensorflow
.edgenet.binomial.loss <- function(y, means, ncol, ...)
{
    obj <- 0
    for (j in seq(1, ncol)) {
        prob <- tfp$distributions$Bernoulli(logits = means[,j])
        obj <- obj + tf$reduce_sum(prob$log_prob(y[,j]))
    }

    -obj
}


#' @noRd
#' @import tensorflow
.lasso <- function(lambda, beta) {
    lambda * tf$reduce_sum(tf$abs(beta))
}
