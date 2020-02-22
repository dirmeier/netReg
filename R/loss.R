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
#' @import tensorflow
edgenet.loss <- function(gx, gy, family)
{
    invlink <- family$linkinv
    loss.function <- family$loss

    loss <- function(alpha, beta, lambda, psigx, psigy, x, y)
    {
        eta <- linear.predictor(alpha, beta, x)
        obj <- loss.function(y, eta, invlink) + .lasso.penalty(lambda, beta)

        if (!is.null(gx)) {
            obj <- obj + psigx * .edgenet.x.penalty(gx, beta)
        }
        if (!is.null(gy)) {
            obj <- obj + psigy * .edgenet.y.penalty(gy, beta)
        }

        obj
    }

    loss
}


#' @noRd
#' @import tensorflow
.edgenet.x.penalty <- function(gx, beta)
{
    tf$linalg$trace(tf$matmul(tf$transpose(beta), tf$matmul(gx, beta)))
}


#' @noRd
#' @import tensorflow
.edgenet.y.penalty <- function(gy, beta)
{
    tf$linalg$trace(tf$matmul(beta, tf$matmul(gy, tf$transpose(beta))))
}


#' @noRd
#' @import tensorflow
group.lasso.loss <- function(grps, family)
{
    invlink <- family$linkinv
    loss.function <- family$loss

    loss <- function(alpha, beta, lambda, x, y)
    {
        eta <- linear.predictor(alpha, beta, x)
        obj <- loss.function(y, eta, invlink) +
            .group.lasso.penalty(lambda, beta, grps)

        obj
    }

    loss
}


#' @noRd
#' @importFrom tensorflow tf
.group.lasso.penalty <- function(lambda, beta, grps) {
    pen <- 0
    iter <- unique(grps[!is.na(grps)])
    for (el in iter) {
        idxs <- which(grps == el)
        grp.pen <-  tf$sqrt(cast_float(length(idxs)))
        for (j in seq(ncol(beta))) {
            pen <- pen + grp.pen * tf$math$reduce_euclidean_norm(beta[idxs, j])
        }
    }
    lambda * pen
}


#' @noRd
#' @importFrom tensorflow tf
.lasso.penalty <- function(lambda, beta)
{
    lambda * tf$reduce_sum(tf$abs(beta))
}


#' @noRd
#' @importFrom tensorflow tf
.ridge.penalty <- function(lambda, beta)
{
    lambda * tf$reduce_sum(tf$square(beta))
}


#' @noRd
.elastic.penalty <- function(alpha, lambda, beta)
{
    lambda * (.ridge.penalty((1 - alpha) / 2, beta) +
              .lasso.penalty(alpha, beta))
}

