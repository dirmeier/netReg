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
edgenet.loss <- function(gx, gy, family) {
  family <- family$family
  loss.function <- switch(
    family,
    "gaussian" = gaussian.loss,
    "binomial" = binomial.loss,
    "poisson" = poisson.loss,
    "Beta regression" = beta.loss, # mgcv::betar()
    "inverse.gaussian" = inverse.gaussian.loss,
    not.supported.yet(family)
  )

  loss <- function(alpha, beta, lambda, psigx, psigy, x, y) {
    eta <- linear.predictor(alpha, beta, x)
    obj <- loss.function(y, eta) + lasso(lambda, beta)

    if (!is.null(gx)) {
      obj <- obj + psigx * .edgenet.x.penalty(gx, beta)
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
