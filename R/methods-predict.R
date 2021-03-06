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
#' @method predict cv.edgenet
#' @importFrom stats predict
predict.cv.edgenet <- function(object, newdata = NULL, ...) {
  predict(object$fit, newdata, ...)
}


#' @export
#' @method predict cv.group.lasso
#' @importFrom stats predict
predict.cv.group.lasso <- function(object, newdata = NULL, ...) {
  predict(object$fit, newdata, ...)
}


#' @export
#' @method predict gaussian.edgenet
predict.gaussian.edgenet <- function(object, newdata = NULL, ...) {
  mean <- function(x) x
  .predict(object, newdata, mean, ...)
}


#' @export
#' @method predict gaussian.group.lasso
predict.gaussian.group.lasso <- function(object, newdata = NULL, ...) {
  predict.gaussian.edgenet(object, newdata, ...)
}


#' @export
#' @method predict binomial.edgenet
predict.binomial.edgenet <- function(object, newdata = NULL, ...) {
  mean <- function(x) 1 / (1 + base::exp(-x))
  .predict(object, newdata, mean, ...)
}


#' @export
#' @method predict binomial.group.lasso
predict.binomial.group.lasso <- function(object, newdata = NULL, ...) {
  predict.binomial.edgenet(object, newdata, ...)
}


#' @export
#' @method predict poisson.edgenet
predict.poisson.edgenet <- function(object, newdata = NULL, ...) {
  mean <- function(x) base::exp(x)
  .predict(object, newdata, mean, ...)
}


#' @export
#' @method predict poisson.group.lasso
predict.poisson.group.lasso <- function(object, newdata = NULL, ...) {
  predict.poisson.edgenet(object, newdata, ...)
}


#' @noRd
#' @importFrom stats coef
.predict <- function(object, newdata, mean, ...) {
  if (is.null(newdata)) stop("newdata is null")
  X <- newdata
  n <- dim(X)[1]
  p <- dim(X)[2]
  coefs <- object$beta
  if (p != dim(coefs)[1]) {
    stop("newdata dimensions do not fit coefficient dimensions!")
  }
  Y.hat <- mean(X %*% coefs + intercept.matrix(n, object$alpha))
  Y.hat
}
