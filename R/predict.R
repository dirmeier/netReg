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
predict.cv.edgenet <- function(object, newdata=NULL, ...)
{
    predict(object$fit, newdata, ...)
}


#' @export
#' @method predict gaussian.edgenet
predict.gaussian.edgenet <- function(object, newdata=NULL, ...)
{
    mean <- function(x) x
   .predict(object, newdata, mean, ...)
}


#' @export
#' @method predict binomial.edgenet
predict.binomial.edgenet <- function(object, newdata=NULL, ...)
{
    mean <- function(x) 1 / (1 + exp(-x))
    .predict(object, newdata, mean, ...)
}


#' @export
#' @method predict poisson.edgenet
predict.poisson.edgenet <- function(object, newdata=NULL, ...)
{
    mean <- function(x) exp(x)
    .predict(object, newdata, mean, ...)
}


#' @export
#' @method predict beta.edgenet
predict.beta.edgenet <- function(object, newdata=NULL, ...)
{
    mean <- function(x) 1 / (1 + exp(-x))
    .predict(object, newdata, mean, ...)
}


#' @export
#' @method predict inverse.gaussian.edgenet
predict.inverse.gaussian.edgenet <- function(object, newdata=NULL, ...)
{
    mean <- function(x) 1 / sqrt(x)
    .predict(object, newdata, mean, ...)
}


#' @noRd
#' @importFrom stats coef
.predict <- function(object, newdata, mean, ...)
{
    if (is.null(newdata)) stop("newdata is null")
    X <- newdata
    n <- dim(X)[1]
    p <- dim(X)[2]
    coefs <- object$beta
    if(p != dim(coefs)[1])
        stop("newdata dimensions do not fit coefficient dimensions!")
    Y.hat <- mean(X %*% coefs + intercept.matrix(n, object$alpha))
    Y.hat
}
