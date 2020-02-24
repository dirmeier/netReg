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
#' @method coef edgenet
coef.edgenet <- function(object, ...)
{
    alpha <- object$alpha
    beta  <- object$beta
    coefs <- rbind(alpha, beta)
    rownames(coefs) <- c("(Intercept)", sprintf("x[%s]", seq(nrow(beta))))
    colnames(coefs) <- sprintf("y[%s]", seq(ncol(beta)))
    coefs
}


#' @export
#' @method coef cv.edgenet
coef.cv.edgenet <- function(object, ...)
{
   coef(object$fit)
}


#' @export
#' @method coef group.lasso
coef.group.lasso <- function(object, ...)
{
    coef.edgenet(object, ...)
}


#' @export
#' @method coef cv.group.lasso
coef.cv.group.lasso <- function(object, ...)
{
    coef(object$fit)
}
