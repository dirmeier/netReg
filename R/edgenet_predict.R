# netReg: graph-regularized linear regression models.
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
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


#' Predict method for gaussian edgenet fits
#' 
#' @export
#' 
#' @description Predicts the estimated Y.hat values for a newdata 
#'  design matrix X similar to the other predict methods, 
#'  e.g. from glm and glmnet
#' 
#' @importFrom stats coef
#' 
#' @param object a fitted object of class \emph{gaussian.edgenet}
#' @param ... further arguments
#' @param newdata a new (\code{m} x \code{p})-dimensional design matrix with a 
#' variable number of observations \code{m}, but a constant number 
#' of co-variables \code{p}
#' 
#' @return A (\code{m} x \code{q})-dimensional matrix
#' 
#' @method predict gaussian.edgenet
#' 
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100*10),100,10)
#' G.X <- matrix(rpois(10*10,1),10)
#' G.X <- t(G.X) + G.X
#' diag(G.X) <- 0
#' 
#' Y <- matrix(rnorm(100*10),100,10)
#' fit <- edgenet(X=X, Y=Y, G.X=G.X, family="gaussian")
#' pred <- predict(fit, X)
#' }
predict.gaussian.edgenet <- function(object, newdata=NULL, ...)
{
  if(is.null(newdata))
    stop("newdata is null")
  X <- newdata
  n <- dim(X)[1]
  p <- dim(X)[2]
  coefs <- stats::coef(object)
  if(p != dim(coefs)[1])
    stop("newdata dimensions do not fit coefficient dimensions!")
  mu <- object$intercept
  Y.hat <- X %*% coefs + intercept.matrix(n=n, mu=mu)
  Y.hat
}