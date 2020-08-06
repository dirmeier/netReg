# netReg: network-regularized linear regression models.
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


#' @title Fit a linear regression model using a group lasso penalty
#'
#' @export
#' @docType methods
#' @rdname grouplasso-methods
#'
#' @importFrom stats gaussian binomial
#'
#' @description  Fit a linear regression model the group LASSO penalty.
#'
#' @param X  input matrix, of dimension (\code{n} x \code{p})
#' where \code{n} is the number of observations and \code{p} is the number
#' of covariables. Each row is an observation vector.
#' @param Y  output matrix, of dimension (\code{n} x \code{q})
#' where \code{n} is the number of observations and \code{q} is the number
#' of response variables. Each row is an observation vector.
#' @param grps  vector of integers or \code{NA_integer_} of length \code{p}
#'  that encodes the grouping of variables, e.g., \code{c(1, 1, 2, 2, NA)}
#' @param lambda  \code{numerical} shrinkage parameter
#' @param thresh  \code{numerical} threshold for optimizer
#' @param maxit  maximum number of iterations for optimizer
#'  (\code{integer})
#' @param learning.rate   step size for Adam optimizer (\code{numerical})
#' @param family  family of response, e.g., \emph{gaussian} or \emph{binomial}
#'
#' @return An object of class \code{edgenet}
#' \item{beta }{ the estimated (\code{p} x \code{q})-dimensional
#'  coefficient matrix B.hat}
#' \item{alpha }{ the estimated (\code{q} x \code{1})-dimensional
#'  vector of intercepts}
#' \item{parameters }{ regularization parameters}
#' \item{lambda }{ regularization parameter lambda)}
#' \item{family }{ a description of the error distribution and link function
#'    to be used. Can be a \code{\link[netReg:family]{netReg::family}} function or a character string
#'    naming a family function, e.g. \code{gaussian} or "gaussian".}
#' \item{call }{ the call that produced the object}
#'
#' @examples
#' X <- matrix(rnorm(100 * 10), 100, 5)
#' b <- rnorm(5)
#' grps <- c(NA_integer_, 1L, 1L, 2L, 2L)
#'
#' # estimate the parameters of a Gaussian model
#' Y <- X %*% b + rnorm(100)
#' fit <- group.lasso(X = X, Y = Y, grps = grps, family = gaussian, maxit = 10)
#' @references
#'  Yuan, Ming and Lin, Yi (2006),
#'  Model selection and estimation in regression with grouped variables. \cr
#'  \emph{Journal of the Royal Statistical Society: Series B}\cr \cr
#'  Meier, Lukas and Van De Geer, Sara and Bühlmann, Peter (2008),
#'  The group lasso for logistic regression. \cr
#'  \emph{Journal of the Royal Statistical Society: Series B}
#'
setGeneric(
  "group.lasso",
  function(X, Y, grps = NULL,
           lambda = 0,
           thresh = 1e-5, maxit = 1e5, learning.rate = 0.01,
           family = gaussian) {
    standardGeneric("group.lasso")
  },
  package = "netReg"
)


#' @rdname grouplasso-methods
setMethod(
  "group.lasso",
  signature = signature(X = "matrix", Y = "numeric"),
  function(X, Y, grps = NULL,
           lambda = 0,
           thresh = 1e-5, maxit = 1e5, learning.rate = 0.01,
           family = gaussian) {
    group.lasso(
      X, as.matrix(Y), grps,
      lambda,
      thresh, maxit, learning.rate,
      family
    )
  }
)


#' @rdname grouplasso-methods
setMethod(
  "group.lasso",
  signature = signature(X = "matrix", Y = "matrix"),
  function(X, Y, grps = NULL,
           lambda = 0,
           thresh = 1e-5, maxit = 1e5, learning.rate = 0.01,
           family = gaussian) {
    stopifnot(
      is.numeric(maxit), is.numeric(thresh),
      is.numeric(learning.rate)
    )

    if (is.null(grps)) {
      grps <- rep(NA_integer_, ncol(X))
    }
    stopifnot(
      all(is.integer(grps)),
      max(grps, na.rm = TRUE) <= ncol(X),
      min(grps, na.rm = TRUE) >= 1
    )

    check.matrices(X, Y)
    check.dimensions(X, Y, nrow(X), ncol(X))
    lambda <- check.param(lambda, 0, `<`, 0)
    maxit <- check.param(maxit, 0, `<`, 1e5)
    thresh <- check.param(thresh, 0, `<`, 1e-5)
    family <- get.family(family)

    # estimate coefficients
    ret <- .group.lasso(
      x = X, y = Y, grps = grps,
      lambda = lambda,
      thresh = thresh, maxit = maxit,
      learning.rate = learning.rate, family = family
    )

    ret$call <- match.call()
    class(ret) <- c(class(ret), "group.lasso")

    ret
  }
)


#' @noRd
.group.lasso <- function(x, y, grps,
                         lambda, thresh, maxit, learning.rate, family) {
  cols.x <- colnames(x)
  cols.y <- colnames(y)
  x <- cast_float(x)
  y <- cast_float(y)

  mod <- model(ncol(x), ncol(y), family)
  loss <- group.lasso.loss(lambda, grps, family)
  res <- fit(mod, loss, x, y, maxit, learning.rate, thresh)

  # finalize output
  beta <- res$beta
  alpha <- res$alpha
  rownames(beta) <- cols.x
  colnames(beta) <- cols.y
  ret <- list(
    beta = beta,
    alpha = alpha,
    parameters = c("lambda" = lambda),
    lambda = lambda
  )

  ret$family <- family
  class(ret) <- paste0(family$family, ".group.lasso")

  ret
}
