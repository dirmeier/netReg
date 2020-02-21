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


#' @title Fit a linear regression model the group lasso penalty
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
#'    to be used. Can a \code{\link{netReg::family}} function or a character string
#'    naming a family function, e.g. \code{gaussian} or "gaussian".}
#' \item{call }{ the call that produced the object}
#'
#' @examples
#' X <- matrix(rnorm(100 * 10), 100, 5)
#' b <- rnorm(5)
#' grps <- c(NA, 1, 1, 2, 2)
#'
#' # estimate the parameters of a Gaussian model
#' Y <- X %*% b + rnorm(100)
#' fit <- group.lasso(X = X, Y = Y, grps = grps, family = gaussian, maxit = 10)
#'
#' # estimation for binomial models
#' eta <- 1 / (1 + exp(-X %*% b))
#' Y <- rbinom(100, 1, eta)
#' fit <- group.lasso(X = X, Y = Y, grps = grps, family = binomial, maxit = 10)
#'
#' # estimation for Poisson models
#' eta <- exp(-X %*% b)
#' Y <- rpois(100, eta)
#' fit <- group.lasso(X = X, Y = Y, G.X = G.X, G.Y, family = poisson, maxit = 1)
setGeneric(
  "group.lasso",
  function(X, Y, grps = NULL,
           lambda = 1,
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
           lambda = 1,
           thresh = 1e-5, maxit = 1e5, learning.rate = 0.01,
           family = gaussian) {
    group.lasso(
      X, as.matrix(Y), G.X, G.Y,
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
           lambda = 1,
           thresh = 1e-5, maxit = 1e5, learning.rate = 0.01,
           family = gaussian) {
    stopifnot(
      is.numeric(maxit), is.numeric(thresh),
      is.numeric(learning.rate)
    )
    stopifnot(all(is.integer(grps)), max(grps) > ncol(X), min(grps) < 1)

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
    class(ret) <- c(class(ret), "edgenet")

    ret
  }
)


#' @noRd
.group.lasso <- function(x, y, grps,
                         lambda,
                         thresh, maxit, learning.rate, family) {
  p <- ncol(x)
  q <- ncol(y)

  reset_graph()

  x <- cast_float(x)
  y <- cast_float(y)

  if (!is.null(gx)) {
    gx <- cast_float(laplacian_(gx))
  }
  if (!is.null(gy)) {
    gy <- cast_float(laplacian_(gy))
  }

  # TODO: think about this
  alpha <- zero_vector(q) + 1
  beta <- zero_matrix(p, q) + 1

  # estimate coefficients
  loss <- edgenet.loss(gx, gy, family)
  objective <- loss(alpha, beta, lambda, psigx, psigy, x, y)
  res <- fit(objective, alpha, beta, maxit, learning.rate, thresh)

  # finalize output
  beta <- res$beta
  alpha <- res$alpha
  rownames(beta) <- colnames(x)
  colnames(beta) <- colnames(y)
  ret <- list(
    beta = beta,
    alpha = alpha,
    parameters = c("lambda" = lambda, "psigx" = psigx, "psigy" = psigy),
    lambda = lambda,
    psigx = psigx,
    psigy = psigy
  )

  ret$family <- family
  class(ret) <- paste0(family$family, ".edgenet")

  ret
}
