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


#' Find the optimal shrinkage parameters for edgenet
#'
#' @export
#' @docType methods
#' @rdname cvedgenet-methods
#'
#' @importFrom stats gaussian binomial
#'
#' @description Finds the optimal regulariztion parameters
#'  using cross-validation for edgenet. We use the BOBYQA algorithm to
#'  find the optimial regularization parameters in a cross-validation framework.
#'
#' @param X  input matrix, of dimension (\code{n} x \code{p})
#'  where \code{n} is the number of observations and \code{p} is the number
#'  of covariables. Each row is an observation vector.
#' @param Y  output matrix, of dimension (\code{n} x \code{q})
#'  where \code{n} is the number of observations and \code{q} is the number
#'  of response variables Each row is an observation vector.
#' @param G.X  non-negativ affinity matrix for \code{X}, of dimensions
#' (\code{p} x \code{p}) where \code{p} is the number of covariables.
#'  Providing a graph \code{G.X} will optimize the regularization
#'  parameter \code{psi.gx}. If this is not desired just set \code{G.X} to
#'  \code{NULL}.
#' @param G.Y  non-negativ affinity matrix for \code{Y}, of dimensions
#'  (\code{q} x \code{q}) where \code{q} is the number of responses \code{Y}.
#'  Providing a graph \code{G.Y} will optimize the regularization
#'  parameter \code{psi.gy}. If this is not desired just set \code{G.Y} to
#'  \code{NULL}.
#' @param lambda  \code{numerical} shrinkage parameter for LASSO. Per default
#' this parameter is
#'  set to \code{NA_real_} which means that \code{lambda} is going to be estimated
#'  using cross-validation. If any \code{numerical} value for \code{lambda}
#'  is set, estimation of the optimal parameter will \emph{not} be conducted.
#' @param psigx  \code{numerical} shrinkage parameter for graph-regularization
#'  of \code{G.X}. Per default this parameter is
#'  set to \code{NA_real_} which means that \code{psigx} is going to be estimated
#'  in the cross-validation. If any \code{numerical} value for \code{psigx} is
#'  set, estimation of the optimal parameter will \emph{not} be conducted.
#' @param psigy  \code{numerical} shrinkage parameter for graph-regularization
#'  of \code{G.Y}. Per default this parameter is
#'  set to \code{NA_real_} which means that \code{psigy} is going to be estimated
#'  in the cross-validation. If any \code{numerical} value for \code{psigy} is
#'  set, estimation of the optimal parameter will \emph{not} be conducted.
#' @param thresh  \code{numerical} threshold for the optimizer
#' @param maxit  maximum number of iterations for the optimizer
#'  (\code{integer})
#' @param learning.rate  step size for Adam optimizer (\code{numerical})
#' @param family  family of response, e.g. \emph{gaussian} or \emph{binomial}
#' @param optim.thresh  \code{numerical} threshold criterion for the
#'  optimization to stop.  Usually 1e-3 is a good choice.
#' @param optim.maxit  the maximum number of iterations for the optimization
#'   (\code{integer}). Usually 1e4 is a good choice.
#' @param nfolds  the number of folds to be used - default is 10
#'
#' @return An object of class \code{cv.edgenet}
#' \item{parameters }{ the estimated, optimal regularization parameters}
#' \item{lambda }{ optimal estimated value for regularization parameter lambda
#'   (or, if provided as argument, the value of the parameter)}
#' \item{psigx }{ optimal estimated value for regularization parameter psigx
#'   (or, if provided as argument, the value of the parameter)}
#' \item{psigy }{ optimal estimated value for regularization parameter psigy
#'   (or, if provided as argument, the value of the parameter)}
#' \item{estimated.parameters }{ names of parameters that were estimated}
#' \item{family }{ family used for estimated}
#' \item{fit }{ an \code{edgenet} object fitted with the optimal, estimated
#'  paramters}
#' \item{call }{ the call that produced the object}
#'
#' @examples
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' b <- matrix(rnorm(100), 10)
#' G.X <- abs(rWishart(1, 10, diag(10))[, , 1])
#' G.Y <- abs(rWishart(1, 10, diag(10))[, , 1])
#' diag(G.X) <- diag(G.Y) <- 0
#'
#' # estimate the parameters of a Gaussian model
#' Y <- X %*% b + matrix(rnorm(100 * 10), 100)
#'
#' ## dont use affinity matrices and estimate lambda
#' fit <- cv.edgenet(
#'   X = X, Y = Y, family = gaussian,
#'   maxit = 1, optim.maxit = 1
#' )
#' ## only provide one matrix and estimate lambda
#' fit <- cv.edgenet(
#'   X = X, Y = Y, G.X = G.X, psigx = 1, family = gaussian,
#'   maxit = 1, optim.maxit = 1
#' )
#' ## estimate only lambda with two matrices
#' fit <- cv.edgenet(
#'   X = X, Y = Y, G.X = G.X, G.Y, psigx = 1, psigy = 1,
#'   family = gaussian, maxit = 1, optim.maxit = 1
#' )
#' ## estimate only psigx
#' fit <- cv.edgenet(
#'   X = X, Y = Y, G.X = G.X, G.Y, lambda = 1, psigy = 1,
#'   family = gaussian, maxit = 1, optim.maxit = 1
#' )
#' ## estimate all parameters
#' fit <- cv.edgenet(
#'   X = X, Y = Y, G.X = G.X, G.Y,
#'   family = gaussian, maxit = 1, optim.maxit = 1
#' )
#' ## if Y is vectorial, we cannot use an affinity matrix for Y
#' fit <- cv.edgenet(
#'   X = X, Y = Y[, 1], G.X = G.X,
#'   family = gaussian, maxit = 1, optim.maxit = 1
#' )
setGeneric(
  "cv.edgenet",
  function(X, Y, G.X = NULL, G.Y = NULL,
           lambda = NA_real_, psigx = NA_real_, psigy = NA_real_,
           thresh = 1e-5, maxit = 1e5, learning.rate = 0.01,
           family = gaussian,
           optim.maxit = 1e2, optim.thresh = 1e-2,
           nfolds = 10) {
    standardGeneric("cv.edgenet")
  },
  package = "netReg"
)


#' @rdname cvedgenet-methods
setMethod(
  "cv.edgenet",
  signature = signature(X = "matrix", Y = "numeric"),
  function(X, Y, G.X = NULL, G.Y = NULL,
           lambda = NA_real_, psigx = NA_real_, psigy = NA_real_,
           thresh = 1e-5, maxit = 1e5, learning.rate = 0.01,
           family = gaussian,
           optim.maxit = 1e2, optim.thresh = 1e-2,
           nfolds = 10) {
    cv.edgenet(
      X, as.matrix(Y), G.X, G.Y,
      lambda, psigx, psigy,
      thresh, maxit, learning.rate,
      family,
      optim.maxit, optim.thresh,
      nfolds
    )
  }
)


#' @rdname cvedgenet-methods
setMethod(
  "cv.edgenet",
  signature = signature(X = "matrix", Y = "matrix"),
  function(X, Y, G.X = NULL, G.Y = NULL,
           lambda = NA_real_, psigx = NA_real_, psigy = NA_real_,
           thresh = 1e-5, maxit = 1e5, learning.rate = 0.01,
           family = gaussian,
           optim.maxit = 1e2, optim.thresh = 1e-2,
           nfolds = 10) {
    stopifnot(
      is.numeric(nfolds), nfolds > 0, is.numeric(learning.rate),
      is.numeric(optim.maxit), is.numeric(optim.thresh),
      is.numeric(maxit), is.numeric(thresh)
    )

    n <- dim(X)[1]
    p <- dim(X)[2]

    if (is.null(G.X)) psigx <- 0
    if (is.null(G.Y)) psigy <- 0

    check.matrices(X, Y)
    check.graphs(X, Y, G.X, G.Y, psigx, psigy)
    check.dimensions(X, Y, n, p)
    lambda <- check.param(lambda, 0, `<`, 0)
    psigx <- check.param(psigx, 0, `<`, 0)
    psigy <- check.param(psigy, 0, `<`, 0)
    maxit <- check.param(maxit, 0, `<`, 1e5)
    optim.maxit <- check.param(optim.maxit, 0, `<`, 1e2)
    thresh <- check.param(thresh, 0, `<`, 1e-5)
    optim.thresh <- check.param(optim.thresh, 0, `<`, 1e-2)
    family <- get.family(family)

    if (ncol(Y) == 1) {
      psigy <- 0
      G.Y <- NULL
    }

    if (n < nfolds) nfolds <- n
    folds <- sample(rep(seq_len(10), length.out = n))

    ret <- .cv.edgenet(
      X, Y, G.X, G.Y,
      lambda, psigx, psigy,
      family,
      thresh, maxit, learning.rate,
      nfolds, folds,
      optim.maxit, optim.thresh
    )

    ret$fit <- edgenet(
      X, Y, G.X, G.Y,
      ret$lambda, ret$psigx, ret$psigy,
      thresh, maxit, learning.rate, family
    )

    ret$call <- match.call()
    class(ret) <- c(class(ret), "cv.edgenet")

    ret
  }
)


#' @noRd
#' @import Rcpp
.cv.edgenet <- function(x, y, gx, gy,
                        lambda, psigx, psigy,
                        family,
                        thresh, maxit, learning.rate,
                        nfolds, folds,
                        optim.maxit, optim.thresh) {
  reg.params <- list(lambda = lambda, psigx = psigx, psigy = psigy)
  estimatable.params <- Filter(is.na, reg.params)
  fixed.params <- Filter(is.finite, reg.params)

  init.params <- rep(0.1, length(estimatable.params))
  if (is.null(gx) & is.null(gy) & !is.na(lambda)) {
    stop("you didn't set graphs and lambda != NA_real_.
             got nothing to estimate", call. = FALSE)
  }
  if (length(init.params) == 0) {
    stop("please set either of lambda/psigx/psigy to NA_real_",
      call. = FALSE
    )
  }

  p <- ncol(x)
  q <- ncol(y)

  if (!is.null(gx)) {
    gx <- cast_float(laplacian_(gx))
  }
  if (!is.null(gy)) {
    gy <- cast_float(laplacian_(gy))
  }


  lambda.tensor <- init_zero_scalar()
  psigx.tensor <- init_zero_scalar()
  psigy.tensor <- init_zero_scalar()

  loss <- edgenet.loss(lambda, psigx, psigy, gx, gy, family)

  optimizer <- keras::optimizer_adam(learning.rate)
  fn <- cross.validate(
    objective, train,
    x, y,
    x.tensor, y.tensor,
    lambda.tensor, psigx.tensor, psigy.tensor,
    nfolds, folds,
    maxit, thresh, learning.rate
  )

  with(session() %as% sess, {
    opt <- optim(fn, init.params,
      var.args = fixed.params,
      sess = sess, alpha = alpha, beta = beta,
      lower = rep(0, length(init.params)),
      upper = rep(100, length(init.params)),
      control = list(
        maxeval = optim.maxit,
        xtol_rel = optim.thresh,
        ftol_rel = optim.thresh,
        ftol_abs = optim.thresh
      )
    )
  })

  ret <- .cv.edgenet.post.process(opt, estimatable.params, fixed.params)
  ret$family <- family
  class(ret) <- paste0(family$family, ".cv.edgenet")

  ret
}


#' @noRd
.cv.edgenet.post.process <- function(opt, estimatable.params, fixed.params) {
  ret <- list(
    parameters = c(),
    "lambda" = NA_real_,
    "psigx" = NA_real_,
    "psigy" = NA_real_
  )

  for (i in seq(estimatable.params)) {
    nm <- names(estimatable.params)[i]
    ret[[nm]] <- opt$par[i]
    ret$estimated.parameters <- c(ret$estimated.parameters, nm)
  }
  for (i in seq(fixed.params)) {
    ret[[names(fixed.params)[i]]] <- as.double(fixed.params[i])
  }

  pars <- c("lambda" = ret$lambda, "psigx" = ret$psigx, "psigy" = ret$psigy)
  for (i in seq(pars))
  {
    if (names(pars)[i] %in% ret$estimated.parameters) {
      names(pars)[i] <- paste(names(pars)[i], "(estimated)")
    } else {
      names(pars)[i] <- paste(names(pars)[i], "(fixed)")
    }
  }

  ret$parameters <- pars
  ret
}
