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


#' Find the optimal shrinkage parameters for edgenet
#'
#' @export
#'
#' @description Finds the optimal shrinkage parameters
#'  using cross-validation for edgenet. We use the BOBYQA algorithm to
#'  find the optimial regularization parameters and coordinate
#'  descent in order to minimize the objective function of the linear model.
#'
#' @param X  input matrix, of dimension (\code{n} x \code{p})
#'  where \code{n} is the number of observations and \code{p} is the number
#'  of covariables. Each row is an observation vector.
#' @param Y  output matrix, of dimension (\code{n} x \code{q})
#'  where \code{n} is the number of observations and \code{q} is the number
#'  of response variables Each row is an observation vector.
#' @param G.X  non-negativ affinity matrix for \code{n}, of dimensions
#'  (\code{p} x \code{p}) where \code{p} is the number of covariables \code{X}.
#'  Providing a graph \code{G.X} will optimize the regularization
#'  parameter \code{psi.gx}. If this is not desired just set \code{G.X} to
#'  \code{NULL}.
#' @param G.Y  non-negativ affinity matrix for \code{n}, of dimensions
#'  (\code{q} x \code{q}) where \code{q} is the number of responses \code{Y}.
#'  Providing a graph \code{G.Y} will optimize the regularization
#'  parameter \code{psi.gy}. If this is not desired just set \code{G.Y} to
#'  \code{NULL}.
#' @param lambda  \code{numerical} shrinkage parameter for LASSO. Per default
#' this parameter is
#'  set to \code{NULL} which means that \code{lambda} is going to be estimated
#'  using cross-validation. If any \code{numerical} value for \code{lambda}
#'  is set, estimation of the optimal parameter will \emph{not} be conducted.
#' @param psigx  \code{numerical} shrinkage parameter for graph-regularization
#'  of \code{G.X}. Per default this parameter is
#'  set to \code{NULL} which means that \code{psigx} is going to be estimated
#'  in the cross-validation. If any \code{numerical} value for \code{psigx} is
#'  set, estimation of the optimal parameter will \emph{not} be conducted.
#' @param psigy  \code{numerical} shrinkage parameter for graph-regularization
#'  of \code{G.Y}. Per default this parameter is
#'  set to \code{NULL} which means that \code{psigy} is going to be estimated
#'  in the cross-validation. If any \code{numerical} value for \code{psigy} is
#'  set, estimation of the optimal parameter will \emph{not} be conducted.
#' @param thresh  \code{numerical} threshold for coordinate descent
#' @param maxit  maximum number of iterations for the coordinate descent
#'  (\code{integer})
#' @param family  family of response, e.g. \emph{gaussian}
#' @param optim.epsilon  \code{numerical} threshold criterion for the
#'  optimization to stop.  Usually 1e-3 is a good choice.
#' @param optim.maxit  the maximum number of iterations for the optimization
#'   (\code{integer}). Usually 1e4 is a good choice.
#' @param nfolds  the number of folds to be used - default is 10
#'  (minimum 3, maximum \code{nrow(X)}).
#'
#' @return An object of class \code{cv.edgenet}
#' \item{call }{ the call that produced the object}
#' \item{lambda }{ the estimated (\code{p} x \code{q})-dimensional
#'  coefficient matrix B.hat}
#' \item{psigx }{ the estimated (\code{q} x \code{1})-dimensional
#'  vector of intercepts}
#' \item{psigy }{ the estimated (\code{q} x \code{1})-dimensional vector
#'  of intercepts}
#'
#' @references
#'  Dirmeier, Simon and Fuchs, Christiane and Mueller, Nikola S and Theis,
#'  Fabian J (2018),
#'  netReg: Network-regularized linear models for biological association
#'  studies. \cr
#'  Friedman J., Hastie T., Hoefling H. and Tibshirani R. (2007),
#'  Pathwise coordinate optimization.\cr
#'  \emph{The Annals of Applied Statistics}\cr \cr
#'  Friedman J., Hastie T. and Tibshirani R. (2010),
#'  Regularization Paths for Generalized Linear Models via
#'   Coordinate Descent. \cr
#'  \emph{Journal of Statistical Software}\cr \cr
#'  Fu W. J. (1998),  Penalized Regression: The Bridge Versus the Lasso.\cr
#'  \emph{Journal of Computational and Graphical Statistics}\cr \cr
#'  Cheng W. and Wang W. (2014), Graph-regularized dual Lasso for
#'   robust eQTL mapping.\cr
#'  \emph{Bioinformatics}\cr \cr
#'  Powell M.J.D. (2009),
#'  The BOBYQA algorithm for bound constrained optimization without
#'   derivatives.\cr
#'  \url{http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf}
#'
#' @examples
#' X <- matrix(rnorm(100*10), 100, 10)
#' b <- rnorm(10)
#' G.X <- matrix(rpois(10*10,1),10)
#' G.X <- t(G.X) + G.X
#' diag(G.X) <- 0
#'
#' # fit a Gaussian model
#' Y <- X%*%b + rnorm(100)
#' cv.edge <- cv.edgenet(X=X, Y=Y, G.X=G.X, family="gaussian")
cv.edgenet <- function(X, Y, G.X=NULL, G.Y=NULL,
                       lambda=NULL, psigx=NULL, psigy=NULL,
                       thresh=1e-5, maxit=1e5,
                       family=c("gaussian"),
                       optim.epsilon=1e-3,
                       optim.maxit=1e4,
                       nfolds=10)
{
    UseMethod("cv.edgenet")
}

#' @export
#' @method cv.edgenet default
cv.edgenet.default <- function(X, Y, G.X=NULL, G.Y=NULL,
                               lambda=NULL, psigx=NULL, psigy=NULL,
                               thresh=1e-5, maxit=1e5,
                               family=c("gaussian"),
                               optim.epsilon=1e-3,
                               optim.maxit=1e4,
                               nfolds=10)
{
    stopifnot(is.numeric(nfolds), nfolds > 0,
              is.numeric(optim.epsilon), is.numeric(10),
              is.numeric(maxit), is.numeric(thresh))
    check.matrices(X, Y)
    n <- dim(X)[1]
    p <- dim(X)[2]
    q <- dim(Y)[2]

    do.lambda <- do.psigx <- do.psigy <- TRUE
    if (is.positive.numeric(lambda)) do.lambda <- FALSE
    if (is.positive.numeric(psigx))  do.psigx  <- FALSE
    if (is.positive.numeric(psigy))  do.psigy  <- FALSE

    if (is.null(G.X)) G.X <- matrix(0, 1, 1)
    if (is.null(G.Y)) G.Y <- matrix(0, 1, 1)
    if (all(G.X == 0)) {
        psigx <- 0
        do.psigx <- FALSE
    }
    if (all(G.Y == 0)) {
        psigy <- 0
        do.psigy <- FALSE
    }

    check.graphs(X, Y, G.X, G.Y, psigx, psigy)
    check.dimensions(X, Y, n, p)
    if (maxit < 0) {
        warning("maxit < 0, setting to 1e5!")
        maxit <- 1e5
    }
    if (thresh < 0) {
        warning("thresh < 0, setting to 1e-5!")
        thresh <- 1e-5
    }
    if (optim.epsilon < 0) {
        warning("epsilon < 0; settint to 1e-3")
        optim.epsilon <- 1e-3
    }
    if (optim.maxit < 0) {
        warning("approx.maxit < 0; settint to 1e4")
        optim.maxit <- 1e4
    }

    foldid <- NULL
    # check if some parameters have values
    if (!is.null(foldid) & is.numeric(foldid)) {
        nfolds <- max(foldid)
        stopifnot(length(foldid) == n)
    }
    if (is.null(foldid)) foldid <- NA_integer_
    if (!is.numeric(foldid))
        stop("Please provide either an integer vector or NULL for foldid")
    if (q == 1)     psigy  <- 0
    if (n < nfolds) nfolds <- n

    # cast nulls to doubles to avoid errors
    if (is.null(lambda)) lambda <- 0
    if (is.null(psigx))  psigx <- 0
    if (is.null(psigy))  psigy <- 0

    # set static to avoid memory overload
    if (n >= 1000 && p >= 500) nfolds <- 5
    family <- match.arg(family)

    # estimate shrinkage parameters
    ret <- .cv.edgenet(
        X=X, Y=Y, G.X=G.X, G.Y=G.Y,
        lambda=lambda, psigx=psigx, psigy=psigy,
        do.lambda=do.lambda, do.psigx=do.psigx, do.psigy=do.psigy,
        family=family, thresh=thresh, maxit=maxit,
        nfolds=nfolds, foldid=foldid,
        optim.epsilon=optim.epsilon, optim.maxit=optim.maxit)

    ret$call   <- match.call()
    class(ret) <- c(class(ret), "cv.edgenet")

    ret
}

#' @noRd
#' @import Rcpp
.cv.edgenet <- function(X, Y, G.X, G.Y,
                        lambda, psigx, psigy,
                        do.lambda, do.psigx, do.psigy,
                        family ,thresh, maxit,
                        nfolds, foldid, optim.maxit, optim.epsilon)
{
    cv <- .Call("cv_edgenet_cpp",
                X, Y, G.X, G.Y,
                as.double(lambda),
                as.double(psigx),
                as.double(psigy),
                as.logical(do.lambda),
                as.logical(do.psigx),
                as.logical(do.psigy),
                as.integer(maxit),
                as.double(thresh),
                as.integer(nfolds),
                as.integer(foldid),
                as.integer(length(foldid)),
                as.character(family),
                as.integer(optim.maxit),
                as.double(optim.epsilon))

    ret <- list(lambda=cv$parameters[1],
                psigx =cv$parameters[2],
                psigy =cv$parameters[3],
                folds =cv$folds+1)
    ret$family <- family
    class(ret) <- paste0(family, ".cv.edgenet")

    ret
}
