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


#' Fit a graph-regularized linear regression model using
#'  edge-based regularization.
#'
#' @export
#' @docType methods
#' @rdname edgenet-methods
#'
#' @description  Fit a graph-regularized linear regression model using
#'  edge-penalization. The coefficients are computed using graph-prior
#'  knowledge in the form of one/two affinity matrices. Graph-regularization is
#'  an extension to previously introduced regularization techniques,
#'  such as the LASSO. See the vignette for details on the objective function of
#'  the model: \href{../doc/edgenet.html}{\code{vignette("edgenet", package="netReg")}}
#'
#' @param X  input matrix, of dimension (\code{n} x \code{p})
#' where \code{n} is the number of observations and \code{p} is the number
#' of covariables. Each row is an observation vector.
#' @param Y  output matrix, of dimension (\code{n} x \code{q})
#' where \code{n} is the number of observations and \code{q} is the number
#' of response variables. Each row is an observation vector.
#' @param G.X  non-negativ affinity matrix for \code{n}, of dimensions
#' (\code{p} x \code{p}) where \code{p} is the number of covariables \code{X}
#' @param G.Y  non-negativ affinity matrix for \code{n}, of dimensions
#' (\code{q} x \code{q}) where \code{q} is the number of responses \code{Y}
#' @param lambda  \code{numerical} shrinkage parameter for LASSO.
#' @param psigx  \code{numerical} shrinkage parameter for graph-regularization
#'  of \code{G.X}
#' @param psigy  \code{numerical} shrinkage parameter for graph-regularization
#'  of \code{G.Y}
#' @param thresh  \code{numerical} threshold for optimizer
#' @param maxit  maximum number of iterations for optimizer
#'  (\code{integer})
#' @param family  family of response, e.g. \emph{gaussian} or \emph{binomial}
#'
#' @return An object of class \code{edgenet}
#' \item{coefficients }{ the estimated (\code{p} x \code{q})-dimensional
#'  coefficient matrix B.hat}
#' \item{intercept }{ the estimated (\code{q} x \code{1})-dimensional
#'  vector of intercepts}
#' \item{call }{ the call that produced the object}
#' \item{family }{ the family of the response}
#'
#' @references
#'  Dirmeier, Simon and Fuchs, Christiane and Mueller, Nikola S and Theis,
#'  Fabian J (2018),
#'  netReg: Network-regularized linear models for biological association
#'  studies. \cr
#'  \emph{Bioinformatics}\cr \cr
#'  Friedman J., Hastie T., Hoefling H. and Tibshirani R. (2007),
#'  Pathwise coordinate optimization.\cr
#'  \emph{The Annals of Applied Statistics}\cr \cr
#'  Friedman J., Hastie T. and Tibshirani R. (2010),
#'  Regularization Paths for Generalized Linear Models via
#'  Coordinate Descent.\cr \emph{Journal of Statistical Software}\cr \cr
#'  Fu W. J. (1998),  Penalized Regression: The Bridge Versus the Lasso.\cr
#'  \emph{Journal of Computational and Graphical Statistics}\cr \cr
#'  Cheng W. and Wang W. (2014), Graph-regularized dual Lasso for robust
#'  eQTL mapping.\cr
#'  \emph{Bioinformatics}
#'
#' @examples
#' X <- matrix(rnorm(100*10), 100, 10)
#' b <- rnorm(10)
#' G.X <- matrix(rpois(100,1), 10)
#' G.X <- t(G.X) + G.X
#' diag(G.X) <- 0
#'
#' # fit a Gaussian model
#' Y <- X %*% b + rnorm(100)
#' fit <- edgenet(X=X, Y=Y, G.X=G.X, family=gaussian)
setGeneric(
    "edgenet",
    function(X, Y, G.X=NULL, G.Y=NULL,
             lambda=1, psigx=1, psigy=1,
             thresh=1e-5, maxit=1e5, learning.rate=0.01,
             family=gaussian)
    {
        standardGeneric("edgenet")
    },
    package = "netReg"
)


#' @noRd
setMethod(
    "edgenet",
    signature = signature(X="matrix", Y="numeric"),
    function(X, Y, G.X=NULL, G.Y=NULL,
             lambda=1, psigx=1, psigy=1,
             thresh=1e-5, maxit=1e5, learning.rate=0.01,
             family=gaussian)
    {
        edgenet(X, as.matrix(Y), G.X, G.Y,
                lambda, psigx, psigy,
                thresh, maxit, learning.rate,
                family)
    }
)

#' @noRd
setMethod(
    "edgenet",
    signature = signature(X="matrix", Y="matrix"),
    function(X, Y, G.X=NULL, G.Y=NULL,
             lambda=1, psigx=1, psigy=1,
             thresh=1e-5, maxit=1e5, learning.rate=0.01,
             family=gaussian)
    {
        if (is.null(G.X)) psigx <- 0
        if (is.null(G.Y)) psigy <- 0

        check.matrices(X, Y)
        check.graphs(X, Y, G.X, G.Y, psigx, psigy)
        check.dimensions(X, Y, nrow(X), ncol(X))
        lambda <- check.param(lambda, 0, `<`, 0)
        psigx <- check.param(psigx, 0, `<`, 0)
        psigy <- check.param(psigy, 0, `<`, 0)
        maxit <- check.param(maxit, 0, `<`, 1e5)
        thresh <- check.param(thresh, 0, `<`, 1e-5)
        family <- get.family(family)

        if (ncol(Y) == 1) {
            psigy <- 0
            G.Y <- NULL
        }

        # estimate coefficients
        ret <- .edgenet(x=X, y=Y, gx=G.X, gy=G.Y,
                        lambda = lambda, psigx = psigx, psigy = psigy,
                        thresh = thresh, maxit = maxit,
                        learning.rate = learning.rate, family = family)

        ret$call   <- match.call()
        class(ret) <- c(class(ret), "edgenet")

        ret
    }
)


#' @noRd
#' @import Rcpp
.edgenet <- function(x, y, gx, gy,
                     lambda, psigx, psigy,
                     thresh, maxit,learning.rate, family)
{
    p <- ncol(x)
    q <- ncol(y)

    x <- cast_float(x)
    y <- cast_float(y)

    if (!is.null(gx))
        gx <- cast_float(netReg:::laplacian_(gx))
    if (!is.null(gy))
        gy <- cast_float(netReg:::laplacian_(gy))

    alpha <- zero_vector(q)
    beta  <- zero_matrix(p, q)

     #estimate coefficients
    loss  <- edgenet.loss(gx, gy, family, q)
    objective <- loss(alpha, beta, lambda, psigx, psigy, x, y)
    res <- fit(objective, alpha, beta, maxit, learning.rate, thresh)

    # finalize output
    beta  <- res$beta
    alpha <- res$alpha
    rownames(beta) <- colnames(x)
    colnames(beta) <- colnames(y)
    ret <- list(beta = beta,
                alpha = alpha,
                lambda = lambda,
                psigx = psigx,
                psigy = psigy)

    ret$family <- family$family
    class(ret) <- paste0(family$family, ".edgenet")

    ret
}
