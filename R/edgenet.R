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


#' Fit a graph-regularized linear regression model using
#'  edge-based regularization.
#'
#' @export
#'
#' @description  Fit a graph-regularized linear regression model using
#'  edge-penalization. The coefficients are computed using graph-prior
#'  knowledge in the form of one/two affinity matrices. Graph-regularization is
#'  an extension to previously introduced regularization techniques,
#'  such as the LASSO. For that reason we are also using coordinate descent
#'  for minimization of the objective function of the linear model.
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
#' @param thresh  \code{numerical} threshold for coordinate descent
#' @param maxit  maximum number of iterations for coordinate descent
#'  (\code{integer})
#' @param family  family of response, e.g. \emph{gaussian}
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
#' Y <- X%*%b + rnorm(100)
#' fit <- edgenet(X=X, Y=Y, G.X=G.X, family="gaussian")

edgenet <- function(X, Y, G.X=NULL, G.Y=NULL,
                    lambda=1, psigx=1, psigy=1,
                    thresh=1e-5, maxit=1e5,
                    family=gaussian)
{
    UseMethod("edgenet")
}


#' @export
#' @method edgenet default
edgenet.default <- function(X, Y, G.X=NULL, G.Y=NULL,
                            lambda=1, psigx=1, psigy=1,
                            thresh=1e-5, maxit=1e5,
                            family=gaussian)
{
    check.matrices(X, Y)
    n <- dim(X)[1]
    p <- dim(X)[2]
    q <- dim(Y)[2]

    if (is.null(G.X)) psigx <- 0
    if (is.null(G.Y)) psigy <- 0

    check.graphs(X, Y, G.X, G.Y, psigx, psigy)
    check.dimensions(X, Y, n, p)
    lambda <- check.param(lambda, 0, `<`, 0)
    psigx <- check.param(psigx, 0, `<`, 0)
    psigy <- check.param(psigy, 0, `<`, 0)
    maxit <- check.param(maxit, 0, `<`, 1e5)
    thresh <- check.param(thresh, 0, `<`, 1e-5)

    if (q == 1) {
        psigy <- 0
        G.Y <- NULL
    }
    family <- get.family(family)

    # estimate coefficients
    ret <- .edgenet(X = X, Y = Y,
                    G.X = G.X, G.Y = G.Y,
                    lambda = lambda, psigx = psigx, psigy = psigy,
                    thresh = thresh, maxit = maxit, family = family)
    ret$call   <- match.call()
    class(ret) <- c(class(ret), "edgenet")

    ret
}


#' @noRd
#' @import Rcpp
.edgenet <- function(X, Y, G.X, G.Y,
                     lambda, psigx, psigy,
                     thresh, maxit, family)
{
    res <- .fit.edgenet(
        X, Y, G.X, G.Y, lambda, psigx, psigy, maxit, thresh, family)

    # finalize output
    beta  <- matrix(res$beta, ncol(X))
    alpha <- res$alpha
    rownames(beta) <- colnames(X)
    colnames(beta) <- colnames(Y)
    ret <- list(beta = beta,
                alpha = alpha,
                lambda = lambda,
                psigx = psigx,
                psigy = psigy)

    ret$family <- family
    class(ret) <- paste0(family, ".edgenet")

    ret
}


#' @noRd
#' @import Rcpp tensorflow
.fit.edgenet <- function(
    X, Y, gx, gy,
    lambda, psigx, psigy,
    maxit, thresh, family)
{

    tf$reset_default_graph()

    X <- tf$cast(X, tf$float32)
    Y <- tf$cast(Y, tf$float32)
    if (!is.null(gx))
        gx <- tf$cast(netReg:::laplacian_(gx), tf$float32)
    if (!is.null(gy))
        gy <- tf$cast(netReg:::laplacian_(gy), tf$float32)

    beta  <- tf$Variable(tf$zeros(shape(ncol(X), ncol(Y))))
    alpha <- tf$Variable(tf$zeros(shape(ncol(Y))))
    ones  <- tf$ones(shape(nrow(Y), 1), tf$float32)

    loss  <- .edgenet.loss(alpha, beta, X, Y, gx, gy, ones,
                           lambda, psigx, psigy, family)
    coefs <- fit(loss, alpha, beta)
    coefs
}

.edgenet.loss <- function(alpha, beta, x, y, gx, gy, ones, lambda, psigx, psigy, family)
{
    loss.function <- switch(family, "gaussian" = .edgenet.gaussian.loss)

    loss <- function(alpha, beta) {
        eta <- tf$matmul(x, beta) + ones * tf$transpose(alpha)
        obj <- loss.function(y, eta, lambda) + .lasso(lambda, beta)

        if (!is.null(gx)) {
            x.penalty <- tf$trace(tf$matmul(tf$transpose(beta), tf$matmul(gx, beta)))
            obj <- obj + psigx * x.penalty
        }
        if (!is.null(gy)) {
            y.penalty <- tf$trace(tf$matmul(beta, tf$matmul(gy, tf$transpose(beta))))
            obj <- obj + psigy * y.penalty
        }

        obj
    }

    loss
}

.edgenet.gaussian.loss <- function(y, mean, lambda)
{
    obj <- tf$reduce_sum(tf$pow(y - mean, 2))
    obj
}

.edgenet.binomial.loss <- function(y, means, lambda, ncol)
{
    obj <-  tf$reduce_sum(
        sapply(seq(ncol), function(j) {
            prob <- tfp$distributions$Bernoulli(logits = means[,j])
            tf$reduce_sum(prob$lob_prob(y[,j]))
        })
    )
    -obj
}

.lasso <- function(lambda, beta) {
    lambda * tf$reduce_sum(tf$abs(beta))
}
