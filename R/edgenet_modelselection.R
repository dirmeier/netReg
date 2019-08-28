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
setGeneric(
    "cv.edgenet",
    function(X, Y, G.X=NULL, G.Y=NULL,
             lambda=NA_real_, psigx=NA_real_, psigy=NA_real_,
             thresh=1e-5, maxit=1e5, learning.rate=0.01,
             family=gaussian,
             optim.maxit=1e2,
             nfolds=10)
    {
        standardGeneric("cv.edgenet")
    },
    package = "netReg"
)


#' @noRd
setMethod(
    "cv.edgenet",
    signature = signature(X="matrix", Y="numeric"),
    function(X, Y, G.X=NULL, G.Y=NULL,
             lambda=NA_real_, psigx=NA_real_, psigy=NA_real_,
             thresh=1e-5, maxit=1e5, learning.rate=0.01,
             family=gaussian,
             optim.maxit=1e2,
             nfolds=10)
    {
        cv.edgenet(X, as.matrix(Y), G.X, G.Y,
                   lambda, psigx, psigy,
                   thresh, maxit, learning.rate,
                   family,
                   optim.maxit,
                   nfolds)
    }
)


#' @noRd
setMethod(
    "cv.edgenet",
    signature = signature(X="matrix", Y="matrix"),
    function(X, Y, G.X=NULL, G.Y=NULL,
             lambda=NA_real_, psigx=NA_real_, psigy=NA_real_,
             thresh=1e-5, maxit=1e5, learning.rate=0.01,
             family=gaussian,
             optim.maxit=1e2,
             nfolds=10)
    {
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
        family <- get.family(family)

        if (ncol(Y) == 1) {
            psigy <- 0
            G.Y <- NULL
        }

        if (n < nfolds) nfolds <- n
        folds <- sample(rep(seq_len(10), length.out=n))

        # estimate shrinkage parameters
        ret <- .cv.edgenet(
            X, Y, G.X, G.Y,
            lambda, psigx, psigy,
            family,
            thresh, maxit, learning.rate,
            nfolds, folds,
            optim.maxit)

        ret$call   <- match.call()
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
                        optim.maxit)
{
    n <- nrow(x)
    p <- ncol(x)
    q <- ncol(y)

    if (!is.null(gx))
        gx <- cast_float(netReg:::laplacian_(gx))
    if (!is.null(gy))
        gy <- cast_float(netReg:::laplacian_(gy))

    alpha <- zero_vector(q)
    beta  <- zero_matrix(p, q)

    x.tensor <- tf$placeholder(tf$float32, shape(NULL, p), name = "x.tensor")
    y.tensor <- tf$placeholder(tf$float32, shape(NULL, q), name = "y.tensor")
    lambda.tensor <- tf$placeholder(tf$float32, shape())
    psigx.tensor <- tf$placeholder(tf$float32, shape())
    psigy.tensor <- tf$placeholder(tf$float32, shape())

    #estimate coefficients
    loss  <- edgenet.loss(gx, gy, family, q)
    objective <- loss(alpha, beta,
                      lambda.tensor, psigx.tensor, psigy.tensor,
                      x.tensor, y.tensor)

    optimizer <- tf$train$AdamOptimizer(learning_rate = learning.rate)
    train <- optimizer$minimize(objective)

    fn <- function(params, ..., sess, alpha, beta)
    {
        params <- .get.params(params, ...)
        print(params)
        losses <- vector(mode = "double", length = nfolds)

        for (fold in seq(nfolds)) {
            x.train <- x[which(folds != fold), ]
            y.train <- y[which(folds != fold),,drop=FALSE]
            x.test  <- x[which(folds == fold), ]
            y.test  <- y[which(folds == fold),,drop=FALSE]

            sess$run(tf$global_variables_initializer())

            target.old <- Inf
            for (step in seq(maxit))
            {
                sess$run(train, feed_dict = dict(x.tensor = x.train,
                                                 y.tensor = y.train,
                                                 lambda.tensor=params[1],
                                                 psigx.tensor=params[2],
                                                 psigy.tensor=params[3]))
                if (step %% 25 == 0) {
                    target <- sess$run(objective, feed_dict = dict(x.tensor = x.test,
                                                                   y.tensor = y.test,
                                                                   lambda.tensor=params[1],
                                                                   psigx.tensor=params[2],
                                                                   psigy.tensor=params[3]))
                    if (sum(abs(target - target.old)) < thresh)
                        break
                    target.old <- target
                }
            }

            losses[fold] <- sess$run(
                objective, feed_dict = dict(x.tensor = x.test,
                                            y.tensor = y.test,
                                            lambda.tensor=params[1],
                                            psigx.tensor=params[2],
                                            psigy.tensor=params[3]))
        }

        sum(losses)
    }

    reg.params <- list(lambda=lambda, psigx=psigx, psigy=psigy)
    params <- rep(0.1, sum(is.na(reg.params)))
    fixed.params <- Filter(is.finite, reg.params)

    print(params)
    print(fixed.params)

    with(tf$Session() %as% sess, {
        opt <- optim(fn, params, var.args=fixed.params,
                     sess=sess, alpha=alpha, beta=beta,
                     lower=rep(0, length(params)), upper=rep(100, length(params)),
                     control=list(maxeval=optim.maxit, xtol_rel = 1e-2))
    })

    ret <- opt
    ret$family <- family$family
    class(ret) <- paste0(family$family, ".cv.edgenet")

    ret
}
