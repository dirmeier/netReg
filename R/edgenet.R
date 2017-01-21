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


#' Fit a graph-regularized linear regression model using edge-based regularization.
#' 
#' @export
#' 
#' @author Simon Dirmeier, \email{mail@@simon-dirmeier.net}
#' 
#' @description  Fit a graph-regularized linear regression model using edge-penalization.
#' The coefficients are computed using graph-prior knowledge in the form of 
#' one/two affinity matrices. Graph-regularization is an extension to previously 
#' introduced regularization techniques, such as the LASSO.
#'
#' @param X  input matrix, of dimension (\code{n} x \code{p}) 
#' where \code{n} is the number of observations and \code{p} is the number 
#' of covariables. Each row is an observation vector.
#' @param Y  output matrix, of dimension (\code{n} x \code{q}) 
#' where \code{n} is the number of observations and \code{q} is the number 
#' of response variables Each row is an observation vector.
#' @param G.X  non-negativ affinity matrix for \code{n}, of dimensions 
#' (\code{p} x \code{p}) where \code{p} is the number of covariables \code{X}
#' @param G.Y  non-negativ affinity matrix for \code{n}, of dimensions 
#' (\code{q} x \code{q}) where \code{q} is the number of covariables \code{Y}
#' @param lambda  shrinkage parameter for LASSO.
#' @param psigx  shrinkage parameter for graph-regularization of \code{G.X}
#' @param psigy  shrinkage parameter for graph-regularization of \code{G.Y}
#' @param thresh  threshold for coordinate descent
#' @param maxit  maximum number of iterations
#' @param family  family of response, e.g. gaussian
#' @param ...  additional params
#' 
#' @return An object of class \code{edgenet}
#' \item{coefficients }{ the estimated (\code{p} x \code{q})-dimensional coefficient matrix B.hat}
#' \item{intercept }{ the estimated (\code{q} x \code{1})-dimensional vector of intercepts}
#' \item{call }{ the call that produced the object}
#' \item{family }{ the family of the response}
#' 
#' @references 
#'  Friedman J., Hastie T., Hoefling H. and Tibshirani R. (2007), 
#'  Pathwise coordinate optimization.\cr
#'  \emph{The Annals of Applied Statistics}\cr \cr
#'  Friedman J., Hastie T. and Tibshirani R. (2010),
#'  Regularization Paths for Generalized Linear Models via Coordinate Descent.\cr
#'  \emph{Journal of Statistical Software}\cr \cr
#'  Fu W. J. (1998),  Penalized Regression: The Bridge Versus the Lasso.\cr
#'  \emph{Journal of Computational and Graphical Statistics}\cr \cr
#'  Cheng W. and Wang W. (2014), Graph-regularized dual Lasso for robust eQTL mapping.\cr
#'  \emph{Bioinformatics}
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100*10), 100, 10)
#' G.X <- matrix(rpois(100,1), 10)
#' G.X <- t(G.X) + G.X
#' diag(G.X) <- 0
#' 
#' # fit a Gaussian model
#' Y <- matrix(rnorm(100),100,1)
#' fit <- edgenet(X=X, Y=Y, G.X=G.X, family="gaussian")
#' 
#' # fit a binomial model
#' Y <- matrix(rbinom(100, 1, .5), 100, 1)
#' fit <- edgenet(X=X, Y=Y, G.X=G.X, family="binomial")
#' }
edgenet <- function(X, Y, G.X=NULL, G.Y=NULL, lambda=1, psigx=1, psigy=1, 
                    thresh=1e-5, maxit=1e5, family=c("gaussian"), ...) 
{
  UseMethod("edgenet")
}

#' @export
#' @method edgenet default
edgenet.default <- function(X, Y, G.X=NULL, G.Y=NULL,
                            lambda=1, psigx=1, psigy=1, 
                            thresh=1e-5, maxit=1e5,
                            family=c("gaussian"), ...)
{
  .check.matrices(X, Y)
  # parse dimensions
  n <- dim(X)[1]                              
  p <- dim(X)[2]     
  q <- dim(Y)[2]  
  if (is.null(G.X)) G.X <- matrix(0, 1, 1)
  if (is.null(G.Y)) G.Y <- matrix(0, 1, 1)
  .check.graphs(X, Y, G.X, G.Y, psigx, psigy)
  .check.dimensions(X, Y, n, p, q)
  # check if graphs are valid
  if (all(G.X == 0)) psigx <- 0
  if (all(G.Y == 0)) psigy <- 0
  if (lambda < 0) 
  {
    warning("lambda < 0, setting to 0!")
    lambda <- 0
  }
  if (psigx < 0)
  {
    warning("psigx < 0, setting to 0!")
    psigx <- 0
  }
  if (psigy < 0)
  {
    warning("psigy < 0, setting to 0!")
    psigy <- 0
  }
  if (maxit < 0)
  {
    warning("maxit < 0, setting to 1e5!")
    maxit <- 1e5
  }
  if (thresh < 0)
  {
    warning("thresh < 0, setting to 1e-5!")
    thresh <- 1e-5
  }
  if (q == 1) psigy <- 0
  family <- match.arg(family)
  # estimate coefficients
  ret <- .fit(X=X, Y=Y, 
              G.X=G.X, G.Y=G.Y,
              lambda=lambda,
              psigx=psigx, psigy=psigy,
              thresh=thresh, maxit=maxit,
              family=family) 
  ret$call   <- match.call()    
  class(ret) <- c(class(ret), "edgenet")
  ret
}
