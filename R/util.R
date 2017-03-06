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


#' @noRd
intercept <-function(Y, X, B, n)
{
    ((t(Y - X %*% B) %*% rep(1, n)) / n)
}

#' @noRd
intercept.matrix <- function(n, mu)
{
    rep(1, n) %*% t(mu)
}

#' @noRd
rss <- function(Y, Y.hat)
{
  sum((Y - Y.hat) ** 2)
}

#' @noRd
#' @importFrom stats runif
cvsets <- function(n, folds = 10, seed = 23)
{
    if (n < 1)
        stop("n<1; need positive integer!")
    if (folds < 0)
        stop("folds<0; need positive integer!")
    n <- as.integer(n)
    folds <- as.integer(folds)
    if (n < folds)
        stop("n<folds; need n>folds!")
    id <- (1:n)[order(stats::runif(n))]
    k <- as.integer(n * seq(1, folds - 1) / folds)
    k <- matrix(c(0, rep(k, each = 2), n), ncol = 2, byrow = TRUE)
    k[, 1] <- k[, 1] + 1
    l <- lapply(seq.int(folds), function(x, k, d)
      list(train = d[!(seq(d) %in% seq(k[x, 1], k[x, 2]))],
           test = d[seq(k[x, 1], k[x, 2])]), k = k, d = id)
    return(l)
}

#' @noRd
check.matrices <- function(X, Y)
{
  stopifnot(is.matrix(X), is.matrix(Y))
}

#' @noRd
check.graphs <- function(X, Y, G.X, G.Y, psigx, psigy)
{
    stopifnot(is.matrix(G.X), is.matrix(G.Y))
    if (psigx != 0 & any(dim(G.X) != dim(X)[2]))
      stop("ncol(X) and dim(G.X) do not fit!")
    if (psigy != 0 & any(dim(G.Y) != dim(Y)[2]))
      stop("ncol(Y) and dim(G.Y) do not fit!")
    if (any(G.X < 0))
      stop("Some elements G.X<0; please use non-negative matrix!")
    if (any(G.Y < 0))
      stop("Some elements G.Y<0; please use non-negative matrix!")
}

#' @noRd
check.dimensions <- function(X, Y, n, p)
{
  if (dim(X)[1] != n)
    stop("X and Y have not same number of observations!")
  if (dim(X)[1] != dim(Y)[1])
    stop("X and Y have not same number of observations!")
  if (n != dim(Y)[1])
    stop("X and Y have not same number of observations!")
  if (p < 2)
    stop("Pls use a X matrix with at least 2 covariables!")
}
