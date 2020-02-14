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


#' @noRd
intercept <- function(Y, X, B, n)
{
    (t(Y - X %*% B) %*% rep(1, n)) / n
}


#' @noRd
intercept.matrix <- function(n, alpha)
{
    rep(1, n) %*% t(alpha)
}


#' @noRd
rss <- function(Y, Y.hat)
{
    sum((Y - Y.hat) ** 2)
}


#' @noRd
check.matrices <- function(X, Y)
{
    stopifnot(is.matrix(X), is.matrix(Y))
}


#' @noRd
check.graphs <- function(X, Y, G.X, G.Y, psigx, psigy)
{
    if (psigx != 0 & any(dim(G.X) != dim(X)[2]))
        stop("ncol(X) and dim(G.X) do not fit!")
    if (psigy != 0 & any(dim(G.Y) != dim(Y)[2]))
        stop("ncol(Y) and dim(G.Y) do not fit!")
    if (is.matrix(G.X) & any(G.X < 0))
        stop("Some elements G.X<0; please use non-negative matrix!")
    if (is.matrix(G.Y) & any(G.Y < 0))
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


#' @noRd
is.positive.numeric <- function(d)
{
    is.numeric(d) && d > 0
}


#' @noRd
check.param <- function(param, comp, op, replace.with)
{
    if (!is.na(param) & op(param, comp)) {
        warning(sprintf("%s < 0, setting to 0!", deparse(substitute(param))))
        param <- replace.with
    }

    param
}


# shamelessly copied from stats::glm
#' @noRd
get.family <- function(family)
{
    if (is.character(family))
        family <- get(family, mode = "function")
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        stop("'family' not recognized", call. = FALSE)
    }

    family
}


#' @noRd
not.supported.yet <- function(family)
{
    err <- sprintf(
        "family '%s' not supported yet. choose 'gaussian'/'binomial' please.",
        family
    )
    stop(err, call. =  FALSE)
}

warn.experimental <- function(family)
    warning(paste("family", family, "is still experimental. enjoy with care."),
            call. = FALSE)
