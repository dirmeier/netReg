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
cross.validate <- function(loss.function, nfolds, folds,
                           x, y, gx, gy, family,
                           maxit=1000, thresh=1e-5, learning.rate=0.01)
{
    fn <- function(params, ..., alpha, beta)
    {
        params <- .get.params(params, ...)
        loss <- loss.function(params[1], params[2], params[3], gx, gy, family, ncol(y))
        losses <- vector(mode = "double", length = nfolds)
        for (fold in seq(nfolds)) {
            x.train <- x[which(folds != fold), ]
            y.train <- y[which(folds != fold),,drop=FALSE]
            x.test  <- x[which(folds == fold), ]
            y.test  <- y[which(folds == fold),,drop=FALSE]

            coefs <- fit(loss, alpha, beta, cast_float(x.train), cast_float(y.train), maxit, thresh, learning.rate)
            losses[fold] <- run(loss(
                constant_float(coefs$alpha), constant_float(coefs$beta), cast_float(x.test), cast_float(y.test)))
        }

        sum(losses)
    }

    fn
}


.get.params <- function(params, ...)
{
    var.args <- list(...)
    par.len <- length(params)
    if (length(var.args) + par.len != 3) stop("you provided wrong params")

    lambda <- if (methods::hasArg(lambda)) {
        var.args$lambda
    } else params[1]

    psigx  <- if (methods::hasArg(psigx))  {
        var.args$psigx
    } else if (methods::hasArg(lambda)) {
        params[1]
    } else params[2]

    psigy  <- if (methods::hasArg(psigy)) {
        var.args$psigy
    } else if (methods::hasArg(lambda) & methods::hasArg(psigx)) {
        params[1]
    } else params[2]

    print(c(lambda, psigx, psigy))
    c(lambda, psigx, psigy)
}
