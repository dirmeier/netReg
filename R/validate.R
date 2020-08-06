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


#' @noRd
cross.validate <- function(mod, loss,
                           x, y,
                           lambda.tensor, psigx.tensor, psigy.tensor,
                           nfolds, folds,
                           maxit = 1000, thresh = 1e-5, learning.rate = 0.01) {
  fn <- function(params, ...) {
    params <- .get.params(params, ...)
    lambda.tensor$assign(params[1])
    psigx.tensor$assign(params[2])
    psigy.tensor$assign(params[3])

    losses <- vector(mode = "double", length = nfolds)

    for (fold in seq(nfolds)) {
      x.train <- x[which(folds != fold), , drop = FALSE]
      y.train <- y[which(folds != fold), , drop = FALSE]
      x.test <- x[which(folds == fold), , drop = FALSE]
      y.test <- y[which(folds == fold), , drop = FALSE]

      mod$reinit()
      invisible(fit(
        mod, loss,
        cast_float(x.train), cast_float(y.train),
        maxit, learning.rate, thresh
      ))
      losses[fold] <- loss(
        mod,
        cast_float(x.test),
        cast_float(y.test)
      )$numpy()
    }

    mean(losses)
  }

  fn
}


.get.params <- function(params, var.args) {
  par.len <- length(params)
  if (length(var.args) + par.len != 3) stop("you provided wrong params")


  has.lambda <- "lambda" %in% names(as.list(var.args))
  has.psigx <- "psigx" %in% names(as.list(var.args))
  has.psigy <- "psigy" %in% names(as.list(var.args))

  lambda <- if (has.lambda) {
    var.args$lambda
  } else {
    params[1]
  }

  psigx <- if (has.psigx) {
    var.args$psigx
  } else if (has.lambda) {
    params[1]
  } else {
    params[2]
  }

  psigy <- if (has.psigy) {
    var.args$psigy
  } else if (has.lambda & has.psigx) {
    params[1]
  } else {
    params[2]
  }

  c(lambda, psigx, psigy)
}
