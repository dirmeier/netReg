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
cross.validate <- function(objective, train,
                           x, y,
                           x.tensor, y.tensor,
                           lambda.tensor, psigx.tensor, psigy.tensor,
                           nfolds, folds,
                           maxit = 1000, thresh = 1e-5, learning.rate = 0.01) {
  fn <- function(params, ..., sess, alpha, beta) {
    params <- .get.params(params, ...)
    losses <- vector(mode = "double", length = nfolds)

    for (fold in seq(nfolds)) {
      x.train <- x[which(folds != fold), , drop = FALSE]
      y.train <- y[which(folds != fold), , drop = FALSE]
      x.test <- x[which(folds == fold), , drop = FALSE]
      y.test <- y[which(folds == fold), , drop = FALSE]

      sess$run(init_variables())

      target.old <- Inf
      for (step in seq(maxit))
      {
        sess$run(train, feed_dict = dict(
          x.tensor = x.train,
          y.tensor = y.train,
          lambda.tensor = params[1],
          psigx.tensor = params[2],
          psigy.tensor = params[3]
        ))
        if (step %% 25 == 0) {
          target <- sess$run(objective,
            feed_dict = dict(
              x.tensor = x.test,
              y.tensor = y.test,
              lambda.tensor = params[1],
              psigx.tensor = params[2],
              psigy.tensor = params[3]
            )
          )
          if (sum(abs(target - target.old)) < thresh) {
            break
          }
          target.old <- target
        }
      }

      losses[fold] <- sess$run(
        objective,
        feed_dict = dict(
          x.tensor = x.test,
          y.tensor = y.test,
          lambda.tensor = params[1],
          psigx.tensor = params[2],
          psigy.tensor = params[3]
        )
      )
    }

    sum(losses)
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
