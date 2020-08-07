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
#' @import tensorflow
#' @importFrom purrr transpose
fit <- function(mod, loss, x, y, maxit = 1000, learning.rate = 0.03, thresh = 1e-4) {
  optimizer <- keras::optimizer_adam(learning.rate)
  lo.old <- Inf

  for (step in seq_len(maxit)) {
    with(tf$GradientTape() %as% t, {
      lo <- loss(mod, x, y)
    })

    gradients <- t$gradient(lo, mod$trainable_variables)
    optimizer$apply_gradients(purrr::transpose(list(
      gradients, mod$trainable_variables
    )))

    if (step %% 25 == 0) {
      if (sum(abs(lo$numpy() - lo.old)) < thresh) {
        break
      }
      lo.old <- lo$numpy()
    }
  }

  list(
    beta = mod$beta$numpy(),
    alpha = mod$alpha$numpy()
  )
}
