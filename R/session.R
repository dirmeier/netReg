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
#' @import tensorflow
fit <- function(loss, alpha, beta, maxit=1000, learning.rate = 0.01, thresh = 1e-4)
{
  optimizer <- tf$train$AdamOptimizer(learning_rate = learning.rate)
  objective <- loss(alpha, beta)
  train <- optimizer$minimize(objective)

  sess <- tf$Session()
  sess$run(tf$global_variables_initializer())

  target.old <- Inf
  for (step in seq(maxit)) {
      sess$run(train)
      if (step %% 25 == 0) {
          target <- sess$run(objective)
          if (sum(abs(target - target.old)) < thresh)
              break
          target.old <- target
      }
  }

  alpha <- sess$run(alpha)
  beta <- sess$run(beta)
  sess$close()

  list(beta=beta, alpha=alpha)
}
