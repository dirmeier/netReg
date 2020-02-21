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
group.lasso.loss <- function(grps, family)
{
    invlink <- family$linkinv
    loss.function <- family$loss

    loss <- function(alpha, beta, lambda, x, y)
    {
        eta <- linear.predictor(alpha, beta, x)
        obj <- loss.function(y, eta, invlink) +
            group.lasso.penalty(lambda, beta, grps)

        obj
    }

    loss
}
