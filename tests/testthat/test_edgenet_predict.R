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


context("netReg edgenet predict")

n   <- 100
p   <- 10
q   <- 10
G.X <- matrix(rbeta(p * p, 1, 1), p, p)
G.X <- t(G.X) + G.X
G.Y <- matrix(rbeta(q * q, 1, 1), q, q)
G.Y <- t(G.Y) + G.Y
diag(G.X) <- diag(G.Y) <- 0
X   <- matrix(rnorm(n * p), n, p)
B   <- matrix(rnorm(p * q), p, q)
Y   <- X %*% B + matrix(rnorm(n * q), n)

e <- edgenet(X, Y, lambda=1, maxit=1000, thresh=1e-5)

test_that("predict throws at NULL", {
    testthat::expect_error(predict(e, NULL))
})

test_that("predict throws at wrong newdata dimension", {
    testthat::expect_error(predict(e, matrix(1, 25)))
})
