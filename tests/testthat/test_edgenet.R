# netReg: graph-regularized linear regression models.
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


context("edgenet")

n <- 100L
p <- 5L
q <- 5L
G.X <- matrix(rbeta(p * p, 1, 1), p, p)
G.X <- t(G.X) + G.X
G.Y <- matrix(rbeta(q * q, 1, 1), q, q)
G.Y <- t(G.Y) + G.Y
diag(G.X) <- diag(G.Y) <- 0
X <- matrix(rnorm(n * p), n, p)
B <- matrix(rnorm(p * q), p, q)
Y <- X %*% B + matrix(rnorm(n * q), n)


test_that("gaussian edgenet returns s3 class", {
  e <- edgenet(X, Y, lambda = 1, maxit = 1, thresh = 1e-5)
  testthat::expect_s3_class(e, "edgenet")
  testthat::expect_s3_class(e, "gaussian.edgenet")
})


test_that("gaussian edgenet throws at wrong X", {
  testthat::expect_error(edgenet(X = "s", Y, G.X, G.Y))
})


test_that("gaussian edgenet throws at wrong Y", {
  testthat::expect_error(edgenet(X, Y = "s", G.X, G.Y))
})


test_that("gaussian edgenet warns at lambda<0", {
  testthat::expect_warning(edgenet(X, Y, G.X, G.Y, lambda = -1, maxit = 1))
})


test_that("gaussian edgenet warns at psigx<0", {
  testthat::expect_warning(edgenet(X, Y, G.X, G.Y, psigx = -1, maxit = 1))
})


test_that("gaussian edgenet warns at psigy<0", {
  testthat::expect_warning(edgenet(X, Y, G.X, G.Y, psigy = -1, maxit = 1))
})


test_that("gaussian edgenet warns at maxit<0", {
  testthat::expect_warning(edgenet(X, Y, G.X, G.Y, maxit = -1))
})


test_that("gaussian edgenet warns at thresh<0", {
  testthat::expect_warning(edgenet(X, Y, G.X, G.Y, thresh = -1))
})


test_that("gaussian edgenet throws at wrong dim[1] G.X", {
  G.X <- matrix(1, p - 1, p)
  testthat::expect_error(edgenet(X, Y, G.X, G.Y))
})


test_that("gaussian edgenet throws at wrong dim[2] G.X", {
  G.X <- matrix(1, p, p - 1)
  testthat::expect_error(edgenet(X, Y, G.X, G.Y))
})


test_that("gaussian edgenet throws at wrong dim[1] G.Y", {
  G.Y <- matrix(1, q - 1, q)
  testthat::expect_error(edgenet(X, Y, G.X, G.Y))
})


test_that("gaussian edgenet throws at wrong dim[2] G.Y", {
  G.Y <- matrix(1, q, q - 1)
  testthat::expect_error(edgenet(X, Y, G.X, G.Y))
})
