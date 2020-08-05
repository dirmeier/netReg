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


context("edgenet modelselection")

n <- 100L
p <- 5L
q <- 2L
G.X <- matrix(rbeta(p * p, 1, 1), p, p)
G.X <- t(G.X) + G.X
G.Y <- matrix(rbeta(q * q, 1, 1), q, q)
G.Y <- t(G.Y) + G.Y
diag(G.X) <- diag(G.Y) <- 0
X <- matrix(rnorm(n * p), n, p)
B <- matrix(rnorm(p * q), p, q)
Y <- X %*% B + matrix(rnorm(n * q), n)

# ecv <- cv.edgenet(X, Y, maxit = 1, thresh = 1, optim.maxit = 1)
#
# test_that("gaussian edgenet modelselection returns s3 class", {
#   testthat::expect_s3_class(ecv, "gaussian.cv.edgenet")
#   testthat::expect_s3_class(ecv, "cv.edgenet")
# })
#
# test_that("gaussian edgenet modelselection throws at wrong X", {
#   testthat::expect_error(cv.edgenet(X = "s", Y, maxit = 1, thresh = 1))
# })
#
# test_that("gaussian edgenet modelselection throws at wrong Y", {
#   testthat::expect_error(cv.edgenet(X, Y = "s", maxit = 1, thresh = 1))
# })
#
# test_that("gaussian edgenet modelselection throws at wrong fold", {
#   testthat::expect_error(cv.edgenet(X, Y, maxit = 1, thresh = 1, nfolds = "s"))
# })
#
# test_that("gaussian edgenet modelselection warngs at maxit < 0 ", {
#   testthat::expect_warning(cv.edgenet(X, Y, thresh = 1, maxit = -1))
# })
#
# test_that("gaussian edgenet modelselection warngs at thresh < 0 ", {
#   testthat::expect_warning(cv.edgenet(X, Y, maxit = 1, thresh = -1))
# })
#
# test_that("gaussian edgenet modelselection warngs at epsilon < 0 ", {
#   testthat::expect_warning(cv.edgenet(X, Y, maxit = 1, thresh = 1, optim.thresh = -1))
# })
#
# test_that("gaussian edgenet modelselection warngs at approx.maxit < 0 ", {
#   testthat::expect_warning(cv.edgenet(X, Y, maxit = 1, thresh = 1, optim.maxit = -1))
# })
