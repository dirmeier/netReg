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


context("group lasso")

n <- 100L
p <- 5L
grps <- rep(1L, p)

X <- matrix(rnorm(n * p), n, p)
B <- rnorm(p)
Y <- rbinom(n, 1, 1 / (1 + exp(-X %*% B)))

test_that("gaussian group.lasso returns s3 class", {
  e <- group.lasso(X, Y, grps, lambda = 1, maxit = 1, thresh = 1e-5)
  testthat::expect_s3_class(e, "group.lasso")
  testthat::expect_s3_class(e, "gaussian.group.lasso")
})


test_that("gaussian group.lasso throws at wrong X", {
  testthat::expect_error(group.lasso(X = "s", Y, grps))
})


test_that("gaussian group.lasso throws at wrong Y", {
  testthat::expect_error(group.lasso(X, Y = "s", grps))
})


test_that("gaussian group.lasso warns at lambda<0", {
  testthat::expect_warning(group.lasso(X, Y, grps, lambda = -1, maxit = 1))
})


test_that("gaussian group.lasso warns at maxit<0", {
  testthat::expect_warning(group.lasso(X, Y, grps, maxit = -1))
})


test_that("gaussian group.lasso warns at thresh<0", {
  testthat::expect_warning(group.lasso(X, Y, grps, thresh = -1, maxit = 1))
})


test_that("gaussian group.lasso throws at wrong grps", {
  grps <- c(rep(1L, p - 1), p + 1L)
  testthat::expect_error(group.lasso(X, Y, grps, thresh = -1, maxit = 1))
  grps <- c(rep(1L, p - 1), 0L)
  testthat::expect_error(group.lasso(X, Y, grps, thresh = -1, maxit = 1))
  grps <- abs(rnorm(p))
  testthat::expect_error(group.lasso(X, Y, grps, thresh = -1, maxit = 1))
})
