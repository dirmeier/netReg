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


context("netReg utility")
 
Y <- matrix(1 + 1 : 10, 10, 1)
X <- matrix(1 : 10, 10, 1)
B <- matrix(1, 1, 1)
n <- 10

test_that("intercept is correct", {
  testthat::expect_equal(as.numeric(intercept(Y, X, B, n)), 1)
})

test_that("intercept matrix returns correctly", {
  m <- intercept.matrix(10, 1:5)
  for (i in 1:10)
  {
    testthat::expect_equal(m[i,], 1:5)
  }
})

test_that("check.matrix works", {
  m <- intercept.matrix(10, 1:5)
  testthat::expect_silent(check.matrices(m, m))
})

test_that("check.matrix throws", {
  testthat::expect_error(check.matrices(1, "s"))
})

test_that("check.graph works", {
  X   <- matrix(1, 10, 5)
  Y   <- matrix(1, 10, 3)
  G.X <- matrix(1, 5, 5)
  G.Y <- matrix(1, 3, 3)
  testthat::expect_silent(check.graphs(X, Y, G.X, G.Y, 1, 1))
})

test_that("check.graph throws at wrong X", {
  X   <- matrix(1, 10, 2)
  Y   <- matrix(1, 10, 3)
  G.X <- matrix(1, 5, 5)
  G.Y <- matrix(1, 3, 3)
  testthat::expect_error(check.graphs(X, Y, G.X, G.Y, 1, 1))
})

test_that("check.graph throws at wrong Y", {
  X   <- matrix(1, 10, 5)
  Y   <- matrix(1, 10, 2)
  G.X <- matrix(1, 5, 5)
  G.Y <- matrix(1, 3, 3)
  testthat::expect_error(check.graphs(X, Y, G.X, G.Y, 1, 1))
})

test_that("check.graph throws at G.X<0", {
  X   <- matrix(1, 10, 5)
  Y   <- matrix(1, 10, 3)
  G.X <- matrix(-1, 5, 5)
  G.Y <- matrix(1, 3, 3)
  testthat::expect_error(check.graphs(X, Y, G.X, G.Y, 1, 1))
})

test_that("check.graph throws at G.Y<0", {
  X   <- matrix(1, 10, 5)
  Y   <- matrix(1, 10, 3)
  G.X <- matrix(1, 5, 5)
  G.Y <- matrix(-1, 3, 3)
  testthat::expect_error(check.graphs(X, Y, G.X, G.Y, 1, 1))
})

test_that("check.dimensions is correct", {
  X <- matrix(1, 10, 5)
  Y <- matrix(1, 10, 3)
  n <- dim(X)[1]
  p <- dim(X)[2]
  testthat::expect_silent(check.dimensions(X, Y, n, p))
})

test_that("check.dimensions throws at wrong X/n", {
  X <- matrix(1, 10, 5)
  Y <- matrix(1, 10, 3)
  n <- dim(X)[1]
  p <- dim(X)[2]
  testthat::expect_error(check.dimensions(X, Y, 3L, p))
})

test_that("check.dimensions throws at wrong dimensions of X/Y", {
  X <- matrix(1, 10, 5)
  Y <- matrix(1, 9, 3)
  n <- dim(X)[1]
  p <- dim(X)[2]
  testthat::expect_error(check.dimensions(X, Y, n, p))
})

test_that("check.dimensions throws at wrong Y/n", {
  X <- matrix(1, 9, 5)
  Y <- matrix(1, 10, 3)
  n <- dim(X)[1]
  p <- dim(X)[2]
  testthat::expect_error(check.dimensions(X, Y, n, p))
})

test_that("check.dimensions throws at wrong Y/n", {
  X <- matrix(1, 9, 5)
  Y <- matrix(1, 10, 3)
  n <- dim(X)[1]
  p <- dim(X)[2]
  testthat::expect_error(check.dimensions(X, Y, n, 1))
})

test_that("rss works correctly", {
  X <- matrix(2, 10, 3)
  Y <- matrix(1, 10, 3)
  testthat::expect_equal(netReg:::rss(X, Y), 1 * 10 * 3)
})
