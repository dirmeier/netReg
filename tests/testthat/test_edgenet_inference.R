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


context("edgenet inference")

set.seed(42)

n <- 100L
p <- 10L
q <- 1L

G.X <- matrix(rbeta(p * p, 1, 1), p, p)
G.X <- t(G.X) + G.X
G.Y <- matrix(rbeta(q * q, 1, 1), q, q)
G.Y <- t(G.Y) + G.Y
diag(G.X) <- diag(G.Y) <- 0
X <- matrix(rnorm(n * p, 0, 0.1), n, p)
B <- rnorm(p, 0, .1)


test_that("binomial dry run", {
  eta <- 1 / (1 + base::exp(-X %*% B))
  Y <- rbinom(n, 1, eta)
  try({
    fit.nr <- edgenet(X, Y,
      lambda = 1, psigx = 0, psigy = 0,
      family = binomial(), maxit = 1
    )
  })
})


test_that("predict throws at NULL", {
  Y <- X %*% B + matrix(rnorm(n * q), n)
  e <- edgenet(X, Y, lambda = 1, maxit = 1, thresh = 1e-5)
  testthat::expect_error(predict(e, NULL))
})


test_that("predict throws at wrong newdata dimension", {
  Y <- X %*% B + matrix(rnorm(n * q), n)
  e <- edgenet(X, Y, lambda = 1, maxit = 1, thresh = 1e-5)
  testthat::expect_error(predict(e, matrix(1, 25)))
})


test_that("gaussian without regularization reproduces stats::glm", {
  Y <- X %*% B + rnorm(n * q, 0, 0.1)
  fit.glm <- glm(Y ~ X, family = stats::gaussian)
  fit.nr <- edgenet(X, Y, lambda = 0, psigx = 0, psigy = 0, family = "gaussian")

  coef.glm <- unname(coef(fit.glm))
  coef.nr <- unname(coef(fit.nr)[, "y[1]"])

  testthat::expect_equal(coef.glm, coef.nr, tolerance = 0.1)
  testthat::expect_visible(predict(fit.nr, X))
})


test_that("binomial without regularization reproduces stats::glm", {
  eta <- 1 / (1 + base::exp(-X %*% B))
  Y <- rbinom(n, 1, eta)
  for (link in c("logit", "probit")) {
    fit.glm <- glm(Y ~ X, family = stats::binomial(link))
    fit.nr <- edgenet(X, Y, lambda = 0, psigx = 0, psigy = 0, family = binomial(link))
    coef.glm <- unname(coef(fit.glm))
    coef.nr <- unname(coef(fit.nr)[, "y[1]"])

    testthat::expect_equal(coef.glm, coef.nr, tolerance = 0.1)
    testthat::expect_visible(predict(fit.nr, X))
  }
})


test_that("poisson without regularization reproduces stats::glm", {
  eta <- base::exp(-X %*% B)
  Y <- rpois(n, eta)
  fit.glm <- glm(Y ~ X, family = stats::poisson)
  fit.nr <- edgenet(X, Y, lambda = 0, psigx = 0, psigy = 0, family = poisson)

  coef.glm <- unname(coef(fit.glm))
  coef.nr <- unname(coef(fit.nr)[, "y[1]"])

  testthat::expect_equal(coef.glm, coef.nr, tolerance = 0.1)
  testthat::expect_visible(predict(fit.nr, X))
})
