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


context("group lasso inference")

options(warn = -1)
set.seed(42)

n <- 100L
p <- 5L
grps <- rep(1L, p)

X <- matrix(rnorm(n * p), n, p)
B <- rnorm(p)


test_that("binomial dry run", {
  eta <- 1 / (1 + base::exp(-X %*% B))
  Y <- rbinom(n, 1, eta)
  try({
    fit.nr <- group.lasso(X, Y, lambda = 1, family = binomial(), maxit = 1)
  })
})


test_that("predict throws at NULL", {
  Y <- X %*% B + rnorm(n)
  e <- group.lasso(X, Y, lambda = 1, maxit = 1, thresh = 1e-5)
  testthat::expect_error(predict(e, NULL))
})


test_that("predict throws at wrong newdata dimension", {
  Y <- X %*% B + rnorm(n)
  e <- edgenet(X, Y, lambda = 1, maxit = 1, thresh = 1e-5)
  testthat::expect_error(predict(e, matrix(1, 25)))
})


test_that("gaussian without regularization reproduces stats::glm", {
  Y <- X %*% B + rnorm(n, 0, 0.1)
  fit.glm <- glm(Y ~ X, family = stats::gaussian)
  fit.nr <- group.lasso(X, Y, lambda = 0, family = "gaussian")

  coef.glm <- unname(coef(fit.glm))
  coef.nr <- unname(coef(fit.nr)[, "y[1]"])

  testthat::expect_equal(coef.glm, coef.nr, tolerance = 0.1)
  testthat::expect_visible(predict(fit.nr, X))
})


test_that("binomial without regularization reproduces stats::glm", {
  eta <- 1 / (1 + base::exp(-X %*% B))
  Y <- rbinom(n, 1, eta)
  for (link in c("logit")) {
    fit.glm <- glm(Y ~ X, family = stats::binomial(link))
    fit.nr <- group.lasso(X, Y, lambda = 0, family = binomial(link))
    coef.glm <- unname(coef(fit.glm))
    coef.nr <- unname(coef(fit.nr)[, "y[1]"])

    testthat::expect_equal(coef.glm, coef.nr, tolerance = 0.1)
    testthat::expect_visible(predict(fit.nr, X))
  }
})


if (requireNamespace("grplasso", quietly = TRUE)) {
  test_that("gaussian group lasso correct output", {
    options(warn = -1)

    Y <- rnorm(n, X %*% B, .1)
    for (lam in c(0)) {
      e0 <- group.lasso(
        X, Y,
        grps = NULL, lambda = lam,
        maxit = 500, thresh = 5 * 10^-8, family = netReg::gaussian
      )
      e1 <- grplasso::grplasso(
        cbind(1, X),
        y = Y, index = c(NA, grps), lambda = lam,
        model = grplasso::LinReg(),
        center = FALSE, standardize = FALSE
      )

      coef.nr <- unname(coef(e0))
      coef.glm <- unname(coef(e1))
      testthat::expect_equal(coef.nr, coef.glm, tolerance = 0.1)
    }
  })
}
