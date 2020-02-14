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
#

context("coefficients")

test_that("gaussian without regularization reproduces stats::glm", {
    # generate data
    set.seed(1)

    n <- 100 # observation count
    p <- 10 # covariate count
    q <- 1 # response count

    X <- matrix(rnorm(n * p), n, p)
    B <- matrix(rnorm(p * q), p, q)

    Y <- X %*% B + matrix(rnorm(n * q, 0, 0.1), n, q)

    # fit glm
    fit.glm <- glm(Y ~ X, family=gaussian())

    # fit edgenet
    fit.nr <- edgenet(
        X, Y,
        lambda=0, psigx=0, psigy=0,
        family="gaussian")

    # compare results
    coef.glm <- unname(coef(fit.glm))
    coef.nr <- unname(coef(fit.nr)[,'y[1]'])

    testthat::expect_equal(c(0, B), coef.glm, tolerance=0.1)
    testthat::expect_equal(coef.glm, coef.nr, tolerance=0.1)
})


test_that("inverse gaussian without regularization reproduces stats::glm", {
    # generate data
    set.seed(1)

    n <- 100 # observation count
    p <- 3 # covariate count
    q <- 1 # response count

    X <- matrix(rnorm(n * p, 10, 4), n, p)
    B <- matrix(rnorm(p * q, 10, 5), p, q)

    mu.prior <- inverse.gaussian()$linkinv(X %*% B)
    Y <- statmod::rinvgauss(n, mean=mu.prior, shape=1)

    # fit glm
    fit.glm <- glm(Y ~ X, family=inverse.gaussian())

    # fit edgenet
    fit.nr <- edgenet(
        X, Y,
        lambda=0, psigx=0, psigy=0,
        family="inverse.gaussian")

    # compare results
    coef.glm <- unname(coef(fit.glm))
    coef.nr <- unname(coef(fit.nr)[,'y[1]'])

    testthat::expect_equal(
        mu.prior,
        inverse.gaussian()$linkinv(X %*% coef.glm[-1]),
        tolerance=0.1
    )
    testthat::expect_equal(coef.glm, coef.nr, tolerance=0.1)
})

test_that("beta without regularization reproduces betareg::betareg", {
    # generate data
    set.seed(1)

    n <- 200
    p <- 5
    q <- 1

    X <- matrix(rnorm(n * p), n, p)
    B <- matrix(rnorm(p * q), p, q)

    mu.prior <- mgcv::betar(theta=1, link="logit")$linkinv(X %*% B)
    phi.prior <- 1

    p.beta <- mu.prior * phi.prior
    q.beta <- (1 - mu.prior) * phi.prior

    Y <- matrix(sapply(1:length(p.beta), function(i) { rbeta(1, p.beta[i], q.beta[i]) }), n)

    # truncate reponse to circumvent stability problems
    eps <- .Machine$double.eps * 1e9 # see ?mgcv::betar  5
    Y[Y > 1-eps] <- 1 - eps
    Y[Y < eps] <- eps

    # fit betareg
    fit.br <- betareg::betareg(Y ~ X, link="logit")

    # fit edgenet
    fit.nr <- edgenet(
        X, Y,
        lambda=0, psigx=0, psigy=0,
        family=mgcv::betar(theta=1, link="logit"))

    # compare results
    coef.br <- unname(coef(fit.br))
    coef.nr <- unname(coef(fit.nr)[,'y[1]'])

    coef.br <- coef.br[-length(coef.br)] # remove phi estimate

    testthat::expect_equal(
        mu.prior,
        mgcv::betar(theta=1, link="logit")$linkinv(X %*% coef.br[-1]),
        tolerance=0.1
    )
    testthat::expect_equal(coef.br, coef.nr, tolerance=0.1)
})
