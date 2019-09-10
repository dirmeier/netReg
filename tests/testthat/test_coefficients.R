context("coefficients")

test_that("gaussian without regularization reproduces stats::glm", {
    # generate data
    set.seed(1)

    n <- 100 # observation count
    p <- 10 # covariate count
    q <- 1 # response count

    X <- matrix(rnorm(n * p), n)
    B <- matrix(rnorm(p * q), p)

    Y <- X %*% B + matrix(rnorm(n * q, 0, 0.1), n)

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

    testthat::expect_equal(coef.glm, coef.nr, tolerance=.1)
})

test_that("inverse gaussian without regularization reproduces stats::glm", {
    # generate data
    set.seed(1)

    n <- 100 # observation count
    p <- 3 # covariate count
    q <- 1 # response count

    X <- matrix(rnorm(n * p, 10), n)
    B <- matrix(rnorm(p * q, 5), p) #matrix(c(1, 2, 3), p)

    mu <- 1 / sqrt(X %*% B)
    Y <- statmod::rinvgauss(n, mean=mu, shape=1)

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

    testthat::expect_equal(coef.glm, coef.nr, tolerance=1)
})
