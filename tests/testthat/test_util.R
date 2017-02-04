context("util")

Y <- matrix(1, 10, 1)
X <- matrix(0, 10, 1)
B <- matrix(1, 1, 1)

test_that("intercept", {
  intr <- intercept(Y, X, B, 10)
  expect_equal(intr[1,1], 1)
})

test_that("intercept matrix", {
  intr <- intercept.matrix(10, 1:5)
  expect_equal(intr[,1], rep(1, 10))
})

test_that("cv sets train", {
  cv <- cvsets(10,10,23)
  for (i in seq(length(cv)))
  {
    for (j in seq(length(cv)))
    {
      if (i == j) expect_true(length(setdiff(cv[[i]]$train, cv[[j]]$train)) == 0) 
      else expect_true(length(setdiff(cv[[i]]$train, cv[[j]]$train)) == 1)
    } 
  }
})

test_that("cv sets test", {
  cv <- cvsets(10,10,23)
  for (i in seq(length(cv)))
  {
    for (j in seq(length(cv)))
    {
      if (i == j) expect_true(cv[[i]]$test == cv[[j]]$test) 
      else expect_false(cv[[i]]$test == cv[[j]]$test) 
    } 
  }
})
