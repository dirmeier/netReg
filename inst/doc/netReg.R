## ------------------------------------------------------------------------
  set.seed(23)
  X <- matrix(rnorm(1000*5), 1000)
  Y <- matrix(rnorm(1000*5), 1000)

## ------------------------------------------------------------------------
  library(netReg)

## ------------------------------------------------------------------------
  aff.mat <- matrix(rbeta(25, 1, 5),  5)
  aff.mat <- (t(aff.mat) + aff.mat) / 2
  diag(aff.mat) <- 0

## ------------------------------------------------------------------------
  fit <- edgenet(X=X, Y=Y, G.X=aff.mat, lambda=1, psigx=1, family="gaussian")
  print(fit)

## ------------------------------------------------------------------------
  X.new <- matrix(rnorm(10*5),10)
  pred  <- predict(fit, X.new)

## ------------------------------------------------------------------------
  cv <- cv.edgenet(X=X, Y=Y, G.X=aff.mat, family="gaussian", maxit=1000)
  print(cv)

## ------------------------------------------------------------------------
  p <- 25
  X <- matrix(rnorm(10*p), 10)
  Y <- matrix(rnorm(10*p), 10)
  aff.mat <- matrix(rgamma(p * p, 5, 1), p)
  aff.mat <- (t(aff.mat) + aff.mat)
  diag(aff.mat) <- 0
  cv <- cv.edgenet(X=X, Y=Y, G.X=aff.mat, family="gaussian", maxit=1000)
  print(cv)

