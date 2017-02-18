source("./edgenet_pure.R")
library(netReg)

n      <- 1000
ps     <- c(10, 50, 100, 200, 500, 1000, 10000)
q      <- 10
lambda <- 10
psigx  <- psigy <- 5
thresh <- 1e-9
maxit  <- 1000000

for (p in 10)
{
  X <- matrix(rnorm(n*p), n, p)
  B <- matrix(rnorm(p*q), p, q)
  E <- matrix(rnorm(n*q), n, q)
  Y <- X %*% B + E
  
  G.X <- matrix(0, p, p) 
  G.X[upper.tri(G.X)] <- rbeta(p * (p-1) / 2, 1, 2)
  G.X <- t(G.X) + G.X
  G.Y <- matrix(0, q, q) 
  G.Y[upper.tri(G.Y)] <- rbeta(q * (q-1) / 2, 1, 2)
  G.Y <- t(G.Y) + G.Y
  
  B.pure <- edgenet.pure.R(X, Y, G.X, G.Y, lambda, psigx, psigy, thresh, maxit)
  B.rcpp <- edgenet(X, Y, G.X, G.Y, lambda, psigx, psigy, thresh, maxit)$coefficients
}