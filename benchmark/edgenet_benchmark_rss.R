library(netReg)
library(glmnet)
library(uuid)

cms    <- commandArgs(trailingOnly=T)
stopifnot(length(cms) != 3)
p      <- as.integer(cms[1])
stopifnot(!is.na(p), is.numeric(p))
q      <- as.integer(cms[2])
stopifnot(!is.na(q), is.numeric(q))
sig      <- as.integer(cms[3])
stopifnot(!is.na(sig), is.numeric(sig))

n      <- 1000
thresh <- 1e-7
maxit  <- 100000

G.X <- matrix(0, p, p) 
G.X[1:floor(p/3), 1:floor(p/3)]                       <- 1
G.X[floor(p/3):floor(2*p/3), floor(p/3):floor(2*p/3)] <- 1
G.X[floor(2*p/3):p, floor(2*p/3):p]                   <- 1
diag(G.X)                                             <- 0

G.Y <- matrix(0, q, q) 
G.Y[1:floor(q/3), 1:floor(q/3)]                       <- 1
G.Y[floor(q/3):floor(2*q/3), floor(q/3):floor(2*q/3)] <- 1
G.Y[floor(2*q/3):q, floor(2*q/3):q]                   <- 1
diag(G.Y)                                             <- 0

X <- matrix(rnorm(n*p), n, p)
B <- matrix(rnorm(p*q), p, q)
E <- matrix(rnorm(n*q, 0, sig), n, q)
Y <- X %*% B + E

for (cv in 1:10)
{
  cve <- list()
  for (i in 1:q)
  {
    cve[[i]] <- cv.glmnet(X, Y[,i], type.measure="mse", thresh=thresh, 
                          standardize=F, maxit=maxit) 
  }
  cv.edge  <- cv.edgenet(X, Y, G.X=G.Y, G.Y=G.Y, thresh=thresh, maxit=maxit, family="gaussian" )
}