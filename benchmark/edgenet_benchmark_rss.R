library(netReg)
library(glmnet)
library(uuid)


n      <- 100
p      <- 10
q      <- 1
thresh <- 1e-7
maxit  <- 10000


G.X <- matrix(1, p, p) 
G.X[1:10, 1:10]   <- 5
G.X[40:50, 40:50] <- 1
G.X[80:90, 80:90] <- 1


X <- matrix(rnorm(n*p), n, p)
B <- matrix(rnorm(p*q), p, q)
E <- matrix(rnorm(n*q), n, q)
Y <- X %*% B + E

G.Y <- matrix(1, q, q) 
G.Y[upper.tri(G.Y)] <- rbeta(q * (q-1) / 2, 1, 2)
G.Y <- t(G.Y) + G.Y

cv.lasso <- cv.glmnet(X, Y, type.measure="mse") 
cv.edge  <- cv.edgenet(X, Y, G.X, thresh=thresh, maxit=maxit, family="gaussian" )
