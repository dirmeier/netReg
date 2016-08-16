library(lassoshooting)

n <- 10
p <- 2
q <- 1
X <- matrix(rnorm(n*p), n)
Y <- matrix(rnorm(n*q), n)

G.X <- matrix(as.numeric(rpois(p * p, 1)), p)
G.X <- t(G.X) + G.X
diag(G.X) <- 0

G.Y <- matrix(as.numeric(rpois(q*q, 1)), q)
G.Y <- t(G.Y) + G.Y
diag(G.X) <- 0

cv.edge <- cv.edgenet(X=X, Y=Y,   family="gaussian")

lambda <- 1000
edgenet(X=X, Y=Y, family="gaussian", lambda=lambda)
lassoshooting(X, Y,lambda=lambda)

##########

sigm <- function(x) 1/(1 + exp(-x))

p <- 1
q <- 1
n <- 100
thresh <- 1e-5
miter <- 1e5
X <- matrix(rnorm(n*p), n)
b.true <- rnorm(p)
y <- rbinom(100,1,.5)
y[y > .5] <- 1
y[y <= .5] <- 0
b <- b.old <- rnorm(p)
iter <- 0
repeat
{
  iter <- iter + 1
  bold <- b
  nu <- as.vector(X %*% b)
  mu <- sigm(nu)
  s <- mu * (1 - mu)
  S <- diag(s)
  z <- nu + solve(S)%*%(y-mu)
  b <- solve(t(X) %*% S %*% X) %*% t(X) %*% S %*% z
  if (sum(abs(b-b.old)) < thresh | iter > miter) break  
}
y.hat <- sigm(X%*%b)
y.hat[y.hat > .5] <- 1
y.hat[y.hat <= .5] <- 0
plot(as.vector(X), y, type="p")
points(as.vector(X), y.hat, col=2)
