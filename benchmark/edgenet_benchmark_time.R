source("./edgenet_pure.R")
library(netReg)
library(microbenchmark)
library(uuid)

p <- 1000

cms    <- commandArgs(trailingOnly=T)
stopifnot(length(cms) != 0)
p      <- as.integer(cms[1])
stopifnot(!is.na(p), is.numeric(p))

n      <- 1000
q      <- 10
lambda <- 10
psigx  <- psigy <- 5
thresh <- 1e-7
maxit  <- 10000

X <- matrix(rnorm(n*p), n, p)
B <- matrix(rnorm(p*q), p, q)
E <- matrix(rnorm(n*q), n, q)
Y <- X %*% B + E

G.X <- matrix(1, p, p) 
G.X[upper.tri(G.X)] <- rbeta(p * (p-1) / 2, 1, 2)
G.X <- t(G.X) + G.X
G.Y <- matrix(1, q, q) 
G.Y[upper.tri(G.Y)] <- rbeta(q * (q-1) / 2, 1, 2)
G.Y <- t(G.Y) + G.Y

m <- microbenchmark(
              Pure=edgenet.pure.R(X, Y, G.X, G.Y, lambda, psigx, psigy, thresh, maxit),
               Cpp=.Call("edgenet_cpp", X, Y, G.X, G.Y,
                     as.double(lambda), as.double(psigx),  as.double(psigy),
                     as.integer(maxit), as.double(thresh),
                     as.character("gaussian"))$coefficients,
               times=10)

uuid <- uuid::UUIDgenerate()
path                         <- "~/PROJECTS/netreg_project/results/"
if (!file.exists(path)) path <- "/cluster/home/simondi/results/netReg/"
path <- paste0(path, "benchmark_", p, "_", uuid, ".rds")
cat(paste("going to:", path, "\n"))
saveRDS(m , path)
