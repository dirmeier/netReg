 source("./edgenet_pure.R")
library(netReg)
library(microbenchmark)
library(uuid)

# args    <- commandArgs(trailingOnly=T)
# stopifnot(length(args) == 4)

# n <- p <- 0
# for (i in seq(1, length(args), by=2))
# {
#   if      (args[i] == "-n") n <- as.integer(args[i + 1])
#   else if (args[i] == "-p") p <- as.integer(args[i + 1])
#   else stop(paste("wrong flag", args[i], "\n"))
# }


.bench <- function(n, p, q, times, with=T)
{
  cat(paste0("Measuring time with n=", n, ", p=", p,  ", q=", q, "\n"))
  lambda <- 10
  psigx  <- psigy <- 5
  thresh <- 1e-7
  maxit  <- 100000

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

  if (with) {
  m <-
    microbenchmark(
      Pure=edgenet.pure.R(X, Y, G.X, G.Y, lambda, psigx, psigy, thresh, maxit),
      Cpp=.Call("edgenet_cpp", X, Y, G.X, G.Y,
                as.double(lambda), as.double(psigx),  as.double(psigy),
                as.integer(maxit), as.double(thresh),
                as.character("gaussian"))$coefficients,
      times=times)
  }
  else
  {
     m <-
      microbenchmark(
          Cpp=.Call("edgenet_cpp", X, Y, G.X, G.Y,
                    as.double(lambda), as.double(psigx),  as.double(psigy),
                    as.integer(maxit), as.double(thresh),
                    as.character("gaussian"))$coefficients,
          times=times)
  }
  m
}


m1 <- .bench(100, 100, 10, 10)
m1
saveRDS(m1, "~/Desktop/m1.rds" )

m2 <- .bench(10, 10, 10, 10)
m2
saveRDS(m2, "~/Desktop/m2.rds" )

m3 <- .bench(1000, 1000, 10, 10)
m3
saveRDS(m3, "~/Desktop/m3.rds" )

m4 <- .bench(10010, 10000, 10, 10, T)
m4
saveRDS(m4, "~/Desktop/m4.rds" )

m5 <- .bench(10010, 10000, 10, 10)
m5
saveRDS(m5, "~/Desktop/m5.rds" )

# uuid                         <- uuid::UUIDgenerate()
# path                         <- "~/PROJECTS/netreg_project/results/"
# if (!file.exists(path)) path <- "/cluster/home/simondi/results/netReg/"
# path                         <- paste0(path, "benchmark_time_n_", n, "_p_", p, "_", uuid, ".rds")
# cat(paste("going to:", path, "\n"))
# saveRDS(m , path)
