library(netReg)
library(uuid)

args    <- commandArgs(trailingOnly=T)
stopifnot(length(args) == 8)

p <- q <- sig <- n <- 0
for (i in seq(1, length(args), by=2))
{
  if      (args[i] == "-s") sig <- as.integer(args[i + 1])
  else if (args[i] == "-n")   n <- as.integer(args[i + 1])
  else if (args[i] == "-p")   p <- as.integer(args[i + 1])
  else if (args[i] == "-q")   q <- as.integer(args[i + 1])
  else stop(paste("wrong flag", args[i], "\n"))
}

cat(paste0("Measuring rss with n=", n, ", p=", p, ", q=", q, ", sig=", sig, "\n"))

thresh <- 1e-12
maxit  <- 100000000

G.X <- matrix(0, p, p)
G.X[1:floor(p/3), 1:floor(p/3)]                       <- 1
G.X[(1+floor(p/3)):floor(2*p/3), (1+floor(p/3)):floor(2*p/3)] <- 1
G.X[(1+floor(2*p/3)):p,(1+floor(2*p/3)):p]                   <- 1

G.Y <- matrix(0, q, q)
G.Y[1:floor(q/3), 1:floor(q/3)]                       <- 1
G.Y[(1+floor(q/3)):floor(2*q/3), (1+floor(q/3)):floor(2*q/3)] <- 1
G.Y[(1+floor(2*q/3)):q, (1+floor(2*q/3)):q]                   <- 1

B <- matrix(0, p, q)
B[1:floor(p/3), 1:floor(q/3)]                           <- rnorm(length(1:floor(p/3)) * length(1:floor(q/3)), 1, 0.25)
B[(1+floor(p/3)):floor(2*p/3), floor(q/3):floor(2*q/3)] <- rnorm(length((1+floor(p/3)):floor(2*p/3)) * length(floor(q/3):floor(2*q/3)), 2, 0.25)
B[(1+floor(2*p/3)):p, (1+floor(2*q/3)):q]               <- rnorm(length((1+floor(2*p/3)):p) * length((1+floor(2*q/3)):q), 3, 0.25)

X <- matrix(rnorm(n*p), n, p)
E <- matrix(rnorm(n*q, 0, sig), n, q)

Y <- X %*% B + E


n.folds <- 5
folds <- sample(cut(seq_along(1:n), n.folds, labels = FALSE))

.rss <- function(obj, X, Y)
{
  Y.hat <- predict(obj, X)
  netReg:::rss(Y, Y.hat)
}

l <- c()
for (cv in 1:n.folds)
{
  X.train  <- X[folds != cv, ]
  Y.train  <- Y[folds != cv, ]
  X.test  <-  X[folds == cv, ]
  Y.test  <-  Y[folds == cv, ]

  cv.lasso <- cv.edgenet(X.train, Y.train,
                         thresh=thresh, maxit=maxit, family="gaussian", nfolds=5)
  las     <-  edgenet(X.train, Y.train,
                     lambda=cv.lasso$lambda,
                     thresh=thresh, maxit=maxit, family="gaussian")

  cv.edge  <- cv.edgenet(X.train, Y.train, G.X=G.X, G.Y=G.Y,
                         thresh=thresh, maxit=maxit, family="gaussian", nfolds=5)
  edge <-  edgenet(X.train, Y.train, G.X=G.X, G.Y=G.Y,
                      lambda=cv.edge$lambda, psigx=cv.edge$psigx, psigy=cv.edge$psigy,
                      thresh=thresh, maxit=maxit, family="gaussian")

  rss.lasso <- .rss(las,  X.test, Y.test)
  rss.edge  <- .rss(edge, X.test, Y.test)
  l <- rbind(l, c(Lasso=rss.lasso, Edgenet=rss.edge))
}

uuid                         <- uuid::UUIDgenerate()
path                         <- "~/PROJECTS/netreg_project/results/"
if (!file.exists(path)) path <- "/cluster/home/simondi/results/netReg/"
path                         <- paste0(path, "benchmark_rss_n_", n, "_p_", p, "_q_", q, "_sig_", sig, "_", uuid, ".rds")
cat(paste("going to:", path, "\n"))
saveRDS(l , path)
