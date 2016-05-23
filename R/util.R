#' @noRd
intercept <- function(Y,X,B,n)  ((t(Y - X %*% B) %*% rep(1,n)) / n )

#' @noRd
intercept.matrix <- function(n, mu)  rep(1, n) %*% t(mu)

#' @noRd
rss <- function(Y, Y.hat) sum((Y - Y.hat) ** 2) 

#' @noRd
cvsets <- 
function
( 
 n,
 folds=10
)
{      
  if (n < 1) stop("n<1; need positive integer!")
  if (folds < 0) stop("folds<0; need positive integer!")
  n <- as.integer(n)
  folds <- as.integer(folds)
  if (n<folds) stop("n<folds; need n>folds!")    
  id <- (1:n)[order(runif(n))]
  k <- as.integer(n * seq(1,folds - 1) / folds)
  k <- matrix(c(0 ,rep(k, each=2), n), ncol=2, byrow=TRUE)
  k[,1] <- k[,1] + 1
  l <- lapply(seq.int(folds),function(x,k,d) 
    list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
         test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}