#' Predict method for gaussian edgenet fits
#' 
#' @export
#' 
#' @author Simon Dirmeier, \email{netreg@@simon-dirmeier.net}
#' @description Predicts the estimated Y.hat values for a newdata design matrix X 
#' similar to the other predict methods, e.g. from glm and glmnet
#' 
#' @param object a fitted object of class \emph{gaussian.edgenet}
#' @param ... further arguments
#' @param newdata a new (\code{m} x \code{p})-dimensional design matrix with a 
#' variable number of observations \code{m}, but a constant number 
#' of co-variables \code{p}
#' 
#' @return A (\code{m} x \code{q})-dimensional matrix
#' 
#' @examples
#' X <- matrix(rnorm(100*10),100,10)
#' Y <- matrix(rnorm(100*10),100,10)
#' G.X <- matrix(rpois(10*10,1),10)
#' G.X <- t(G.X) + G.X
#' diag(G.X) <- 0
#' fit <- edgenet(X=X, Y=Y, G.X=G.X, family="gaussian")
#' pred <- predict(fit, X)
predict.gaussian.edgenet <- 
function
(
 object, 
 newdata=NULL,
 ...
)
{
  if(is.null(newdata))
    stop("newdata is null")
  X <- newdata
  n <- dim(X)[1]
  p <- dim(X)[2]
  coefs <- coef(object)   
  if(p != dim(coefs)[1]) stop("newdata dimensions do not fit coefficient dimensions!")
  mu <- object$intercept
  Y.hat <- X %*% coefs + intercept.matrix(n=n, mu=mu)
  Y.hat
}