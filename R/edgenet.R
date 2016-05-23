#' Fit a graph-regularized linear regression model using edge-based regularization.
#' 
#' @export
#' @useDynLib netReg
#' @author Simon Dirmeier, \email{netreg@@simon-dirmeier.net}
#' @description  Fit a graph-regularized linear regression model using edge-penalization.
#' The coefficients are computed using graph-prior knowledge in the form of 
#' one/two affinity matrices. Graph-regularization is an extension to previously 
#' introduced regularization techniques, such as the LASSO.
#'
#' @param X  input matrix, of dimension (\code{n} x \code{p}) 
#' where \code{n} is the number of observations and \code{p} is the number 
#' of covariables. Each row is an observation vector.
#' @param Y  output matrix, of dimension (\code{n} x \code{q}) 
#' where \code{n} is the number of observations and \code{q} is the number 
#' of response variables Each row is an observation vector.
#' @param G.X  non-negativ affinity matrix for \code{n}, of dimensions 
#' (\code{p} x \code{p}) where \code{p} is the number of covariables \code{X}
#' @param G.Y  non-negativ affinity matrix for \code{n}, of dimensions 
#' (\code{q} x \code{q}) where \code{q} is the number of covariables \code{Y}
#' @param lambda  shrinkage parameter for LASSO.
#' @param psigx  shrinkage parameter for graph-regularization of \code{G.X}
#' @param psigy  shrinkage parameter for graph-regularization of \code{G.Y}
#' @param thresh  threshold for coordinate descent
#' @param maxit  maximum number of iterations
#' @param family  family of response, e.g. gaussian
#' 
#' @return An object of class \code{edgenet}
#' \item{coefficients }{ the estimated (\code{p} x \code{q})-dimensional coefficient matrix B.hat}
#' \item{intercept }{ the estimated (\code{q} x \code{1})-dimensional vector of intercepts}
#' \item{call }{ the call that produced the object}
#' \item{residuals }{ the residuals, that is reponse minus fitted values}
#' \item{fitted.values }{ the estimated response values}
#' 
#' @references 
#'  Friedman J., Hastie T., Hoefling H. and Tibshirani R. (2007), 
#'  Pathwise coordinate optimization. \cr
#'  \emph{The Annals of Applied Statistics}\cr \cr
#'  Friedman J., Hastie T. and Tibshirani R. (2010),
#'  Regularization Paths for Generalized Linear Models via Coordinate Descent. \cr
#'  \emph{Journal of statistical software}\cr \cr
#'  Fu W. J. (1998),  Penalized Regression: The Bridge Versus the Lasso. \cr
#'  \emph{Journal of computational and graphical statistics}\cr \cr
#'  Cheng W. and Wang W. (2014), Graph-regularized dual Lasso for robust eQTL mapping. \cr
#'  \emph{Bioinformatics}
#'
#' @examples
#' X <- matrix(rnorm(100*10),100,10)
#' Y <- matrix(rnorm(100),100,1)
#' G.X <- matrix(rpois(10*10,1),10)
#' G.X <- t(G.X) + G.X
#' diag(G.X) <- 0
#' fit <- edgenet(X=X, Y=Y, G.X=G.X, family="gaussian")
edgenet <-
function
(
  X, Y,
  G.X=NULL, 
  G.Y=NULL, 
  lambda=1, psigx=1, psigy=1, 
  thresh=1e-5,
  maxit=1e5,
  family=c("gaussian")
)
UseMethod("edgenet")

#' @export
#' @noRd
edgenet.default <-
function
(
 X, Y, 
 G.X=NULL, 
 G.Y=NULL,
 lambda=1, psigx=1, psigy=1, 
 thresh=1e-5,
 maxit=1e5, 
 family=c("gaussian")
) 
{
  if (!is.matrix(X)) stop ("X is no matrix!")      
  if (!is.matrix(Y)) stop ("Y is no matrix!")
  # parse dimensions
  n <- dim(X)[1]                              
  p <- dim(X)[2]     
  q <- dim(Y)[2]  
  if (is.null(G.X)) 
  {
    G.X <- matrix(0, 1, 1)
    psigx <- 0
  }
  if (!is.matrix(G.X)) stop("GX is no matrix!")
  if (is.null(G.Y)) 
  {
    G.Y <- matrix(0, 1, 1)
    psigy <- 0 
  }
  if (!is.matrix(G.Y)) stop("GY is no matrix!")
  # check if X and Y are valid
  if (n != dim(Y)[1]) stop("X and Y have not same number of observations!")        
  if (p < 2) stop("Pls use a X matrix with at least 2 covariables!")  
  # check if graphs are valid
  if (psigx != 0 & any(dim(G.X)!=dim(X)[2])) stop("ncol(X) and dim(G.X) do not fit!")
  if (psigy != 0 & any(dim(G.Y)!=dim(Y)[2])) stop("ncol(Y) and dim(G.Y) do not fit!")
  if (lambda < 0) 
  {
    warning("lambda < 0, setting to 0!")
    lambda <- 0
  }
  # TODO
  if (psigx < 0)
  {
    warning("psigx < 0, setting to 0!")
    psigx <- 0
  }
  if (psigy < 0)
  {
    warning("psigy < 0, setting to 0!")
    psigy <- 0
  }
  if (maxit < 0)
  {
    warning("maxit < 0, setting to 1e5!")
    maxit <- 1e5
  }
  if (thresh < 0)
  {
    warning("thresh < 0, setting to 1e-5!")
    thresh <- 1e-5
  }
  if (any(G.X < 0))  stop("Some elements G.X<0; please use non-negative matrix!")
  if (any(G.Y < 0))  stop("Some elements G.Y<0; please use non-negative matrix!")  
  if (all(G.X == 0)) psigx <- 0
  if (all(G.Y == 0)) psigy <- 0
  if (q == 1)        psigy <- 0
  # estimate coefficients
  obj <- fit(X=X, Y=Y, 
             G.X=G.X, G.Y=G.Y,
             lambda=lambda,
             psigx=psigx, psigy=psigy,
             thresh=thresh, maxit=maxit,
             family=family)    
  obj$call <- match.call()    
  class(obj) <- "edgenet"
  obj
}

#' @noRd
fit <-
function
(
  X, Y, 
  G.X, G.Y, 
  lambda, psigx, psigy, 
  thresh, maxit, family
)
{
  # parse dimensions
  n <- dim(X)[1]                              
  p <- dim(X)[2]     
  q <- dim(Y)[2]
  # TODO use family
  family = match.arg(family)
  # make C call to estimate coefficients and posterior networks
  res <- .Call("edgenet", 
               X, Y,
               G.X, G.Y, 
               as.integer(n), as.integer(p), as.integer(q),
               as.double(lambda), 
               as.double(psigx),  as.double(psigy),
               as.integer(maxit), as.double(thresh), 
               PACKAGE="netReg") 
  # coefficients
  B  <- matrix(1, p, q)#res$B
  # intercepts
  mu <- rep(1, q)#res$MU 
  rownames(B) <- colnames(X)
  colnames(B) <- colnames(Y)
  # TODO: calculate the dfs. how to do that in such a model??
  # calc degrees of freedom
  # df <- n - p  
  # if (df < 0) df <- 1
  fitted.values <- as.matrix(X%*%B + intercept.matrix(n=n, mu=mu))
  residuals <- Y - fitted.values
  # standard deviation of residuals would be great to have :)
  # sigma2 <- sum((residuals)^2) / df  
  list(coefficients=B,
       intercept=mu,
       residuals=residuals,
       fitted.values=fitted.values)
}