#' Find the optimal shrinkage parameters for edgenet
#' 
#' @export
#' @useDynLib netReg
#' 
#' @author Simon Dirmeier, \email{netreg@@simon-dirmeier.net}
#' 
#' @description Finds the optimal shrinkage parameters using cross-validation for edgenet. 
#' We use the BOBYQA algorithm to minimize the sum of squared residuals objective function.
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
#' @param thresh  threshold for coordinate descent
#' @param maxit  maximum number of iterations
#' @param family  family of response, e.g. gaussian
#' @param nfolds  the number of folds to be used - default is 10 (minimum 3, maximum nrow(X)). 
#' @param foldid  an optional vector of length \code{nrow(X)} of values between 1 and \code{nfold} identifying what fold each observation is in.
#' @param ...  additional parameters

#' @return An object of class \code{cv.edgenet}
#' \item{call }{ the call that produced the object}
#' \item{lambda }{ the estimated (\code{p} x \code{q})-dimensional coefficient matrix B.hat}
#' \item{psigx }{ the estimated (\code{q} x \code{1})-dimensional vector of intercepts}
#' \item{psigy }{ the estimated (\code{q} x \code{1})-dimensional vector of intercepts}
#' \item{foldid }{ the vector of fold assignments yoused}
#' 
#' @references 
#'  Friedman J., Hastie T., Hoefling H. and Tibshirani R. (2007), 
#'  Pathwise coordinate optimization.\cr
#'  \emph{The Annals of Applied Statistics}\cr \cr
#'  Friedman J., Hastie T. and Tibshirani R. (2010),
#'  Regularization Paths for Generalized Linear Models via Coordinate Descent. \cr
#'  \emph{Journal of Statistical Software}\cr \cr
#'  Fu W. J. (1998),  Penalized Regression: The Bridge Versus the Lasso.\cr
#'  \emph{Journal of Computational and Graphical Statistics}\cr \cr
#'  Cheng W. and Wang W. (2014), Graph-regularized dual Lasso for robust eQTL mapping.\cr
#'  \emph{Bioinformatics}\cr \cr
#'  Powell M.J.D. (2009), 
#'  The BOBYQA algorithm for bound constrained optimization without derivatives.\cr
#'  \url{http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf}
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100*10),100,10)
#' Y <- matrix(rnorm(100),100,1)
#' G.X <- matrix(rpois(10*10,1),10)
#' G.X <- t(G.X) + G.X
#' diag(G.X) <- 0
#' cv.edge <- cv.edgenet(X=X, Y=Y, G.X=G.X, family="gaussian")
#' }
cv.edgenet <- function (X, Y, G.X=NULL, G.Y=NULL, 
                        thresh=1e-5, maxit=1e5, family=c("gaussian"),
                        nfolds=10, foldid=NULL, ...)
  UseMethod("cv.edgenet")

#' @export
#' @noRd
cv.edgenet.default <- function (X, Y, G.X=NULL, G.Y=NULL, 
                                thresh=1e-5, maxit=1e5, family=c("gaussian"),
                                nfolds=10, foldid=NULL, ...) 
{
  if (!is.matrix(X)) stop ("X is no matrix!")      
  if (!is.matrix(Y)) stop ("Y is no matrix!")
  # parse dimensions
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  if (is.null(G.X)) G.X <- matrix(0, 1, 1)
  if (!is.matrix(G.X)) stop("GX is no matrix!")
  if (is.null(G.Y)) G.Y <- matrix(0, 1, 1)
  if (!is.matrix(G.Y)) stop("GY is no matrix!")
  # check if X and Y are valid
  if (n != dim(Y)[1]) stop("X and Y have not same number of observations!")        
  if (p < 2) stop("Pls use a X matrix with at least 2 covariables!")  
  psigx <- psigy <- -1
  # check if graphs are valid
  if (all(G.X == 0)) psigx <- 0
  if (all(G.Y == 0)) psigy <- 0
  if (psigx != 0 & any(dim(G.X)!=dim(X)[2])) stop("ncol(X) and dim(G.X) do not fit!")
  if (psigy != 0 & any(dim(G.Y)!=dim(Y)[2])) stop("ncol(Y) and dim(G.Y) do not fit!")
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
  # check if some parameters have values
  if (!is.null(foldid) & is.numeric(foldid)) 
  {
    nfolds <- max(foldid)
    if (length(foldid) != n) stop("length(foldid) != sample size (n)!")
  }
  if (is.null(foldid)) foldid <- NA_integer_
  if (!is.numeric(foldid)) stop("Please provide either an integer vector or NULL for foldid")
  psigx <- psigy <- -1
  if (all(G.X == 0)) psigx <- 0
  if (all(G.Y == 0)) psigy <- 0
  if (q == 1)        psigy <- 0
  class(G.X) <- class(G.Y) <- "numeric"
  if (n < nfolds) nfolds <- n
  family <-  match.arg(family)
  # estimate shrinkage parameters
  ret <- .cv(X=X, Y=Y, 
             G.X=G.X, G.Y=G.Y,
             psigx=psigx, psigy=psigy,
             thresh=thresh, maxit=maxit,
             family=family,
             nfolds=nfolds,
             foldid=foldid)    
  ret$call <- match.call()
  class(ret) <- c(class(ret), "cv.edgenet")
  ret
}
