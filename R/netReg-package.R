# netReg: graph-regularized linear regression models.
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of netReg.
#
# netReg is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# netReg is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with netReg. If not, see <http://www.gnu.org/licenses/>.


#' netReg
#'
#' \emph{netReg} is a package for generalized linear regression that includes 
#' prior graphs in the models objective function.
#' 
#' \emph{netReg} uses \emph{Armadillo}, \emph{OpenBLAS}, 
#' \emph{BLAS} and \emph{LAPACK} for fast matrix computations and
#' \emph{Dlib} for constrained derivate-free optimization.
#' 
#' @name netReg-package
#' @author Simon Dirmeier | \email{mail@@simon-dirmeier.net}
#' @docType package
#' @keywords package
#' 
#' @useDynLib netReg
#' 
#' @references 
#'  Friedman J., Hastie T., Hoefling H. and Tibshirani R. (2007), 
#'  Pathwise coordinate optimization.\cr
#'  \emph{The Annals of Applied Statistics}\cr \cr
#'  Friedman J., Hastie T. and Tibshirani R. (2010),
#'  Regularization Paths for Generalized Linear Models 
#'   via Coordinate Descent. \cr
#'  \emph{Journal of Statistical Software}\cr \cr
#'  Fu W. J. (1998),  Penalized Regression: The Bridge Versus the Lasso.\cr
#'  \emph{Journal of Computational and Graphical Statistics}\cr \cr
#'  Cheng W. and Wang W. (2014), Graph-regularized dual Lasso 
#'   for robust eQTL mapping.\cr
#'  \emph{Bioinformatics}\cr \cr
#'  Powell M.J.D. (2009), 
#'  The BOBYQA algorithm for bound constrained optimization
#'   without derivatives.\cr
#'  \url{http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf}
NULL
