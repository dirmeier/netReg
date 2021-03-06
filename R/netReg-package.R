# netReg: network-regularized linear regression models.
#
# Copyright (C) 2015 - 2020 Simon Dirmeier
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


#' Network-Regularized Regression Models
#'
#' \emph{netReg} fits linear regression models using network-regularization.
#' Network-regularization used graph prior as penality for the coefficients of linear models.
#' The graph can represent any relationship between the covariables/responses of the
#' model, for instance, some quantifiable biological relationship such as coexpression.
#'
#' \emph{netReg} uses \emph{Armadillo} and \emph{TensorFlow} for
#' fast matrix computations and optimization.
#'
#' @name netReg-package
#' @author Simon Dirmeier | \email{simon.dirmeier@@web.de}
#' @docType package
#' @keywords package
#'
#' @useDynLib netReg, .registration = TRUE
#'
#' @import methods
#' @import tensorflow
#' @import tfprobability
#' @importFrom tensorflow tf
#'
#' @references
#'  Dirmeier, Simon and Fuchs, Christiane and Mueller, Nikola S and Theis,
#'  Fabian J (2018),
#'  netReg: Network-regularized linear models for biological association
#'  studies. \cr
#'  \emph{Bioinformatics}\cr \cr
#'  Abadi, Martín et al. (2016),
#'  Tensorflow: A system for large-scale machine learning. \cr
#'  \emph{12th USENIX Symposium on Operating Systems Design and Implementation (OSDI 16)}\cr \cr
#'  Powell M.J.D. (2009),
#'  The BOBYQA algorithm for bound constrained optimization
#'   without derivatives.\cr
#'  \url{http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf}
#'  Eddelbuettel, Dirk and Sanderson, Conrad (2014),
#'  RcppArmadillo: Accelerating R with high-performance C++ linear algebra.
#'  \emph{Computational Statistics & Data Analysis}
NULL
