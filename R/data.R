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


#' A sample yeast data set for regression
#'
#' The yeast data set is a \code{list} containing three matrices that
#' can be used as an example for using \code{netReg}. The data have been taken
#' from the references listed below.
#'
#' \itemize{
#'     \item \code{X}  (112 x 500)-dimensional binary matrix of 500
#'       genetic markers for 112 yeast samples
#'     \item \code{Y}  (112 x 231)-dimensional double matrix of 231
#'       gene expression values for 112 yeast samples
#'     \item \code{GY} (231 x 231)-dimensional adjaceny matrix representing
#'       protein-protein interactions for 231 yeast genes
#'  }
#'
#' @name yeast
#' @references
#'  Brem, Rachel B., et al. (2005),
#'  Genetic interactions between polymorphisms that affect gene expression in
#'  yeast. \cr
#'  \emph{Nature} \cr \cr
#'  Storey, John D., Joshua M. Akey, and Leonid Kruglyak (2005),
#'  Multiple locus linkage analysis of genomewide expression in yeast. \cr
#'  \emph{PLoS Biology} \cr \cr
#'
#' @docType data
#' @keywords datasets data
#' @usage data(yeast)
#' @format A \code{list} containing three matrices
NA
