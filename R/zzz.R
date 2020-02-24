# netReg: graph-regularized linear regression models.
#
# Copyright (C) 2015 - 2019 Simon Dirmeier
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


#' @importFrom reticulate py_module_available
#' @importFrom tensorflow install_tensorflow
.onLoad <- function(libname, pkgname) {
  if (!reticulate::py_module_available("tensorflow")) {
    warning("Could not find tensorflow installation")
  }
  if (!reticulate::py_module_available("tensorflow_probability")) {
    warning("Could not find tensorflow-probability installation")
  }
  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 2)
  options(tensorflow.one_based_extract = TRUE)

  library(tensorflow)
  library(tfprobability)
}


.onUnload <- function(libpath) {
  library.dynam.unload("netReg", libpath)
}
