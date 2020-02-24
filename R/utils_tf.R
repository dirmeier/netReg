# netReg: graph-regularized linear regression models.
#
# Copyright (C) 2015 - 202 0Simon Dirmeier
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


#' @noRd
#' @import tensorflow
init_variables <- function() {
  tf$compat$v1$global_variables_initializer()
}


#' @noRd
#' @import tensorflow
adam <- function(learning.rate) {
  tf$compat$v1$train$AdamOptimizer(learning_rate = learning.rate)
}


#' @noRd
#' @import tensorflow
session <- function() {
  tf$compat$v1$Session()
}


#' @noRd
#' @import tensorflow
reset_graph <- function() {
  tensorflow::tf$compat$v1$reset_default_graph()
}


#' @noRd
#' @import tensorflow
cast_float <- function(x) {
  tensorflow::tf$cast(x, tensorflow::tf$float32)
}


#' @noRd
#' @import tensorflow
constant_float <- function(x) {
  tensorflow::tf$constant(x, tensorflow::tf$float32)
}


#' @noRd
#' @import tensorflow
placeholder <- function(shape, name = NULL) {
  if (!is.null(name)) {
    tensorflow::tf$compat$v1$placeholder(
      tensorflow::tf$float32, shape,
      name = name
    )
  } else {
    tensorflow::tf$compat$v1$placeholder(
      tensorflow::tf$float32, shape
    )
  }
}


#' @noRd
#' @import tensorflow
zero_matrix <- function(m, n) {
  tensorflow::tf$Variable(tensorflow::tf$zeros(shape(m, n)))
}


#' @noRd
#' @import tensorflow
zero_vector <- function(m) {
  tensorflow::tf$Variable(tensorflow::tf$zeros(shape(m)))
}
