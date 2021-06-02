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


#' @noRd
#' @importFrom tensorflow tf
gaussian.loss <- function(y, mu.hat, ...) {
  obj <- tf$reduce_mean(tf$square(y - mu.hat))
  obj
}


#' @noRd
#' @importFrom tensorflow tf
#' @importFrom tfprobability tfd_bernoulli
binomial.loss <- function(y, mu.hat, ...) {
  obj <- 0
  for (j in seq(ncol(y))) {
    prob <- tfprobability::tfd_bernoulli(probs = mu.hat[, j])
    obj <- obj + tf$reduce_sum(prob$log_prob(y[, j]))
  }

  -obj
}


#' @noRd
#' @importFrom tensorflow tf
#' @importFrom tfprobability tfd_poisson
poisson.loss <- function(y, mu.hat, ...) {
  obj <- 0
  for (j in seq(ncol(y))) {
    prob <- tfprobability::tfd_poisson(rate = mu.hat[, j])
    obj <- obj + tf$reduce_sum(prob$log_prob(y[, j]))
  }

  -obj
}


#' @noRd
#' @importFrom tensorflow tf
#' @importFrom tfprobability tfd_gamma
gamma.loss <- function(y, mu.hat, ...) {
  obj <- 0
  for (j in seq(ncol(y))) {
    prob <- tfprobability::tfd_gamma(rate = mu.hat[, j], concentration = 1)
    obj <- obj + tf$reduce_sum(prob$log_prob(y[, j]))
  }

  -obj
}


#' @noRd
#' @importFrom tensorflow tf
#' @importFrom tfprobability tfd_beta
beta.loss <- function(y, mu.hat, ...) {
  obj <- 0
  eps <- .Machine$double.eps * 1e9
  for (j in seq(ncol(y))) {
    mu <- mu.hat[, j]
    phi <- 1 # TODO: replace this with tf$variable

    # reparametrize
    # concentration1 := alpha = mu * phi
    p <- mu * phi
    # concentration0 := beta = (1. - mu) * phi
    q <- (1 - mu) * phi

    # deal with numerical instabilities
    p.trans <- tf$math$maximum(p, eps)
    q.trans <- tf$math$maximum(q, eps)

    prob <- tfprobability::tfd_beta(
      concentration1 = p.trans, concentration0 = q.trans
    )
    obj <- obj + tf$reduce_sum(prob$log_prob(y[, j]))
  }

  -obj
}


#' @noRd
#' @importFrom tensorflow tf
#' @importFrom tfprobability tfd_beta tfd_uniform tfd_mixture tfd_categorical
bum.loss <- function(y, mu.hat, ...) {
  obj <- 0
  eps <- .Machine$double.eps * 1e9
  N <- mu.hat$shape[[1]]
  for (j in seq(ncol(y))) {
    mu <- mu.hat[, j]
    mixture_weight <- tf$constant(matrix(0.75, nrow = N, ncol = 2), dtype = tf$float32) # TODO: learn this parameter
    phi <- 1 # TODO: replace this with tf$variable

    # reparametrize
    # concentration1 := alpha = mu * phi
    p <- mu * phi
    # concentration0 := beta = (1. - mu) * phi
    q <- (1 - mu) * phi

    # deal with numerical instabilities
    p.trans <- tf$math$maximum(p, eps)
    q.trans <- tf$math$maximum(q, eps)

    # need correct batch dimensions for mixture
    p.trans <- tf$stack(list(p.trans, p.trans), 0L)
    q.trans <- tf$stack(list(q.trans, q.trans), 0L)

    prob <- tfprobability::tfd_mixture(
      cat = tfprobability::tfd_categorical(
        probs = c(mixture_weight, 1 - mixture_weight)
      ),
      components = c(
        tfprobability::tfd_beta(
          concentration1 = p.trans, concentration0 = q.trans
        ),
        tfprobability::tfd_uniform(
          low = tf$constant(tf$Variable(matrix(0, nrow = 2, ncol = N), dtype=tf$float32)),
          high = tf$constant(tf$Variable(matrix(1, nrow = 2, ncol = N), dtype=tf$float32))
        )
      )
    )

    obj <- obj + tf$reduce_sum(prob$log_prob(y[, j]))
  }

  -obj
}


#' @noRd
#' @importFrom tensorflow tf
#' @importFrom tfprobability tfd_inverse_gaussian
inverse.gaussian.loss <- function(y, mu.hat, ...) {
  obj <- 0
  for (j in seq(ncol(y))) {
    # loc := mu
    # concentration := lambda (shape)
    prob <- tfprobability::tfd_inverse_gaussian(loc = mu.hat[, j], concentration = 1)
    obj <- obj + tf$reduce_sum(prob$log_prob(y[, j]))
  }

  -obj
}
