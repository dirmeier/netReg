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


#' @title Family objects for models
#'
#' @export
#' @docType methods
#' @rdname family-methods
#'
#' @description Family objects provide a convenient way to specify the details
#'  of the models used by \code{netReg}.
#'  See also \code{\link[stats:family]{stats::family}} for more details.
#'
#' @param link  name of a link function
#' @param object  a object for which the family shoulr be retured
#'  (e.g. \code{edgenet})
#' @param ... further arguments passed to methods
#'
#' @return An object of class \code{netReg.family}
#'  \item{family }{ name of the family}
#'  \item{link }{ name of the link function}
#'  \item{linkinv }{ inverse link function}
#'  \item{loss }{ loss function}
#' @examples
#' gaussian()
#' binomial("probit")$link
#' poisson()$linkinv
#' gamma()$linkinv
#' beta()$loss
#' inverse.gaussian()$loss
family <- function(object, ...) UseMethod("family")


#' @export
#' @noRd
#' @method family edgenet
family.edgenet <- function(object, ...) object$family


#' @export
#' @noRd
#' @method family cv.edgenet
family.cv.edgenet <- function(object, ...) family.edgenet(object, ...)


#' @export
#' @rdname family-methods
gaussian <- function(link = c("identity")) {
  link <- match.arg(link)
  linkinv <- switch(
    link,
    "identity" = identity,
    stop("did not recognize link function", call. = FALSE)
  )

  .as.family(
    "gaussian",
    link,
    linkinv,
    gaussian.loss
  )
}


#' @export
#' @rdname family-methods
binomial <- function(link = c("logit", "probit", "log")) {
  link <- match.arg(link)
  linkinv <- switch(
    link,
    "logit" = logistic,
    "log" = exp,
    "probit" = gcdf,
    stop("did not recognize link function", call. = FALSE)
  )

  .as.family(
    "binomial",
    link,
    linkinv,
    binomial.loss
  )
}


#' @export
#' @rdname family-methods
poisson <- function(link = c("log")) {
  link <- match.arg(link)
  linkinv <- switch(
    link,
    "log" = exp,
    stop("did not recognize link function", call. = FALSE)
  )

  .as.family(
    "poisson",
    link,
    linkinv,
    poisson.loss
  )
}


#' @export
#' @rdname family-methods
gamma <- function(link = c("inverse")) {
  warn.experimental("gamma")
  link <- match.arg(link)
  linkinv <- switch(
    link,
    "inverse" = inverse,
    stop("did not recognize link function", call. = FALSE)
  )

  .as.family(
    "gamma",
    link,
    linkinv,
    gamma.loss
  )
}


#' @export
#' @rdname family-methods
beta <- function(link = c("logit", "probit", "log")) {
  warn.experimental("beta")
  link <- match.arg(link)
  linkinv <- switch(
    link,
    "logit" = logistic,
    "log" = exp,
    "probit" = gcdf,
    stop("did not recognize link function", call. = FALSE)
  )

  .as.family(
    "beta",
    link,
    linkinv,
    beta.loss
  )
}


#' @export
#' @rdname family-methods
bum <- function(link = c("logit", "probit", "log")) {
  warn.experimental("bum")
  link <- match.arg(link)
  linkinv <- switch(
    link,
    "logit" = logistic,
    "log" = exp,
    "probit" = gcdf,
    stop("did not recognize link function", call. = FALSE)
  )

  .as.family(
    "bum",
    link,
    linkinv,
    bum.loss
  )
}


#' @export
#' @rdname family-methods
inverse.gaussian <- function(link = c("1/mu^2")) {
  warn.experimental("inverse.gaussian")
  link <- match.arg(link)
  linkinv <- switch(
    link,
    "1/mu^2" = inverse.sqrt,
    stop("did not recognize link function", call. = FALSE)
  )

  .as.family(
    "inverse.gaussian",
    link,
    linkinv,
    inverse.gaussian.loss
  )
}


#' @noRd
.as.family <- function(family, link, linkinv, loss) {
  structure(
    list(
      family = family,
      link = link,
      linkinv = linkinv,
      loss = loss
    ),
    class = "netReg.family"
  )
}
