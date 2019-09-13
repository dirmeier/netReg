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


<<<<<<< HEAD
<<<<<<< HEAD
#' @method family cv.edgenet
family.edgenet <- function(object, ...) object$family


#' @method family cv.edgenet
family.cv.edgenet <- function(object, ...) family.edgenet(object, ...)


=======
>>>>>>> 2ced4bf... Add family and link functions
=======


>>>>>>> bc8280c... Add binomial
#' @title Family objects for models
#'
#' @export
#' @docType methods
#' @rdname family-methods
#'
<<<<<<< HEAD
#' @description Family objects provide a convenient way to specify the details
#'  of the models used by \code{netReg}.
=======
#' @export
#'
#' @description Family objects provide a convenient way to specify the details
#'  of the models used by \code{netReg}}.
>>>>>>> 2ced4bf... Add family and link functions
#'  See also \code{\link[stats:family]{stats::family}} for more details.
#'
#' @param link  name of the link function
#' @param object  the function \code{family} accesses the family objects
#'  which are stored within objects created by modelling functions
#'   (e.g., \code{\link{edgenet}} or \code{\link{cv.edgenet}})
#'
#' @return An object of class \code{netReg.family}
#'  \item{family }{ name of the family}
#'  \item{link }{ name of the link function}
#'  \item{linkinv }{ inverse link function}
<<<<<<< HEAD
#'  \item{loss }{ loss function}
gaussian <- function(link = c("identity"))
{
    link <- match.arg(link)
    linkinv <- switch(
        link,
        "identity"=identity,
        stop("did not recognize link function", call. = FALSE)
    )

    .as.family("gaussian",
               link,
               linkinv,
               gaussian.loss)
}


#' @export
#' @rdname family-methods
binomial <- function(link=c("logit", "probit", "log"))
{
    link <- match.arg(link)
    linkinv <- switch(
        link,
        "logit"=logistic,
        "log"=exp,
        "probit"=gcdf,
        stop("did not recognize link function", call. = FALSE)
    )

    .as.family("binomial",
               link,
               linkinv,
               binomial.loss)
}


#' @export
#' @rdname family-methods
poisson <- function(link=c("log"))
{
    link <- match.arg(link)
    linkinv <- switch(
        link,
        "log"=exp,
        stop("did not recognize link function", call. = FALSE)
    )

    .as.family("poisson",
               link,
               linkinv,
               poisson.loss)
}

inverse.gaussian <- function() stop("not implemented")
beta <- function()  stop("not implemented")
gamma <- function()  stop("not implemented")


#' @noRd
.as.family <- function(family, link, linkinv, loss)
{
    structure(
        list(family=family,
             link=link,
             linkinv=linkinv,
             loss=loss),
        class="netReg.family"
    )
}


=======
#'
family.edgenet <- function(object, ...)

#' @rdname family-methods
gaussian <- function(link = c("identity", "log", "inverse"))
{
    link <- match.arg(link)
    linkinv <- switch(
        link,
        "identity"=identity,
        "log"=exp,
        "inverse"=inverse,
        stop("did not recognize link function", call. = FALSE)
    )

    .as.family(as.character(match.call()[[1]]), link, linkinv)
}


#' @rdname family-methods
binomial <- function(link=c("logit", "probit", "log"))
{
    link <- match.arg(link)
    linkinv <- switch(
        link,
        "logit"=logistic,
        "log"=exp,
        "probit"=gcdf,
        stop("did not recognize link function", call. = FALSE)
    )

    .as.family(as.character(match.call()[[1]]), link, linkinv)
}


inverse.gaussian <- function() stop("not implemented")
beta <- function()  stop("not implemented")
poisson <- function()  stop("not implemented")
gamma <- function()  stop("not implemented")
<<<<<<< HEAD
>>>>>>> 2ced4bf... Add family and link functions
=======


#' @noRd
.as.family <- function(family, link, linkinv)
{
    structure(
        list(family=family,
             link=link,
             linkinv=linkinv),
        class="netReg.family"
    )
}


>>>>>>> bc8280c... Add binomial
