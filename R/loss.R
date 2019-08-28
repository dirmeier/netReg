#' @noRd
#' @import tensorflow
gaussian.loss <- function(y, eta, ...)
{
    obj <- tf$reduce_sum(tf$square(y - eta))
    obj
}


#' @noRd
#' @import tensorflow
binomial.loss <- function(y, eta, ...)
{
    obj <- 0
    for (j in seq(ncol(y))) {
        prob <- tfp$distributions$Bernoulli(logits = eta[,j])
        obj <- obj + tf$reduce_sum(prob$log_prob(y[,j]))
    }

    -obj
}


#' @noRd
#' @import tensorflow
lasso <- function(lambda, beta)
{
    lambda * tf$reduce_sum(tf$abs(beta))
}


#' @noRd
#' @import tensorflow
ridge<- function(lambda, beta)
{
    lambda * tf$reduce_sum(tf$square(beta))
}


#' @noRd
#' @import tensorflow
elastic <- function(alpha, lambda, beta)
{
    lambda * (ridge((1 - alpha) / 2, beta) +
              lasso(alpha, beta))
}



