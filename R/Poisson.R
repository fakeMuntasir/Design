#' @title Poisson Estimate
#' @description
#' Estimation of the parameter & confidence interval from the sample
#' of a Poisson Distribution
#'
#' @param data Data set/sample
#' @param cL confidence interval
#'
#' @return estimated parameter lambda, confidence interval
#' @export
#'
#' @examples
#' poisson(rpois(100,.5),.95)
#' poisson(rpois(1000,.45),.95)
Poisson <- function(data,cL){

  n <- length(data)
  est_lemda <- sum(data)/n

  z_value <-qnorm(cL/2,lower.tail = F)
  margin_error <- z_value * sqrt(est_lemda/n)
  lower_bound <- est_lemda - margin_error
  upper_bound <- est_lemda + margin_error
  confidence_interval = c(lower_bound,upper_bound)

  result <- list(
    Lemda = est_lemda,
    ci = confidence_interval
  )
  return(result)

}
