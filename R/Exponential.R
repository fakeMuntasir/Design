#' @title Exponential Estimate
#' @description
#' Estimation of the parameter & confidence interval from an exponential distribution
#' @param data A Data set/Sample
#' @param cL  Confidence interval
#' @return estimated parameter lamda,confidence interval
#' @examples
#' exponential(rexp(1000,5),.95)
#' exponential(rexp(10000,34),.98)
exponential <- function(data,cL){

  n <- length(data)

  est_lemda <- 1/(sum(data)/n)

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
