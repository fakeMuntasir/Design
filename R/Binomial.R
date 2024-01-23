#' @title Binomial Estimate
#' @description
#' Estimates of parameter p of a binomial distribution & its confidence interval
#' n always 1.
#' @param data Data set/sample from distribution
#' @param CL confidence interval
#' @param trials number of trials
#'
#' @return estimate of probability, confidence interval
#' @export
#'
#' @examples
#' Binomial(rbinom(1,10,0.49),.95,10)
#' Binomial(rbinom(1,20,0.7),.95,20)
Binomial <- function(data,cL,trials){


  #n <- length(data)

  est_p <- data/trials#(sum(data)/n)/trials

  z_value <-qnorm(cL/2,lower.tail = F)
  margin_error <- z_value * sqrt(est_p*(1-est_p)/trials)
  lower_bound <- est_p - margin_error
  upper_bound <- est_p + margin_error
  confidence_interval = c(lower_bound,upper_bound)

  result <- c(p = est_p,
              ci = confidence_interval
  )

  return(result)

}
