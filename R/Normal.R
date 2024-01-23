#' @title Normal Estimate
#' @description  Estimation of Mean, Standard Deviation and confidence interval
#' of Normal Distribution
#' @param data Data set/sample
#' @param cL   confidence interval
#'
#' @return estimated mean, estimated standard deviation & confidence interval
#' @export
#'
#' @examples Normal(rnorm(2000,25,4),.95)
Normal <- function(data,cL){

  n <- length(data)

  est_mean <- sum(data)/n
  est_sd <- sqrt(sum((data-est_mean)^2)/n)

  z_value <-qnorm(cL/2,lower.tail = F)
  margin_error <- z_value * est_sd/sqrt(n)

  lower_bound <- est_mean - margin_error
  upper_bound <- est_mean + margin_error
  confidence_interval = c(lower_bound,upper_bound)

  result <- c(Mean = est_mean,
              SD = est_sd,
              ci = confidence_interval)

  return(result)
}
