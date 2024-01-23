
#' This function performs Complete Randomize Design on a given data set.
#'
#' @param a Parameter a is the number of treatment levels for the given experiment.
#' @param n Parameter n is the total number of observation/response for the given experiment.
#' @param response Response is a vector of length n to perform the CRD on.
#' @param alp Parameter alp is the level of significance. This is required for confidence intervals and drawing conclusions.
#'
#' @return This function returns:
#' 1. The resulting ANOVA table for CRD,
#' 2. Decision i.e. whether to reject or accept Null hypothesis,
#' 3. Result of the tukey's test,
#' 4. Confidence interval for every treatment mean,
#' 5. Limits of the Bonferroni's test.
#'
#' @export
#' @names crd
#'
#' @examples



library(tidyverse)
library(dplyr)


crd <- function(a, n, response, alp){



  treatments <- LETTERS[1:a] # Assuming there are 5 treatments



  N <- length(response)

  dat <- as.tibble(matrix(response,nrow = a,ncol = b))
  dat1 <-dat %>%
    mutate(treat = treatments,
           yi_dot = rowSums(dat),
           yi_dot_bar = yi_dot / n ) %>%
    select(treat, everything())
  dat1

  ggplot(data = dat)+
    geom_boxplot(
      mapping = aes ( x = rep(treatments, times = n),
                      y = response)
    )

  cf <- sum(response)^2 / N
  cf

  sst <- sum( response^2 )-cf
  sst

  sstreat <- sum( (dat1$yi_dot)^2 ) /b -cf
  sstreat
  mstreat <- sstreat / (a-1)

  sse <- sst - sstreat
  sse
  mse <- sse / (N-a)

  f_cal <- mstreat / mse
  f_tab <-  qf(0.05, (a-1), (N-a), lower.tail = FALSE)


  dec_crd<-ifelse(abs(f_cal) > f_tab , "Reject H0.", "Fail to reject H0. ")


  dat2<-matrix(c(sstreat, sse, sst, a-1, N-a, N-1, mstreat, mse, NA, f_cal, NA, NA), ncol=4, byrow=F)
  colnames(dat2) = c('sumsq','df','mean sq','F0')
  rownames(dat2) = c('treat','error','total')
  ano_crd<-as.table(dat2)

  u <- dat1$yi_dot_bar
  k<-1:6
  m<-1

  for(i in 1:3){
    for(j in (i+1):4){
      k[m]<-u[i]-u[j]
      m <-m+1
    }
  }
  k

  q <- qtukey(0.95,4,16)

  dif <- q*(sqrt(mse/b))
  result_tukey <- LETTERS[1:6]

  for (i in 1:6) {
    result_tukey[i] <-  ifelse(abs(k[i]) > dif, "differ", "don't differ")
  }

  ##ci

  se <- qt((1-(alp/2)),(N-a))*sqrt(mse/n)
  ll <- u-se
  ul <- u + se
  ci <- tibble(lower_limit=ll, upper_limit=ul)

  ##bonferonni
  alp1<-alp/(2*n)
  se_bon <- qt((1-(alp1/2)),(N-a))*sqrt(mse/n)
  ll_bon <- u - se_bon
  ul_bon <- u + se_bon
  ci_bon <- tibble(ll_bon=ll_bon, ul_bon=ul_bon)

  return(list(dat1, ano_crd, dec_crd, result_tukey, ci, ci_bon))
}
