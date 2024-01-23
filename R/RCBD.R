
#' This function performs Randomized complete block design on a given data set.

#' @param data Data is the data frame to perform RCBD on. Input should be a data frame that contains 3 columns named - raw_data, treat & block

#' @param alp Alpha is the level of significance. This is required to compute the confidence level and drawing conclusions.
#'
#' @return This function returns the ANOVA table, Decision i.e. to accept or reject the null hypothesis, and two plots : comparative box-plot and  Normal Q-Q plot.

#' @export
#'
#' @examples rcbd(data = #Dataset as tibble / dataframe, alp = # significance level)
#' names rcbd

library(tidyverse)

rcbd <- function(data, alp){

  N <- length(data$raw_data)
  a <- length(levels(treat))
  b <- length(levels(block))
  #sum(dat1$treat)
  #sum(raw_data)

  cf <- (sum(data$raw_data))^2/N

  sst <- sum(data$raw_data^2)-cf
  sst
  sstreat <- sum((tapply(data$raw_data, data$treat, sum))^2/b)-cf
  #tapply(dat1, treat, sum )
  ssblocks <- sum((tapply(data$raw_data, data$block, sum))^2/a)-cf
  sse <- sst-sstreat-ssblocks

  mstreat <- sstreat/(a-1)
  msblocks <- ssblocks/(b-1)
  mse <- sse/((a-1)*(b-1))

  f_cal <- mstreat/mse
  f_tab <- qf(1-(alp/2), a-1, N-a)
  dec <- ifelse(f_cal>f_tab, 'H0 REJECT','H0 CANT BE REJECTED')

  dat2 <- matrix(c(sstreat,ssblocks,sse, sst, a-1,b-1,(a-1)*(b-1),N-1, mstreat, msblocks,mse, NA, f_cal, NA, NA, NA),
                 ncol=4, byrow=F)

  colnames(dat2) = c('sumsq','df','mean sq','F0')
  rownames(dat2) <- c('treat','blocks','error','total')
  ano1 <- as.table(dat2)

  plot1<-ggplot(data=data) +
    geom_boxplot(
      mapping = aes(x = as.factor(treat),
                    y = raw_data))+
    labs(
      x = "Treatments",
      y = "Response",
      title = "Comparative Box plot")

  error<-function(x)
  {
    return(x-mean(x))
  }

  errors<-as.vector(sapply(split(data$raw_data,data$treat),error))
  plot2<-qqnorm(errors)

  return(list(ano1, dec, plot1, plot2))

}
