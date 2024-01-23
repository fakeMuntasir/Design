library(tidyverse)

#' Test of hypothesis which determines the p-value,confidence interval and test
#' statistic given the arguments
#'
#' @param parameter argument to determine whether the parameter is known
#' @param independent For matched pair data set, this parameter is false, and otherwise true.
#' @param x A random variable from a normal distribution
#' @param y Another random variable from a distribution
#' @param mu_0 mean of a specific distribution
#' @param sigma standard of a specific distribution
#' @param n1  number of rows matrix
#' @param n2 number of columns matrix
#' @param B number of replication
#' @param alpha level of significance
#' @param H1 alternate hypothesis
#'
#' @return p-value,confidence-interval,test-statistic
#'
#' @examples test_of_hypothesis(parameter = "known", independent = FALSE, x = x,
#' y = y, mu_0 = mu_0, alpha = alpha, H1 = H1 )



test_of_hypothesis<-function(parameter,independent=NULL,x,y=NULL,mu_0=NULL,sigma=NULL,
                             n1=NULL,n2=NULL,B=NULL,alpha,H1){

  if(parameter == "known"){



    if( is.vector(x) & !is.character(x) & !is.null(mu_0)){

      mu=mean(x)
      n=length(x)

      if(!is.null(sigma) & n>30){

        z_cal=((mu-mu_0)*sqrt(n))/sigma
        z_tab<-qnorm(alpha,0,1,F)
        z_tab2<-qnorm(alpha/2,0,1,F)


        b<-function(H1){
          if(H1=="mu>mu_0"){
            p=pnorm(z_cal,lower.tail = F)
            lo=NULL
            up=mu+((z_tab*sigma)/sqrt(n))

          }
          else if(H1=="mu<mu_0"){
            p=pnorm(z_cal,lower.tail = T)
            lo=mu-((z_tab*sigma)/sqrt(n))
            up=NULL


          }
          else if(H1== "mu !=mu_0"){
            p=2*(1-pnorm(abs(z_cal)))
            lo=mu-((z_tab2*sigma)/sqrt(n))
            up=mu+((z_tab2*sigma)/sqrt(n))


          }
          H0<-"mu=mu_0"

          return(c(p_value=p,alpha=alpha,lower_bound=lo,upper_bound=up ,Test_statistics=z_cal,H0=H0,H1=H1))
        }
      }
      else if(n<30){
        s=sqrt((sum((x-mu)^2))/(n-1))

        t_cal=((mu-mu_0)*sqrt(n))/s
        t_tab<-qt((alpha),(n-1),lower.tail = F)
        t_tab_2<-qt((alpha/2),(n-1),lower.tail = F)

        b<-function(H1){
          if(H1=="mu>mu_0"){
            p=pt(t_cal,(n-1),lower.tail = F)
            lo=NULL
            up=mu+((t_tab*s)/sqrt(n))


          }
          else if(H1=="mu<mu_0"){
            p=pt(t_cal,(n-1),lower.tail = T)
            lo=mu-((t_tab*s)/sqrt(n))
            up=NULL


          }
          else if(H1== "mu !=mu_0"){
            p=2*pt(-abs(t_cal),(n-1))
            lo=mu-((t_tab_2*s)/sqrt(n))
            up=mu+((t_tab_2*s)/sqrt(n))


          }
          H0<-"mu=mu_0"
          return(c(p_value=p,alpha=alpha,lower_bound=lo,upper_bound=up,Test_statistics=t_cal,H0=H0,H1=H1))
        }
      }
      return(b(H1))

    }


    if(!is.null(independent)  &  !is.null(y)){


      if(independent == "TRUE" ){
        xbar <- mean(x)
        ybar <- mean(y)
        n <- length(x)
        m <- length(y)
        sx <- sd(x)
        sy <- sd(y)
        sp <- sqrt(((n-1)*sx^2+(m-1)*sy^2)/(n+m-2))
        t_statistic <- (xbar - ybar) / sp*(sqrt(1/n+1/m))
        ttab <- qt((alpha/2), (n+m-2),lower.tail = F)

        p_value<- 2*(pt(-abs(t_statistic),(n+m-2)))
        lo<-(xbar-ybar)-(ttab*sp*(sqrt(1/n+1/m)))
        up<-(xbar-ybar)+(ttab*sp*(sqrt(1/n+1/m)))


        result <- c(t_statistic = t_statistic, p_value = p_value,alpha=alpha,lower_bound=lo,upper_bound=up)
      }
      else if(independent == "FALSE")
      {
        differences <- (y -x)
        n <- length(differences)

        mean_diff <- mean(differences)
        sd_diff <- sqrt(sum((differences - mean_diff)^2) / (n - 1))

        t_statistic <- mean_diff / (sd_diff / sqrt(n))

        df <- n - 1

        p_value <- 2 * pt(-abs(t_statistic), df)

        t_tab<- qt((alpha/2),df,lower.tail = F)

        lo=mean_diff-((t_tab*sd_diff)/sqrt(n))
        up=mean_diff+((t_tab*sd_diff)/sqrt(n))


        H0<-"mu=mu_0"

        result <- c(t_statistic = t_statistic, p_value = p_value,alpha=alpha,lower_bound=lo,upper_bound=up,H0=H0,H1=H1)
      }
      return(result)
    }



    if(is.character(x) & !is.null(y) & !is.null(n1) & !is.null(n2)){



      datt <- data.frame(x,y )

      N <- n1*n2

      constant <- sum(y)^2 / N

      sst <- sum(y^2) - constant

      yi. <- sapply(split(y,x), sum )
      sstreat <- (sum(yi.^2)/n2) - constant

      sse <- sst - sstreat

      mstreat <- sstreat / (n1-1)
      mse <- sse / (N-n1)

      f_cal <- mstreat / mse

      f_tab <- qf((alpha), n1-1, N-n1, lower.tail = F)



      H0 <- "Treatment means are equal"

      pvalue <-pf(f_cal,n1-1, N-n1,lower.tail = F)
      t_tab<-qt((alpha/2),(N-n1),lower.tail = F)
      t_tab
      val<-t_tab*(sqrt(mse/n2))
      val
      a<-datt %>%
        group_by(x) %>%
        reframe(CI_lower=(mean(y)-val),CI_upper=(mean(y)+val))

      e<-list(a)


      return(c(H0= H0, H1 = H1, F_statistic =  f_cal,p_value=pvalue,alpha=alpha,e))
    }


    if(is.matrix(x) & !is.null(n1)  & !is.null(n2)){



      df=(n1-1)*(n2-1)



      a<-NULL
      b<-NULL
      for(i in 1:n1){
        a[i]<-sum(x[i,])
      }

      for(j in 1:n2){
        b[j]<-sum(x[,j])
      }


      eij <- matrix(0, nrow = n1, ncol = n2)
      for(i in 1:n1){
        for(j in 1:n2){
          eij[i,j]<-a[i]*b[j]/sum(x)
        }
      }

      q<-0

      for(i in 1:n1){
        for(j in 1:n2){
          q=q+sum((x[i,j]-eij[i,j])^2/eij[i,j])
        }
      }
      pvalue<-pchisq(q,df,lower.tail = F)

      H0<-"There is no association"

      return(c(Test_statistic=q,p_value=pvalue,alpha=alpha,H0=H0,H1=H1))

    }
  }



  if(parameter == "unknown"){

    if( !is.null(y) & !is.null(B)){

      boot_stat_vec <- c()

      for(i in 1:B){

        dat <- data.frame(y, x)

        index <- sample(1:length(x), size = length(x), replace = T)

        dat_boot <- dat[index,]

        boot_stat_vec[i] <-  dat_boot %>%
          group_by(x) %>%
          summarise(mean = mean(y)) %>%
          summarise(test_stat = mean[x == 1] - mean[x == 0]) %>%
          pull()
      }
      test_stat_main <- dat %>%
        group_by(x) %>%
        summarise(mean = mean(y)) %>%
        summarise(test_stat = mean[x == 1] - mean[x == 0]) %>%
        pull()

      p_val = (sum(boot_stat_vec >= test_stat_main))/B

      H0<-"There is no treatment effect"

      return(c(p_val=p_val,alpha=alpha,H0=H0,H1=H1))

    }
  }

}
