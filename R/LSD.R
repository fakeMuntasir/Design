#' Calculates the sum of the squares of the elements of a provided vector.
#'
#' @param x A vector
#'
#' @return Outputs the sum of the squares of the elements of the supplied vector.
#' @export
#'
#' @examples sumsq(c(1:5))


sumsq <- function(x) {
  res <- sum(x^2)
  return(res)
}


#' Performs Latin Square design on a given data set.
#'
#' @param x The provided data set to perform Latin square design on.
#' It is to be remembered that in Latin square design, the number of columns, rows and treatments should be equal, otherwise it'll theoretically be inappropriate.
#'
#' @return Returns the Analysis of variance table with the F-values with respect to the treatment, row factor and column factor.
#' @export
#'
#' @examples LSD(x = #Dataset as tibble / dataframe)



LSD <- function(x) {
  data <- x

  num_of_col <- ncol(data)
  num_of_row <- num_of_col
  N <- num_of_col * num_of_row

  SST <- sum(data^2) - (sum(data))^2/N

  row.names(data) <- letters[seq( from = 1, to = num_of_row )]

  yi..2 <-  tapply (data,
                    rownames <- rownames(data),
                    sum)


  SSrow <- (sumsq(yi..2) / num_of_col) - ((sum(data))^2/N)

  y..k2 <- c(0)
  for(i in 1:5) {
    y..k2[i] <- sum(data[,i])
  }


  SScol <- (sumsq(y..k2) / num_of_row) - ((sum(data))^2/N)

  random_greek_letters <- matrix(data = rep(NA, N), nrow = num_of_col, byrow = T)


  for(j in 1 : num_of_col) {
    random_greek_letters[,j] <- as.matrix(x = sample(1 : num_of_col, replace = F),
                                          ncol = num_of_col,
                                          nrow = num_of_row)

  }


  sum_of_greek_letters <- c(0)

  for (p in 1 : num_of_col){
    for(q in 1 : num_of_col) {
      sum_of_greek_letters[p] <- sum(data[random_greek_letters[p, ], q])
    }

  }

  sum_of_greek_letters

  y.j.2 <- sumsq(sum_of_greek_letters)


  SStreatment <- (y.j.2 / num_of_col) - ((sum(data))^3/N)

  SSE <- SST - SSBlkrow - SSBlkcol - SStreatment


  df_col <- num_of_col - 1
  df_row <- num_of_row - 1
  df_treatmeant <- num_of_col - 1
  df_error <- (num_of_col - 2) * (num_of_row - 1)
  df_total <- num_of_col^2 - 1

  ms_col <- SSBlkcol / df_col
  ms_row <- SSBlkrow / df_row
  ms_treatment <- SStreatment / df_treatmeant
  ms_error <- SSE / df_error
  ms_total <- SST / df_total


  F_treat <- ms_treatment / ms_error
  F_row <- ms_row / ms_error
  F_col <- ms_col / ms_error

  output <- matrix(
    c(SStreatment, SSBlkrow, SSBlkcol, SSE, SST, df_treatmeant, df_row, df_col, df_error, df_total, ms_treatment, ms_row, ms_col, ms_error,ms_total, F_treat, F_row, F_col, NA, NA),
    ncol = 4, byrow = F
  )

  colnames(output) <- c("Sum of Squares", "Degrees of freedom", "Mean Square", "F_cal values")

  rownames(output) <- c("Treatment", "Row Factor", "Column_factor", "Error", "Total")

  output <- as.table(output)

  return(output)
}




