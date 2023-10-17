#' Generate Data that is Missing At Random (MAR)
#'
#' Simulates missing value at random as NA for a given matrix.
#'
#' @param x a matrix to be used to fill in missing values as NA.
#' @param miss.rate a value of missing rate within the range (0, 1) for variables that contain missing values.
#' @param miss.var proportion of variables (columns) that contain missing values.
#'
#' @return x a matrix with missing values in "NA".
#' @seealso \code{\link{misspi}}
#'
#' @docType package
#'
#' @examples
#'
#' \donttest{
#' set.seed(0)
#' data(toxicity, package = "misspi")
#' toxicity.miss <- missar(toxicity, 0.4, 1)
#' toxicity.miss[1:5, 1:5]}
#'
#' @author Zhongli Jiang \email{jiang548@purdue.edu}
#'
#' @export
missar <- function(x, miss.rate = 0.2, miss.var = 1){
 # browser()
  if(miss.rate<=0 | miss.rate>=1){
    stop("Please control the missing rate to be within interval (0, 1)")
  }
  if(miss.var<=0 | miss.var>1){
    stop("Please control the proportion of variables missing to be within interval (0, 1]")
  }
   if (!is.matrix(x)) {
    x <- as.matrix(x)
   }

  n <- nrow(x)
  p <- ncol(x)

  x.miss.col <- sample(1:p, round(p*miss.var))
  miss.idx <- sample(1:(n*length(x.miss.col)), round((n*length(x.miss.col)*miss.rate)))
  x[, x.miss.col][miss.idx] <- NA
  return(x)
}

