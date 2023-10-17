#' Evaluate the Imputation Quality
#'
#' Calculates Root Mean Squared Error (RMSE), Mean Absolute Error (MAE) and Normalized Root Mean Squared Error (NRMSE). It also performs visualization for imputation quality evaluation.
#'
#' @param x.true a vector with true values.
#' @param x.impute a vector with estimated values.
#' @param plot a Boolean that indicates whether to plot or not.
#' @param interactive a Boolean that indicates whether to use interactive plot when the plot option is invoked (plot = "TRUE").
#'
#' @import ggplot2
#' @importFrom plotly ggplotly
#'
#' @return rmse root mean squared error.
#' @return mae mean absolute error.
#' @return nrmse normalized root mean squared error.
#'
#' @seealso
#' \code{\link{misspi}}, \code{\link{missar}}
#'
#' @examples
#' # A very quick example
#' n <- 100
#' x.true <- rnorm(n)
#' x.est <- x.true
#' na.idx <- sample(1:n, 20)
#' x.est[na.idx] <- x.est[na.idx] + rnorm(length(na.idx), sd = 0.1)
#'
#' # Default plot
#' er.eval <- evaliq(x.true[na.idx], x.est[na.idx])
#'
#' # Interactive plot
#' er.eval <- evaliq(x.true[na.idx], x.est[na.idx], interactive = TRUE)
#'
#' # Turn off plot
#' # All of the three case will return the value of error
#' er.eval <- evaliq(x.true[na.idx], x.est[na.idx], plot = FALSE)
#' er.eval
#'
#'
#' \donttest{
#' # Real data example
#' set.seed(0)
#' data(toxicity, package = "misspi")
#' toxicity.miss <- missar(toxicity, 0.4, 0.2)
#' impute.res <- misspi(toxicity.miss)
#' x.imputed <- impute.res$x.imputed
#'
#' na.idx <- which(is.na(toxicity.miss))
#' evaliq(toxicity[na.idx], x.imputed[na.idx])
#' evaliq(toxicity[na.idx], x.imputed[na.idx], interactive = TRUE)}
#'
#' @author Zhongli Jiang \email{jiang548@purdue.edu}

#' @importFrom stats as.formula coef lm median na.omit predict var
#' @export


# Please be aware that the order matters !!!
evaliq <- function(x.true, x.impute, plot = TRUE, interactive = FALSE){

  if (!is.vector(x.true)) {
    x.true <- as.vector(x.true)
  }
  if (!is.vector(x.impute)) {
    x.impute <- as.vector(x.impute)
  }
  mse <- mean((x.true - x.impute)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(x.true - x.impute))
  nrmse <- sqrt(mse / var(x.true))

  # begin plot
  if(plot){
    data <- as.data.frame(cbind(x.impute, x.true))
    colnames(data) <- c("x.true", "x.impute")
    model <- lm(x.impute ~ x.true)
    eq <- substitute(Y == a + b %.% X*"," ~~ R^2~"="~r2,
                     list(a = round(unname(coef(model)[1]), 2),
                          b = round(unname(coef(model)[2]), 2),
                          r2 = round(summary(model)$r.squared, 2)
                     ))

    eq <- as.character(as.expression(eq))

    # Create a scatterplot of X vs. Y with a white background and centered title
    p <- ggplot(data, aes(x = x.true, y = x.impute)) +
      geom_point(size=1/log10(length(x.impute)), alpha=0.7) +
      geom_smooth(method = "lm", se = FALSE,
                  linewidth=1/log10(length(x.impute))) +
      labs(x = "True values", y = "Imputed values", title="True values vs Imputed values") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))

    if(interactive){
      print(ggplotly(p))
 #     plot.res <- p
    }else{
      p.static <- p +
        annotate(geom="text", x = min(data$x.true), y = max(data$x.impute), label=eq, parse=T, hjust=0, vjust=1)
      plot(p.static)
 #     plot.res <- p.static
    }
  }else{
 #     plot.res <- NULL
  }

#  res <- list(rmse=rmse, mae=mae, nrmse=nrmse, plot.res=plot.res)
   res <- list(rmse=rmse, mae=mae, nrmse=nrmse)
  return(res)
}
