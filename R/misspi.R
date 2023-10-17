#' Missing Value Imputation in Parallel
#'
#' Enables embarrassingly parallel computing for imputation. Some of the advantages include
#' \itemize{
#' \item Provides fast implementation especially for high dimensional datasets.
#' \item Accepts a variety of machine learning models as methods with friendly user portal.
#' \item Supports multiple initializations.
#' \item Supports early stopping that prohibits unnecessary iterations.
#' }
#'
#' @param x a matrix of numerical values for imputation, missing value should all be "NA".
#' @param ncore number of cores to use, will be set to the cores detected as default.
#' @param init.method initializing method to fill in the missing value before imputation. Support "rf" for random forest imputation as default, "mean" for mean imputation, "median" for median imputation.
#' @param method method name for the imputation, support "rf" for random forest, "lgb" for lightgbm, "lasso" for LASSO, or "customize" if you want to use your own method.
#' @param earlystopping a Boolean which indicates whether to stop the algorithm if the relative difference stop decreasing, with TRUE as default.
#' @param ntree number of trees to use for imputation when method is "rf" or "gbm".
#' @param init.ntree number of trees to use for initialization when method is "rf"
#' @param viselect the number of variables with highest variable importance calculated from random forest initialization to work on if the value is not NULL. This would only work when init.method is "rf", and method is "rf" or "gbm".
#' @param model.train machine learning model to be invoked for customizing the imputation. Only invoked when parameter method = "customize". The input model should be able to take y~x for fitting process where y, and x are matrices, also make sure that it could be called using method "predict" for model prediction. You could pass the parameters for the model through the additional arguments ...
#' @param lgb.params parameters to customize for lightgbm models, could be invoked when method is "rf" or "gbm".
#' @param lgb.params0 parameters to customize for initialization using random forest, could be invoked when init.method is "rf".
#' @param pmm a Boolean which indicated whether to use predictive mean matching.
#' @param nn number of neighbors to use for prediction if predictive mean matching is invoked (pmm is "TRUE").
#' @param intcol a vector of indices of columns that are know to be integer, and will be round to integer in every iteration.
#' @param maxiter maximum number of iterations for imputation.
#' @param rdiff.thre relative difference threshold for determining the imputation convergence.
#' @param verbose a Boolean that indicates whether to print out the intermediate steps verbally.
#' @param progress a Boolean that indicates whether to show the progress bar.
#' @param nlassofold number of folds for cross validation when the method is "lasso".
#' @param isis a Boolean that indicates whether to use isis if the method is "lasso", recommended to use for ultra high dimension.
#' @param char a character to use which also accept unicode for progress bar. For example, u03c, u213c for pi, u2694 for swords, u2605 for star, u2654 for king, u26a1 for thunder, u2708 for plane.
#' @param iteration a Boolean that indicates whether use iterative algorithm.
#' @param ndecimal number of decimals to round for the result, with NULL meaning no intervention.
#' @param ... other arguments to be passed to the method.
#'
#' @returns a list that contains the imputed values, time consumed and number of iterations.
#' @return x.imputed the imputed matrix.
#' @return time.elapsed time consumed for the algorithm.
#' @return niter number of iterations used in the algorithm.
#'
#'
#' @references Azur, M. J., Stuart, E. A., Frangakis, C., & Leaf, P. J. (2011). Multiple imputation by chained equations: what is it and how does it work?. International journal of methods in psychiatric research, 20(1), 40-49.
#' @references Ke, G., Meng, Q., Finley, T., Wang, T., Chen, W., Ma, W., ... & Liu, T. Y. (2017). Lightgbm: A highly efficient gradient boosting decision tree. Advances in neural information processing systems, 30.
#' @references Stekhoven, D. J., & Bühlmann, P. (2012). MissForest—non-parametric missing value imputation for mixed-type data. Bioinformatics, 28(1), 112-118.
#' @references Fan, J., & Lv, J. (2008). Sure independence screening for ultrahigh dimensional feature space. Journal of the Royal Statistical Society Series B: Statistical Methodology, 70(5), 849-911.
#' @references Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society Series B: Statistical Methodology, 58(1), 267-288.
#' @references Breiman, L. (2001). Random forests. Machine learning, 45, 5-32.
#'
#'
#'
#' @seealso \code{\link{missar}}
#' @examples
#'
#' \donttest{
#' # Quick example 1
#' # Load a small data
#' data(iris)
#' # Keep numerical columns
#' num.col <- which(sapply(iris, is.numeric))
#' iris.numeric <- as.matrix(iris[, num.col])
#' set.seed(0)
#' iris.miss <- missar(iris.numeric, 0.3, 1)
#' iris.impute <- misspi(iris.miss)
#' iris.impute
#'
#' # Quick example 2
#' # Load a high dimensional data
#' data(toxicity, package = "misspi")
#' set.seed(0)
#' toxicity.miss <- missar(toxicity, 0.4, 0.2)
#' toxicity.impute <- misspi(toxicity.miss)
#' toxicity.impute
#'
#' # Change cores
#' iris.impute.5core <- misspi(iris.miss, ncore = 5)
#'
#' # Change initialization and maximum iterations (no iteration in the example)
#' iris.impute.mean.5iter <- misspi(iris.miss, init.method = "mean", maxiter = 0)
#'
#' # Change fun shapes for progress bar
#' iris.impute.king <- misspi(iris.miss, char = " \u2654")
#'
#'
#' # Use variable selection
#' toxicity.impute.vi <- misspi(toxicity.miss, viselect = 128)
#'
#'
#' # Use different machine learning algorithms as method
#' # linear model
#' iris.impute.lm <- misspi(iris.miss, model.train = lm)
#'
#' # From external packages
#' # Support Vector Machine (SVM)
#'
#' library(e1071)
#' iris.impute.svm.radial <- misspi(iris.miss, model.train = svm)
#'
#'
#' # Neural Networks
#'
#' library(neuralnet)
#' iris.impute.nn <- misspi(iris.miss, model.train = neuralnet)}
#'
#'
#' @docType package
#'
#' @author Zhongli Jiang \email{jiang548@purdue.edu}
#'
#' @import doParallel doSNOW glmnet lightgbm SIS foreach
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom parallel detectCores makeCluster stopCluster
#'
#'
#' @export

misspi <- function(x, ncore = NULL, init.method = "rf", method = "rf", earlystopping = TRUE,
                    ntree = 100, init.ntree = 100, viselect = NULL, lgb.params = NULL, lgb.params0 = NULL,
                    model.train = NULL, pmm = TRUE, nn=3, intcol = NULL, maxiter = 10, rdiff.thre = 0.01,
                    verbose = TRUE, progress = TRUE, nlassofold = 5, isis = FALSE, char = " * ",
                    iteration = TRUE, ndecimal = NULL, ...){

  s.time <- Sys.time()

  if(!is.numeric(x)){
    stop("We currently only support numeric matrices ...\n")
  }

  if(is.null(colnames(x))){
    colnames(x) <- paste0("x", 1:ncol(x))
  }

  if(!is.null(model.train)){
    method <- "customize"
    message("Since your model.train is not null, we have set method to customize ... \n")
  }

  if(method!="rf" & method!="gbm" & method!="lasso" & method!="customize"){
    stop("method not supported ... ")
  }

  if(method == "customize"){
    if(is.null(model.train)){
      stop("model.train must be provided in the customize mode")
    }
  }

  if(ncol(x) > 1000 & is.null(viselect) & method!="customize" & init.method=="rf" & method!="lasso"){
    message("We highly recommend activate viselect since your data is in high dimension. This may highly improve the speed and imputation accuracy ...\n")
  }

#  if(parallel){
    if(is.null(ncore)){
      ncore <- detectCores()
    }
    if(ncore > ncol(x)){
      ncore <- ncol(x)
      if(verbose){
      cat(paste0("Set number of cores used to ", ncore, " necessary ... \n"))
      }
    }
 # }else{
 #    if(is.null(ncore)){
 #   ncore <- 1
 #    }
 # }

  if(ncore > detectCores()){
    message(paste0("number of cores overloaded, reset the number of cores to the number of cores detected: ", detectCores(), "... \n"))
    ncore <- detectCores()
  }

  ncore.limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(ncore.limit) && ncore.limit == "TRUE") {
    ncore <- 2
    message("Set cores to 1 due to the limit ... \n")
  }


  #if(parallel){

  cl <- makeCluster(ncore)
  #  registerDoParallel(cl)
  registerDoSNOW(cl)
  if(verbose){
    cat(paste0("Parallel computing using ", ncore, " cores ...\n"))
  }
  #  }

  # Seperate variables types

  if("tbl" %in% x) {
    x <- as.data.frame(x)
  }

 # numerical.col <- which(sapply(x, is.numeric))
#  categorical.col <- c(which(sapply(x, is.factor)), which(sapply(x, is.character)), which(sapply(x, is.logical)))

  if (!is.matrix(x)) {
  x <- as.matrix(x)
  }

  n <- nrow(x)
  p <- ncol(x)

  if(verbose){
    cat("Imputing a matrix with ", n, " rows and ", p, " columns ... \n")
  }

  # Neighbor mean matching.  For a number x, find closest k value in the vector y, and return the corresponding mean in y2
  nmm <- function(x, y, y2, k){
    diff <- abs(y-x)
    if(length(diff)<k){ k <- length(diff) }

    res <- mean(y2[order(diff)[1:k]])
    return(res)
  }

  # Initialize
  # corresponds to missing value low cation
  miss.cord <- which(is.na(x), arr.ind=T)
  miss.table <- table(miss.cord[, 2])
  # missing value index
  miss.ind <- as.vector(as.numeric(names(miss.table)))
  nmiss <- length(miss.ind)
  miss.rate <- as.vector(miss.table/n)


  if(verbose){
    cat(paste0("Highest missing rate is ", max(miss.rate), " ... \n"))
  }

  x.imputed <- x

  # initialize
  if(verbose){
    cat("Initializing ... \n")
  }

  if (progress) {
  pb <- txtProgressBar(max=nmiss, style = 3, char=char)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  }else{
    opts <- NULL
  }

  if(init.method=="mean"){
    x.imputed.tmp <- foreach(i = miss.ind, .combine = "cbind", .options.snow=opts) %dopar% {
      miss.row.id <- miss.cord[which(miss.cord[, 2]==i), 1]
      predict.miss <- mean(x[-miss.row.id, i])
      x.imputed[miss.row.id, i] <- predict.miss
      return(x.imputed[, i])
    }
    x.imputed[, miss.ind] <- x.imputed.tmp

  }else if(init.method=="median"){
    x.imputed.tmp <- foreach(i = miss.ind, .combine = "cbind", .options.snow=opts) %dopar% {
      miss.row.id <- miss.cord[which(miss.cord[, 2]==i), 1]
      predict.miss <- median(x[-miss.row.id, i])
      x.imputed[miss.row.id, i] <- predict.miss
      return(x.imputed[, i])
    }
    x.imputed[, miss.ind] <- x.imputed.tmp

  }else if(init.method=="rf"){

    if(!is.null(viselect)){
      if(viselect %% 1 !=0 | viselect <= 0){
        stop("viselelct must be a positive integer !!!")
      }else{
        viselect.list <- lapply(1:length(miss.ind), function(x) NULL)
      }
    }

    if(is.null(lgb.params0)){
      params <- list(
        objective = "regression",
        boosting_type = "rf",        # Use Random Forest boosting type
        bagging_freq = 1,            # Frequency for bagging (1 means full bagging)
        bagging_fraction = 0.8,      # Fraction of data to be used for bagging
        feature_fraction = 0.33       # Fraction of features to be used for building trees
      )
      params$n_estimators <- init.ntree
      #  if(n<100){
      #  params$bagging_fraction = 1
      #   }
    }else{
      params <- lgb.params0
      params$n_estimators <- init.ntree
    }

    # initialize global variable
    i <- 1

    if(is.null(viselect)){

      x.imputed.tmp <- foreach(i = miss.ind, .combine = "cbind", .options.snow=opts) %dopar% {
        miss.row.id <- miss.cord[which(miss.cord[, 2]==i), 1]
        data.lgb <- lightgbm::lgb.Dataset(data = x[-miss.row.id, -i], label = x[-miss.row.id, i])
        model <- lightgbm::lgb.train(params=params, data=data.lgb, verbose=-1, ...)

        predict.miss <- predict(model, x[miss.row.id, -i, drop=F])

        if(pmm){
          predict.obs <- predict(model, x.imputed[-miss.row.id, -i, drop=F])
          predict.mm <- sapply(predict.miss, nmm, y=predict.obs, y2=x.imputed[-miss.row.id, i], k=nn)
          x.imputed[miss.row.id, i] <- predict.mm
          res <- x.imputed[, i]
        }else{
          x.imputed[miss.row.id, i] <- predict.miss
          res <- x.imputed[, i]
        }
        return(res)
      }
      x.imputed[, miss.ind] <- x.imputed.tmp

    }else{
      # integrating imporatance
      x.imputed.tmp <- foreach(i = miss.ind, .combine = "rbind", .options.snow=opts) %dopar% {

        res <- list()
        miss.row.id <- miss.cord[which(miss.cord[, 2]==i), 1]
        data.lgb <- lightgbm::lgb.Dataset(data = x[-miss.row.id, -i], label = x[-miss.row.id, i])
        model <- lightgbm::lgb.train(params=params, data=data.lgb, verbose=-1, ...)

        # Get feature importance scores from the full-feature model
        importance.scores <- lightgbm::lgb.importance(model)
        top.features <- na.omit(unlist(importance.scores[order(-importance.scores$Gain), "Feature"][1:viselect]))

        predict.miss <- predict(model, x[miss.row.id, -i, drop=F])

        if(pmm){
          predict.obs <- predict(model, x.imputed[-miss.row.id, -i, drop=F])
          predict.mm <- sapply(predict.miss, nmm, y=predict.obs, y2=x.imputed[-miss.row.id, i], k=nn)
          x.imputed[miss.row.id, i] <- predict.mm
          res$imputed <- x.imputed[, i]
        }else{
          x.imputed[miss.row.id, i] <- predict.miss
          res$imputed <- x.imputed[, i]
        }
        res$top.features <- top.features
        return(res)
      }
      x.imputed[, miss.ind] <- do.call(cbind, x.imputed.tmp[, "imputed"])
      viselect.list[miss.ind] <- x.imputed.tmp[, "top.features"]

    }

  }else{
    stop("Initialization method is not supported ...\n")
  }

  if(!is.null(intcol)){
    x.imputed[, intcol] <- apply(x.imputed[, intcol], 2, round)
  }

  # sort the variables with missing rate from lowest to highest
  miss.ind.sort <- miss.ind[order(miss.rate, decreasing=F)]

  niter <- 0
  rdiff <- 1
  w0 <- rep(1, p-1)

  if(as.numeric(maxiter)<1){
    iteration <- F
    if(verbose){
      cat("\n No iterations, imputation completed ... \n")
    }
  }

  while(1){
    if(iteration==F){
      break
    }
    rdiff.pre <- rdiff
    niter <- niter + 1
    if(verbose){
    cat(paste0("\nIteration ", niter, " ... \n"))
    }

    x.imputed.new <- x.imputed

    # count the loop inside miss.ind.sort
    count <- 0

    if(method=="rf"){

      x.imputed.new.tmp <- foreach(i = miss.ind.sort, .combine = "cbind", .options.snow=opts) %dopar% {
        # count <- count + 1
        miss.row.id <- miss.cord[which(miss.cord[, 2]==i), 1]

        if(is.null(lgb.params)){
          params <- list(
            objective = "regression",
            boosting_type = "rf",        # Use Random Forest boosting type
            bagging_freq = 1,            # Frequency for bagging (1 means full bagging)
            bagging_fraction = 0.8,      # Fraction of data to be used for bagging
            feature_fraction = 0.33       # Fraction of features to be used for building trees
          )
          params$n_estimators <- ntree
          #  if(n<100){
          #  params$bagging_fraction = 1
          #   }
        }else{
          params <- lgb.params
          params$n_estimators <- ntree
        }

        if(!is.null(viselect)){
          active.set <- viselect.list[[i]]
        }else{
          active.set <- -i
        }

        obs.data.lgb <- lightgbm::lgb.Dataset(data = x.imputed[-miss.row.id, active.set], label = x.imputed[-miss.row.id, i])
        model <- lightgbm::lgb.train(params=params, data=obs.data.lgb, verbose=-1, ...)

        predict.miss <- predict(model, x.imputed[miss.row.id, active.set, drop=F])

        if(pmm){
          predict.obs <- predict(model, x.imputed[-miss.row.id, active.set, drop=F])
          predict.mm <- sapply(predict.miss, nmm, y=predict.obs, y2=x.imputed[-miss.row.id, i], k=nn)
          x.imputed.new[miss.row.id, i] <- predict.mm
          res <- x.imputed.new[, i]
        }else{
          x.imputed.new[miss.row.id, i] <- predict.miss
          res <- x.imputed.new[, i]
        }
        return(res)

        #  if(progress){progress(count, nmiss)}
      }
      x.imputed.new[, miss.ind.sort] <-  x.imputed.new.tmp
    }

    if(method=="gbm"){

      x.imputed.new.tmp <- foreach(i = miss.ind.sort, .combine = "cbind", .options.snow=opts) %dopar% {
        # count <- count + 1
        miss.row.id <- miss.cord[which(miss.cord[, 2]==i), 1]

        if(is.null(lgb.params)){
          params <- list(
            objective = "regression"
            , boosting_type = "gbdt"
            , metric = "l2"
            , min_data = min(20, round(n*0.2))
            , learning_rate = 0.1
            , lambda_l2 = 0.1
            , feature_fraction = 1
            , bagging_fraction = 1
          )
          params$n_estimators <- ntree
        }else{
          params <- lgb.params
          params$n_estimators <- ntree
        }

        if(!is.null(viselect)){
          active.set <- viselect.list[[i]]
        }else{
          active.set <- -i
        }

        obs.data.lgb <- lightgbm::lgb.Dataset(data = x.imputed[-miss.row.id, active.set], label = x.imputed[-miss.row.id, i])

        model <- lightgbm::lgb.train(params=params, data=obs.data.lgb, verbose=-1, ...)

        predict.miss <- predict(model, x.imputed[miss.row.id, active.set, drop=F])

        if(pmm){
          predict.obs <- predict(model, x.imputed[-miss.row.id, active.set, drop=F])
          predict.mm <- sapply(predict.miss, nmm, y=predict.obs, y2=x.imputed[-miss.row.id, i], k=nn)
          x.imputed.new[miss.row.id, i] <- predict.mm
          res <- x.imputed.new[, i]
        }else{
          x.imputed.new[miss.row.id, i] <- predict.miss
          res <- x.imputed.new[, i]
        }

        return(res)

        #  if(progress){progress(count, nmiss)}
      }
      x.imputed.new[, miss.ind.sort] <-  x.imputed.new.tmp
    }

    if(method=="lasso"){
      x.imputed.new.tmp <- foreach(i = miss.ind.sort, .combine = "cbind", .options.snow=opts) %dopar% {
        # for(i in miss.ind.sort){
        # count <- count + 1
        miss.row.id <- miss.cord[which(miss.cord[, 2]==i), 1]
        if(isis){
          sis.res <- SIS::SIS(x = x.imputed[-miss.row.id, -i], y = x.imputed[-miss.row.id, i], iter = T, nfolds=nlassofold)
          w0[-sis.res$ix] <- Inf
        }
        cv.res <- glmnet::cv.glmnet(x = x.imputed[-miss.row.id, -i],  y = x.imputed[-miss.row.id, i], alpha = 1, nfolds=nlassofold, penalty.factor=w0, ...)
        model <- glmnet::glmnet(x = x.imputed[-miss.row.id, -i], y = x.imputed[-miss.row.id, i],
                                alpha = 1, lambda = cv.res$lambda.min, ...)
        predict.miss <- predict(model, x.imputed[miss.row.id, -i, drop=F])

        if(pmm){
          predict.obs <- predict(model, x.imputed[-miss.row.id, -i, drop=F])
          predict.mm <- sapply(predict.miss, nmm, y=predict.obs, y2=x.imputed[-miss.row.id, i], k=nn)
          x.imputed.new[miss.row.id, i] <- predict.mm
          res <- x.imputed.new[, i]
        }else{
          x.imputed.new[miss.row.id, i] <- predict.miss
          res <- x.imputed.new[miss.row.id, i]
        }
        return(res)
        #    if(progress){progress(count, nmiss)}
      }
      x.imputed.new[, miss.ind.sort] <- x.imputed.new.tmp
    }

    #browser()
    if(method=="customize"){
      x.imputed.new.tmp <- foreach(i = miss.ind.sort, .combine = "cbind", .options.snow=opts) %dopar% {
        # count <- count + 1
        miss.row.id <- miss.cord[which(miss.cord[, 2]==i), 1]

        x.imputed <- as.data.frame(x.imputed)

        train.formula <- as.formula(paste0(colnames(x.imputed)[i], " ~ ."))

        model <- model.train(train.formula, data=x.imputed, ...)

        predict.miss <- predict(model, as.data.frame(x.imputed[miss.row.id, -i, drop=F]))

        if(pmm){
          predict.obs <- predict(model, as.data.frame(x.imputed[-miss.row.id, -i, drop=F]))
          predict.mm <- sapply(predict.miss, nmm, y=predict.obs, y2=x.imputed[-miss.row.id, i], k=nn)
          x.imputed.new[miss.row.id, i] <- predict.mm
          res <- x.imputed.new[, i]
        }else{
          x.imputed.new[miss.row.id, i] <- predict.miss
          res <- x.imputed.new[, i]
        }
        return(res)

        #  if(progress){progress(count, nmiss)}
      }
      x.imputed.new[, miss.ind.sort] <-  x.imputed.new.tmp
    }

    if(!is.null(intcol)){
      x.imputed.new[, intcol] <- apply(x.imputed.new[, intcol], 2, round)
    }

    rdiff <- sum((x.imputed.new[miss.cord] - x.imputed[miss.cord])^2)/sum(x.imputed.new[miss.cord]^2)

    x.imputed <- x.imputed.new

    if(verbose){
      cat(paste0("\nRelative squared difference is ", rdiff, " ... \n"))
    }

    if( (rdiff < rdiff.thre) | (niter==maxiter)){
      if(niter==maxiter){
        break
      }
      if(verbose){
        cat(paste0("\nAlgorithm converged at iter ", niter, " ... \n"))
      }
      break
    }

    if( (niter>=2) * (rdiff>=rdiff.pre) * (earlystopping)){
      if(verbose){
      cat("Early stopping invoked ... \n")
      }
      break
    }

  }

  #  if(parallel){
  stopCluster(cl)
  #  }

  if(!is.null(ndecimal)){
  x.imputed <- round(x.imputed, ndecimal)
  }

  e.time <- Sys.time()
  time.elapsed <- e.time - s.time



  return(list(x.imputed=x.imputed, time.elapsed=time.elapsed, niter=niter))

}


# # Easter egg: Use Ramanujan-Sato series to approximate the name of the package
# options(digits = 15)
# ramanujan.pi <- function(x){
# res <- 1/(2 * sqrt(2) / 9801 * sum(factorial(4 * x) * (1103 + 26390 * x ) / (factorial(x)^4 * 396^(4 * x))))
# return(res)
# }
# ramanujan.pi(seq(0, 1))
# options(NULL)
# # The approximation is amazing even using only 2 terms !!!!!!!!!!!!!
