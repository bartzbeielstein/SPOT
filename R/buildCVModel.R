#' buildCVModel
#' 
#' Build a set of models trained on different folds of cross-validated data.
#' Can be used to estimate the uncertainty of a given model type at any point.
#' 
#' @param x design matrix (sample locations)
#' @param y vector of observations at \code{x}
# @param modellingFunction the model that shall be fitted to each data fold
#' @param control (list), with the options for the model building procedure:\cr
#' \code{types} a character vector giving the data type of each variable. All but "factor" will be handled as numeric, "factor" (categorical) variables will be subject to the hamming distance.\cr
#' \code{target} target values of the prediction, a vector of strings. Each string specifies a value to be predicted, e.g., "y" for mean, "s" for standard deviation.
#' This can also be changed after the model has been built, by manipulating the respective \code{object$target} value.\cr
#' \code{uncertaintyEstimator} a character vector specifying which uncertaintyEstimator should be used.
#' "s" or the linearlyAdapted uncertrainty "sLinear". Default is "sLinear"
#'
#' @return set of models (class cvModel)
#' @export
buildCVModel <- function(x, y, control=list()){
    ## Load Control list
    con<-list(nFolds = 10,
              modellingFunction = buildKriging,
              target = c("y","s"),
              uncertaintyEstimator = "sLinear")
    con[names(control)] <- control
    control<-con
    
    control$nFolds <- min(control$nFolds, nrow(x))
    
    modellingFunction <- control$modellingFunction
    
    ## Empty List for the final model
    cvModel <- list()
    
    ## Save initial x and y in model
    cvModel$x <- x
    cvModel$y <- y
    cvModel[names(control)] <- control
    
    #Randomly shuffle the data
    shuffleIndexes <- sample(nrow(x))
    x <- x[shuffleIndexes,,drop = F]
    y <- y[shuffleIndexes, drop = F]
    
    #Create nFolds equally sized folds
    folds <- cut(seq(1,nrow(x)),breaks=control$nFolds,labels=FALSE)
    
    createSingleModel <- function(i){
        leaveOutIndex <- which(folds==i,arr.ind=TRUE)
        trainX <- x[-leaveOutIndex,, drop=F]
        trainY <- as.matrix(y[-leaveOutIndex])
        model <- modellingFunction(trainX, trainY, control = control)
        return(model)
    }
    
    cvModel$models <- lapply(1:control$nFolds, createSingleModel)
    class(cvModel)<- "cvModel"
    return(cvModel)
}

#' maxNearestNeighbourDistance
#'
#' Find the maximum distance between 2 nearest neighbours in a data set
#'
#' @param x matrix with candidate solutions
#'
#' @return maximum euclidean distance between two nearest neighbours
maxNearestNeighbourDistance <- function(x){
    minDists <- NULL
    for(i in 1:nrow(x)){
        currentDists <- abs(t(t(x[-i,, drop=F])-x[i,]))
        minDists <- c(minDists,sqrt(min(apply(currentDists,1,function(x){sum(x^2)}))))
    }
    return(sqrt(max(minDists)))
}

#' linearAdaptedSE
#' 
#' Linearly adapt the uncertainty estimation of a CV model regarding its distance to known neighbours
#'
#' @param sOld numeric vector, old uncertainty values
#' @param newdata matrix, new data points for which the uncertainty is estimated
#' @param x matrix, already evaluated data points
#'
#' @return numeric vector, adapted uncertainty values
linearAdaptedSE <- function(sOld, newdata, x){
    ifelse(is.null(nrow(newdata)),nr <- 1,nr <- nrow(newdata))
    newdata <- matrix(newdata, nrow = nr)
    if(nr <= 1){
        for(i in 1:nrow(newdata)){
            minDist <- min(abs(x-newdata[i]))
            sOld[i] <- sOld[i] * minDist/max(diff(sort(x)))
        }
    }else{
        for(i in 1:nrow(newdata)){
            dists <- abs(t(t(x)-newdata[i,]))
            minDist <- sqrt(min(apply(dists,1,function(x){sum(x^2)})))
            sOld[i] <- sOld[i] * minDist/maxNearestNeighbourDistance(x)
        }
    }
    sOld * 2 # ?
}

#' predict.cvModel
#'
#' Predict with the cross validated model
#'
#' @param object Kriging model (settings and parameters) of class \code{kriging}.
#' @param newdata design matrix to be predicted
#' @param ... Additional parameters passed to the model
#'
#' @return prediction results: list with predicted mean ('y'), estimated uncertainty ('y'), linearly adapted uncertainty ('sLinear')
#' @export
predict.cvModel <- function(object,newdata,...){
    predictSingle <- function(model){
        return(predict(model,as.matrix(newdata),...)$y)
    }
    
    if(is.null(object$uncertaintyEstimator)){
        object$uncertaintyEstimator <- "s"
    }
    
    results <- list()
    results$all <- sapply(object$models,predictSingle)
    
    ifelse(is.null(nrow(results$all)),nr <- 1,nr <- nrow(results$all))
    if(nr > 1){
        results$y <- apply(results$all,1,mean)
        funSE <- function(x){
            sd(x)/sqrt(length(x))
        }
        results$s <- apply(results$all,1 , funSE)
        
        if(tolower(object$uncertaintyEstimator) == "slinear"){
            results$s <- linearAdaptedSE(results$s, newdata, object$x)
        }else if(!(tolower(object$uncertaintyEstimator) %in% c("s", "slinear"))){
            stop("unrecognized option for modelControl$uncertaintyEstimator")
        }
    }else{
        results$y <- mean(results$all)
        results$s <- sd(results$all)/sqrt(length(results$all))
        
        if(tolower(object$uncertaintyEstimator) == "slinear"){
            results$s <- linearAdaptedSE(results$s, newdata, object$x)
        }else if(!(tolower(object$uncertaintyEstimator) %in% c("s", "slinear"))){
            stop("unrecognized option for modelControl$uncertaintyEstimator")
        }
    }
    
    return(results)
}