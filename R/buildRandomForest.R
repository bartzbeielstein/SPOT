
###################################################################################################
#' Random Forest Interface
#'
#' This is a simple wrapper for the randomForest function from the randomForest package.
#' The purpose of this function is to provide an interface as required by SPOT, to enable
#' modeling and model-based optimization with random forest.
#'
#' @param x matrix of input parameters. Rows for each point, columns for each parameter.
#' @param y one column matrix of observations to be modeled.
#' @param control list of control parameters, currently not used.
#' 
#' @return an object of class \code{"spotRandomForest"}, with a \code{predict} method and a \code{print} method.
#'
#' @export
#'
#' @examples
#'\dontrun{
#' ## Test-function:
#' braninFunction <- function (x) {	
#' 	(x[2]  - 5.1/(4 * pi^2) * (x[1] ^2) + 5/pi * x[1]  - 6)^2 + 
#'	10 * (1 - 1/(8 * pi)) * cos(x[1] ) + 10
#' }
#' ## Create design points
#' set.seed(1)
#' x <- cbind(runif(20)*15-5,runif(20)*15)
#' ## Compute observations at design points (for Branin function)
#' y <- as.matrix(apply(x,1,braninFunction))
#' ## Create model
#' fit <- buildRandomForest(x,y,control = list(algTheta=optimLHD))
#' ## Print model parameters
#' print(fit)
#' ## Predict at new location
#' predict(fit,cbind(1,2))
#' ## True value at location
#' braninFunction(c(1,2))
#'}
###################################################################################################
#todo: example
buildRandomForest <- function(x, y, control=list()){ 
  fit <- list()
  fit$rfFit <- randomForest(x, y)
	fit$pNames <- colnames(x)
	fit$x <- x
	fit$y <- y
  class(fit) <- "spotRandomForest"
  fit
}

#update.spotRandomForest <- function(x,y,fit){
#  
#}

###################################################################################################
#' Prediction method for random forest
#'
#' Wrapper for \code{predict.randomForest}.
#'
#' @param object fit of the model, an object of class \code{"spotRandomForest"}, produced by \code{\link{buildRandomForest}}.
#' @param newdata matrix of new data.
#' @param ... not used
#' 
#' @export
#' @keywords internal
###################################################################################################
predict.spotRandomForest <- function(object,newdata,...){
  if(!all(colnames(newdata) %in% object$pNames))
    colnames(newdata) <- object$pNames
  res <- predict(object$rfFit,newdata,...)  
  list(y=res)
}

###################################################################################################
#' Print method for random forest
#' 
#' Wrapper for \code{print.randomForest}.
#'
#' @param object fit of the model, an object of class \code{"spotRandomForest"}, produced by \code{\link{buildRandomForest}}.
#' @param ... not used
#' 
#' @export
#' @keywords internal
###################################################################################################
print.spotRandomForest <- function(x,...){
  print(x$rfFit)
}