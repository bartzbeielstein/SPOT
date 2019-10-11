
###################################################################################################
#' ranger Interface
#'
#' This is a simple wrapper for the \code{ranger} function from the \code{ranger} package.
#' The purpose of this function is to provide an interface as required by SPOT, to enable
#' modeling and model-based optimization with \code{ranger}.
#'
#' @param x matrix of input parameters. Rows for each point, columns for each parameter.
#' @param y one column matrix of observations to be modeled.
#' @param control list of control parameters. These are all configuration parameters 
#' of the \code{ranger} function, and will be passed on to it.
#' 
#' @return an object of class \code{"spotRanger"}, with a \code{predict} method and a \code{print} method.
#'
#' @export
#'
#' @examples
#'\dontrun{
#' ## Create a simple training data set
#' testfun <- function (x) x[1]^2
#' x <- cbind(sort(runif(30)*2-1))
#' y <- as.matrix(apply(x,1,testfun))
#' ## test data:
#' xt <- cbind(sort(runif(3000)*2-1))
#' ## Example with default model (standard randomforest)
#' fit <- buildRanger(x,y)
#' yt <- predict(fit,data.frame(x=xt))
#' plot(xt,yt$y,type="l")
#' points(x,y,col="red",pch=20)
#' ## Example with extratrees, an interpolating model
#' fit <- buildRanger(x,y,
#'                    control=list(rangerArguments = 
#'                                 list(replace = F,
#'                                    sample.fraction=1,
#'                                    min.node.size = 1,
#'                                    splitrule = "extratrees")))
#' yt <- predict(fit,data.frame(x=xt))
#' plot(xt,yt$y,type="l")
#' points(x,y,col="red",pch=20)
#'}
#' @importFrom ranger ranger
###################################################################################################
buildRanger <- function(x, y, control=list()){ 

	## to data frame
	x <- as.data.frame(x)
	y <- as.data.frame(y)
	colnames(y) <- "y"
	df <- cbind(y,x)
	
	## store parameter names
	fit <- list()	
	fit$pNames <- colnames(x)
	
	## if not give, initialize empty argument list for ranger
	if(is.null(control$rangerArguments))
		control$rangerArguments <- list()
		
	## formula for control		
	control$rangerArguments$formula <- "y ~ . "
	
	## data into control
	control$rangerArguments$data <- df
	
	## call ranger, with parameters taken from control
  fit$rangerFit <- do.call(ranger,control$rangerArguments)
	fit$x <- x
	fit$y <- y
  class(fit) <- "spotRanger"
  fit
}

###################################################################################################
#' Predictor for spotExtraTrees
#'
#' Wrapper for \code{predict.ranger}.
#'
#' @param object fit of the model, an object of class \code{"spotRandomForest"}, produced by \code{\link{buildRandomForest}}.
#' @param newdata matrix of new data.
#' @param ... not used
#' 
#' @export
#' @keywords internal
###################################################################################################
predict.spotRanger <- function(object,newdata,...){
  if(!all(colnames(newdata) %in% object$pNames))
    colnames(newdata) <- object$pNames
  res <- predict(object$rangerFit,newdata,...)$predictions
  list(y=res)
}

###################################################################################################
#' Print method for random forest
#' 
#' Wrapper for \code{print.ranger}.
#'
#' @param object fit of the model, an object of class \code{"spotRandomForest"}, produced by \code{\link{buildRandomForest}}.
#' @param ... not used
#' 
#' @export
#' @keywords internal
###################################################################################################
print.spotRanger <- function(x,...){
  print(x$rangerFit)
}