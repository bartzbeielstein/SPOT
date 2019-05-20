
###################################################################################################
#' Linear Model Interface
#'
#' This is a simple wrapper for the lm function, which fits linear models.
#' The purpose of this function is to provide an interface as required by SPOT, 
#' to enable modeling and model-based optimization with linear models. 
#' The linear model is build with main effects. 
#' Optionally, the model is also
#' subject to the AIC-based stepwise algorithm, 
#' using the \code{step} function from the \code{stats} package.
#'
#' @param x matrix of input parameters. Rows for each point, columns for each parameter. 
#' @param y one column matrix of observations to be modeled.
#' @param control list of control parameters, currently only with 
#' parameters \code{useStep} and \code{formula}. 
#' The \code{useStep} boolean specifies whether the \code{step} function is used.
#' The \code{formula} is passed to the lm function.
#' Without a formula, a second order model will be built.
#' 
#' @return an object of class \code{"spotLinearModel"}, 
#' with a \code{predict} method and a \code{print} method.
#'
#' @export
#'
#' @examples
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
#' fit <- buildLM(x,y,control = list(algTheta=optimLHD))
#' ## Print model parameters
#' print(fit)
#' ## Predict at new location
#' predict(fit,cbind(1,2))
#' ## True value at location
#' braninFunction(c(1,2))
###################################################################################################
buildLM<-function(x,y,control=list()){
  
  ## Control settings
  con<-list(useStep=FALSE)
  con[names(control)] <- control
  control<-con
	
	control$x <- x
	control$y <- y
	
	## Convert inputs to combined data frame
	x <- as.data.frame(x)
	y <- as.data.frame(y)
	colnames(y) <- "y"
	df <- cbind(y,x)
	
	## Extract parameter names (iputs and output)
  pNames<-colnames(x)
  yName<-"y"
	 
  ## build formula
  if(is.null(control$formula)){
    Mainterms<-"." #Main effects 
    Fullformula <- as.formula(paste0(yName, "~",Mainterms) ) #combine to full formula
  }else{
    Fullformula <- control$formula
  }
    
	## Fit linear model	
	fit <- lm(Fullformula, data=df)
	
	## Improve stepwise
  if (control$useStep==TRUE){
    fit <- step(fit, trace=0, k=2)
  }
	
	## Prepare return
  control$fit <- fit  
  control$pNames <- pNames
  control$yName <- yName	
  class(control)<-"spotLinearModel"
  
	## Return
  return(control)
}


###################################################################################################
#' Prediction method for linear model
#'
#' Wrapper for \code{predict.lm}.
#'
#' @param object fit of the model, an object of class \code{"spotLinearModel"}, produced by \code{\link{buildLM}}.
#' @param newdata matrix of new data.
#' @param ... not used
#' 
#' @export
#' @keywords internal
###################################################################################################
predict.spotLinearModel<-function(object, newdata, ...){
	newdata <- as.data.frame(newdata)
  if(!all(colnames(newdata) %in% object$pNames))
    colnames(newdata) <- object$pNames
  res<-predict(object$fit,newdata)
  list(y=res)
}

###################################################################################################
#' Print method for linear model
#' 
#' Wrapper for \code{print.lm}.
#'
#' @param object fit of the model, an object of class \code{"spotLinearModel"}, produced by \code{\link{buildLM}}.
#' @param ... not used
#' 
#' @export
#' @keywords internal
###################################################################################################
print.spotLinearModel <- function(x,...){
  print(summary(x$fit))
}