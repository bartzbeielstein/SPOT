###################################################################################
#' Build LOESS Model
#'
#' Build an interpolation model using the \code{loess} function. Essentially a SPOT-style
#' interface to that function.
#'
#' @param x design matrix (sample locations), rows for each sample, columns for each variable.
#' @param y vector of observations at \code{x}
#' @param control named list, with the options for the model building procedure loess. These
#' will be passed to loess as arguments. Please refrain from setting the formula or data arguments
#' as these will be supplied by the interface, based on \code{x} and \code{y}.
#'
#' @return returns an object of class \code{spotLOESS}.
#'
#' @seealso \code{\link{predict.spotLOESS}}
#'
#' @examples
#' ## Create a test function: branin
#' braninFunction <- function (x) {	
#' 	(x[2]  - 5.1/(4 * pi^2) * (x[1] ^2) + 5/pi * x[1]  - 6)^2 + 
#'	10 * (1 - 1/(8 * pi)) * cos(x[1] ) + 10
#' }
#' ## Create design points
#' set.seed(1)
#' x <- cbind(runif(40)*15-5,runif(40)*15)
#' ## Compute observations at design points
#' y <- as.matrix(apply(x,1,braninFunction))
#' ## Create model with default settings
#' fit <- buildLOESS(x,y)
#' fit
#' ## Predict new point
#' predict(fit,cbind(1,2))
#' ## True value at location
#' braninFunction(c(1,2))
#' ## Change model control
#' fit <- buildLOESS(x,y,control=list(parametric=c(TRUE,FALSE)))
#' fit
#' 
#' @export
###################################################################################
buildLOESS <- function(x, y, control=list()){ #nugget -1 means that the nugget will be optimized in lme
	con <- list(degree = 1, span = 0.2, normalize=FALSE)
	con[names(control)] <- control
	control<-con
			
	## number of variables 
	nParam <- ncol(x)
	
	## initialize fit list
	fit <- list(x=x,y=y)
	
	## to data frame
	x <- as.data.frame(x)
	y <- as.data.frame(y)
	colnames(y) <- "y"
	df <- cbind(y,x)
	
	## Extract parameter names (iputs and output)
  pNames <- colnames(x)
		
	## create a formula for variable coding (rescaling)
	fmla <- paste("y ~", pNames[1])
	if(nParam > 1){
		for (i in 2:nParam) 
			fmla <- paste(fmla, "*", pNames[i])
	}
	fmla <- as.formula(fmla)
	control$formula <- fmla
	control$data <- df
	
	fit <- do.call(loess,control)
	
  fit$loessfit=fit
	fit$fmla=fmla
	fit$target=control$target
	fit$min=min(fit$y)
	class(fit)<- "spotLOESS"
	fit
}

###################################################################################
#' Predict loess model
#' 
#' Predict with model produced by \code{\link{buildLOESS}}.
#'
#' @param object loess model (settings and parameters) of class \code{spotLOESS}.
#' @param newdata design matrix to be predicted
#' @param ... not used
#'
#' @return list with predicted value \code{y}, standard error 
#' \code{s} (if object$target contains \code{"s"}) and
#' \code{ei} (if object$target contains \code{"ei"}) and
#'
#' @seealso \code{\link{buildLOESS}}
#'
#' @export
#' @keywords internal
###################################################################################
predict.spotLOESS <- function(object,newdata,...){  
	x <- as.data.frame(newdata)
	colnames(x) <- colnames(object$x)	
	
	if (any(object$target %in% c("s","ei"))){
		res <- predict(object$loessfit,x,se=TRUE)	
		ret <- list(
			y=matrix(res$fit,nrow(x),1),			
			s=matrix(res$se.fit,nrow(x),1)	
		)
		if(any(object$target == "ei")){
      ret$ei <- expectedImprovement(ret$y,ret$s,object$min)
    }  
	}
	else{
		res <- predict(object$loessfit,x)	
		ret <- list(y=matrix(res,nrow(x),1))
	}
	ret	
}


###################################################################################################
#' Print method for loess model
#' 
#' Wrapper for \code{summary.loess}.
#'
#' @param object fit of the model, an object of class \code{"spotLOESS"}, produced by \code{\link{buildLOESS}}.
#' @param ... not used
#'
#' @seealso \code{\link{buildLOESS}}
#'
#' @export
#' @keywords internal
###################################################################################################
print.spotLOESS <- function(x,digits = max(3L, getOption("digits") - 3L), ...){
	x$loessfit$call <- NULL
  print(summary(x$loessfit))
}