###################################################################################
#' Build Response Surface Model
#'
#' Using the \code{rsm} package, this function builds a linear response surface model.
#'
#' @param x design matrix (sample locations), rows for each sample, columns for each variable.
#' @param y vector of observations at \code{x}
#' @param control (list), with the options for the model building procedure:
#' \describe{
#'		\item{\code{mainEffectsOnly}}{Logical, defaults to FALSE. Set to TRUE if a model with main effects only is desired (no interactions, second order effects).}
#'		\item{\code{canonical}}{Logical, defaults to FALSE. If this is TRUE, use the canonical path to descent from saddle points. Else, simply use steepest descent}
#' }
#'
#' @return returns an object of class \code{spotRSM}.
#'
#' @seealso \code{\link{predict.spotRSM}}
#'
#' @examples
#' ## Create a test function: branin
#' braninFunction <- function (x) {	
#' 	(x[2]  - 5.1/(4 * pi^2) * (x[1] ^2) + 5/pi * x[1]  - 6)^2 + 
#'	10 * (1 - 1/(8 * pi)) * cos(x[1] ) + 10
#' }
#' ## Create design points
#' x <- cbind(runif(20)*15-5,runif(20)*15)
#' ## Compute observations at design points
#' y <- as.matrix(apply(x,1,braninFunction))
#' ## Create model with default settings
#' fit <- buildRSM(x,y)
#' ## Predict new point
#' predict(fit,cbind(1,2))
#' ## True value at location
#' braninFunction(c(1,2))
#' ## plots
#' plot(fit)
#' ## path of steepest descent
#' descentSpotRSM(fit)
#' 
#' @export
###################################################################################
buildRSM <- function(x, y, control=list()){ #nugget -1 means that the nugget will be optimized in lme
    con <- list(canonical = FALSE, mainEffectsOnly= FALSE)
    con[names(control)] <- control
    control<-con
    
    
    ## number of data samples
    nExp <- nrow(x)
    ## number of variables 
    nParam <- ncol(x)
    
    ## to data frame
    xx <- x
    yy <- y
    x <- as.data.frame(x)
    y <- as.data.frame(y)
    colnames(y) <- "y"
    df <- cbind(y,x)
    
    ## Extract parameter names (iputs and output)
    pNames <- colnames(x)
    
    ## get data bounds
    lower <- apply(x,2,min)
    upper <- apply(x,2,max)
    
    ## create a formula for variable coding (rescaling)
    fmla <- NULL				
    for (i in 1:nParam) {
        a <- lower[i]
        b <- upper[i]
        v1 <- mean(c(a,b))
        v2 <- (b-a)/2	
        fmla <- c(fmla,
                  as.formula(paste(paste("x",i,sep=""),
                                   "~ (", pNames[i], " - ", v1, ")/", v2
                  )
                  )
        )
    }
    
    ## code data
    codedData <- coded.data(df, formulas = fmla)
    
    #### TODO: the following seems to be redundant?
    ## check feasibility of data points in codedData:
    #df2x <- codedData[,1:nParam]	
    #codedData <- codedData[apply(codedData, 1, function(x) all(x >= -1 & x <= 1)), ]
    
    
    ## Determine number of required samples for a ...
    ## ... Linear model without interactions:
    nRequired1 <- 1 + nParam
    ## ... Linear model with interactions:
    nRequired2 <- 1 + nParam * (nParam + 1) / 2 
    ## ... Full quadratic model:
    nRequired3 <- 1 + nParam + nParam * (nParam + 1) / 2
    
    ## Create the most powerful model possible given the current number of parameters...
    makeNNames <- function(n) { Map(function(i) paste("x", i, sep = ""), 1:n) }
    makeNParameters <- function(n) { Reduce(function(n1, n2) paste(n1, n2, sep=","), makeNNames(n)) }
    paramString <- makeNParameters(nParam) # the string "x1,x2,...,xnParams"
    if ((nExp >= nRequired1 && nExp < nRequired2) | control$mainEffectsOnly) {
        rsmFormula <- as.formula(sprintf("y ~ FO(%s)", paramString))
    }	else if (nExp >= nRequired2 && nExp < nRequired3) {
        rsmFormula <- as.formula(sprintf("y ~ FO(%s) + TWI(%s)", paramString, paramString))
    }	else if (nExp >= nRequired3) {
        rsmFormula <- as.formula(sprintf("y ~ FO(%s) + TWI(%s) + PQ(%s)", paramString, paramString, paramString))
    }
    fit <- rsm(formula = rsmFormula, data = codedData)
    
    fit <- list(rsmfit=fit,fmla=fmla,nParam=nParam,codeddata=codedData,canonical=control$canonical) 
    fit$x <- xx
    fit$y <- yy
    class(fit)<- "spotRSM"
    fit$pNames <- pNames
    fit
}

###################################################################################
#' Predict RSM model
#' 
#' Predict with model produced by \code{\link{buildRSM}}.
#'
#' @param object RSM model (settings and parameters) of class \code{spotRSM}.
#' @param newdata design matrix to be predicted
#' @param ... not used
#'
#' @return list with predicted value \code{y}
#'
#' @seealso \code{\link{buildRSM}}
#'
#' @export
#' @keywords internal
###################################################################################
predict.spotRSM <- function(object,newdata,...){  
    if(!all(colnames(newdata) %in% object$pNames))
        colnames(newdata) <- object$pNames
    
    x <- as.data.frame(newdata)
    x <- coded.data(x, formulas = object$fmla)
    
    res <- predict(object$rsmfit,x)
    list(y=matrix(res,nrow(x),1))
}

###################################################################################
#' Plot RSM model
#' 
#' Plot model produced by \code{\link{buildRSM}}.
#'
#' @param x RSM model (settings and parameters) of class \code{spotRSM}.
#' @param ... parameters passed to plotting function (\code{contour})
#'
#' @export
#' @keywords internal
###################################################################################
plot.spotRSM <- function(x,...){
    nCol <- ceiling(sqrt(sum(1:(x$nParam-1))))	
    makeNNames <- function(n) { Map(function(i) paste("x", i, sep = ""), 1:n) }
    makeNParametersSum <- function(n) { Reduce(function(n1, n2) paste(n1, n2, sep="+"), makeNNames(n)) }
    mf <-  par("mfrow")
    par(mfrow=c(nCol,nCol) )
    contour(x$rsmfit, as.formula(paste("~",makeNParametersSum(x$nParam))), image=TRUE,...)
    par(mfrow=mf)
}

###################################################################################
#' Descent RSM model
#' 
#' Generate steps along the path of steepest descent for a RSM model.
#' This is only intended as a manual tool to use together with 
#' \code{\link{buildRSM}}.
#'
#' @param object RSM model (settings and parameters) of class \code{spotRSM}.
#' @return list with
#' \describe{
#'		\item{\code{x}}{list of points along the path of steepest descent}
#'		\item{\code{y}}{corresponding predicted values}
#' }
#'
#' @seealso \code{\link{buildRSM}}
#'
#' @export
###################################################################################
descentSpotRSM <- function(object){
    if(object$canonical){
        steepestDesc <- as.data.frame(canonical.path(object$rsmfit, descent=TRUE, dist = seq(-0.2,0.2, by = 0.1))[,2:eval(object$nParam+1)])
    }else{ # start at origin in one direction
        steepestDesc <- as.data.frame(steepest(object$rsmfit, descent=TRUE, dist = seq(0.1,1, by = 0.1))[,2:eval(object$nParam+1)])
    }
    ## ensure feasibility			
    #steepestDesc<- steepestDesc[apply(steepestDesc, 1, function(x) all(x > -1 & x < 1)), ]
    yHat <- predict(object$rsmfit,steepestDesc)
    yHat <- matrix(yHat,length(yHat),1)
    xNew <- code2val(steepestDesc, codings = codings(object$codeddata))	
    list(x=xNew, y=yHat)
}


###################################################################################################
#' Print method for RSM model
#' 
#' Wrapper for \code{summary.rsm}.
#'
#' @param object fit of the model, an object of class \code{"spotRSM"}, produced by \code{\link{buildRSM}}.
#' @param ... not used
#'
#' @seealso \code{\link{buildRSM}}
#'
#' @export
#' @keywords internal
###################################################################################################
print.spotRSM <- function(x,...){
    print(summary(x$rsmfit))
}