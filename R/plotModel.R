
###################################################################################################
#' Surface plot of a model
#'
#' A (filled) contour or perspective plot of a fitted model.
#'
#' @param object fit created by a modeling function, e.g., \code{\link{buildRandomForest}}.
#' @param which a vector with two elements, each an integer giving the two independent variables of the plot 
#' (the integers are indices of the respective data set). 
#' @param constant a numeric vector that states for each variable a constant value that it will take on
#' if it is not varied in the plot. This affects the parameters not selected by the \code{which} parameter.
#' By default, this will be fixed to the best known solution, i.e., the one with minimal y-value, according
#' to \code{which.min(object$y)}. The length of this numeric vector should be the same as the number of columns in \code{object$x}
#' @param xlab a vector of characters, giving the labels for each of the two independent variables.
#' @param ylab character, the value of the dependent variable predicted by the corresponding model.
#' @param type string describing the type of the plot:  \code{"filled.contour"} (default), \code{"contour"}, 
#' \code{"persp"} (perspective), or \code{"persp3d"} plot.
#' Note that "persp3d" is based on the plotly package and will work in RStudio, but not in the standard RGui.
#' @param ... additional parameters passed to the \code{contour} or \code{filled.contour} function.
#'
#' @examples
#' ## generate random test data
#' testfun <- function (x) sum(x^2)
#' set.seed(1)
#' k <- 30
#' x <- cbind(runif(k)*15-5,runif(k)*15,runif(k)*2-7,runif(k)*5+22)
#' y <- as.matrix(apply(x,1,testfun))
#' fit <- buildLM(x,y)
#' plotModel(fit)
#' plotModel(fit,type="contour")
#' plotModel(fit,type="persp")
#' plotModel(fit,which=c(1,4))
#' plotModel(fit,which=2:3)
#'
#' @seealso \code{\link{plotFunction}}, \code{\link{plotData}}
#' 
#' @export
###################################################################################################
plotModel <- function(object,which=if(ncol(object$x)>1 & tolower(type) != "singledim"){1:2}else{1},
						constant=object$x[which.min(object$y),], #best known solution. default.
						xlab= paste("x",which,sep=""),ylab="y",type="filled.contour",...){
    xlab <- xlab[order(which)]
    which <- sort(which)
    #number of variables
    nvar <- ncol(object$x)
    #bounds
    if(length(which) == 1){
        lower <- min(object$x)
        upper <- max(object$x)
    }else{
        lower <- apply(object$x[,which],2,min)
        upper <- apply(object$x[,which],2,max)
    }
    
    #varied variables
    vary <-  (1:nvar) %in% which
    force(object)
    force(nvar)
    force(vary)
    force(constant)
    
    ## Pre Checkup
    if(nvar < 2 & tolower(type) != "singledim"){
        stop("The specified plot type is only available for 2 or more dimensions")
    }
    
    ## Preparation offunction for 'plotFunction()'
    if(nvar == 1){
        plotfun <- evaluateModel(object, infillCriterion = infillGetFullPrediction)
    }else if(nvar == 2){
        if(tolower(type) == "singledim"){
            plotfun2 <- evaluateModel(object, infillCriterion = infillGetFullPrediction)
            plotfun2
            plotfun <- function(xx){ #fix constants
                z2 <- matrix(constant,length(xx),nvar,byrow=TRUE)
                z2[,which(vary)[1]] <- xx
                plotfun2(z2)
            }  
        }else{
            plotfun <- evaluateModel(object, infillCriterion = infillGetFullPrediction) 
        }
    }else if(nvar > 2){
        plotfun2 <- evaluateModel(object, infillCriterion = infillGetFullPrediction)
        plotfun2
        if(tolower(type) == "singledim"){
            plotfun <- function(xx){ #fix constants
                z2 <- matrix(constant,length(xx),nvar,byrow=TRUE)
                z2[,which(vary)[1]] <- xx
                plotfun2(z2)
            }  
        }else{
            plotfun <- function(xx){ #fix constants
                z2 <- matrix(constant,nrow(xx),nvar,byrow=TRUE)
                z2[,which(vary)[1]] <- xx[,1]
                z2[,which(vary)[2]] <- xx[,2]
                plotfun2(z2)
            }	
        }
    }else{
        stop("Dimensionality does not meet plot type")
    }
    
    ##Wrapper for plotFun$y
    plotfuny <- function(xx){
        res <- plotfun(xx)
        if(is.list(res)){
            return(res$y)
        }
        res
    }
    
    ## Call to plotFunction()
    if(type=="persp3d"){
        plotFunction(f=plotfuny,lower=lower,upper=upper,
                     type=type,
                     xlab=xlab[1],ylab=xlab[2],zlab=ylab,points1=cbind(object$x[,which],object$y),...)	
    }else if(tolower(type) == "singledim"){
        plotSingleDimFunction(plotfun, lower = lower, upper = upper, object$target)
    }else{
        plotFunction(f=plotfuny,lower=lower,upper=upper,
                     type=type,
                     xlab=xlab[1],ylab=xlab[2],zlab=ylab,points1=object$x[,which],...)	
    }
}