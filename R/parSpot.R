
#' Parallel Sequential Parameter Optimization
#' 
#' Parallel adaptation of the standard SPOT method.
#' This adaptation enables the user to specify multiple control lists to be evaluated in parallel.
#' In each control list, a model and model optimizer etc. is specified.
#' They are evaluated on a parallel cluster. The results are joined in a synchronization phase.
#'
#' @param x is an optional start point (or set of start points), specified as a matrix. One row for each point, and one column for each optimized parameter.
#' @param fun is the objective function. It should receive a matrix x and return a matrix y. In case the function uses external code and is noisy, an additional seed parameter may be used, see the \code{control$seedFun} argument below for details.
#' @param lower is a vector that defines the lower boundary of search space. This determines also the dimensionality of the problem.
#' @param upper is a vector that defines the upper boundary of search space.
#' @param sequentialControlList is a list with control settings for spot. See \code{\link{spotControl}}.
#' @param parallelControlLists is a list of lists. Each list should contain a model + optimizer pair to be evaluated in parallel.
#' See \code{\link{spotControl}}.
#' @param ... additional parameters passed to \code{fun}.
#'
#' @return This function returns a list with:
#' \describe{
#'		\item{\code{xbest}}{Parameters of the best found solution (matrix).}
#'		\item{\code{ybest}}{Objective function value of the best found solution (matrix).}
#'		\item{\code{x}}{Archive of all evaluation parameters (matrix).}
#'		\item{\code{y}}{Archive of the respective objective function values (matrix).}
#'		\item{\code{count}}{Number of performed objective function evaluations.}
#'		\item{\code{msg}}{Message specifying the reason of termination.}
#'		\item{\code{modelFit}}{The fit of the last build model, i.e., an object returned by the last call to the function specified by \code{control$model}.}
#' }
#' 
#' @export
#'
#' @examples
#' res <- parSpot(x= NULL, funSphere, lower = c(-5,-5,-5), upper = c(5,5,5), 
#'     sequentialControlList = list(),
#'     parallelControlLists = list(list(model = buildKriging), 
#'                                 list(model = buildRandomForest),
#'                                 list(model = buildLM)))
parSpot <- function(x=NULL, fun, 
                    lower, upper, sequentialControlList = list(), 
                    parallelControlLists = list(list()), cl = NULL, nCores = NULL, ...){
    
    #Initial Input Checking
    initialInputCheck(x,lower,upper,sequentialControlList)
    
    ## Set up a parallel cluster if none is provided previously
    if(is.null(cl)){
        cl <- parallel::getDefaultCluster()
        if(is.null(cl)){
            if(is.null(nCores)){
                nCores <- parallel::detectCores()
            }
            parallel::setDefaultCluster(parallel::makeCluster(nCores))
        }
    }
    
    ## Filling sequentialControlList with standard spotControl
    ## 
    dimension <- length(lower)
    sequentialControlList <- spotFillControlList(sequentialControlList, dimension)
    
    ## Initial design generation
    set.seed(sequentialControlList$seedSPOT)
    x <- sequentialControlList$design(x=x,lower=lower,upper=upper,control=sequentialControlList$designControl)
    
    ## Rounding values produced by the design function to integers, etc.
    x <- repairNonNumeric(x,sequentialControlList$types)
    
    ## Evaluate initial design with objective function
    y <- objectiveFunctionEvaluation(x=NULL,xnew=x,fun=fun,seedFun=sequentialControlList$seedFun,noise=sequentialControlList$noise,...)
    
    #######
    ## 
    ## After the evaluation of the initial design each spot control should be run iteration by iteration in parallel
    ## 
    #######
    
    if(length(parallelControlLists) == 0){
        parallelControlLists <- list(spotFillControlList(list(), dimension))
    }else if(any(sapply(parallelControlLists, typeof) != "list")){
        error("parallelControlLists should be a list of lists.")
    }else{
        parallelControlLists <- lapply(parallelControlLists, function(cList){spotFillControlList(cList, dimension)})
    }
    
    totalBudget <- sequentialControlList$funEvals
    
    parallelSpotResults <- NULL
    parallel::clusterCall(cl = cl ,fun = function(){
        require(SPOT)
    })
    while((nrow(x) + length(parallelControlLists)) <= totalBudget){
        doParallelSpotIter <- function(cList){
            cList$funEvals <- nrow(x) + 1
            spotLoop(x=x,y=y,fun=function(x){x},lower=lower,upper=upper,control=cList,parallelCall = T, ...)
        }
        parallelSpotResults <- parallel::parLapply(cl = cl, parallelControlLists, doParallelSpotIter)
        newX <- do.call(rbind.data.frame, parallelSpotResults)
        newY <- objectiveFunctionEvaluation(x=x,xnew=newX,fun=fun,seedFun=sequentialControlList$seedFun,
                                            noise=sequentialControlList$noise,...)
        x <- rbind(x, newX)
        y <- rbind(y, newY)
    }
    
    indexBest <- which.min(y)
    parallelSpotResults$x <- x
    parallelSpotResults$y <- y
    parallelSpotResults$xbest <- x[indexBest,,drop=FALSE]
    parallelSpotResults$ybest <- y[indexBest,,drop=FALSE]
    parallelSpotResults$msg <- "budget exhausted"
    parallelSpotResults
}