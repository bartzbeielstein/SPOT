
###################################################################################################
#' Function Evaluation Wrapper
#'
#' This is a simple wrapper that turns a function of type \code{y=f(x)}, where x is a vector and y is a scalar,
#' into a function that accepts and returns matrices, as required by \code{\link{spot}}.
#' Note that the wrapper essentially makes use of the apply function. This is effective, but not necessarily
#' efficient. The wrapper is intended to make the use of spot easier, but it could be faster
#' if the user spends some time on a more efficient vectorization of the target function.
#'
#' @param fun the function \code{y=f(x)} to be wrapped, with x a vector and y a numeric
#' @return a function in the style of \code{y=f(x)}, accepting and returning a matrix
#'
#' @examples
#' ## example function
#' branin <- function (x) {  
#'   y <- (x[2]  - 5.1/(4 * pi^2) * (x[1] ^2) + 5/pi * x[1]  - 6)^2 + 
#'     10 * (1 - 1/(8 * pi)) * cos(x[1] ) + 10
#'   y  
#' }
#' ## vectorize / wrap
#' braninWrapped <-wrapFunction(branin)
#' ## test original
#' branin(c(1,2))
#' branin(c(2,2))
#' branin(c(2,1))
#' ## test wrapped
#' braninWrapped(matrix(c(1,2,2,2,2,1),3,2,byrow=TRUE))
#'
#' @export
###################################################################################################
wrapFunction <- function(fun){
  return(function(x){matrix(apply(x,1,fun),,1)})
}

###################################################################################################
#' Parallelized Function Evaluation Wrapper
#'
#' This is a simple wrapper that turns a function of type \code{y=f(x)}, where x is a vector and y is a scalar,
#' into a function that accepts and returns matrices, as required by \code{\link{spot}}.
#' While doing so, the wrapper will use the parallel package in order to parallelize the execution of each function
#' evaluation. This function will create a computation cluster if no cluster is specified and there is no
#' default cluster setup!
#'
#' @param fun the function that shall be evaluated in parallel
#' @param cl Optional, an existing computation cluster
#' @param nCores Optional, amount of cores to use for creating a new computation cluster. Default is all cores.
#'
#' @return numeric vector, result of the parallelized evaluation
#' @export
wrapFunctionParallel <- function(fun, cl = NULL, nCores = NULL){
    if(is.null(cl)){
        cl <- parallel::getDefaultCluster()
        if(is.null(cl)){
            if(is.null(nCores)){
                nCores <- parallel::detectCores()
            }
            parallel::setDefaultCluster(parallel::makeCluster(nCores))
        }
    }
    return(
        function(x){
            res <- parallel::parApply(cl,x,1,fun)
            matrix(res,,1)
        }
    )
}


###################################################################################################
#' wrapBatchTools
#' 
#' Wrap a given objective function to be evaluated via the batchtools package and make it accessible
#' for SPOT.
#'
#' @param fun function to wrap
#' @param reg batchtools registry, if none is provided, then one will be created automatically
#' @param clusterFunction batchtools clusterFunction, default: makeClusterFunctionsInteractive()
#' @param resources resource list that is passed to batchtools, default NULL
#'
#' @return callable function for SPOT
#' @export
###################################################################################################
wrapBatchTools <- function(fun, reg = NULL,
                           clusterFunction = batchtools::makeClusterFunctionsInteractive(), resources = NULL){
    
    if(is.null(reg)){
        reg <- tryCatch(batchtools::getDefaultRegistry(),error = {
            batchtools::makeRegistry(NA)
        })
    }
    
    reg$cluster.functions <- clusterFunction
    batchtools::setDefaultRegistry(reg)
    
    function(x){
        batchtools::clearRegistry(reg)
        
        if(!is.null(nrow(x))){
            x <- split(x, rep(1:nrow(x), each = ncol(x)))
            results <- as.vector(unlist(batchtools::btlapply(x, fun, resources = resources, reg = reg)))
        }else{
            results <- batchtools::btlapply(x, fun, resources = resources, reg = reg)
        }
        
        results
    }
}

#' wrapSystem_parseMatrixToString
#' 
#' Create a String that can be passed via the command line from an R-Matrix
#'
#' @param m a matrix
#'
#' @return parsed string
#' @keywords internal
wrapSystem_parseMatrixToString <- function(m){
    parseVecToString <- function(x){
        return(paste(x, sep =",", collapse = ","))
    }
    return(apply(m,1,parseVecToString))
}

#' wrapSystemCommand
#' 
#' Optimize parameters for a script that is accessible via Command Line
#'
#' @param systemCall String that calls the command line script. 
#'
#' @return callable function for SPOT
#' @export
#' @examples
#' \donttest{
#' exampleScriptLocation <- system.file("consoleCallTrialScript.R",package = "SPOT")
#' f <- wrapSystemCommand(paste("Rscript", exampleScriptLocation))
#' spot(,f,c(1,1),c(100,100))
#' }
wrapSystemCommand <- function(systemCall){
    return(
        function(x){
            paramStrings <- wrapSystem_parseMatrixToString(x)
            doSysCall <- function(p){
                system(paste(systemCall, p), intern=T)
            }
            res <- lapply(paramStrings, doSysCall)
            res <- matrix(as.numeric(unlist(res)), ncol=1)
            return(res)
        }
    )
}
