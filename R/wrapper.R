
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
#'
#' @examples
parallelWrapFunction <- function(fun, cl = NULL, nCores = NULL){
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
