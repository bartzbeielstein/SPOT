
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

