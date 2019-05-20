
###################################################################################################
#' Minimization by L-BFGS-B
#' 
#' For minimization, this function uses the \code{"L-BFGS-B"} method from the \code{optim} function, 
#' which is part of the code{stats} package. It is basically a wrapper, to enable L-BFGS-B for usage
#' in SPOT.
#'
#' @param x optional matrix of points. Only first point (row) is used as startpoint.
#' @param fun objective function, which receives a matrix x and returns observations y
#' @param lower boundary of the search space
#' @param upper boundary of the search space
#' @param control list of control parameters
#' \describe{
#'   \item{\code{funEvals}}{Budget, number of function evaluations allowed. Default is 100.}
#' }
#' All other \code{control} parameters accepted by the \code{optim} function can be used, too, and are passed to \code{optim}.
#' @param ... passed to \code{fun}
#'
#' @return list, with elements
#' \describe{
#'   \item{\code{x}}{NA, not used}
#'   \item{\code{y}}{NA, not used}
#'   \item{\code{xbest}}{best solution}
#'   \item{\code{ybest}}{best observation}
#'   \item{\code{count}}{number of evaluations of \code{fun} 
#'        (estimated from the more complicated \code{"counts"} variable returned by \code{optim})}
#'   \item{\code{message}}{termination message returned by \code{optim}}
#' }
#'
#' @examples
#' res <- optimLBFGSB(,fun = funSphere,lower = c(-10,-20),upper=c(20,8))
#' res$ybest
#' @export
###################################################################################################
optimLBFGSB<-function(x=NULL,fun,lower,upper,control=list(),...){
  ## generate random start point
	if(is.null(x)) 
		x <- lower + runif(length(lower)) * (upper-lower)
	else
		x <- x[1,] #optim requires vector start point
	## handle default controls
	con<-list(funEvals=100)
	con[names(control)] <- control
	control<-con
	control$maxit <- control$funEvals
	control$types <- NULL
	
	## vectorization not properly handled by optim, so build a wrapper
	fn <- function(x,...)fun(matrix(x,1),...) 
	
	## remove other settings to avoid warnings
	control$funEvals <- NULL
	
	## start optim
  res <- optim(par=x,fn=fn,lower=lower,upper=upper,control=control,method="L-BFGS-B",...)
  list(x=NA,y=NA,xbest=matrix(res$par,1),ybest=res$value,count= res$counts[[1]] +res$counts[[2]] * 2 * length(res$par),msg=res$message)
}