
###################################################################################################
#' Minimization by Latin Hypercube Sampling
#' 
#' This uses Latin Hypercube Sampling (LHS) to optimize a specified target function.
#' A Latin Hypercube Design (LHD) is created with \code{\link{designLHD}}, then evaluated
#' by the objective function. All results are reported, including the best (minimal)
#' objective value, and corresponding design point.
#'
#' @param x optional matrix of points to be included in the evaluation
#' @param fun objective function, which receives a matrix x and returns observations y
#' @param lower boundary of the search space
#' @param upper boundary of the search space
#' @param control list of control parameters
#' \describe{
#'   \item{\code{funEvals}}{Budget, number of function evaluations allowed. Default: 100.}
#'   \item{\code{retries}}{Number of retries for design generation, used by \code{\link{designLHD}}. Default: 100.}
#' }
#' @param ... passed to \code{fun}
#'
#' @return list, with elements
#' \describe{
#'   \item{\code{x}}{archive of evaluated solutions}
#'   \item{\code{y}}{archive of observations}
#'   \item{\code{xbest}}{best solution}
#'   \item{\code{ybest}}{best observation}
#'   \item{\code{count}}{number of evaluations of \code{fun}}
#'   \item{\code{message}}{success message}
#' }
#'
#' @examples
#' res <- optimLHD(,fun = funSphere,lower = c(-10,-20),upper=c(20,8))
#' res$ybest
#' @export
###################################################################################################
optimLHD<-function(x=NULL,fun,lower,upper,control=list(),...){
	#if (length(par)==0) stop("dimension of par is null")
	con<-list(funEvals=100,retries=100,types= rep("numeric",length(lower)))
	con[names(control)] <- control
	control<-con
	  
  if(is.null(x)){
    k=0
  }else{
    k=nrow(x)
    if(k>=control$funEvals){
      stop("Design size in optimLHD is zero or negative, due to large number of rows of user-specified x matrix.")
    }
  }
  
  x <- designLHD(x, lower, upper, control=list(size=control$funEvals-k,retries=control$retries,types=control$types)) 
  y <- matrix(fun(x,...),,1)
  indexBest <- which.min(y)
  list(x=x,y=y,xbest=x[indexBest,,drop=FALSE],ybest=y[indexBest,,drop=FALSE],count=control$funEvals,msg="success")
}