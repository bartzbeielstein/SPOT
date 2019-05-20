
###################################################################################################
#' Minimization by Differential Evolution
#' 
#' For minimization, this function uses the \code{"DEoptim"} method from the
#' code{DEoptim} package. It is basically a wrapper, to enable DEoptim for usage
#' in SPOT.
#'
#' @param x optional start point, not used in DEoptim
#' @param fun objective function, which receives a matrix x and returns observations y
#' @param lower boundary of the search space
#' @param upper boundary of the search space
#' @param control list of control parameters
#' \describe{
#'   \item{\code{funEvals}}{Budget, number of function evaluations allowed. Default is 200.}
#'   \item{\code{populationSize}}{Population size or number of particles in the population. Default is 10*dimension.}
#' }
#' @param ... passed to \code{fun}
#'
#' @return list, with elements
#' \describe{
#'   \item{\code{x}}{archive of the best member at each iteration}
#'   \item{\code{y}}{archive of the best value of fn at each iteration}
#'   \item{\code{xbest}}{best solution}
#'   \item{\code{ybest}}{best observation}
#'   \item{\code{count}}{number of evaluations of \code{fun}}
#' }
#'
#' @examples
#' res <- optimDE(,lower = c(-10,-20),upper=c(20,8),fun = funSphere)
#' res$ybest
#' @export
###################################################################################################
optimDE<-function(x=NULL,fun,lower,upper,control=list(),...){
  con<-list(funEvals=200,populationSize= (10 * length(lower)),trace=0)
  con[names(control)] <- control
  control<-con
  funEvals <- control$funEvals
  NP <- control$populationSize
	
  ## recalculate funEvals to DE iteration on basis of population size
  itermax <- floor((funEvals - NP) / NP)
  if(itermax < 1) itermax= 1
	
  control$NP <- control$populationSize
  control$itermax <- itermax
	
	## Delete unused settings
  control$populationSize <- NULL
  control$funEvals <- NULL
  control$types <- NULL
	
  ## wrapper for matrix inputs to fun
  fn <- function(x,...)fun(matrix(x,1),...) 

  ## start optim
  res <- DEoptim(fn=fn,lower=lower,upper=upper,control=control,...)
  list(x=res$member$bestmemit,y=res$member$bestvalit,xbest=matrix(res$optim$bestmem,1),ybest=matrix(res$optim$bestval,1),count= res$optim$nfeval*nrow(res$member$pop))
}