
###################################################################################################
#' Minimization by GENetic Optimization Using Derivatives
#' 
#' For minimization, this function uses the \code{"genoud"} method from the
#' code{rgenoud} package. It is basically a wrapper, to enable genoud for usage
#' in SPOT.
#'
#' @param x optional start point, not used
#' @param fun objective function, which receives a matrix x and returns observations y
#' @param lower boundary of the search space
#' @param upper boundary of the search space
#' @param control list of control parameters
#' \describe{
#'   \item{\code{funEvals}}{Budget, number of function evaluations allowed. Default is 100.}
#'   \item{\code{populationSize}}{Population size, number of individuals in the population. Default is 10*dimension.}
#' }
#' @param ... passed to \code{fun}
#'
#' @return list, with elements
#' \describe{
#'   \item{\code{x}}{NULL, currently not used}
#'   \item{\code{y}}{NULL, currently not used}
#'   \item{\code{xbest}}{best solution}
#'   \item{\code{ybest}}{best observation}
#'   \item{\code{count}}{number of evaluations of \code{fun}}
#' }
#'
#' @examples
#' res <- optimGenoud(,fun = funSphere,lower = c(-10,-20),upper=c(20,8))
#' res$ybest
#' @export
###################################################################################################
optimGenoud<-function(x=NULL,fun,lower,upper,control=list(),...){
  con<-list(funEvals=1000,populationSize= (10 * length(lower))) #TODO: funEvals may not be the correct upper bound for function evaluations
  con[names(control)] <- control
  control <- con
  funEvals <- control$funEvals
  NP <- control$populationSize
  itermax <- floor(funEvals/NP) #TODO: most probably not correct limit.
  if(itermax < 1) itermax <- 1
  #umrechnung maxevals to DE iteration on basis op population size
  pop.size <- control$populationSize
  int.seed <- runif(1)
  unif.seed <- runif(1)
  max.generations <- itermax
  nvars <- length(lower)
  Domains <- cbind(as.numeric(lower),as.numeric(upper))
  control$populationSize <- NULL
  control$funEvals <- NULL
  control$types <- NULL
  ## vectorization wrapper
  fn <- function(x,...)fun(matrix(x,1),...) 
  ## start genoud
  res <- genoud(fn=fn,nvars=nvars,pop.size=pop.size,max.generations=max.generations,
								Domains=Domains,boundary.enforcement=2,
								int.seed=int.seed,unif.seed=unif.seed,control=control,...)
  list(x=NULL,y=NULL,xbest=matrix(res$par,1),ybest=matrix(res$value,1),
				count= res$generations * pop.size) 
				#TODO: this count value is not correct.
}