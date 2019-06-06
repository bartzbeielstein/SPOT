
###################################################################################################
#' Sequential Parameter Optimization
#'
#' This is one of the main interfaces for using the SPOT package. Based on a user-given objective function
#' and configuration, \code{spot} finds the parameter setting that yields the lowest objective value (minimization).
#' To that end, it uses methods from the fields of design of experiment, statistical modeling / machine learning
#' and optimization.
#' 
#' @param x is an optional start point (or set of start points), specified as a matrix. One row for each point, and one column for each optimized parameter.
#' @param fun is the objective function. It should receive a matrix x and return a matrix y. In case the function uses external code and is noisy, an additional seed parameter may be used, see the \code{control$seedFun} argument below for details.
#' @param lower is a vector that defines the lower boundary of search space. This determines also the dimensionality of the problem.
#' @param upper is a vector that defines the upper boundary of search space.
#' @param control is a list with control settings for spot. See \code{\link{spotControl}}.
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
#' @examples
#' ## Most simple example: Kriging + LHS + predicted 
#' ## mean optimization (not expected improvement)
#' res <- spot(,funSphere,c(-2,-3),c(1,2))
#' res$xbest
#' ## With expected improvement
#' res <- spot(,funSphere,c(-2,-3),c(1,2),
#'    control=list(modelControl=list(target="ei"))) 
#' res$xbest
#' ### With additional start point:
#' #res <- spot(matrix(c(0.05,0.1),1,2),funSphere,c(-2,-3),c(1,2))
#' #res$xbest
## #larger budget:
#' #res <- spot(,funSphere,c(-2,-3),c(1,2),
#' #    control=list(funEvals=50)) 
#' #res$xbest
#' ### Use a local optimizer instead of LHS
#' #res <- spot(,funSphere,c(-2,-3),c(1,2),
#' #    control=list(optimizer=optimLBFGSB)) 
#' #res$xbest
#' ### Random Forest instead of Kriging
#' #res <- spot(,funSphere,c(-2,-3),c(1,2),
#' #    control=list(model=buildRandomForest)) 
#' #res$xbest    
#' ### LM instead of Kriging
#' #res <- spot(,funSphere,c(-2,-3),c(1,2),
#' #    control=list(model=buildLM)) #lm as surrogate
#' #res$xbest
#' ### LM and local optimizer (which for this simple example is perfect)
#' #res <- spot(,funSphere,c(-2,-3),c(1,2),
#' #    control=list(model=buildLM, optimizer=optimLBFGSB)) 
#' #res$xbest
#' ### Or a different Kriging model:
#' #res <- spot(,funSphere,c(-2,-3),c(1,2),
#' #    control=list(model=buildKrigingDACE, optimizer=optimLBFGSB)) 
#' #res$xbest
#' ## With noise: (this takes some time)
#' #res1 <- spot(,function(x)funSphere(x)+rnorm(nrow(x)),c(-2,-3),c(1,2),
#' #		control=list(funEvals=100,noise=TRUE)) #noisy objective
#' #res2 <- spot(,function(x)funSphere(x)+rnorm(nrow(x)),c(-2,-3),c(1,2),
#' #		control=list(funEvals=100,noise=TRUE,replicates=2,
#' #		designControl=list(replicates=2))) #noise with replicated evaluations
#' #res3 <- spot(,function(x)funSphere(x)+rnorm(nrow(x)),c(-2,-3),c(1,2),
#' #		control=list(funEvals=100,noise=TRUE,replicates=2,OCBA=T,OCBABudget=1,
#' #		designControl=list(replicates=2))) #and with OCBA
#' ### Check results with non-noisy function:
#' #funSphere(res1$xbest)
#' #funSphere(res2$xbest)
#' #funSphere(res3$xbest)
#' ## The following is for demonstration only, to be used for random number 
#' ## seed handling in case of external noisy target functions.
#' #res3 <- spot(,function(x,seed){set.seed(seed);funSphere(x)+rnorm(nrow(x))},
#' #    c(-2,-3),c(1,2),control=list(funEvals=100,noise=TRUE,seedFun=1))
#' ## 
#' ## Next Example: Handling factor variables
#' ## Note: factors should be coded as integer values, i.e., 1,2,...,n
#' ## create a test function:
#' braninFunctionFactor <- function (x) {  
#'   y <- (x[2]  - 5.1/(4 * pi^2) * (x[1] ^2) + 5/pi * x[1]  - 6)^2 + 
#'     10 * (1 - 1/(8 * pi)) * cos(x[1] ) + 10
#'   if(x[3]==1)
#'     y <- y +1
#'   else if(x[3]==2)
#'     y <- y -1
#'   y  
#' }
#' ## vectorize
#' objFun <- function(x){apply(x,1,braninFunctionFactor)}
#' set.seed(1)
#' res <- spot(fun=objFun,lower=c(-5,0,1),upper=c(10,15,3),     
#'     control=list(model=buildKriging, 
#'       types= c("numeric","numeric","factor"), 
#'       optimizer=optimLHD)) 
#' res$xbest
#' res$ybest
#' @export
###################################################################################################
spot <- function(x=NULL,fun, #mostly, fun must have format y=f(x,...). 
		## If a noisy function requires some specific seed handling (e.g., in some other non-R code) 
		## a seed can be passed to fun. For that purpose, the user must specify noise=TRUE and fun should 
		## be fun(x,seed,...)
		lower,upper,control=list(),...){
    
    #Initial Input Checking
    initialInputCheck(x,fun,lower,upper,control)
    
	## default settings
	dimension <- length(lower)
	con <- spotControl(dimension)
	con[names(control)] <- control
	control <- con
	rm(con)

	## All functions should deal with the same data types 
	control$designControl$types <- control$types
	control$modelControl$types <- control$types
	control$optimizerControl$types <- control$types
	
	## Initial design generation
	set.seed(control$seedSPOT)
  x <- control$design(x=x,lower=lower,upper=upper,control=control$designControl)

  ## Rounding values produced by the design function to integers, etc.
  x <- repairNonNumeric(x,control$types)
    
	## Evaluate initial design with objective function
  y <- objectiveFunctionEvaluation(x=NULL,xnew=x,fun=fun,seedFun=control$seedFun,noise=control$noise,...)
	
	result <- spotLoop(x=x,y=y,fun=fun,lower=lower,upper=upper,control=control, ...)
	result
} 

###################################################################################################
#' Default Control list for spot 
#'
#' This function returns the default controls for the functions \code{\link{spot}} and \code{\link{spotLoop}}.
#' Control is a list of the settings: 
#' \describe{
#'   \item{\code{funEvals}}{ This is the budget of function evaluations (spot uses no more than funEvals evaluations of fun), defaults to 20.}
#'   \item{\code{types}}{ Vector of data type of each variable as a string, defaults \code{"numeric"} for all variables.}
#'   \item{\code{design}}{A function that creates an initial design of experiment. Functions that accept the same parameters, 
#'				and return a matrix like \code{\link{designLHD}} or \code{\link{designUniformRandom}} can be used. Default is \code{\link{designLHD}}.}
#'   \item{\code{designControl}}{A list of controls passed to the \code{control} list of the \code{design} function. See help 
#'				of the respective function for details. Default is an empty \code{list}.}
#'   \item{\code{model}}{A function that builds a statistical model of the observed data. Functions that accept the same 
#'				parameters, and return a matrix like \code{\link{buildKriging}} or \code{\link{buildRandomForest}} 
#'				can be used. Default is \code{\link{buildKriging}}.}
#'   \item{\code{modelControl}}{A list of controls passed to the \code{control} list of the \code{model} function. 
#'				See help of the respective function for details.Default is an empty \code{list}.}
#'   \item{\code{optimizer}}{A function that is used to optimize based on \code{model}, finding the most promising 
#'				candidate solutions. Functions that accept the same parameters, and return a matrix like \code{\link{optimLHD}} 
#'				or \code{\link{optimLBFGSB}} can be used. Default is \code{\link{optimLHD}}.}
#'   \item{\code{optimizerControl}}{A list of controls passed to the \code{control} list of the \code{optimizer} function. 
#'				See help of the respective function for details. Default is an empty \code{list}.}
#'   \item{\code{noise}}{Boolean, whether the objective function has noise or not. Default is non-noisy, that is, \code{FALSE}.}
#'   \item{\code{OCBA}}{Boolean, indicating whether Optimal Computing Budget Allocation (OCBA) should be used in case of a noisy 
#'				objective function or not. OCBA controls the number of replications for each candidate solution. 
#' 				Note, that \code{replicates} should be larger than one in that case, and that the initial experimental design 
#'				(see \code{design}) should also have replicates larger one. Default is \code{FALSE}.}
#'   \item{\code{OCBAbudget}}{The number of objective function evaluations that OCBA can distribute in each iteration. Default is 3.}
#'   \item{\code{replicates}}{The number of times a candidate solution is initially evaluated, that is, in the initial design, 
#'				or when created by the optimizer. Default is \code{1}.}
#'   \item{\code{seedFun}}{An initial seed for the objective function in case of noise, by default \code{NA}. The default means that no seed is set.
#'				The user should be very careful with this setting. It is intended to generate reproducible experiments for each objective
#'				function evaluation, e.g., when tuning non-deterministic algorithms. If the objective function uses a constant number
#'				of random number generations, this may be undesirable. Note, that this seed is by default set prior to each evaluation. A replicated
#'				evaluation will receive an incremented value of the seed.
#'				Sometimes, the user may want to call external code using random numbers. To allow for that case, the user can specify an objective function (\code{fun}),
#' 				which has a second parameter \code{seed}, in addition to first parameter (matrix \code{x}). This seed can then be passed
#'				to the external code, for random number generator initialization. See end of examples section for a demonstration.}
#'   \item{\code{seedSPOT}}{This value is used to initialize the random number generator. It ensures that experiments are reproducible. Default is \code{1}.}
#'   \item{\code{duplicate}}{In case of a deterministic (non-noisy) objective function, this handles duplicated candidate solutions.
#'				By default (\code{duplicate = "EXPLORE"}), duplicates are replaced by new candidate solutions, generated by random 
#'				sampling with uniform distribution. If desired, the user can set this to "STOP", which means that the optimization
#'				stops and results are returned to the user (with a warning). This may be desirable, as duplicates can be a indicator
#'				for convergence, or for a problem with the configuration.
#'				In case of noise, duplicates are allowed.}
#'   \item{\code{plots}}{Whether progress should be tracked by a line plot, default is false}
#' }
#' @param dimension dimensionality of the problem, that is, the number of optimzed parameters.
#' @return a list
#' @export
#' @keywords internal
#' spotDefaultControls()
###################################################################################################
spotControl <- function(dimension){
list(
    funEvals=20,
		types= rep("numeric",dimension),
    design = designLHD,
    designControl = list(),
    model = buildKriging,
    modelControl = list(),
    optimizer = optimLHD,
    optimizerControl = list(),
		plots = FALSE,
    OCBA=FALSE,
    OCBABudget=1, #the budget available to OCBA, to be distributed to replications of "old" solutions
    replicates=1, #the number of replications for all "new" solutions (unless generated by the initial design, which handles them separately)
		noise= FALSE, #whether or not the target function is non-deterministic. 
		seedFun = NA, #start RNG seed for the target function (only important if non-deterministic, i.e., if noise==TRUE). NA means that seed is not set before running fun, this is important for cases where an initial seed for the target function is undesirable. (e.g., functions with constant and identical number of calls to random number generator)
		seedSPOT =1 #start RNG seed for the SPOT loop 
	  ## Note: in case of optimization algorithm tuning, where the optimization algorithm again solves 
		## a noisy objective function, the user should handle proper random number generator seeds for the 
		## noisy objective function.
    ## SPOT only sets/handles seeds for the main process as well as the tuned algorithm. 
    ## Check objectiveFunctionEvaluation to see how to properly disentangle different "random streams"
		## in R.
    ## Note: all evaluations are passed to the model, including repeated entries for replicated samples 
		## (e.g., based on OCBA). The model itself should handle repeated samples.
    ## The exception from this rule should be if noise==FALSE, in that case, duplicates are removed.
		## (?) TODO -> what if only one duplicate suggested -> will be suggested again in next iteration?
  )
}

###################################################################################################
#' Sequential Parameter Optimization Main Loop
#'
#' SPOT is usually started via the function \code{\link{spot}}. However, SPOT runs can be continued
#' (i.e., with a larger budget specified in \code{control$funEvals}) by using \code{spotLoop}.
#' This is the main loop of SPOT iterations. It requires the user to give the same inputs as
#' specified for \code{\link{spot}}.
#' 
#' @param x are the known candidate solutions that the SPOT loop is started with, specified as a matrix. One row for each point, and one column for each optimized parameter.
#' @param y are the corresponding observations for each solution in \code{x}, specified as a matrix. One row for each point.
#' @param fun is the objective function. It should receive a matrix x and return a matrix y. In case the function uses external code and is noisy, an additional seed parameter may be used, see the \code{control$seedFun} argument below for details.
#' @param lower is a vector that defines the lower boundary of search space. This determines also the dimensionality of the problem.
#' @param upper is a vector that defines the upper boundary of search space.
#' @param control is a list with control settings for spot. See \code{\link{spotControl}}.
#' @param ... additional parameters passed to \code{fun}.
#'
#' @return This function returns a list with:
#' \describe{
#'		\item{{xbest}}{Parameters of the best found solution (matrix).}
#'		\item{\code{ybest}}{Objective function value of the best found solution (matrix).}
#'		\item{\code{x}}{Archive of all evaluation parameters (matrix).}
#'		\item{\code{y}}{Archive of the respective objective function values (matrix).}
#'		\item{\code{count}}{Number of performed objective function evaluations.}
#'		\item{\code{msg}}{Message specifying the reason of termination.}
#'		\item{\code{modelFit}}{The fit of the last build model, i.e., an object returned by the last call to the function specified by \code{control$model}.}
#' }
#'
#' @examples
#' ## Most simple example: Kriging + LHS + predicted 
#' ## mean optimization (not expected improvement)
#' control <- list(funEvals=20)
#' res <- spot(,funSphere,c(-2,-3),c(1,2),control)
#' ## now continue with larger budget
#' control$funEvals <- 25
#' res2 <- spotLoop(res$x,res$y,funSphere,c(-2,-3),c(1,2),control)
#' res2$xbest
#' res2$ybest
#' @export
###################################################################################################
spotLoop <- function(x,y,fun,lower,upper,control,...){
    
    #Initial Input Checking
    initialInputCheck(x,fun,lower,upper,control, inSpotLoop = T)
    
	## default settings
	dimension <- length(lower)
	con <- spotControl(dimension)
	con[names(control)] <- control
	control <- con
	rm(con)
	
	## Initialize evaluation counter
  count <- nrow(y)
  
	## Main Loop
	modelFit <- NA
	while(count < control$funEvals){ 
		## Model building 
    modelFit <- control$model(x=x,y=y,control=control$modelControl) #todo return modelControl to allow memory?
		
		## Generate a surrogate target function from the model
		funSurrogate <- evaluateModel(modelFit)
		
		## Model optimization
    optimRes <- control$optimizer(,funSurrogate,lower,upper,control$optimizerControl) #todo return optimizerControl to allow memory?
    xnew <- optimRes$xbest    

		## Handling of duplicates 
		xnew <- duplicateAndReplicateHandling(xnew,x,lower,upper,control)

		## Rounding non-numeric values produced by the optimizer
		xnew <- repairNonNumeric(xnew,control$types)

		## If desired, use OCBA to handle replications of old solutions
    if(control$noise & control$OCBA){
      xnew <- rbind(xnew,repeatsOCBA(x,y,control$OCBABudget))
		}
    
    ## Prevent exceeding the budget:
    xnew <- xnew[1:min(max(control$funEvals-count,1),nrow(xnew)),,drop=FALSE]
		
		## Evaluation with objective function
		ynew <- objectiveFunctionEvaluation(x=x,xnew=xnew,fun=fun,seedFun=control$seedFun,noise=control$noise,...)
    
		## Plots, output, etc 
		if(control$plots){
			plot(y,type="l")
			abline(a=0,b=0)
		}
		
    ## 
    colnames(xnew) <- colnames(x)
    x <- rbind(x,xnew)
    y <- rbind(y,ynew)
    count <- count + nrow(ynew)
  }
    
  indexBest <- which.min(y)
  list(xbest=x[indexBest,,drop=FALSE],ybest=y[indexBest,,drop=FALSE],x=x, y=y, count=count,msg="budget exhausted",modelFit=modelFit)
}	