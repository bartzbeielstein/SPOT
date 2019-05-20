
###################################################################################################
#' Objective Function Evaluation
#'
#' This function handles the evaluation of the objective function in \code{\link{spot}}.
#' This includes handling of the random number generator stream as well as the actual evaluation.
#'
#' @param x matrix of already known solutions, to determine whether RNG seeds for new solutions need to be incremented.
#' @param xnew matrix of new solutions.
#' @param fun objective function to evaluate the solutions in \code{xnew}.
#' @param seedFun initial seed to be used for the random number generator seed. Set to NA to avoid using a fixed seed.
#' @param noise parameter specifying whether the target function is noisy.
#' @param ... parameters passed to \code{fun}.
#' 
#' @return the matrix ynew, which are the observations for fun(xnew)
#'
#' @seealso \code{\link{spot}} for more details on the parameters, e.g., \code{fun}
#'
#' @export
#' @keywords internal
###################################################################################################
objectiveFunctionEvaluation <- function(x=NULL,xnew,fun,seedFun=NA,noise=FALSE,...){ # TODO: vectorization
	if(!is.null(x))
		x <- data.matrix(x) 
	xnew <- data.matrix(xnew) #TODO: same as in ocba. else, problems with identical() due to names

  ## if xnew is empty, return.
	if(nrow(xnew)==0)
		return(numeric(0))
	
	## save seed status (to avoid disturbing the main loops RNG stream)
  ## note: theoretically only needed in case of noise==TRUE, but always done to avoid mistakes.
	if(exists(as.character(substitute(.Random.seed))))
    SAVESEED<-.Random.seed
  else
    SAVESEED=NULL
      
	if(noise & !is.na(seedFun)){	      
    ## calculate seeds for each evaluation
		seed <- numeric(nrow(xnew))
		for(i in 1:nrow(xnew)){
				xnewi <- xnew[i,]
				x <- rbind(x,xnewi)
				repetitions <- sum(apply(x,1,identical,xnewi)) -1
				seed[i] <- seedFun + repetitions
		}
		
		## either pass seed to objective function fun (if there is a parameter to receive it)
		nms <- names(formals(fun))		
		passSeed <- FALSE
		if(length(nms)>1){
			passSeed <- names(formals(fun)[2])=="seed"
		}			
		if(passSeed){
			ynew <- fun(xnew,seed,...) 
		}else{ ## or set seed directly here
			ynew <- NULL
			for(i in 1:nrow(xnew)){
				set.seed(seed[i])
				ynew <- rbind(ynew,fun(xnew[i,,drop=FALSE])) 
			}
		}
	}else{
		ynew <- fun(xnew,...)
	}
  
  ## load seed status (to avoid disturbing the main loops RNG stream), see save at start of function.
  if(!is.null(SAVESEED))
    assign(".Random.seed", SAVESEED, envir=globalenv())
	
	if(is.numeric(ynew))
		ynew <- matrix(ynew,,1) #convert to column matrix
		
	ynew
}
