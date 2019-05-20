####################################################################################
#' Uniform Design Generator
#'
#' Create a simple experimental design based on uniform random sampling.
#'
#' @param x optional data.frame x to be part of the design
#' @param lower vector with lower boundary of the design variables (in case of categorical parameters, please map the respective factor to a set of contiguous integers, e.g., with lower = 1 and upper = number of levels)
#' @param upper vector with upper boundary of the design variables (in case of categorical parameters, please map the respective factor to a set of contiguous integers, e.g., with lower = 1 and upper = number of levels)
#' @param control list of controls:\cr
#'  \code{size} number of design points\cr
#'  \code{types} this specifies the data type for each design parameter, as a vector of either "numeric","integer","factor". (here, this only affects rounding)\cr
#'  \code{replicates} integer for replications of each design point. E.g., if replications is two, every design point will occur twice in the resulting matrix.
#'
#' @return matrix \code{design} \cr
#' - \code{design} has \code{length(lower)} columns and \code{(size + nrow(x))*control$replicates} rows.
#' All values should be within \code{lower <= design <= upper}
#' @export
#' @examples
#' set.seed(1) #set RNG seed to make examples reproducible
#' design <- designUniformRandom(,1,2) #simple, 1-D case
#' design
#' design <- designUniformRandom(,1,2,control=list(replicates=3)) #with replications
#' design
#' design <- designUniformRandom(,c(-1,-2,1,0),c(1,4,9,1),
#' 		control=list(size=5, types=c("numeric","integer","factor","factor")))
#' design
#' x <- designUniformRandom(,c(1,-10),c(2,10),control=list(size=5))
#' x2 <- designUniformRandom(x,c(1,-10),c(2,10),control=list(size=5))
#' plot(x2)
#' points(x, pch=19)
####################################################################################
designUniformRandom <- function(x=NULL, lower, upper, control=list()) {
  ## number of parameters
  n <- length(lower)
	
  ## defaults:
  con<-list(size=10, #number of design points
			replicates=1, #replications for each design point. leads to replicates*size rows.
      types=rep("numeric",n))#data type of each column
	con[names(control)] <- control
	control<-con
    
  design <- matrix(runif(n*control$size),control$size,n) #t(replicate(control$size,lower+runif(n)*(upper-lower)))
	
	design <- rbind(x,design)
	for (i in 1:n){
		lowerBound <-  lower[i]
		upperBound <-  upper[i]    
    if(control$types[i] != "numeric"){
      lowerBound <- lowerBound - 0.5
			upperBound <- upperBound + 0.4999999999999
    }
		design[,i] <- lowerBound + design[,i] * (upperBound-lowerBound)
    if(control$types[i] != "numeric") #use rounding if not numeric. note that categorical parameters are mapped to integers.
      design[,i] <- floor(design[,i]+0.5)
	}
  
  if(control$replicates>1){
    m <- nrow(design)
    design <- design[rep(1:m,control$replicates),]
  }
  
	as.matrix(design)
}
#Tests: check nrow, ncol of design, check class of design, check class of design[,i] ...
#check lower <= design <= upper