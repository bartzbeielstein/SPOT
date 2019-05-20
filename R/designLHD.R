####################################################################################
#' Normalized LHD Design
#'
#' Produces a normalized design and calculates the minimal distance 
#' if required.
#' A design  is a matrix with \code{dim} columns and \code{size} rows.
#' Distance can be calculated with respect to a fixed, nested design.
#' 
#' @param dim number, dimension of the problem (will be no. of columns of the result matrix)
#' @param size number of points with that dimension needed. (will be no. of rows of the result matrix).
#' @param calcMinDistance Boolean to indicate whether a minimal distance should be calculated.
#' @param nested nested design to be considered during distance calculation.
#' @param inequalityConstraint inequality constraint function, smaller zero for infeasible points. Used to replace infeasible points with random points. Has to evaluate points in interval [0;1].
#'
#' @seealso This function is used as a basis for \code{\link{designLHD}}.
#'
#' @return list \code{L}  \cr
#' - \code{L} consists of a matrix \code{L$design} and (if required) a minimal distance \code{L$minDistance}
#' @keywords internal
#' @author Original code by Christian Lasarczyk, adaptations by Martin Zaefferer
####################################################################################
designLHDNorm <- function(dim,size, calcMinDistance=FALSE, nested=NULL, inequalityConstraint=NULL){
	step <- 1/size;
	design <- replicate(dim, sample(0:(size-1),size) * step + runif(size) * step);

	if(!is.null(inequalityConstraint)){ #TODO: this may be inefficient if the feasible space is small.
		feasible <- apply(design,1,inequalityConstraint) <= 0
		if(any(feasible))
			design <- design[feasible,,drop=FALSE]
		else
			design <- matrix(NA,0,ncol(design))
		while(nrow(design)<size){
			newP <- runif(dim)
			if(inequalityConstraint(newP)<=0){
				design<-rbind(design,as.numeric(newP))
			}
		}
	}
	
	des <- rbind(design, nested) #adds nested design points, to be considered in distance calculation.
	
	if (calcMinDistance)
		minDistance <- min(dist(des))
	else
		minDistance <- NA;
	
	list(  design=design
			, minDistance=minDistance); 
}


####################################################################################
#' Latin Hypercube Design Generator
#'
#' Creates a latin Hypercube Design (LHD) with user-specified dimension and number of design points.
#' LHDs are created repeatedly created at random. For each each LHD, the minimal pair-wise distance between design points is computed.
#' The design with the maximum of that minimal value is chosen.
#'
#' @param x optional matrix x, rows for points, columns for dimensions. This can contain one or more points which are part of the design, but specified by the user. These points are added to the design, 
#' and are taken into account when calculating the pair-wise distances. They do not count for the design size. E.g., if \code{x} has two rows, \code{control$replicates} is one and \code{control$size} is ten, the returned design will have
#' 12 points (12 rows). The first two rows will be identical to \code{x}. Only the remaining ten rows are guaranteed to be a valid LHD. 
#' @param lower vector with lower boundary of the design variables (in case of categorical parameters, please map the respective factor to a set of contiguous integers, e.g., with lower = 1 and upper = number of levels)
#' @param upper vector with upper boundary of the design variables (in case of categorical parameters, please map the respective factor to a set of contiguous integers, e.g., with lower = 1 and upper = number of levels)
#' @param control list of controls:
#' \describe{
#'  \item{\code{size}}{number of design points}
#'  \item{\code{retries}}{number of retries during design creation}
#'  \item{\code{types}}{this specifies the data type for each design parameter, as a vector of either "numeric","integer","factor". (here, this only affects rounding)}
#'  \item{\code{inequalityConstraint}}{inequality constraint function, smaller zero for infeasible points. Used to replace infeasible points with random points.}
#'  \item{\code{replicates}}{integer for replications of each design point. E.g., if replications is two, every design point will occur twice in the resulting matrix.}
#' }
#'
#' @return matrix \code{design} \cr
#' - \code{design} has \code{length(lower)} columns and \code{(size + nrow(x))*control$replicates} rows.
#' All values should be within \code{lower <= design <= upper}
#' @export
#' @author Original code by Christian Lasarczyk, adaptations by Martin Zaefferer
#' @examples
#' set.seed(1) #set RNG seed to make examples reproducible 
#' design <- designLHD(,1,2) #simple, 1-D case
#' design
#' design <- designLHD(,1,2,control=list(replicates=3)) #with replications
#' design
#' design <- designLHD(,c(-1,-2,1,0),c(1,4,9,1),
#'		control=list(size=5, retries=100, types=c("numeric","integer","factor","factor")))
#' design
#' x <- designLHD(,c(1,-10),c(2,10),control=list(size=5,retries=100))
#' x2 <- designLHD(x,c(1,-10),c(2,10),control=list(size=5,retries=100))
#' plot(x2)
#' points(x, pch=19)
####################################################################################
designLHD <- function(x=NULL, lower, upper, control=list()) {
  ## TODO as matrix x?
  ## number of parameters
  n <- length(lower)
	
  ## defaults:
  con<-list(size=10, #number of design points
			retries=10, #number of randomly created designs, best selected based on max min distance
			replicates=1, #replications for each design point. leads to replicates*size rows.
			inequalityConstraint=NULL,
      types=rep("numeric",n)#data type of each column. possib
     )
	con[names(control)] <- control
	control<-con

  ## normalize x
  if(!is.null(x)){
    for (i in 1:length(lower)){
      lowerBound <-  lower[i]
      upperBound <-  upper[i]
      x[,i] <- (x[,i] - lowerBound)/ (upperBound-lowerBound)
    }
	}
	
	#constraint function: needs to be scaled in [0,1]
	if(!is.null(control$inequalityConstraint)){
		ineqConstraint <- control$inequalityConstraint
		force(ineqConstraint)
		force(lower)
		force(upper)
		ineqConstraint01 <- function(xx){
			xx <- lower + xx * (upper - lower)
			ineqConstraint(xx)
		}
	}else{
		ineqConstraint01 <- NULL
	}
  
	## Bei einer Wiederholung muss die Distanz nicht berechnet werden
	best <- designLHDNorm(length(lower),control$size,calcMinDistance=control$retries>1,nested=x,inequalityConstraint=ineqConstraint01)
	
	if (control$retries>1) {
		for (i in 1:(control$retries-1)) {
			tmpDes <- designLHDNorm(length(lower),control$size,calcMinDistance=TRUE,nested=x,inequalityConstraint=ineqConstraint01)
			## maximize minimal distance
			if (tmpDes$minDistance > best$minDistance)
				best <- tmpDes
		}
	}
	
	design <- rbind(x,best$design)
	for (i in 1:n){
		lowerBound <-  lower[i]
		upperBound <-  upper[i]    
    if(control$types[i] != "numeric"){
      lowerBound <- lowerBound - 0.5 #offset required so that lower bound has equal probability to be drawn randomly
			upperBound <- upperBound + 0.4999999999999  #offset required so that upper bound has equal probability to be drawn randomly
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
