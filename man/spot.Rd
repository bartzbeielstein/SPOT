% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spot.R
\name{spot}
\alias{spot}
\title{Sequential Parameter Optimization}
\usage{
spot(x = NULL, fun, lower, upper, control = list(), ...)
}
\arguments{
\item{x}{is an optional start point (or set of start points), specified as a matrix. One row for each point, and one column for each optimized parameter.}

\item{fun}{is the objective function. It should receive a matrix x and return a matrix y. In case the function uses external code and is noisy, an additional seed parameter may be used, see the \code{control$seedFun} argument below for details.}

\item{lower}{is a vector that defines the lower boundary of search space. This determines also the dimensionality of the problem.}

\item{upper}{is a vector that defines the upper boundary of search space.}

\item{control}{is a list with control settings for spot. See \code{\link{spotControl}}.}

\item{...}{additional parameters passed to \code{fun}.}
}
\value{
This function returns a list with:
\describe{
	\item{\code{xbest}}{Parameters of the best found solution (matrix).}
	\item{\code{ybest}}{Objective function value of the best found solution (matrix).}
	\item{\code{x}}{Archive of all evaluation parameters (matrix).}
	\item{\code{y}}{Archive of the respective objective function values (matrix).}
	\item{\code{count}}{Number of performed objective function evaluations.}
	\item{\code{msg}}{Message specifying the reason of termination.}
	\item{\code{modelFit}}{The fit of the last build model, i.e., an object returned by the last call to the function specified by \code{control$model}.}
}
}
\description{
This is one of the main interfaces for using the SPOT package. Based on a user-given objective function
and configuration, \code{spot} finds the parameter setting that yields the lowest objective value (minimization).
To that end, it uses methods from the fields of design of experiment, statistical modeling / machine learning
and optimization.
}
\examples{
## Most simple example: Kriging + LHS + predicted 
## mean optimization (not expected improvement)
res <- spot(,funSphere,c(-2,-3),c(1,2),control=list(funEvals=15))
res$xbest
## With expected improvement
res <- spot(,funSphere,c(-2,-3),c(1,2),
   control=list(funEvals=15,modelControl=list(target="ei"))) 
res$xbest
### With additional start point:
#res <- spot(matrix(c(0.05,0.1),1,2),funSphere,c(-2,-3),c(1,2))
#res$xbest
#res <- spot(,funSphere,c(-2,-3),c(1,2),
#    control=list(funEvals=50)) 
#res$xbest
### Use local optimization instead of LHS
#res <- spot(,funSphere,c(-2,-3),c(1,2),
#    control=list(optimizer=optimLBFGSB)) 
#res$xbest
### Random Forest instead of Kriging
#res <- spot(,funSphere,c(-2,-3),c(1,2),
#    control=list(model=buildRandomForest)) 
#res$xbest    
### LM instead of Kriging
#res <- spot(,funSphere,c(-2,-3),c(1,2),
#    control=list(model=buildLM)) #lm as surrogate
#res$xbest
### LM and local optimizer (which for this simple example is perfect)
#res <- spot(,funSphere,c(-2,-3),c(1,2),
#    control=list(model=buildLM, optimizer=optimLBFGSB)) 
#res$xbest
### Or a different Kriging model:
#res <- spot(,funSphere,c(-2,-3),c(1,2),
#    control=list(model=buildKrigingDACE, optimizer=optimLBFGSB)) 
#res$xbest
## With noise: (this takes some time)
#res1 <- spot(,function(x)funSphere(x)+rnorm(nrow(x)),c(-2,-3),c(1,2),
#		control=list(funEvals=100,noise=TRUE)) #noisy objective
#res2 <- spot(,function(x)funSphere(x)+rnorm(nrow(x)),c(-2,-3),c(1,2),
#		control=list(funEvals=100,noise=TRUE,replicates=2,
#		designControl=list(replicates=2))) #noise with replicated evaluations
#res3 <- spot(,function(x)funSphere(x)+rnorm(nrow(x)),c(-2,-3),c(1,2),
#		control=list(funEvals=100,noise=TRUE,replicates=2,OCBA=T,OCBABudget=1,
#		designControl=list(replicates=2))) #and with OCBA
### Check results with non-noisy function:
#funSphere(res1$xbest)
#funSphere(res2$xbest)
#funSphere(res3$xbest)
## The following is for demonstration only, to be used for random number 
## seed handling in case of external noisy target functions.
#res3 <- spot(,function(x,seed){set.seed(seed);funSphere(x)+rnorm(nrow(x))},
#    c(-2,-3),c(1,2),control=list(funEvals=100,noise=TRUE,seedFun=1))
## 
## Next Example: Handling factor variables
## Note: factors should be coded as integer values, i.e., 1,2,...,n
## create a test function:
braninFunctionFactor <- function (x) {  
  y <- (x[2]  - 5.1/(4 * pi^2) * (x[1] ^2) + 5/pi * x[1]  - 6)^2 + 
    10 * (1 - 1/(8 * pi)) * cos(x[1] ) + 10
  if(x[3]==1)
    y <- y +1
  else if(x[3]==2)
    y <- y -1
  y  
}
## vectorize
objFun <- function(x){apply(x,1,braninFunctionFactor)}
set.seed(1)
res <- spot(fun=objFun,lower=c(-5,0,1),upper=c(10,15,3),     
    control=list(model=buildKriging, 
      types= c("numeric","numeric","factor"), 
      optimizer=optimLHD)) 
res$xbest
res$ybest
}
