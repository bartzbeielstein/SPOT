
###################################################################################################
#' Evolution Strategy 
#'
#' This is an implementation of an Evolution Strategy.
#' 
#' @param x optional start point, not used
#' @param fun objective function, which receives a matrix x and returns observations y
#' @param lower is a vector that defines the lower boundary of search space (this also defines the dimensionality of the problem)
#' @param upper is a vector that defines the upper boundary of search space (same length as lower)
#' @param control list of control parameters. The \code{control} list can contain the following settings:
#' \describe{
#' 			\item{funEvals}{number of function evaluations, stopping criterion, default is \code{500}}
#' 			\item{mue}{number of parents, default is \code{10}}
#' 			\item{nu}{selection pressure. That means, number of offspring (lambda) is mue multiplied with nu. Default is \code{10}}
#' 			\item{mutation}{string of mutation type, default is \code{1}}
#' 			\item{sigmaInit}{initial sigma value (step size), default is \code{1.0}}
#' 			\item{nSigma}{number of different sigmas, default is \code{1}}
#' 			\item{tau0}{number, default is \code{0.0}. tau0 is the general multiplier.}
#' 			\item{tau}{number, learning parameter for self adaption, i.e. the local multiplier for step sizes (for each dimension).default is \code{1.0}}
#' 			\item{rho}{number of parents involved in the procreation of an offspring (mixing number), default is \code{"bi"}}
#' 			\item{sel}{number of selected individuals, default is \code{1}}
#' 			\item{stratReco}{Recombination operator for strategy variables. \code{1}: none. \code{2}: dominant/discrete (default). \code{3}: intermediate. \code{4}: variation of intermediate recombination. }
#' 			\item{objReco}{Recombination operator for object variables. \code{1}: none. \code{2}: dominant/discrete (default). \code{3}: intermediate. \code{4}: variation of intermediate recombination. }
#' 			\item{maxGen}{number of generations, stopping criterion, default is \code{Inf}}
#' 			\item{seed}{number, random seed, default is \code{1}}
#' 			\item{noise}{number, value of noise added to fitness values, default is \code{0.0}}
#' 			\item{verbosity}{defines output verbosity of the ES, default is \code{0}}
#' 			\item{plotResult}{boolean, specifies if results are plotted, default is \code{FALSE}}
#' 			\item{logPlotResult}{boolean, defines if plot results should be logarithmic, default is \code{FALSE}}
#' 			\item{sigmaRestart}{number, value of sigma on restart, default is \code{0.1}}
#' 			\item{preScanMult}{initial population size is multiplied by this number for a pre-scan, default is \code{1}}
#' 			\item{globalOpt}{termination criterion on reaching a desired optimum value, default is \code{rep(0,dimension)}}
#'}
#' @param ... additional parameters to be passed on to \code{fun}
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
#' cont <- list(funEvals=100)
#' optimES(fun=funSphere,lower=rep(0,2), upper=rep(1,2), control= cont) 
#'
#' @export
###################################################################################################
optimES <- function(x=NULL
                      , fun
                      , lower
                      , upper
                      , control=list()
											, ...){					  
	con <- list(
		funEvals=100,
		seed=1,
		mue = 10,
		nu = 10,
		mutation = 2,
		sigmaInit = 1.0,
		nSigma = 1,
		tau0 = 0.0,
		tau = 1.0,
		rho = "bi",
		sel = -1,
		stratReco = 1,
		objReco = 2,
		maxGen = Inf,
		noise = 0.0, 
		verbosity=0,
		plotResult=FALSE,
		logPlotResult=FALSE,
		sigmaRestart = 0.1,
		preScanMult= 1);
	con[names(control)] <- control
	control <- con;
	
	force(fun)
	## vectorization wrapper
  fn <- function(x,...)fun(matrix(x,1),...) 
	
	dimension <- length(lower)
	
### call SPOT	
	esResult <- spotAlgEs(
		mue = control$mue,
		nu = control$nu,
		dimension = dimension,
		mutation = control$mutation,
		sigmaInit = control$sigmaInit,
		nSigma = control$nSigma,
		tau0 = control$tau0,
		tau = control$tau,
		rho = control$rho,
		sel = control$sel,
		stratReco = control$stratReco,
		objReco = control$objReco,
		maxGen = control$maxGen,
		maxIter = control$funEvals,
		seed = control$seed,
		noise = control$noise,
		fName = fn,
		lowerLimit = lower,
		upperLimit = upper,
		verbosity = control$verbosity,
		plotResult = control$plotResult,
		logPlotResult = control$logPlotResult,
		sigmaRestart = control$sigmaRestart,
		preScanMult = control$preScanMult,
		globalOpt = control$globalOpt,		
		...)
	
	### return Results	
	list(x=NULL,y=NULL,xbest=matrix(esResult$X,1),ybest=matrix(esResult$Y,1),
				count= esResult$counts) 
	
}




### Main Loop #######################################################  

###################################################################################################
#' Evolution Strategy Implementation
#'
#' This function is used by \code{\link{optimES}} as a main loop for running
#' the Evolution Strategy with the given parameter set specified by SPOT.
#'
#' @param mue number of parents, default is \code{10}
#' @param nu selection pressure. That means, number of offspring (lambda) is mue multiplied with nu. Default is \code{10}
#' @param dimension dimension number of the target function, default is \code{2}
#' @param mutation mutation type, either \code{1} or \code{2}, default is \code{1}
#' @param sigmaInit initial sigma value (step size), default is \code{1.0}
#' @param nSigma number of different sigmas, default is \code{1}
#' @param tau0 number, default is \code{0.0}. tau0 is the general multiplier.
#' @param tau number, learning parameter for self adaption, default is \code{1.0}. tau is the local multiplier for step sizes (for each dimension).
#' @param rho number of parents involved in the procreation of an offspring (mixing number), default is \code{"bi"}
#' @param sel number of selected individuals, default is \code{1}
#' @param stratReco Recombination operator for strategy variables. \code{1}: none. \code{2}: dominant/discrete (default). \code{3}: intermediate. \code{4}: variation of intermediate recombination. 
#' @param objReco Recombination operator for object variables. \code{1}: none. \code{2}: dominant/discrete (default). \code{3}: intermediate. \code{4}: variation of intermediate recombination. 
#' @param maxGen number of generations, stopping criterion, default is \code{Inf}
#' @param maxIter number of iterations (function evaluations), stopping criterion, default is \code{100}
#' @param seed number, random seed, default is \code{1}
#' @param noise number, value of noise added to fitness values, default is \code{0.0}
# @param thrs threshold string, default is \code{"no"}
# @param thrsConstant number, default is \code{0.0}
#' @param fName function, fitness function, default is \code{\link{funSphere}}
#' @param lowerLimit number, lower limit for search space, default is \code{-1.0}
#' @param upperLimit number, upper limit for search space, default is \code{1.0}
#' @param verbosity defines output verbosity of the ES, default is \code{0}
#' @param plotResult boolean, asks if results are plotted, default is \code{FALSE}
#' @param logPlotResult boolean, asks if plot results should be logarithmic, default is \code{FALSE}
#' @param sigmaRestart number, value of sigma on restart, default is \code{0.1}
#' @param preScanMult initial population size is multiplied by this number for a pre-scan, default is \code{1}
#' @param globalOpt termination criterion on reaching a desired optimum value, should be a vector of length dimension (LOCATION of the optimum). Default to NULL, which means it is ignored.
# @param conf config number passed to the result file, default is \code{-1}
#' @param ... additional parameters to be passed on to \code{fName}
#'
#' @export
###################################################################################################
spotAlgEs <- function(mue = 10,
											nu = 10,
											dimension = 2,
											mutation = 2,
											sigmaInit = 1.0,
											nSigma = 1,
											tau0 = 0.0,
											tau = 1.0,
											rho = "bi",
											sel = -1,
											stratReco = 1,
											objReco = 2,
											maxGen = Inf,
											maxIter = Inf,
											seed = 1,
											noise = 0.0,
											#thrs = "no",
											#thrsConstant = 0.0,
											fName = funSphere,
											lowerLimit = -1.0,
											upperLimit = 1.0,
											verbosity=0,
											plotResult=FALSE,
											logPlotResult=FALSE,
											sigmaRestart = 0.1,
											preScanMult= 1,
											globalOpt=NULL,
											...){  
	### Parameter corrections
	lambda <- round(mue*nu)
	mue <- round(mue)
	nSigma <- round(nSigma)
	### Currently, there are 4 recombination operators:
	stratReco <- max(round(stratReco) %% 5,1) 
	objReco <- max(round(objReco) %% 5,1) 	
	if(nSigma>dimension){
		nSigma <- dimension # correction: if nSigma > dimension
		warning("In Evolution Strategy (spotAlgEs): nSigma should not be larger than dimension.")
	}
	lambda <- max(mue,lambda) # correction: if lambda < mue
	if(plotResult==TRUE){
		logPlotResult=FALSE
	}
	###
	if(nSigma==1){
		tau0=0.0
		###  tau <- tau/sqrt(dimension)
	}
	#else{
	#  tau0 <- tau/sqrt(2*dimension)
	#  tau <- tau/sqrt(2*sqrt(dimension))
	#}
	###  necessary for (1,+\mue):
	if(rho=="bi") rhoVal=min(2,mue) else rhoVal=mue   
	set.seed(seed)
	parentPop <- NULL
	gen <- 0
	#succ <- 1
	if(verbosity>=1){
		sigmaAvgList <- rep(sigmaInit,nSigma)
		sigmaMedList <- rep(sigmaInit,nSigma)
	}
	bestFitness <- NULL
	realBest <-  NULL
	realBestPar <- NULL
	allTimeBest <- NULL
	alg.currentBest <- 0.0
	###
	## Perform pre-scan:
	## the initial population size is multiplied by preScanMult, say 10. Then the mue best ind are
	## selected from 10*mue individuals
	#if (gen == 0)
	###
	parentPop <- spotAlgEsInitParentPop(sigmaInit, dimension, nSigma, noise, fName, gen, lowerLimit, upperLimit, round(mue*preScanMult))
	iter <- nrow(parentPop) #number of initial function evaluations
	parentPop <- data.frame(parentPop)
	parentPop <- parentPop[order(parentPop$fitness),]
	###
	if(verbosity==2) print(parentPop)	
	bestInd <- parentPop[1,]
	bestFitness <- parentPop$fitness[[1]]
	allTimeBest <- bestFitness
	realBest <- parentPop$realFitness[[1]]
	parentPop<-parentPop[1:mue,]
	### End pre-scan
	
	#browser()
	##while(iter < maxIter && gen < maxGen ){
	while(alg.currentBest < Inf & spotAlgEsTermination(it=iter, maxIt=maxIter, ge=gen, maxGe=maxGen, xk=bestInd[1:dimension], xOpt=globalOpt)){
		gen <- gen +1
		###
		if(verbosity==2)   print(paste("Gen:", gen))
		
		offspringPop <- NULL
		for(i in 1:lambda){
			marriagePop <- NULL
			marriagePop <- spotAlgEsMarriage(parentPop, rhoVal)
			###
			if(verbosity==2)   print(paste("MarriagePop for Offspring:", i))
			if(verbosity==2)   print(marriagePop)
			# Recombination
			stratRecombinant <- switch(stratReco,
																 as.matrix(marriagePop[1,I(dimension+1):I(dimension+nSigma)]), ### 1 = perform no reco
																 spotAlgEsDominantReco(marriagePop, rhoVal, dimension, nSigma, objType="strat"),
																 spotAlgEsInterReco(marriagePop, rhoVal, dimension, nSigma, objType="strat"),
																 spotAlgEsInterRecoBeSw02(marriagePop, dimension, nSigma, objType="strat")
			)
			
			objRecombinant <- switch(objReco,
															 as.matrix(marriagePop[1,1:dimension]), ### 1 = perform no reco
															 spotAlgEsDominantReco(marriagePop, rhoVal, dimension, nSigma, objType="obj"),
															 spotAlgEsInterReco(marriagePop, rhoVal, dimension, nSigma, objType="obj"),					
															 spotAlgEsInterRecoBeSw02(marriagePop, dimension, nSigma, objType="obj")
			)
			###
			if(verbosity==2)   print(paste("stratRecombinant:", i))
			if(verbosity==2)   print(stratRecombinant)
			if(verbosity==2)   print(paste("objRecombinant:", i))
			if(verbosity==2)   print(objRecombinant)
			
			sigmaNew <- switch(mutation,
												 stratRecombinant, ### 1 = perform no mutation
												 spotAlgEsStratMutation(stratRecombinant, tau0, tau, sigmaRestart, sigmaInit))
			####
			xNew <- switch(mutation,
										 objRecombinant, ### 1 = perform no mutation
										 spotAlgEsObjMutation(objRecombinant,sigmaNew))
			###
			if(verbosity==2)   print(paste("sigmaNew:", i))
			if(verbosity==2)   print(sigmaNew)
			if(verbosity==2)   print(paste("xNew:", i))
			if(verbosity==2)   print(xNew)
			
			### repair lower / upper limits
			xNew <- pmax(xNew,lowerLimit)
			xNew <- pmin(xNew,upperLimit)
			
			### evaluate
			result <- as.numeric(fName(xNew,...))
			realFitness <- result
			fitness <- result
			
			iter <- iter +1
			offspring <- NULL
			offspring <- c(x= xNew,
										 sigma= sigmaNew ,
										 realFitness=realFitness <-  realFitness,
										 fitness=fitness <- fitness,
										 generation= gen)
			###
			if(verbosity==2)   print(offspring)
			offspringPop <- rbind(offspringPop, offspring)
			if(!spotAlgEsTermination(it=iter, maxIt=maxIter, ge=gen, maxGe=maxGen, xk=bestInd[1:dimension], xOpt=globalOpt)) break;
		}
		row.names(offspringPop) <- 1:nrow(offspringPop)
		offspringPop <- data.frame(offspringPop)
		###
		if(verbosity==2)   print("OffspringPop")
		if(verbosity==2)   print(offspringPop)
		parentPop <- spotAlgEsSelection(parentPop,offspringPop, sel, gen, iter, maxIter)
		
		###
		if(verbosity==2)   print("parentPop")
		if(verbosity==2)   print(parentPop)
		
		#succ <- spotAlgEsGetSuccessRate(gen,parentPop)
		alg.currentBest <- parentPop$fitness[[1]]
		currentReal <- parentPop$realFitness[[1]]
		currentPar <- as.numeric(parentPop[setdiff(names(parentPop),c("sigma","realFitness","fitness","generation"))][1,])
		
		if(verbosity>=1) {
			#      currentSigma<- parentPop[1,I(dimension+1)]
			#      currentSigma<- parentPop[1,I(dimension+1):I(dimension+nSigma)]
			currentSigma<- parentPop[1,I(dimension+1):I(dimension+nSigma)]
			if(verbosity==2){
				print(c(alg.currentBest, currentReal, currentSigma))
			}
			sigmaAvgList <- c(sigmaAvgList, mean(unlist(currentSigma)))
			if(verbosity==2){
				cat("sigmaAvgList:",sigmaAvgList,"\n")
			}
			sigmaMedList <- c(sigmaMedList, median(unlist(currentSigma)))
			if(verbosity==2){
				cat("sigmaMedList:",sigmaMedList,"\n")
			}
		}
		# store only the doe values
		if(verbosity==0){
			realBest <- currentReal
			realBestPar <- currentPar # TODO also log par for alltime best
			bestFitness <- alg.currentBest 
			if (alg.currentBest < allTimeBest){
				allTimeBest <- alg.currentBest
				bestInd <- parentPop[1,]
			}
		}
		else{
			## log every generation:
			##
			realBest <-  c(realBest, currentReal)
			realBestPar <- rbind(realBestPar, currentPar)  
			allBest <- allTimeBest[length(allTimeBest)]
			bestFitness <- c(bestFitness, alg.currentBest)
			####
			####cat(gen, "fitness: ", alg.currentBest, "\n")
			if (alg.currentBest < allBest) {
				allTimeBest <- c(allTimeBest, alg.currentBest)
				bestInd <- parentPop[1,]
			}
			else{
				allTimeBest <- c(allTimeBest, allBest)
			}
		}
		
		#### Plotting ##################################
		if(sel==-1){stratName="Plus"} else{	stratName="Kappa"}
		
		if(verbosity>=1 && plotResult==TRUE){
			par(mfrow=c(2,1))
			#xlim=c(1,maxIter/lambda)
			plot(bestFitness, ylim=c(min(min(bestFitness), min(realBest)), max(max(bestFitness), max(realBest))), type="l",
					 xlab = paste( "Generation: " , as.character(gen), "Fitness: ", as.character(alg.currentBest)),
					 main = paste(as.character(mue), stratName, as.character(lambda), ", nSgm: ",as.character(nSigma), ", tau0: ", as.character(tau0), ", tau: ", as.character(tau)  )
			)
			lines(realBest, col="red", type="b")
			plot(sigmaAvgList, type="l", col="green", xlab = paste( "Generation: " , as.character(gen), "AvgSgm: ",
																															as.character(sigmaAvgList[length(sigmaAvgList)])))
			lines(sigmaMedList, col="blue", type="l")
			## TBB 25 2 2009:
			##
		}
		####
		if(verbosity>=1 && logPlotResult==TRUE){
			par(mfrow=c(2,1))
			##xlim=c(1,maxIter/lambda)
			plot(log(bestFitness), ylim=c(min(min(log(bestFitness)), min(log(realBest))), max(max(log(bestFitness)), max(log(realBest)))), type="b",
					 xlab = paste( "Generation: " , as.character(gen), "Fitness: ", as.character(alg.currentBest)),
					 main = paste(as.character(mue), stratName, as.character(lambda), ", nSgm: ",as.character(nSigma), ", tau0: ", as.character(tau0), ", tau: ", as.character(tau)  )           
			)
			lines(log(realBest), col="red", type="b")
			plot(log(sigmaAvgList), type="l", col="green")
			lines(log(sigmaMedList), col="blue", type="l")
		}
	}
	if(verbosity>=1){
		print(list(real=realBest,
							 realPar=realBestPar,
							 best=bestFitness,
							 allTime=allTimeBest,
							 bestInd=bestInd))
	}
	##print(  realBest[[length(realBest)]])
	##
	##data frame written to the DOE file:
	##(1) bestInd is the allTimeBestInd, and *not* the bestInd in the last generation.
	##(2) realBest is the fitness of the bestInd (that is not necessarily a member of the last generation)
	##(3) allTimeBest is the best fitness found during the whole optimization run, the (noisy) fitness of the bestInd. 
	##(4) NoisyFitness is the best Fitness in the last generation.
	
	## special treatment for factors. Their indices should be written as strings.
	#recoType <- c("no", "disc","inter", "interRecoBeSw02")
	#OBJRECO = paste('"',toString((1:length(recoType))[recoType==objReco]),'"', sep="")
	#STRATRECO = paste('"',toString((1:length(recoType))[recoType==stratReco]),'"', sep="")
	## TODO: use which for the following selection:
	#OBJRECO = toString((1:length(recoType))[recoType==objReco])
	#STRATRECO = toString((1:length(recoType))[recoType==stratReco])
	if(is.matrix(realBestPar)){
		X=realBestPar[nrow(realBestPar),]
	}else{
		X=realBestPar
	}	
	#return:
	list(Y=realBest[[length(realBest)]], # last value     #TODO: should be alltime best?
			 X=X,
			 NoisyFitness=bestFitness[[length(bestFitness)]],
			 AllTimeBest=allTimeBest[[length(allTimeBest)]],
			 Generation=bestInd$generation,
			 Percentage=bestInd$generation/gen*100,
			 counts=iter)
}

### Termination #################################################################
###
#' Termination
#'
#' Handles the termination functions for the ES. 
#'
#' @param it iteration (function evaluations)
#' @param maxIt Maximum number of iterations (function evaluations)
#' @param ge generation
#' @param maxGe Maximum number of generations
#' @param xk current best value of the optimization run
#' @param xOpt target value of the optimization run
#'
#' @return \code{boolean} \cr
#' - TRUE as long as the current value has not yet reached its limit. Once the
#' given termination criterion is reached the function returns FALSE.
#'
#' @seealso  \code{\link{spotAlgEs}} \code{\link{spotAlgEsHps}}
#' @keywords internal
###
spotAlgEsTermination <- function(it, maxIt, ge, maxGe, xk, xOpt=NULL){
	continue <- (it < maxIt) & (ge < maxGe)
	if(!is.null(xOpt))
		continue <- continue & spotAlgEsHps(xk,xOpt)
	continue		
}

###
#' Termination hps
#'
#' Termination function for the ES. Terminates the ES at a given target value.
#'
#' @param xk current best value of the optimization run
#' @param xOpt target value of the optimization run
#'
#' @return \code{boolean} \cr
#' - TRUE as long as the current value has not yet reached its limit. Once the
#' given termination criterion is reached the function returns FALSE.
#' @keywords internal
###  
spotAlgEsHps <- function(xk, xOpt){
	n <- length(xk)
	tfVec <- 10*sqrt(n)*(xk-xOpt)<=1
	is.element(FALSE, tfVec)
}

###
#' get Success Rate
#'
#' This function determines the success rate.
#'
#' @param gen Generation number
#' @param pop Population to be evaluated (only the generation \code{gen} will be evaluated)
#'
#' @return number \cr
#' - the success rate of the given population
#'
#' @seealso  \code{\link{spotAlgEs}} 
#' @keywords internal
###
spotAlgEsGetSuccessRate <- function(gen,pop){
	succ <- length(subset(pop, pop$generation==gen)[,1])
	len <- length(pop[,1])
	succ/len
}  

### Initialization ########################################################

###
#' Individual Initialization
#'
#' Creates a new Individual for the Evolution Strategy implemented in SPOT.
#'
#' @param s sigma, step size
#' @param n n, number of diff. step sizes
#' @param dimension number of target function dimension
#' @param noise noise to be added to fitness value
#' @param fName target function
#' @param gen generation
#' @param low lower limit
#' @param high upper limit
#' @param des des scaling for placement between low and high
#' @param ... additional parameters to be passed on to \code{fName}
#'
#' @return numeric vector \cr
#' - contains x value, sigma value, real fitness value, fitness with noise, and generation number
#'
#' @seealso  \code{\link{spotAlgEs}}
#' @keywords internal
###
spotAlgEsIndividualInitial <- function(s,dimension,n,noise=0,fName,gen,low=-1.0,high=1.0,	des,...){
	x <- low + (high-low)*des
	sigma <- rep(s,n)
	realFitness <- as.numeric(fName(x,...))
	fitness <- realFitness + rnorm(1,0,noise)
	c(x = x,     #  *runif(dimension),
		sigma = sigma,
		realFitness =  realFitness,#spotAlgEsF(x,fName),
		fitness = fitness,
		generation = gen
	)
}

###
#' Initialize Parent Population
#'
#' Creates initial parent population
#'
#' @param sigmaInit initial sigma value (standard deviation)
#' @param dimension number of target function dimension
#' @param nSigma number of standard deviations
#' @param noise noise to be added to fitness value
#' @param fName target function
#' @param gen generation
#' @param low lower limit
#' @param high upper limit
#' @param mue number of parents in the ES
#' @param ... additional parameters to be passed on to \code{fName}
#'
#' @return matrix \cr
#' - holds the parent population created by this function
#'
#' @seealso  \code{\link{spotAlgEs}} \code{\link{spotAlgEsIndividualInitial}}
#' @keywords internal
###
spotAlgEsInitParentPop <- function(sigmaInit, dimension, nSigma, noise, fName, gen, low, high, mue,...)
{
	parentPop <- NULL	
	ld <- designLHDNorm(dimension,mue,FALSE)$design
	###ld <- optimumLHS(mue, dimension, 2, 0.1)
	for(i in 1:mue){      
		parentPop <- rbind(parentPop,
											 spotAlgEsIndividualInitial(s=sigmaInit,
											 													 dimension=dimension,
											 													 n=nSigma,
											 													 noise=noise,
											 													 fName=fName,
											 													 gen=gen,
											 													 low=low,
											 													 high=high,
											 													 des = ld[i,],...))
	}                                        
	parentPop
}

### Recombination ######################################################
###
#' Marriage with replace
#'
#' Recombination function for the Evolution Strategy.
#'
#' @param pop Population
#' @param rhoVal number of parents involved in the procreation of an offspring
#'
#' @return \code{pop} \cr
#' - \code{pop} Population
#'
#' @seealso  \code{\link{spotAlgEs}} \code{\link{spotAlgEsMarriageWithReplace}} 
#' @keywords internal
###
spotAlgEsMarriageWithReplace <- function(pop,rhoVal){
	pop[sample(nrow(pop), rhoVal, replace=TRUE),]
}

###
#' Marriage
#'
#' Recombination function for the Evolution Strategy.
#'
#' @param pop Population
#' @param rhoVal number of parents involved in the procreation of an offspring
#'
#' @return \code{pop} \cr
#' - \code{pop} Population
#'
#' @seealso  \code{\link{spotAlgEs}} \code{\link{spotAlgEsMarriage}} 
#' @keywords internal
###
spotAlgEsMarriage <- function(pop,rhoVal){
	pop[sample(nrow(pop), rhoVal, replace=FALSE),]
}

###
#' spotAlgEsInterRecoBeSw02
#'
#' Recombination function for the Evolution Strategy.
#'
#' @param parents Parent individuals
#' @param dimension number of dimensions
#' @param nSigma number of standard deviations
#' @param objType string, default is "obj"
#'
#' @seealso \code{\link{spotAlgEs}} \code{\link{spotAlgEsInterReco}} \code{\link{spotAlgEsDominantReco}} 
#' @keywords internal
###
spotAlgEsInterRecoBeSw02 <- function(parents, dimension, nSigma, objType="obj"){
	rObject <- NULL
	if(objType=="obj"){
		for(i in 1:dimension)  rObject <- c(rObject, mean(parents[,i]))
	}
	else{
		for(i in 1:nSigma)  rObject <- c(rObject, mean(parents[,I(dimension+i)]))
	}
	rObject
}

###
#' spotAlgEsInterReco
#'
#' Recombination function for the Evolution Strategy.
#'
#' @param parents Parent individuals
#' @param rhoVal number of parents involved in the procreation of an offspring
#' @param dimension number of dimensions
#' @param nSigma number of standard deviations
#' @param objType string, default is "obj"
#'
#' @seealso \code{\link{spotAlgEs}} \code{\link{spotAlgEsInterRecoBeSw02}} \code{\link{spotAlgEsDominantReco}} 
#' @keywords internal
###
spotAlgEsInterReco <- function(parents, rhoVal, dimension, nSigma, objType="obj"){
	rObject <- NULL
	if(objType=="obj"){
		for(i in 1:dimension){
			select2 <- sample(rhoVal,min(2,rhoVal)) # necessary for (1,+\mue)
			# if(verbosity==2)   print(paste("Selected for XReco: ", select2))
			rObject <- c(rObject, mean(parents[select2,i]))
		}
	}
	else{
		for(i in 1:nSigma){
			select2 <- sample(rhoVal,min(2,rhoVal)) # necessary for (1,+\mue)
			# if(verbosity==2)   print(paste("Selected for SReco: ", select2))
			rObject <- c(rObject, mean(parents[select2,I(dimension+i)]))
		}
	}
	rObject
}

###
#' spotAlgEsDominantReco
#'
#' Recombination function for the Evolution Strategy.
#'
#' @param parents Parent individuals
#' @param rhoVal number of parents involved in the procreation of an offspring
#' @param dimension number of dimensions
#' @param nSigma number of standard deviations
#' @param objType string, default is "obj"
#'
#' @seealso \code{\link{spotAlgEs}} \code{\link{spotAlgEsInterRecoBeSw02}}  \code{\link{spotAlgEsInterReco}} 
#' @keywords internal
###
spotAlgEsDominantReco <- function(parents, rhoVal, dimension, nSigma, objType="obj"){
	rObject <- NULL
	if(objType=="obj"){
		for (i in 1:dimension) rObject <- c(rObject, parents[sample(rhoVal,1),i])
	}
	else
	{
		for (i in (1:nSigma))
			rObject <- c(rObject, parents[sample(rhoVal,1),I(dimension+i)])
	}
	rObject
}


### Mutation #############################################################

###
#' spotAlgEsStratMutation
#'
#' Mutation function for the ES strategy parameter sigma. 
#'
#' @param strat Strategy parameter (sigma), can be a vector or a single number
#' @param tau0 the global/general step size multiplier for self adaption
#' @param tau the step size multiplier for self adaption in each dimension
#' @param sigmaRestart sigma values are reset to initial values when a uniformly distributed random number is smaller than sigmaRestart
#' @param sigmaInit initial sigma value
#'
#' @return number \code{s} \cr
#' - \code{s} is the new sigma value
#' @keywords internal
###
spotAlgEsStratMutation <- function(strat, tau0, tau, sigmaRestart, sigmaInit){
	if (runif(1) < sigmaRestart){
		s <- rep(sigmaInit, length(strat))
	}
	else{
		s <- exp(tau0*rnorm(1,0,1))*as.numeric(strat)*exp(tau*rnorm(length(strat),0,1))
	}
	s
}
###
#' spotAlgEsObjMutation
#'
#' Mutation function for the ES individual parameters
#'
#' @param obj Object to be mutated
#' @param strat Strategy parameter sigma to mutate with
#' @keywords internal
###
spotAlgEsObjMutation <- function(obj, strat){
	obj + as.numeric(strat)*rnorm(length(obj),0,1)
}

### Selection ##################################################################
###
#' spotAlgEsSelection
#'
#' Selection function for the ES
#'
#' @param parentPop Parent population
#' @param offspringPop Offspring population
#' @param sel number of individuals to be selected
#' @param gen generation number
#' @param iter current iteration number
#' @param maxIter maximum iteration number
#' @keywords internal
###      
spotAlgEsSelection <- function(parentPop, offspringPop, sel, gen, iter, maxIter){
	mue <- nrow(parentPop)
	names(offspringPop) <- names(parentPop)
	if(sel==-1)
	{ # plus selection
		parentPop <- rbind(parentPop,offspringPop)
	}
	else
	{ ## kappa selection 
		##parentPop <- parentPop[parentPop$generation>I(gen-sel)]
		##parentPop <- rbind(parentPop,offspringPop)
		## print(c(gen, sel))
		## print(parentPop)
		parentPop <- parentPop[parentPop[,"generation"] > gen - sel,]
		## cat("reduced parentPop: \n")
		## print(parentPop)
		parentPop <- rbind(parentPop,offspringPop)
		## cat("new parentPop: \n")
		## print(parentPop)
		
	}
	parentPop <- parentPop[order(parentPop$fitness),]
	row.names(parentPop) <- 1:nrow(parentPop)
	parentPop[1:mue,]
}


