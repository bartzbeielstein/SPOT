
parSpot <- function(x=NULL, fun, 
                    lower, upper, sequentialControlList = list(), 
                    parallelControlLists = list(list()), cl = NULL, nCores = NULL, ...){
    
    #Initial Input Checking
    initialInputCheck(x,lower,upper,sequentialControlList)
    
    ## Set up a parallel cluster if none is provided previously
    if(is.null(cl)){
        cl <- parallel::getDefaultCluster()
        if(is.null(cl)){
            if(is.null(nCores)){
                nCores <- parallel::detectCores()
            }
            parallel::setDefaultCluster(parallel::makeCluster(nCores))
        }
    }
    
    ## Filling sequentialControlList with standard spotControl
    ## 
    dimension <- length(lower)
    sequentialControlList <- spotFillControlList(sequentialControlList, dimension)
    
    ## Initial design generation
    set.seed(sequentialControlList$seedSPOT)
    x <- sequentialControlList$design(x=x,lower=lower,upper=upper,control=sequentialControlList$designControl)
    
    ## Rounding values produced by the design function to integers, etc.
    x <- repairNonNumeric(x,sequentialControlList$types)
    
    ## Evaluate initial design with objective function
    y <- objectiveFunctionEvaluation(x=NULL,xnew=x,fun=fun,seedFun=sequentialControlList$seedFun,noise=sequentialControlList$noise,...)
    
    #######
    ## 
    ## After the evaluation of the initial design each spot control should be run iteration by iteration in parallel
    ## 
    #######
    
    if(length(parallelControlLists) == 0){
        parallelControlLists <- list(spotFillControlList(list(), dimension))
    }else if(any(sapply(parallelControlLists, typeof) != "list")){
        error("parallelControlLists should be a list of lists.")
    }else{
        parallelControlLists <- lapply(parallelControlLists, function(cList){spotFillControlList(cList, dimension)})
    }
    
    totalBudget <- sequentialControlList$funEvals
    
    parallelSpotResults <- NULL
    parallel::clusterCall(cl = cl ,fun = function(){
        require(SPOT)
    })
    #vars <- ls()
    #vars <- vars[vars!="cl"]
    #parallel::clusterExport(cl = cl,varlist = ls())
    while((nrow(x) + length(parallelControlLists)) <= totalBudget){
        doParallelSpotIter <- function(cList){
            cList$funEvals <- nrow(x) + 1
            spotLoop(x=x,y=y,fun=function(x){x},lower=lower,upper=upper,control=cList,parallelCall = T, ...)
        }
        parallelSpotResults <- parallel::parLapply(cl = cl, parallelControlLists, doParallelSpotIter)
        newX <- do.call(rbind.data.frame, parallelSpotResults)
        newY <- objectiveFunctionEvaluation(x=x,xnew=newX,fun=fun,seedFun=sequentialControlList$seedFun,
                                            noise=sequentialControlList$noise,...)
        x <- rbind(x, newX)
        y <- rbind(y, newY)
    }
    
    indexBest <- which.min(y)
    parallelSpotResults$x <- x
    parallelSpotResults$y <- y
    parallelSpotResults$xbest <- x[indexBest,,drop=FALSE]
    parallelSpotResults$ybest <- y[indexBest,,drop=FALSE]
    parallelSpotResults$msg <- "budget exhausted"
    parallelSpotResults
}