
parSpot <- function(x=NULL, fun, 
                    lower, upper, sequentialControlList = list(), parallelControlLists = list(), cl = NULL, nCores = NULL, ...){
    
    #Initial Input Checking
    initialInputCheck(x,fun,lower,upper,sequentialControlList)
    
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
    ## After the evaluation of the initial design each spot control should be ran iteration by iteration in parallel
    ## 
    #######
    
    parallelControlLists <- lapply(parallelControlLists, function(cList){spotFillControlList(cList, dimension)})
    
    totalBudget <- sequentialControlList$funEvals
    
    parallelSpotResults <- NULL
    while((nrow(x) + length(parallelControlLists)) <= totalBudget){
        doParallelSpotIter <- function(cList){
            cList$funEvals <- nrow(x) + 1
            spotLoop(x=x,y=y,fun=fun,lower=lower,upper=upper,control=cList, ...)
        }
        parallelSpotResults <- parallel::parLapply(cl = cl, parallelControlLists, doParallelSpotIter)
        
        amntRows <- nrow(x)
        for(res in parallelSpotResults){
            x <- rbind(x, res$x[-c(1:amntRows), , drop=F])
            y <- rbind(y, res$y[-c(1:amntRows), , drop=F])
        }
    }
    
    indexBest <- which.min(y)
    parallelSpotResults$x <- x
    parallelSpotResults$y <- y
    parallelSpotResults$xbest <- x[indexBest,,drop=FALSE]
    parallelSpotResults$ybest <- y[indexBest,,drop=FALSE]
    parallelSpotResults$msg <- "budget exhausted"
    parallelSpotResults
}