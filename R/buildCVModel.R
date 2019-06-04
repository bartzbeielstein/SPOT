
buildCVModel <- function(x, y, modellingFunction, control=list()){
    ## Load Control list
    con<-list(nFolds = 10)
    con[names(control)] <- control
    control<-con
    
    control$nFolds <- min(control$nFolds, nrow(x))
    
    ## Empty List for the final model
    cvModel <- list()
    
    ## Save initial x and y in model
    cvModel$x <- x
    cvModel$y <- y
    
    #Randomly shuffle the data
    shuffleIndexes <- sample(nrow(x))
    x <- x[shuffleIndexes,,drop = F]
    y <- y[shuffleIndexes, drop = F]
    
    #Create nFolds equally sized folds
    folds <- cut(seq(1,nrow(x)),breaks=control$nFolds,labels=FALSE)
    
    createSingleModel <- function(i){
        leaveOutIndex <- which(folds==i,arr.ind=TRUE)
        trainX <- x[-leaveOutIndex,, drop=F]
        trainY <- as.matrix(y[-leaveOutIndex])
        model <- modellingFunction(trainX, trainY, control = control)
        return(model)
    }
    
    cvModel$models <- lapply(1:control$nFolds, createSingleModel)
    class(cvModel)<- "cvModel"
    return(cvModel)
}

maxNearestNeighbourDistance <- function(x){
    minDists <- NULL
    for(i in 1:nrow(x)){
        currentDists <- abs(t(t(x[-i,, drop=F])-x[i,]))
        minDists <- c(minDists,sqrt(min(apply(currentDists,1,function(x){sum(x^2)}))))
    }
    return(sqrt(max(minDists)))
}

linearAdaptedSE <- function(sOld, newdata, x){
    ifelse(is.null(nrow(newdata)),nr <- 1,nr <- nrow(newdata))
    if(nr <= 1){
        for(i in 1:length(newdata)){
            minDist <- min(abs(x-newdata[i]))
            sOld[i] <- sOld[i] * minDist/max(diff(sort(x)))
        }
    }else{
        for(i in 1:nrow(newdata)){
            dists <- abs(t(t(x)-newdata[i,]))
            minDist <- sqrt(min(apply(dists,1,function(x){sum(x^2)})))
            sOld[i] <- sOld[i] * minDist/maxNearestNeighbourDistance(x)
        }
    }
    sOld *2 # ?
}

predict.cvModel <- function(object,newdata,...){
    predictSingle <- function(model){
        return(predict(model,newdata,...)$y)
    }
    
    results <- list()
    results$all <- sapply(object$models,predictSingle)
    
    ifelse(is.null(nrow(results$all)),nr <- 1,nr <- nrow(results$all))
    if(nr > 1){
        results$y <- apply(results$all,1,mean)
        
        funSE <- function(x){
            sd(x)/sqrt(length(x))
        }
        results$s <- apply(results$all,1 , funSE)
        results$sLinear <- linearAdaptedSE(results$s, newdata, object$x)
    }else{
        results$y <- mean(results$all)
        results$s <- sd(results$all)/sqrt(length(results$all))
        
        results$sLinear <- linearAdaptedSE(results$s, newdata, object$x)
    }
    
    return(results)
}

singleDimPlot <- function(f, lower, upper, withSE = T, useAdapted = F){
    x <- seq(lower, upper, length.out = 200)
    y <- f(x)$y
    if(!useAdapted){
        s <- f(x)$s 
    }else{
        s <- f(x)$sLinear
    }
    
    loadNamespace("ggplot2")
    if(!withSE){
        ggplot(data=data.frame(x,y), aes(x=x, y=y, group=1)) +
            geom_line()
    }else{
        ggplot(data=data.frame(x,s), aes(x=x)) +
            geom_line(aes(y = y)) + 
            geom_ribbon(aes(ymax=y + s, ymin=y - s), fill="pink", alpha=.5)
    }
}

simpleTest <- function(){
    x <- matrix(runif(20), ncol = 2)
    y <- funSphere(x)
    
    krig <- buildKriging(x,y)
    predict(krig, matrix(c(0.5,1,2,3,3,4), ncol = 2))
    plotModel(krig)
    
    cvMod <- buildCVModel(x,y,buildKriging)
    pred <- predict(cvMod, matrix(c(0.5,1), ncol = 2))
    plotModel(cvMod)
    
    ## Single Dimension SE Plot:
    ## 
    x <- matrix(c(1,2,3,5,7,8,9), ncol = 1)
    y <- matrix(sin(x))
    krig <- buildKriging(x,y)
    krig$target <- c("y","s","ei")
    
    f <- function(x){
        predict(krig,x)
    }
    singleDimPlot(f, min(x), max(x))
    
    cvMod <- buildCVModel(x,y,buildKriging)
    pred <- predict(cvMod, c(0.5,0.6))
    f <- function(x){
        predict(cvMod,x)
    }
    singleDimPlot(f, min(x), max(x))
    
    singleDimPlot(f, min(x), max(x), useAdapted = T)
    
    ## Single dim plot on multi dim fun
    ## 
    x <- matrix(runif(20), ncol = 2)
    y <- spotGUI:::funRastrigin(x)
    
    par(mfrow=c(2,2))
    
    krig <- buildKriging(x,y)
    plotModel(krig)
    krig$target <- c("y","s","ei")
    f <- function(p){
        predict(krig,matrix(c(p,rep(min(x[,2]),length(p))),ncol=2))
    }
    singleDimPlot(f, min(x), max(x))
    
    cvMod <- buildCVModel(x,y,buildKriging)
    plotModel(cvMod)
    f <- function(p){
        predict(cvMod,matrix(c(p,rep(min(x[,2]),length(p))),ncol=2))
    }
    singleDimPlot(f, min(x), max(x), useAdapted = T)
}