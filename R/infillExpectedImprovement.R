
infillExpectedImprovement <- function(predictionList, model){
    mean <- predictionList$y
    if("slinear" %in% tolower(model$target)){
        sd <- predictionList$sLinear
    }else{
        sd <- predictionList$s
    }
    modelMin <- min(model$y)
    
    return(expectedImprovement(mean,sd,modelMin))
}