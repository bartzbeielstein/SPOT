
#' infillExpectedImprovement
#' 
#' Expected Improvement infill criterion that can be passed to control$modelControl$infillCriterion in order to be used during the optimization in SPOT.
#' 
#' @param predictionList The results of a predict.model call
#' @param model The surrogate model which was used for the prediction
#'
#' @return numeric vector, expected improvement results
#' @export
#'
#' @examples
#' spot(,funSphere,c(-2,-3),c(1,2), control = list(infillCriterion = infillExpectedImprovement, modelControl = list(target = c("y","s"))))
infillExpectedImprovement <- function(predictionList, model){
    budgetScaling <- 1
    if(!is.null(model$globalBudget)){
        if((model$globalBudget-model$globalInitialBudget) == 1 | (nrow(model$x) + 1) >= model$globalBudget){
            return(predictionList$y)
        }
        budgetScaling <- min((nrow(model$x)-model$globalInitialBudget)/(model$globalBudget-model$globalInitialBudget-1),1)
    }
    
    mean <- predictionList$y
    sd <- predictionList$s
    modelMin <- min(model$y)
    
    if(model$eiUseWeightedBudgetSum && (!is.null(model$globalBudget))){
        EITermOne=(modelMin-mean)*pnorm((modelMin-mean)/sd)
        EITermTwo=sd*(1/sqrt(2*pi))*exp(-(1/2)*((modelMin-mean)^2/(sd^2)))
        eiTerm <- EITermOne + EITermTwo
        
        mean <- (mean-min(model$y))/(max(model$y)-min(model$y))
        maxEI <- sqrt(model$ssq) * (1/sqrt(2*pi))
        eiTerm <- (eiTerm)/(maxEI)
        
        return(mean * budgetScaling - eiTerm * (1 - budgetScaling))
    }
    
    sd <- sd * budgetScaling
    
    return(expectedImprovement(mean,sd,modelMin))
}