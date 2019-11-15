
#' infillLCB
#' 
#' Lower Confidence Bound infill criterion that can be passed to control$modelControl$infillCriterion in order to be used during the optimization in SPOT.
#' 
#' @param predictionList The results of a predict.model call
#' @param model The surrogate model which was used for the prediction
#'
#' @return numeric vector, LCB results
#' @export
#'
#' @examples
#' spot(,funSphere,c(-2,-3),c(1,2), control = list(infillCriterion = infillLCB, modelControl = list(target = c("y","s"))))
infillLCB <- function(predictionList, model){
    mean <- predictionList$y
    sd <- predictionList$s
    return(mean-sd)
}