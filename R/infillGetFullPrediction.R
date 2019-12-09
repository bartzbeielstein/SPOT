#' infillGetFullPrediction
#'
#' Wrapper to get the full model Prediction instead of just one target for plotting purposes.
#'
#' @param predictionList The results of a predict.model call
#' @param model The surrogate model which was used for the prediction
#'
#' @return list, full model prediction
#' @keywords internal
infillGetFullPrediction <- function(predictionList, model){
    return(predictionList)
} 