
###################################################################################################
#' Evaluate Model 
#'
#' This function produces an objective function with y=f(x) from a provided model fit.
#' Important note: this function expects \code{predict(object,newdata)} to return
#' a list. The \code{object$target} parameter is a string that determins which list item
#' is returned by the created function. If not set (NULL), \code{object$target} is set to \code{"y"}.
#'
#' @param object fit created by a modeling function, e.g., \code{\link{buildRandomForest}}
#' @param infillCriterion optional parameter, a function that accepts prediction results and a model object. The function should use these
#' to alter the prediction result in a user desired way. For example turning the prediction results of a kriging model (mean and sd) into the expected 
#' improvement criterion
#' @return a function in the style of \code{y=f(x)}, which uses the fitted object to predict \code{y} for sample \code{x}. 
#' @export
#' @keywords internal
###################################################################################################
evaluateModel <- function(object, infillCriterion = NULL){
    if(is.null(object$target)) 
        object$target <- "y"
    
    evalModelFun <- function(x){  
        res <- predict(object=object,newdata=x)[object$target]
        if(length(res) == 1){
            res <- res[[1]]
        }
        return(res)
    }
    
    if(is.null(infillCriterion)){
        return(
            #evalModelFun
            function(x){
                res <- evalModelFun(x)
                if(is.list(res)){
                    return(res[[1]])
                }
                return(res)
            }
        )
    }
    
    return(
        function(x){
            return(infillCriterion(evalModelFun(x), object))
        }
    )
}