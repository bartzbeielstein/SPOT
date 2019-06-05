#' plotSingleDimFunction
#' 
#' Plot a single dimensional Function
#' 
#' @param evalFun 
#' @param lower 
#' @param upper 
#' @param target 
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_ribbon
plotSingleDimFunction <- function(evalFun, lower, upper, target){
    force(evalFun)
    
    x <- seq(lower, upper, length.out = 200)
    
    withSE = F
    if(is.null(target) | length(target) == 1){
        y <- evalFun(x)
    }else{
        y <- evalFun(x)$y
        if("s" %in% target){
            withSE = T
            s <- evalFun(x)$s 
        }
        if("sLinear" %in% target){
            withSE = T
            s <- evalFun(x)$sLinear
        }
    }
    
    if(!withSE){
        ggplot(data=data.frame(x,y), aes(x=x, y=y, group=1)) +
            geom_line()
    }else{
        ggplot(data=data.frame(x,s), aes(x=x)) +
            geom_line(aes(y = y)) + 
            geom_ribbon(aes(ymax=y + s, ymin=y - s), fill="pink", alpha=.5)
    }
}