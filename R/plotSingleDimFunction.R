#' plotSingleDimFunction
#' 
#' Plot a single dimensional Function
#' 
#' @param evalFun function to be plotted. The function should either be able to take two vectors or one matrix specifying sample locations. i.e. \code{z=f(X)} or \code{z=f(x2,x1)} where Z is a two column matrix containing the sample locations \code{x1} and \code{x2}.
#' @param lower boundary for x1 and x2 (defaults to \code{c(0,0)}).
#' @param upper boundary (defaults to \code{c(1,1)}).
#' @param target String, which type of uncertainty estimation should be plotted? default: NULL - no estimation plotted. 's' standard estimation. 'sLinear' linearly adapted estimation.
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_ribbon
#' 
#' @export
#' @keywords internal
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