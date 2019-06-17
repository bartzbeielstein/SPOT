#' funOptimLecture
#'
#' A testfunction used in the optimizaton lecture of the AIT Masters course at TH Koeln
#'
#' @param vec input vector or matrix of candidate solution
#'
#' @return vector of objective function values
#' @export
funOptimLecture <- function(vec){
    a <- 1
    b <- 5.1/(4*pi^2)
    c <- 5/pi
    r <- 6
    s <- 10
    t <- 1/(8*pi)
    
    evalTime = 0.25
    if(length(dim(vec))<=1){
        x1 <- vec[1] * 1.2 + 47.5
        x2 <- vec[2] * 0.6 + 123
        term1 <- a * (x2 - b*x1^2 + c*x1 - r)^2
        term2 <- s*(1-t)*cos(x1)
        y <- term1 + term2 + s - 0.397887  + (77.7777)
        Sys.sleep(evalTime)
        return(matrix(y, 1, 1))
    }else{
        return(res = matrix(apply(vec,1,"funOptimLecture"), 1, 1))
    }
}