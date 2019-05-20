###################################################################################################
#' Sphere Test Function
#'
#' @param x matrix of points to evaluate with the sphere function. Rows for points and columns for dimension.
#'
#' @return 1-column matrix with resulting function values
#' @examples
#' funSphere(matrix(runif(18),,3))
#' @export
###################################################################################################
funSphere <- function(x)matrix(apply(x,1,function(x)sum(x^2)),,1)
