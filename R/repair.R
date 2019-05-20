###################################################################################################
#' Repair Non-numeric Values
#'
#' Round non-numeric columns of a matrix, specified by a vector of data given data types.
#'
#' @param x matrix to be rounded
#' @param types data types of the respective columns, numeric columns are specified by \code{"numeric"}.
#' @export
#' @keywords internal
#' @examples
#' x <- matrix(10*runif(12),4,3)
#' types <- c("numeric","factor","factor")
#' repairNonNumeric(x,types)
###################################################################################################
repairNonNumeric <- function(x,types){
	for (i in 1:length(types)){
    if(types[i] != "numeric") #use rounding if not numeric. note that categorical parameters are mapped to integers.
			x[,i] <- round(x[,i])
	}
	x
}