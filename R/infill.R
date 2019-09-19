
###################################################################################################
#' Expected Improvement
#'
#' Compute the negative logarithm of the Expected Improvement of a set of candidate solutions.
#' Based on mean and standard deviation of a candidate solution,
#' this estimates the expectation of improvement. Improvement
#' considers the amount by which the best known value (best observed value)
#' is exceeded by the candidates.
#'
#' @param mean vector of predicted means of the candidate solutions.
#' @param sd vector of estimated uncertainties / standard deviations of the candidate solutions.
#' @param min minimal observed value.
#' 
#' @return a vector with the negative logarithm of the expected improvement values, -log10(EI).
#'
#' @export
#' @examples
#' mean <- 1:10 #mean of the candidates
#' sd <- 10:1 #st. deviation of the candidates
#' min <- 5 #best known value
#' EI <- expectedImprovement(mean,sd,min)
#' EI
###################################################################################################
expectedImprovement <- function(mean,sd,min){ #NegLogExpImp 
	EITermOne=(min-mean)*pnorm((min-mean)/sd)
	EITermTwo=sd*(1/sqrt(2*pi))*exp(-(1/2)*((min-mean)^2/(sd^2)))
	-log10(EITermOne+EITermTwo+(.Machine$double.xmin))
}
