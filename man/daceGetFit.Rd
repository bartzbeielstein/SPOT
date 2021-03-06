% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildKrigingDACE.R
\name{daceGetFit}
\alias{daceGetFit}
\title{Get DACE fit}
\usage{
daceGetFit(theta, pars)
}
\arguments{
\item{theta}{model parameters to be evaluated}

\item{pars}{model option list, as created with \code{\link{dacePrepareFit}}}
}
\value{
list of model variables, with the following elements: \cr
		\code{regr} regression function used \cr
		\code{corr} correlation function used \cr
		\code{theta} model parameters\cr
		\code{beta} generalized least squares estimate\cr
		\code{gamma} correlation factors\cr
		\code{sigma2} maximum Likelihood estimate of the process variance\cr
		\code{S} scaled design points\cr
		\code{Ssc} scaling factors for design arguments\cr
		\code{Y} scaled observations\cr
		\code{Ysc} scaling factors for design ordinates\cr
		\code{C}  Cholesky factor of correlation matrix\cr
		\code{Ft} Decorrelated regression matrix\cr
		\code{G} From QR factorization: Ft = Q*t(G)\cr
}
\description{
Evaluate the fit of a certain set of model parameters (\code{theta}), and
get all relevant variables of the model.
}
\seealso{
\code{\link{buildKrigingDACE}} \code{\link{daceLikelihood}} \code{\link{daceEvalFit}}
}
\author{
The authors of the original DACE Matlab code are Hans Bruun Nielsen, Soren Nymand Lophaven and Jacob Sondergaard. \cr
Extension of the Matlab code by Tobias Wagner \email{wagner@isf.de}. \cr 
Porting and adaptation to R and further extensions by Martin Zaefferer \email{martin.zaefferer@fh-koeln.de}.
}
\keyword{internal}
