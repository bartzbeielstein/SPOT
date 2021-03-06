% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimGenoud.R
\name{optimGenoud}
\alias{optimGenoud}
\title{Minimization by GENetic Optimization Using Derivatives}
\usage{
optimGenoud(x = NULL, fun, lower, upper, control = list(), ...)
}
\arguments{
\item{x}{optional start point, not used}

\item{fun}{objective function, which receives a matrix x and returns observations y}

\item{lower}{boundary of the search space}

\item{upper}{boundary of the search space}

\item{control}{list of control parameters
\describe{
  \item{\code{funEvals}}{Budget, number of function evaluations allowed. Default is 100.}
  \item{\code{populationSize}}{Population size, number of individuals in the population. Default is 10*dimension.}
}}

\item{...}{passed to \code{fun}}
}
\value{
list, with elements
\describe{
  \item{\code{x}}{NULL, currently not used}
  \item{\code{y}}{NULL, currently not used}
  \item{\code{xbest}}{best solution}
  \item{\code{ybest}}{best observation}
  \item{\code{count}}{number of evaluations of \code{fun}}
}
}
\description{
For minimization, this function uses the \code{"genoud"} method from the
code{rgenoud} package. It is basically a wrapper, to enable genoud for usage
in SPOT.
}
\examples{
res <- optimGenoud(,fun = funSphere,lower = c(-10,-20),upper=c(20,8))
res$ybest
}
