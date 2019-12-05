# Package Description for Roxygen:
#' Sequential Parameter Optimization Toolbox
#'
#' SPOT uses a combination statistic models and optimization algorithms 
#' for the purpose of parameter optimization. Design of Experiment methods
#' are employed to generate an initial set of candidate solutions, which
#' are evaluated with a user-provided objective function.
#' The resulting data is used to fit a model, which in turn is subject
#' to an optimization algorithm, to find the most promising candidate solution(s).
#' These are again evaluated, after which the model is updated with the new results.
#' This sequential procedure of modeling, optimization and evaluation is iterated
#' until the evaluation budget is exhausted. 
#'
#' Note, that versions >= 2.0.1 of the package are a complete rewrite of the interfaces and conventions in SPOT.
#' The rewritten SPOT package aims to improve the following issues of the older package:\cr
#' - A more modular architecture is provided, that allows the user to easily customize parts of the SPO procedure.\cr
#' - Core functions for modeling and optimization use interfaces more similar to algorithms from other packages / core-R, 
#'   hence making them easier accessible for new users. Also, these can now be more easily used separately from the main SPO 
#'   approach, e.g., only for modeling.\cr
#' - Reducing the unnecessarily large number of choices and parameters.\cr
#' - Removal of extremely rarely used / un-used features, to reduce overall complexity of the package.\cr
#' - Improving documentation and accessibility in general.\cr
#' - Speed-up of frequently used procedures.\cr
#' We appreciate feedback about any bugs or other issues with the package. 
#' Feel free to send feedback by mail to the maintainer.
#'
#' \tabular{ll}{
#' Package: \tab SPOT\cr
#' Type: \tab Package\cr
#' Version: \tab 2.0.5\cr
#' Date: \tab 2019-12-05\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @name SPOT-package
#' @aliases SPOT
#' @docType package
#' @title Sequential Parameter Optimization Toolbox
#' @author Thomas Bartz-Beielstein \email{thomas.bartz-beielstein@@th-koeln.de}, Joerg Stork, Martin Zaefferer (Maintainer, \email{martin.zaefferer@@gmx.de})
#' with contributions from: C. Lasarczyk, M. Rebolledo, J. Ziegenhirt, W. Konen, O. Flasch, P. Koch, M. Friese, L. Gentile, F. Rehbach.
#' @keywords package
#' @seealso Main interface function is \code{\link{spot}}.
#' @import randomForest
#' @import stats
#' @import MASS
#' @import utils
#' @import graphics
#' @import grDevices
#' @import DEoptim
#' @import rgenoud
#' @import rsm
#' @import nloptr
#' 
#' @section Acknowledgments: 
#' This work has been partially supported by the Federal Ministry of Education
#' and Research (BMBF) under the grants CIMO (FKZ 17002X11) and
#' MCIOP (FKZ 17N0311).
#'
#' @section Maintainer:
#' Martin Zaefferer \email{martin.zaefferer@@gmx.de} 
#End of Package Description
NA 
