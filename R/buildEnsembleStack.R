###################################################################################
#' Ensemble: Stacking
#'
#' Generates an ensemble of surrogate models with stacking (stacked generalization).
#'
#' @param x design matrix (sample locations), rows for each sample, columns for each variable.
#' @param y vector of observations at \code{x}
#' @param control (list), with the options for the model building procedure:\cr
#' \code{modelL1} Function for fitting the L1 model (default: \code{buildLM}) which combines the results of the L0 models. \cr
#' \code{modelL1Control} List of control parameters for the L1 model (default: \code{list()}).\cr
#' \code{modelL0} A list of functions for fitting the L0 models (default: \code{list(buildLM,buildRandomForest,buildKriging)}). \cr
#' \code{modelL0Control} List of control lists for each L0 model (default: \code{list(list(),list(),list())}).\cr
#'
#' @return returns an object of class \code{ensembleStack}.
#'
#' @seealso \code{\link{predict.ensembleStack}}
#'
#' @references Bartz-Beielstein, Thomas. Stacked Generalization of Surrogate Models-A Practical Approach. Technical Report 5/2016, TH Koeln, Koeln, 2016.
#' @references David H Wolpert. Stacked generalization. Neural Networks, 5(2):241-259, January 1992.
#'
#' @note Loosely based on the code by Emanuele Olivetti https://github.com/emanuele/kaggle_pbr/blob/master/blend.py
#'
#' @examples
#' ## Create a test function: branin
#' braninFunction <- function (x) {	
#' 	(x[2]  - 5.1/(4 * pi^2) * (x[1] ^2) + 5/pi * x[1]  - 6)^2 + 
#'	10 * (1 - 1/(8 * pi)) * cos(x[1] ) + 10
#' }
#' ## Create design points
#' x <- cbind(runif(20)*15-5,runif(20)*15)
#' ## Compute observations at design points
#' y <- as.matrix(apply(x,1,braninFunction))
#' ## Create model with default settings
#' fit <- buildEnsembleStack(x,y)
#' ## Predict new point
#' predict(fit,cbind(1,2))
#' ## True value at location
#' braninFunction(c(1,2))
#' 
#' @export
###################################################################################
buildEnsembleStack <- function(x, y, control=list()){ #nugget -1 means that the nugget will be optimized in lme
	con <- list(modelL1 = buildLM  #L1 model: Linear Model
          ,modelL1Control = list() #default: empty list, use default parameters for L1 model
				  ,modelL0 = list(buildLM
				                  ,buildKriging
				                  ) #L0 models: Linear model and Kriging
          ,modelL0Control = list(list()
                                 ,list()
                                  ) #default: empty lists, use default parameters for all L0 models
          ,k = 10 #number of folds (k-fold CV)
				)
	con[names(control)] <- control
	control<-con
		
	
	## number of data samples
	n <- nrow(x)
	## number of folds
	k <- control$k
  ## number of ensemble models on level 0
  p <- length(control$modelL0)
    ## generate k-fold data split
  folds <- sample(cut(seq(1,n),breaks=k,labels=FALSE))
  ## for each model
  fit0 <- list()
  ythat <- matrix(NA,n,p)
  for(j in 1:p){
    fit0[[j]] <- list()
    ## do k-fold cross-validation
    for(i in 1:k){
      xt <- x[folds!=i,,drop=FALSE]
      yt <- y[folds!=i]
      fit0[[j]][[i]] <- control$modelL0[[j]](xt
                                             ,as.matrix(yt)
                                             ,control$modelL0Control[[j]])
      ythat[folds==i,j] <- predict(fit0[[j]][[i]]
                                   ,x[folds==i,,drop=FALSE])$y
      }
  }
  ## fit level 1 model
  fit1 <- control$modelL1(ythat, y, control$modelL1Control)
  fit <- list(fit0=fit0, fit1=fit1, ythat=ythat,y=y,x=x,k=k,p=p) 
	fit$x <- x
	fit$y <- y
	class(fit)<- "ensembleStack"
	fit
}

###################################################################################
#' Predict Stacked Ensemble
#' 
#' Predict with ensemble model produced by \code{\link{buildEnsembleStack}}.
#'
#' @param object Kriging model (settings and parameters) of class \code{kriging}.
#' @param newdata design matrix to be predicted
#' @param ... not used
#'
#' @return list with predicted value \code{y}
#'
#' @seealso \code{\link{buildKriging}}
#' @export
#' @keywords internal
###################################################################################
predict.ensembleStack <- function(object,newdata,...){  
  ## number of folds
  k <- object$k
  ## number of level 0 models
  p <- object$p  
  #for each level 0 model
  ythat <- NULL
  for(j in 1:p){ 
    # TODO this may be different than as described in Bart16j: 
    # averaging on ythat, not on pred(l1,ythat)
    ## for each of k models (k-fold cross-validation)
    ## predict, then average.
    ytemp <- NULL
    for(i in 1:k){
       ytemp <- cbind(ytemp,predict(object$fit0[[j]][[i]],newdata)$y)
    }
    ythat <- cbind(ythat,rowMeans(ytemp))
  }
  predict(object$fit1,ythat)
}