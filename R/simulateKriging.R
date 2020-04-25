
###################################################################################
#' Kriging Simulation: Spectral Method
#' 
#' (Conditional) Simulation via spectral method.
#'
#' @param object fit of the Kriging model (settings and parameters), of class \code{kriging}.
#' @param conditionalSimulation  logical, if set to TRUE (default), the simulation is conditioned with the training data of the Kriging model.
#' Else, the simulation is non-conditional.
#' @param Ncos number of cosine functions used to construct the simulation
#' @param ... further arguments, not used
#'
#' @return Returned value depends on the setting of \code{object$simulationReturnAll}
#'
#' @references N. A. Cressie. Statistics for Spatial Data. JOHN WILEY & SONS INC, 1993.
#' @references C. Lantuejoul. Geostatistical Simulation - Models and Algorithms. Springer-Verlag Berlin Heidelberg, 2002.
#'
#' @seealso \code{\link{buildKriging}}, \code{\link{simulationDecompose}}
#' @keywords internal
###################################################################################
simulationSpectral <-function(object,conditionalSimulation=FALSE,Ncos=100){
  theta <- object$dmodeltheta
  lsq <- 1/theta
  mu <- object$mu
  sigmas  <- object$ssq
  dimension <- length(theta)
  omegavar <- 2/lsq
  omega <- matrix(rnorm(Ncos*dimension,0,sqrt(omegavar)),dimension)
  phi <- runif(Ncos,-pi,pi)
  multiplier <- sqrt(sigmas) * sqrt(2/Ncos)
  omega 
  phi
  object
  multiplier
  mu
  fun <- function(xx){
    xx <- normalizeMatrix2(data.matrix(xx),0,1,object$normalizexmin,object$normalizexmax)
    y <-  multiplier * colSums( cos(t(xx%*%omega)+ phi))
    y <- y+mu
    return(y)
  }  
  force(fun)
  if(conditionalSimulation){
    objectSim <- object
    xsim <- objectSim$x
    #
    #ysim <- fun(xsim)
    ysim <- object$y-fun(xsim)
    #
    objectSim$y <- ysim
    #
    objectSim$yonemu <- ysim #  - objectSim$mu   
    objectSim$mu <- 0
    #
    force(objectSim)
    funConditional <- function(xx){
      #
      #y <- predict(object,xx)$y + fun(xx) - predict(objectSim,xx)$y
      y <- fun(xx) + predict(objectSim,xx)$y
      #
      return(y)
    }
    return(funConditional)
  }else{
    return(fun)
  }
}

###################################################################################
#' Kriging Simulation: Decomposition
#' 
#' (Conditional) Simulation via decomposition approach.
#'
#' @param object fit of the Kriging model (settings and parameters), of class \code{kriging}.
#' @param xsim list of samples in input space, to be simulated
#' @param nsim number of simulations
#' @param conditionalSimulation  logical, if set to TRUE (default), the simulation is conditioned with the training data of the Kriging model.
#' Else, the simulation is non-conditional.
#' @param returnAll if set to TRUE, a list with the simulated values (y) and the corresponding covariance matrix (covar)
#' of the simulated samples is returned. 
#' @param ... further arguments, not used
#'
#' @return Returned value depends on the setting of \code{object$simulationReturnAll}
#'
#' @references N. A. Cressie. Statistics for Spatial Data. JOHN WILEY & SONS INC, 1993.
#' @references C. Lantuejoul. Geostatistical Simulation - Models and Algorithms. Springer-Verlag Berlin Heidelberg, 2002.
#'
#' @seealso \code{\link{buildKriging}}, \code{\link{simulationSpectral}}
#' @keywords internal
###################################################################################
simulationDecompose <- function(object,nsim=1,xsim,conditionalSimulation=TRUE,returnAll=FALSE,...){
  #
  len <- nrow(xsim) #number of simulated samples (points  
  noise <- matrix(rnorm(len*nsim),len, nsim)
  #
  covar <- getCorrelationMatrix(object,xsim)
  #covar <- res$psi
  #
  
  if(conditionalSimulation){
    object$returnCrossCor <- TRUE
    ret <- predict(object,xsim)
    y <- ret$y
    psi <- ret$psi
    
    covarDifference <- covar - psi %*% object$Psinv %*% t(psi)
    eigv <- eigen(object$ssq *covarDifference,symmetric=T) #eigen decomposition
    covarDecomposed <- eigv$vectors %*% diag(sqrt(abs(eigv$values))) %*% eigv$vectors
    ysim <- covarDecomposed %*% noise
    
    #and the following adds the simulation part to the predictor
    y <- matrix(y,len,nsim) + ysim	
  }else{
    eigv <- eigen(object$ssq * covar,symmetric=T) #eigen decomposition
    covarDecomposed <- eigv$vectors %*% diag(sqrt(abs(eigv$values))) %*% eigv$vectors
    y <- object$mu + covarDecomposed %*% noise
  }
  res <- list(y=y,psi=covar)
  if(returnAll)
    return(res)	
  else
    return(y)
}

###################################################################################
#' Kriging Simulation
#' 
#' (Conditional) Simulation at given locations, with a model fit resulting from \code{\link{buildKriging}}.
#' In contrast to prediction or estimation, the goal is to reproduce the covariance 
#' structure, rather than the data itself. Note, that the conditional simulation 
#' also reproduces the training data, but
#' has a two times larger error than the Kriging predictor.
#'
#' @param object fit of the Kriging model (settings and parameters), of class \code{kriging}.
#' @param xsim list of samples in input space, to be simulated at
#' @param method \code{"decompose"} (default) or \code{"spectral"}, specifying the method used for simulation. 
#' Note that \code{"decompose"} is can be preferable, since it is exact but may be computationally infeasible for high-dimensional xsim.
#' On the other hand, \code{"spectral"} yields a function that can be evaluated at arbitrary sample locations.
#' @param nsim number of simulations
#' @param seed random number generator seed. Defaults to NA, in which case no seed is set
#' @param conditionalSimulation  logical, if set to TRUE (default), the simulation is conditioned with the training data of the Kriging model.
#' Else, the simulation is non-conditional.
#' @param Ncos number of cosine functions (used with \code{method="spectral"} only)
#' @param returnAll if set to TRUE, a list with the simulated values (y) and the corresponding covariance matrix (covar)
#' of the simulated samples is returned. 
#' @param ... further arguments, not used
#'
#' @return Returned value depends on the setting of \code{object$simulationReturnAll}
#'
#' @references N. A. Cressie. Statistics for Spatial Data. JOHN WILEY & SONS INC, 1993.
#' @references C. Lantuejoul. Geostatistical Simulation - Models and Algorithms. Springer-Verlag Berlin Heidelberg, 2002.
#'
#' @seealso \code{\link{buildKriging}}, \code{\link{predict.kriging}}
#' @export
###################################################################################
simulate.kriging <- function(object,nsim=1,seed=NA,xsim,method="decompose",conditionalSimulation=TRUE,Ncos=10,returnAll=FALSE,...){
  if (!is.na(seed)){
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
      runif(1)
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }	
  if(method=="decompose"){
    simresult <- simulationDecompose(object=object,xsim=xsim,nsim=nsim,conditionalSimulation=conditionalSimulation,returnAll=returnAll)
    return(simresult)
  }else if(method=="spectral"){
    simresult <- NULL
    simfun <- NULL
    for(i in 1:nsim){
      res <- simulationSpectral(object=object,conditionalSimulation=conditionalSimulation,Ncos=Ncos)
      simresult <- cbind(simresult,res(xsim))
      simfun <- c(simfun,res)
    }
    if(returnAll){
      return(list(y=simresult,simfun=simfun))
    }else{
      return(simresult)
    }
  }else{
    stop("The specified method used in simulate.kriging does not exist. Use 'decompose' or 'spectal'")
  }
}

###################################################################################
#' Compute Correlation Matrix
#' 
#' Compute the correlation matrix of samples x, given the model object.
#'
#' @param object fit of the Kriging model (settings and parameters), of class \code{kriging}.
#' @param x list of samples / data
#'
#' @return the correlation matrix
#'
#' @seealso \code{\link{simulate.kriging}}
#' @seealso \code{\link{predict.kriging}}
#' @keywords internal
###################################################################################
getCorrelationMatrix <- function(object,x){
  x <- normalizeMatrix2(data.matrix(x),0,1,object$normalizexmin,object$normalizexmax)
  k <- ncol(x)
  n <- nrow(x)
  
  A<-matrix(0,k,n*n)
  
  for(i in 1:k){
    if(object$types[i]!="factor"){
      A[i,]<-as.numeric(as.matrix(dist(x[,i]))) #euclidean distance
    }else	{
      tmp <- outer(x[,i],x[,i],'!=') #hamming distance
      class(tmp) <- "numeric"
      A[i,]<-tmp
    }
  }
  
  theta <- object$dmodeltheta
  
  if(object$optimizeP)
    p <- object$P
  else
    p <- rep(2,k)	
  
  A <- abs(A)^p
  
  
  Psi <- exp(-matrix(colSums(theta*A),n,n)) 
  
  if(object$useLambda){
    lambda <- object$dmodellambda
    Psi <- Psi+diag(lambda,n)
  }
  
  Psi  
}

###################################################################################
#' Simulation-based Function Generator
#' 
#' Generate functions via simulation of Kriging models, e.g.,
#' for assessment of optimization algorithms with
#' non-conditional or conditional simulation, based on real-world data.
#'
#' @param object an object generated by \code{\link{buildKriging}}
#' @param nsim the number of simulations, or test functions, to be created
#' @param seed a random number generator seed. Defaults to NA; which means no seed is set. For sake of reproducibility, set this to some integer value.
#' @param method \code{"decompose"} (default) or \code{"spectral"}, specifying the method used for simulation. 
#' Note that \code{"decompose"} is can be preferable, since it is exact but may be computationally infeasible for high-dimensional xsim.
#' On the other hand, \code{"spectral"} yields a function that can be evaluated at arbitrary sample locations.
#' @param xsim list of samples in input space, for simulation (only used for decomposition-based simulation, not for spectral method)
#' @param Ncos number of cosine functions (used with \code{method="spectral"} only)
#' @param conditionalSimulation whether (TRUE) or not (FALSE) to use conditional simulation
#'
#' @return a list of functions, where each function is the interpolation of one simulation realization. The length of the list depends on the nsim parameter.
#' 
#' @seealso \code{\link{buildKriging}}, \code{\link{simulate.kriging}}
#' 
#' @references N. A. Cressie. Statistics for Spatial Data. JOHN WILEY & SONS INC, 1993.
#' @references C. Lantuejoul. Geostatistical Simulation - Models and Algorithms. Springer-Verlag Berlin Heidelberg, 2002.
#'
#' @export
###################################################################################
simulateFunction <- function(object,nsim=1, seed=NA, method="spectral", xsim=NA, Ncos=10,
                             conditionalSimulation=TRUE){
  if(any(is.na(xsim)) & method=="decompose"){
    stop("simulatedBenchmarkFunction can not create a benchmark function via simulation, if method=='decompose' and xsim is not provided (is NA).")
  }
  force(object)
  #todo seed handling?
  if (!is.na(seed)){
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
      runif(1)
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }	
  
  if(method=="decompose"){
    simfit <- simulate(object=object,nsim=nsim,seed=NA,xsim=xsim,conditionalSimulation=conditionalSimulation,returnAll=TRUE)
    ynew <- simfit$y
    
    object$Psi <- simfit$psi
    object$Psinv <- MASS::ginv(object$Psi) 
    object$x <- xsim
    object$scaledx <- normalizeMatrix2(data.matrix(xsim),0,1,object$normalizexmin,object$normalizexmax) 
    PsinvSaved  <- object$Psinv
    
    n <- length(xsim)
    
    fun <- list()
    for(i in 1:nsim){
      object$y <- ynew[,i,drop=FALSE]
      object$yonemu <- ynew[,i,drop=FALSE] - object$mu 
      testFun <- NULL
      assign("testFun", eval(substitute(
        function(x){
          predict(object,x)$y
        }, 
        list(object=object)
      )
      ),
      envir=environment())	      
      fun[[i]] <- testFun
    }	
  }else if(method=="spectral"){
    fun <- list()
    for(i in 1:nsim){
      fun[[i]] <- simulationSpectral(object,conditionalSimulation=conditionalSimulation,Ncos=Ncos)
    }
  }else{
    stop("The specified method in simulatedBenchmarkFunction does not exist. Use 'decompose' or 'spectal'")
  }
  return(fun)
}

