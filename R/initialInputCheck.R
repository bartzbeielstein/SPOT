
###################################################################################################
#' Initial Input Check of Spot Configuration 
#'
#' This function takes the same inputs as given to the \code{spot} call and checks for
#' possible problems in the users configuration
#'
#' @param x is an optional start point (or set of start points), specified as a matrix. One row for each point, and one column for each optimized parameter.
#' @param fun is the objective function. It should receive a matrix x and return a matrix y. In case the function uses external code and is noisy, an additional seed parameter may be used, see the \code{control$seedFun} argument below for details.
#' @param lower is a vector that defines the lower boundary of search space. This determines also the dimensionality of the problem.
#' @param upper is a vector that defines the upper boundary of search space.
#' @param control is a list with control settings for spot. See \code{\link{spotControl}}.
#' 
#' @keywords internal
###################################################################################################
initialInputCheck <- function(x=NULL,fun, 
                              lower,upper,control=list(), inSpotLoop = F){
    checkTypesOfInput(x,lower,upper,control)
    checkInputDimensionsionalityCorrect(x,lower,upper)
    checkLowerNotEqualsUpper(lower,upper)
    checkLowerSmallerThanUpper(lower,upper)
    checkForNAs(x,lower,upper)
    
    checkInputTypesInControl(control)
    checkFunEvalsDesignSize(x, lower, control, inSpotLoop)
    
    checkControlVarTypes(lower, control)
    
    return(TRUE)
}

###################################################################################################
#' Check for NAs in x lower and upper
#'
#' Creates an error message if there are any nas given in x lower or upper
#'
#' @param x is an optional start point (or set of start points), specified as a matrix. One row for each point, and one column for each optimized parameter.
#' @param lower is a vector that defines the lower boundary of search space. This determines also the dimensionality of the problem.
#' @param upper is a vector that defines the upper boundary of search space.
#' 
#' @keywords internal
###################################################################################################
checkForNAs <- function(x, lower, upper){
    if(any(is.na(x))){
        stop("SPOT Configuration Error: spotInput 'x' contains NAs")
    }
    if(any(is.na(lower))){
        stop("SPOT Configuration Error: spotInput 'lower' contains NAs")
    }
    if(any(is.na(upper))){
        stop("SPOT Configuration Error: spotInput 'upper' contains NAs")
    }
}

###################################################################################################
#' Check Input Types
#'
#' Creates an error message if any of the input types are not numeric even though they are configured to be numeric
#'
#' @param x is an optional start point (or set of start points), specified as a matrix. One row for each point, and one column for each optimized parameter.
#' @param lower is a vector that defines the lower boundary of search space. This determines also the dimensionality of the problem.
#' @param upper is a vector that defines the upper boundary of search space.
#' @param control is a list with control settings for spot. See \code{\link{spotControl}}.
#' 
#' @keywords internal
###################################################################################################
checkTypesOfInput <- function(x,lower,upper, control){
    #If control$types is null then no types are specified, everything shoud be numeric
    if(is.null(control$types)){
        if(!is.null(x)){
            if(!is.numeric(x)){
                stop("SPOT Configuration Error: spotInput 'x' contains non-numerics")
            }
        }
        if(!is.numeric(lower)){
            stop("SPOT Configuration Error: spotInput 'lower' contains non-numerics")
        }
        if(!is.numeric(upper)){
            stop("SPOT Configuration Error: spotInput 'upper' contains non-numerics")
        }
    }
}


###################################################################################################
#' Check Dimensions of spotInputs
#'
#' The dimensionality of x, lower and upper should match. If this is not the case, errors are generated
#'
#' @param x is an optional start point (or set of start points), specified as a matrix. One row for each point, and one column for each optimized parameter.
#' @param lower is a vector that defines the lower boundary of search space. This determines also the dimensionality of the problem.
#' @param upper is a vector that defines the upper boundary of search space.
#' 
#' @keywords internal
###################################################################################################
checkInputDimensionsionalityCorrect <- function(x,lower,upper){
    if(is.null(x)){
        if(!(length(lower) == length(upper))){
            stop("SPOT Configuration Error: lengths of spotInputs 'lower' and 'upper' are not equal")
        }
    }else{
        equals <- all(sapply(list(length(lower), length(upper)), function(arg){arg == ncol(x)}))
        if(!equals){
            stop("SPOT Configuration Error: lengths of spotInputs 'lower', 'upper' and 'x' are not equal")
        }
    }
}

###################################################################################################
#' Check That Lower and Upper are not Equal
#'
#' If any entries in lower and upper are equal, the parameter has no range and cant be optimized. 
#' In that case, an error is generated.
#'
#' @param lower is a vector that defines the lower boundary of search space. This determines also the dimensionality of the problem.
#' @param upper is a vector that defines the upper boundary of search space.
#' 
#' @keywords internal
###################################################################################################
checkLowerNotEqualsUpper <- function(lower,upper){
    if(any(lower == upper)){
        stop("SPOT Configuration Error: Entries in 'lower' and 'upper' should not be equal")
    }
}

###################################################################################################
#' Check That Lower is smaller than Upper
#'
#' Check if lower actually contains smaller values than upper. Otherwise a warning is generated.
#'
#' @param lower is a vector that defines the lower boundary of search space. This determines also the dimensionality of the problem.
#' @param upper is a vector that defines the upper boundary of search space.
#' 
#' @keywords internal
###################################################################################################
checkLowerSmallerThanUpper <- function(lower,upper){
    if(any(lower > upper)){
        warning("SPOT Configuration Warning: Entries in 'lower' are higher than entries in 'upper'")
    }
}

###################################################################################################
#' Check funEvals Setting against designSize
#'
#' Checks if the designSize will result in a larger value than funEvals. If so, return an error.
#'
#' @param x is an optional start point (or set of start points), specified as a matrix. One row for each point, and one column for each optimized parameter.
#' @param fun is the objective function. It should receive a matrix x and return a matrix y. In case the function uses external code and is noisy, an additional seed parameter may be used, see the \code{control$seedFun} argument below for details.
#' @param lower is a vector that defines the lower boundary of search space. This determines also the dimensionality of the problem.
#' @param upper is a vector that defines the upper boundary of search space.
#' @param inSpotLoop Boolean indicating whether the check is called from within spotLoop or not
#' 
#' @keywords internal
###################################################################################################
checkFunEvalsDesignSize <- function(x, lower, control, inSpotLoop){
    # Rules for checking funEvals:
    # funEvals should be equal to length(y) passed back by spot. Thus it is the total amount of candidate solutions
    # This also includes any optional points in x that are specified by the user. So if nrow(x) = 50 any funEvals below 50 cant work
    
    if(is.null(x)){
        lenX <- 0
    }else{
        lenX <- nrow(x)
    }
    
    con <- spotControl(length(lower))
    con[names(control)] <- control
    control <- con
    
    funEvals <- control$funEvals
    designSize <- control$designControl$size
    replicates <- control$designControl$replicates
    
    if(is.null(designSize)){designSize <- 10}
    if(is.null(replicates)){replicates <- 1}
    if(inSpotLoop){
        designSize <- 0
        replicates <- 1
    }
    
    if(funEvals == (designSize * replicates + lenX * replicates)){
        warning("SPOT Configuration Warning: The intial design will be as large as funEvals. 
                    SPOT is not run, you are only building a design!
                    Increase funEvals in the spot control list.")
    }
    if(funEvals < (designSize * replicates + lenX * replicates)){
        stop("SPOT Configuration Error: (designControl$size+nrow(x))*designControl$replicates exceeds control$funEvals. 
                    Your design is larger than SPOTs allowed Budget!
                    Increase funEvals in the spot control list.")
    }
}


###################################################################################################
#' Check input types in the spotControl list.
#'
#' specified variables are tested for their type in the control list. If a type mismatch is found, an error is thrown.
#'
#' @param control is a list with control settings for spot. See \code{\link{spotControl}}.
#' 
#' @keywords internal
###################################################################################################
checkInputTypesInControl <- function(control){
    #List of all variables that should be tested for type numeric
    typesNumeric <- list("funEvals", "designControl$size", "designControl$replicates", 
                         "designControl$retries", "OCBAbudget", "replicates")
    
    for(ele in typesNumeric){
        splits <- strsplit(ele,"$", fixed = T)[[1]]
        val <- control[[splits[1]]]
        if(length(splits)>1){
            for(i in 2:length(splits)){
                val <- val[[splits[i]]]
            }
        }
        if(!is.null(val)){
            if(!is.numeric(val)){
                stop(paste0("SPOT Configuration Error: element ", ele, " should be type 'numeric' but isn't"))
            }
        }
    }
}

checkControlVarTypes <- function(lower, control){
    #If control$types is null, test can be skipped as it will be set to default
    if(!is.null(control$types)){
        standardSpotTypes <- c("numeric", "integer", "factor")
        if(!(length(lower) == length(control$types))){
            stop(paste0("SPOT Configuration Error: length of control$types does not match input dimensionality."))
        }
    
        #If at least one non standard type is configured give a warning
        if(!all(control$types %in% standardSpotTypes)){
            warning(paste("SPOT Configuration Warning: Not all configured types in control$types are known to spot.
                      Unknown types will be mapped to integers!
                      Known types are:", paste0("'",standardSpotTypes,"'", collapse=" ")))
        }
    }
}