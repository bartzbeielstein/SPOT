library(SPOT)

args = commandArgs(trailingOnly=TRUE)

paramMatrix <- SPOT::wrapSystem_parseMatrixFromString(args[1])
set.seed(7)

runSANN <- function(params){
    optim(runif(10), function(x)sum(x^2), method = "SANN", control = list(params[1], params[2]))$value
}

res <- matrix(apply(paramMatrix,1,runSANN),,1)
SPOT::wrapSystem_parseMatrixToString(res)
