library(SPOT)

set.seed(7)

args = commandArgs(trailingOnly=TRUE)
params <- as.numeric(strsplit(args[1],",", fixed = T)[[1]])
res <- sum(params^2)
cat(res)
