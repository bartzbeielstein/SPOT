library(testthat)
library(SPOT)

options(spot.run.full.test = F)

test_check("SPOT")

options(spot.run.full.test = NULL)