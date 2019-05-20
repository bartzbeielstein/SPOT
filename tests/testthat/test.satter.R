
context("misc")

test_that("check Satterthwaite function", {
  res <- satter(MScoeff= c(1/4, -1/4)
                , MSi = c(394.9, 73.3)
                , dfi = c(4,3)
                , alpha = 0.1)
  #real <- c(80.4,   3,  30.8648, 685.5266) 
  real <- 3
	expect_equal(res[2], real)
  }
  )
