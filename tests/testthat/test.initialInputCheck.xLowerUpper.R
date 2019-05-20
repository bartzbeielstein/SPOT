
context("InitialInputChecking - xLowerUpper")

test_that("Check conditions with wrong input dimensions", {
    xCorrect <- matrix(1:12,ncol = 3)
    lCorrect <- c(-5,-5,-5)
    uCorrect <- c(15,15,15)
    
    warning <- 1
    error <- 2
    
    #Define a list of wrong scenarios for the lower bounds
    lWrongs <- list(c(-5,-5),c(20,20,20),c(20,20),c(-5,-5,-5,-5))
    lWrongs.expected <- c(error,warning,error,error)
    
    #Define a list of wrong scenarios for the upper bounds
    uWrongs <- list(c(15,15),c(-6,-6,-6),c(-5,-5),c(15,15,15,15))
    uWrongs.expected <- c(error,warning,error,error)
    
    #Define a list of wrong scenarios for x
    xWrongs <- list(matrix(1:12,ncol = 2),
                    matrix(c(1:11,NA),ncol = 3),
                    matrix(1:12,ncol = 4))
    xWrongs.expected <- c(error,error,error)
    
    #Test all wrong lowerbounds
    for(l in seq_along(lWrongs)){
        expected <- lWrongs.expected[l]
        if(expected == warning){
            expect_warning({
                spot(x = xCorrect, fun = funSphere, lower = lWrongs[[l]], upper = uCorrect)
            }, regexp = "SPOT Configuration Warning:")
        }else{
            expect_error({
                spot(x = xCorrect, fun = funSphere, lower = lWrongs[[l]], upper = uCorrect)
            }, regexp = "SPOT Configuration Error:")
        }
    }
    
    #Test all wrong upper bounds
    for(u in seq_along(uWrongs)){
        expected <- uWrongs.expected[u]
        if(expected == warning){
            expect_warning({
                spot(x = xCorrect, fun = funSphere, lower = lCorrect, upper = uWrongs[[u]])
            }, regexp = "SPOT Configuration Warning:")
        }else{
            expect_error({
                spot(x = xCorrect, fun = funSphere, lower = lCorrect, upper = uWrongs[[u]])
            }, regexp = "SPOT Configuration Error:")
        }
    }
    
    #Test all combinations of upper and lower bounds together with correct input for x
    for(l in seq_along(lWrongs)){
        for(u in seq_along(uWrongs)){
            expected <- max(lWrongs.expected[[l]], uWrongs.expected[[u]])
            if(expected == warning){
                expect_warning({
                    spot(x = xCorrect, fun = funSphere, lower = lWrongs[[l]], upper = uWrongs[[u]])
                })
            }else{
                expect_error({
                    spot(x = xCorrect, fun = funSphere, lower = lWrongs[[l]], upper = uWrongs[[u]])
                }, regexp = "SPOT Configuration Error:")
            }
        }
    }
    
    #Test wrong input in x
    for(xW in seq_along(xWrongs)){
        expected <- xWrongs.expected[[xW]]
        if(expected == warning){
            expect_warning({
                spot(x = xWrongs[[xW]], fun = funSphere, lower = lCorrect, upper = uCorrect)
            }, regexp = "SPOT Configuration Warning:")
        }else{
            expect_error({
                spot(x = xWrongs[[xW]], fun = funSphere, lower = lCorrect, upper = uCorrect)
            }, regexp = "SPOT Configuration Error:")
        }
    }
    
    #Test for error if lower and upper are same
    expect_error({
        spot(x = xCorrect, fun = funSphere, lower = c(5,5,5), upper = c(5,5,5))
    }, regexp = "SPOT Configuration Error:")
    expect_error({
        spot(x = xCorrect, fun = funSphere, lower = c(1,8,2.6), upper = c(11,8,17.4))
    }, regexp = "SPOT Configuration Error:")
    expect_error({
        spot(x = xCorrect, fun = funSphere, lower = c(12,1,2), upper = c(12,4,5))
    }, regexp = "SPOT Configuration Error:")
    
    #Test for Error when NAs are present in bounds
    expect_error({
        spot(x = xCorrect, fun = funSphere, lower = c(5,NA,5), upper = c(5,5,5))
    }, regexp = "SPOT Configuration Error:")
    expect_error({
        spot(x = xCorrect, fun = funSphere, lower = c(5,5,5), upper = c(NA,NA,NA))
    }, regexp = "SPOT Configuration Error:")
    
    #Test for Error in strange cases with wrong inputs
    expect_error({
        spot(x = xCorrect, fun = funSphere, lower = c(5,"someString",5), upper = c(5,5,5))
    }, regexp = "SPOT Configuration Error:")
    expect_error({
        spot(x = xCorrect, fun = funSphere, lower = c(5,5,5), upper = c("Cause","an","Error"))
    }, regexp = "SPOT Configuration Error:")
})

test_that("InputChecking Works with correct inputs",{
    xCorrect <- matrix(1:12,ncol = 3)
    lCorrect <- c(-5,-5,-5)
    uCorrect <- c(15,15,15)
    
    expect_equal(initialInputCheck(xCorrect, funSphere, lCorrect, uCorrect),TRUE)
})