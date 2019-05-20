context("InitialInputChecking - ControlLists")

test_that("Wrong Control Settings are blocked",{
    xCorrect <- matrix(1:12,ncol = 3)
    lCorrect <- c(-5,-5,-5)
    uCorrect <- c(15,15,15)
    
    errorControlLists <- list(
        list(funEvals = "50",
             designControl = list(size = 10)),
        list(funEvals = "50",
             designControl = list(size = 10, replicates = "B")),
        list(funEvals = 50,
             designControl = list(size = 10, replicates = "B")),
        list(funEvals = 50,
             designControl = list(size = 10, replicates = 10)),
        list(funEvals = 8),
        
        #This should still crash as there is a 4row input matrix in x
        list(funEvals = 13),
        
        list(funEvals = 3,designControl = list(size=1))
    )
    
    warningControlLists <- list(
        list(funEvals = 50 + 5*4,
             designControl = list(size = 10, replicates = 5)),
        list(funEvals = 50 + 4,
             designControl = list(size = 50)),
        list(funEvals = 14)
    )
    
    goodControlLists <- list(
        list(funEvals <- 50,
             designControl <- list(size = 10, replicates = 3)),
        list(funEvals <- 8+4, ## +4 because of the 4 row input in x
             designControl <- list(size = 7, replicates = 1))
    )
    
    for(c in goodControlLists){
        expect_error({
            spot(x = xCorrect, fun = funSphere, lower = lCorrect, upper = uCorrect, control = c)
        },regexp = NA)
    }
    
    for(c in warningControlLists){
        expect_warning({
            spot(x = xCorrect, fun = funSphere, lower = lCorrect, upper = uCorrect, control = c)
        },regexp = "SPOT Configuration Warning:")
    }
    
    for(c in errorControlLists){
        expect_error({
            spot(x = xCorrect, fun = funSphere, lower = lCorrect, upper = uCorrect, control = c)
        },regexp = "SPOT Configuration Error:")
    }
})

test_that("SPOT throws an error if control$types is configured wrong",{
    expect_error({
        spot(NULL, fun = funSphere, lower = c(0,0), upper = c(5,5), control = list(types = c("numeric")))
    },regexp = "SPOT Configuration Error:")
    expect_error({
        spot(NULL, fun = funSphere, lower = c(0,0), upper = c(5,5), control = list(types = c("numeric","numeric","numeric")))
    },regexp = "SPOT Configuration Error:")
    
    expect_warning({
        spot(NULL,fun = funSphere, lower = c(0,0), upper = c(5,5), control = list(types = c("numeric","num")))
    },regexp = "SPOT Configuration Warning:")
    expect_warning({
        spot(NULL,fun = funSphere, lower = c(0,0), upper = c(5,5), control = list(types = c("float","numeric")))
    },regexp = "SPOT Configuration Warning:")
    expect_warning({
        spot(NULL,fun = funSphere, lower = c(0,0), upper = c(5,5), control = list(types = c("numeric",88)))
    },regexp = "SPOT Configuration Warning:")
})

test_that("SPOT uses the correct amount of function evaluations",{
    xCorrect <- matrix(1:12,ncol = 3)
    lCorrect <- c(-5,-5,-5)
    uCorrect <- c(15,15,15)
    
    #Some sample controlLists
    controlLists <- list(NULL,list(),
                         list(funEvals = 23),
                         list(funEvals = 15,designControl = list(size=10)),
                         list(funEvals = 15),
                         list(funEvals = 9,designControl = list(size=4)),
                         list(funEvals = 6,designControl = list(size=1)))
    
    #Expected result lengths of those lists
    expectedLengths <- c(20,20,23,15,15,9,6)
    
    for(i in seq_along(controlLists)){
        expect_equal({
            length(spot(x = xCorrect, fun = funSphere, 
                        lower = lCorrect, upper = uCorrect, 
                        control = controlLists[[i]])$y)},expectedLengths[i])
    }
})