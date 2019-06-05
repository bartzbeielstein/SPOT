context("Plots")

test_that("check plots are generated without errors", {
    ## 
    ## Single Dimension SE Plot:
    ## 
    x <- matrix(c(1,2,3,5,7,8,9), ncol = 1)
    y <- matrix(sin(x))
    krig <- buildKriging(x,y)
    krig$target <- c("y","s")
    
    expect_error(plotModel(krig, type = "singleDim"), NA)
    
    cvMod <- buildCVModel(x,y,buildKriging)
    expect_error(plotModel(cvMod, type = "singleDim"), NA)
    cvMod$target <- c("y","s")
    expect_error(plotModel(cvMod, type = "singleDim"), NA)
    cvMod$target <- c("y","sLinear")
    expect_error(plotModel(cvMod, type = "singleDim"), NA)
    
    ## This should fail saying only 2 dimensional models can be plotted like that
    expect_error(plotModel(krig),"The specified plot type is only available for 2 or more dimensions")
    
    
    ## 
    ## 2 Dimensional Plots
    ## 
    
    ## Fun Rastrigin from spotGUI package, ignore this
    funRast <- function (vec) {
        if (length(dim(vec)) <= 1) {
            sum = 0
            for (i in vec) {
                sum = sum + (i^2 - 10 * cos(2 * 3.1415926 * i))
            }
            return(matrix((10 * length(vec) + sum), , 1))
        }
        res = matrix(apply(vec, 1, "funRast"), , 1)
        return(res)
    }
    
    x <- matrix(runif(10), ncol = 2)
    y <- funRast(x)
    
    krig <- buildKriging(x,y)
    krig$target <- c("y","s","ei")
    expect_error(plotModel(krig), NA)
    expect_error(plotModel(krig, type = "singleDim"), NA)
    expect_error(plotModel(krig, type = "persp3d"), NA)
    
    cvMod <- buildCVModel(x,y,buildKriging)
    expect_error(plotModel(cvMod), NA)
    expect_error(plotModel(cvMod, type = "singleDim"), NA)
    cvMod$target <- c("y","sLinear")
    expect_error(plotModel(cvMod, type = "singleDim"), NA)
    expect_error(plotModel(cvMod, type = "persp3d"), NA)
    
    ## 
    ## n Dimensional Plots
    ## 
    x <- matrix(runif(4*30), ncol = 4)
    y <- funRast(x)
    
    krig <- buildKriging(x,y)
    krig$target <- c("y","s","ei")
    expect_error(plotModel(krig), NA)
    expect_error(plotModel(krig, type = "singleDim"), NA)
    expect_error(plotModel(krig, type = "persp3d"), NA)
    
    cvMod <- buildCVModel(x,y,buildKriging)
    expect_error(plotModel(cvMod), NA)
    expect_error(plotModel(cvMod, type = "singleDim"), NA)
    cvMod$target <- c("y","sLinear")
    expect_error(plotModel(cvMod, type = "singleDim"), NA)
    expect_error(plotModel(cvMod, type = "persp3d"), NA)
})