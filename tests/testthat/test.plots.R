context("Plots")

manualTestPlots <- F

test_that("check plots are generated without errors", {
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
    
    if(manualTestPlots){
        x <- matrix(seq(-7,7,1.3), ncol = 1)
        y <- matrix(funRast(x))
        
        ## Build Models
        krig <- buildKriging(x,y)
        plotModel(krig, type = "singleDim")
        cvKrig <- buildCVModel(x,y,control = list(modellingFunction = buildKriging))
        plotModel(cvKrig, type = "singleDim")
        cvRF <- buildCVModel(x,y,control = list(modellingFunction = buildRandomForest))
        
        spotResEI <- spot(x, funRast, -5,5, control = list(funEvals = 25, 
                                                           infillCriterion = infillExpectedImprovement,
                                                           model = buildCVModel,
                                                           modelControl = 
                                                               list(modellingFunction = buildKriging),
                                                           optimizer = optimDE,
                                                           designControl = list(size = 0)))
        spotResBP <- spot(x, funRast, -5,5, control = list(funEvals = 25,
                                                           model = buildCVModel,
                                                           modelControl = 
                                                               list(modellingFunction = buildKriging),
                                                           optimizer = optimDE,
                                                           designControl = list(size = 0)))
        
        krig <- buildKriging(spotResBP$x[1:12,,drop=F],spotResBP$y[1:12,,drop=F])
        plotModel(krig, type = "singleDim")
        
        ind <- 14
        cvKrig <- buildCVModel(spotResBP$x[1:ind,,drop=F],spotResBP$y[1:ind,,drop=F],
                               control = list(modellingFunction = buildKriging))
        plotModel(cvKrig, type = "singleDim")
        cvKrig <- buildCVModel(spotResEI$x[1:ind,,drop=F],spotResBP$y[1:ind,,drop=F],
                               control = list(modellingFunction = buildKriging))
        plotModel(cvKrig, type = "singleDim")
        
        
        ## Kriging standard Plot
        krig$target <- c("y")
        plotModel(krig, type = "singleDim")
        
        ## Plot with uncertainty
        krig$target <- c("y", "s")
        plotModel(krig, type = "singleDim")
        
        ## Kriging CV Model and Plot
        cvKrig$target <- c("y")
        plotModel(cvKrig, type = "singleDim")
        
        cvKrig$target <- c("y","s")
        cvKrig$uncertaintyEstimator <- "s"
        plotModel(cvKrig, type = "singleDim")
        
        cvKrig$target <- c("y","s")
        cvKrig$uncertaintyEstimator <- "sLinear"
        plotModel(cvKrig, type = "singleDim")
        
        ## RF CV Model and Plot
        cvRF$target <- c("y")
        plotModel(cvRF, type = "singleDim")
        
        cvRF$target <- c("y","s")
        cvRF$uncertaintyEstimator <- "s"
        plotModel(cvRF, type = "singleDim")
        
        cvRF$target <- c("y","s")
        cvRF$uncertaintyEstimator <- "sLinear"
        plotModel(cvRF, type = "singleDim")
    }
    
    ## 
    ## Single Dimension SE Plot:
    ## 
    x <- matrix(c(1,2,3,5,7,8,9), ncol = 1)
    y <- matrix(sin(x))
    krig <- buildKriging(x,y)
    krig$target <- c("y")
    
    expect_error(plotModel(krig, type = "singleDim"), NA)
    
    cvMod <- buildCVModel(x,y,control = list(modellingFunction = buildKriging))
    expect_error(plotModel(cvMod, type = "singleDim"), NA)
    cvMod$target <- c("y")
    expect_error(plotModel(cvMod, type = "singleDim"), NA)
    cvMod$target <- c("y","s")
    cvMod$uncertaintyEstimator <- "s"
    expect_error(plotModel(cvMod, type = "singleDim"), NA)
    
    cvMod <- suppressWarnings(buildCVModel(x,y,control = list(modellingFunction = buildRandomForest)))
    cvMod$target <- c("y","s")
    cvMod$uncertaintyEstimator <- "sLinear"
    expect_error(plotModel(cvMod, type = "singleDim"), NA)
    
    if(getOption("spot.run.full.test")){
        cvMod <- buildCVModel(x,y,control = list(modellingFunction = buildLM))
        cvMod$target <- c("y","sLinear")
        expect_error(plotModel(cvMod, type = "singleDim"), NA)
        
        cvMod <- buildCVModel(x,y,control = list(modellingFunction = buildRSM))
        cvMod$target <- c("y","sLinear")
        expect_error(plotModel(cvMod, type = "singleDim"), NA)
    }
    
    ## This should fail saying only 2 dimensional models can be plotted like that
    expect_error(plotModel(krig),"The specified plot type is only available for 2 or more dimensions")
    
    
    ## 
    ## 2 Dimensional Plots
    ## 
    
    x <- matrix(runif(20), ncol = 2)
    y <- funRast(x)
    
    krig <- buildKriging(x,y)
    krig$target <- c("y","s","ei")
    expect_error(plotModel(krig), NA)
    if(getOption("spot.run.full.test")){
        expect_error(plotModel(krig, type = "singleDim"), NA)
        expect_error(plotModel(krig, type = "persp3d"), NA)
        
        cvMod <- buildCVModel(x,y,control = list(modellingFunction = buildKriging))
        expect_error(plotModel(cvMod), NA)
    
        expect_error(plotModel(cvMod, type = "singleDim"), NA)
        cvMod$target <- c("y","sLinear")
        expect_error(plotModel(cvMod, type = "singleDim"), NA)
        expect_error(plotModel(cvMod, type = "persp3d"), NA)
        
        cvMod <- suppressWarnings(buildCVModel(x,y,control = list(modellingFunction = buildRandomForest)))
        cvMod$target <- c("y","sLinear")
        expect_error(plotModel(cvMod, type = "singleDim"), NA)
        expect_error(plotModel(cvMod), NA)
    }
    
    ## 
    ## n Dimensional Plots
    ## 
    if(getOption("spot.run.full.test")){
        x <- matrix(runif(4*30), ncol = 4)
        y <- funRast(x)
        
        krig <- buildKriging(x,y)
        krig$target <- c("y","s","ei")
        expect_error(plotModel(krig), NA)
        expect_error(plotModel(krig, type = "singleDim"), NA)
        expect_error(plotModel(krig, type = "persp3d"), NA)
        
        cvMod <- buildCVModel(x,y,control = list(modellingFunction = buildKriging))
        expect_error(plotModel(cvMod), NA)
        expect_error(plotModel(cvMod, type = "singleDim"), NA)
        cvMod$target <- c("y","sLinear")
        
        expect_error(plotModel(cvMod, type = "singleDim"), NA)
        expect_error(plotModel(cvMod, type = "persp3d"), NA)
        
        
        cvMod <- suppressWarnings(buildCVModel(x,y,control = list(modellingFunction = buildRandomForest)))
        cvMod$target <- c("y","sLinear")
        expect_error(plotModel(cvMod, type = "singleDim"), NA)
        expect_error(plotModel(cvMod), NA)
    }
})