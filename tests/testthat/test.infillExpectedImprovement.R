context("infillExpectedImprovement")

test_that("check that infillExpectedImprovement has same result as target=ei", {
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
    
    ## InfillExpectedImprovement should produce exactly the same results as target=ei. Target=y should not!
    set.seed(1)
    resA <- spot(,funSphere,c(-2,-3),c(1,2), control = list(infillCriterion = infillExpectedImprovement, modelControl = list(target = c("y","s"))))
    set.seed(1)
    resB <- spot(,funSphere,c(-2,-3),c(1,2), control = list(modelControl = list(target = c("ei"))))
    set.seed(1)
    resC <- spot(,funSphere,c(-2,-3),c(1,2), control = list(modelControl = list(target = c("y"))))
    expect_true(all(resA$x == resB$x))
    expect_true(!all(resA$x == resC$x))
    
    if(getOption("spot.run.full.test")){
        set.seed(1)
        resA <- spot(,funRast,c(-2,-3),c(1,2), control = list(infillCriterion = infillExpectedImprovement, modelControl = list(target = c("y","s"))))
        set.seed(1)
        resB <- spot(,funRast,c(-2,-3),c(1,2), control = list(modelControl = list(target = c("ei"))))
         set.seed(1)
        resC <- spot(,funRast,c(-2,-3),c(1,2), control = list(modelControl = list(target = c("y"))))
        expect_true(all(resA$x == resB$x))
        expect_true(!all(resA$x == resC$x))
    }
})

test_that("infillExpectedImprovement works with cvModels and sLinear", {
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
    
    ## Specifying target y should fail because expected improvement needs an uncertainty
    expect_error(spot(,funSphere,c(-2,-3),c(1,2), control = list(infillCriterion = infillExpectedImprovement, model = buildCVModel, modelControl = 
                                                                     list(target = c("y"),
                                                                          modellingFunction = buildKriging))))
    
    ## Target y,sLinear should be the same as target NULL but should be different from target s
    if(getOption("spot.run.full.test", FALSE)){
        models <- c(buildKriging, buildRandomForest, buildLM)
        functions <- c(funSphere, funRast)
    }else{
        models <- c(buildKriging)
        functions <- c(funRast)
    }
    
    for(m in models){
        for(f in functions){
            set.seed(1)
            resA <- spot(,f,c(-2,-3),c(1,2), control = list(infillCriterion = infillExpectedImprovement, model = buildCVModel, modelControl = 
                                                                        list(target = c("y","s"),
                                                                             modellingFunction = m,
                                                                             uncertaintyEstimator = "s")))
            set.seed(1)
            resB <- spot(,f,c(-2,-3),c(1,2), control = list(infillCriterion = infillExpectedImprovement, model = buildCVModel, modelControl = 
                                                                        list(modellingFunction = m)))
            set.seed(1)
            resC <- spot(,f,c(-2,-3),c(1,2), control = list(infillCriterion = infillExpectedImprovement, model = buildCVModel, modelControl = 
                                                                        list(target = c("y","s"),
                                                                             modellingFunction = m,
                                                                             uncertaintyEstimator = "sLinear")))
            expect_true(all(resB$x == resC$x))
            expect_true(!all(resA$x == resC$x))
        }
    }
})

