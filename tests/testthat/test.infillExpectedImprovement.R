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
    resA <- spot(,funSphere,c(-2,-3),c(1,2), control = list(modelControl = list(infillCriterion = infillExpectedImprovement, target = c("y","s"))))
    set.seed(1)
    resB <- spot(,funSphere,c(-2,-3),c(1,2), control = list(modelControl = list(target = c("ei"))))
    set.seed(1)
    resC <- spot(,funSphere,c(-2,-3),c(1,2), control = list(modelControl = list(target = c("y"))))
    expect_true(all(resA$x == resB$x))
    expect_true(!all(resA$x == resC$x))
    
    set.seed(1)
    resA <- spot(,funRast,c(-2,-3),c(1,2), control = list(modelControl = list(infillCriterion = infillExpectedImprovement, target = c("y","s"))))
    set.seed(1)
    resB <- spot(,funRast,c(-2,-3),c(1,2), control = list(modelControl = list(target = c("ei"))))
    set.seed(1)
    resC <- spot(,funRast,c(-2,-3),c(1,2), control = list(modelControl = list(target = c("y"))))
    expect_true(all(resA$x == resB$x))
    expect_true(!all(resA$x == resC$x))
})