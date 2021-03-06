
context("SPOT general")

test_that("check that spot function works without problems", {
  res <- spot(,funSphere,c(-2,-3),c(1,2))
	expect_equal(c(nrow(res$x),length(res$y)), c(20,20))
  #res <- spot(,funSphere,c(-2,-3),c(1,2),control=list(model=buildRandomForest))	
  res <- spot(,funSphere,c(-2,-3),c(1,2),control=list(model=buildLM))	
	expect_equal(c(nrow(res$x),length(res$y)), c(20,20))
  res <- spot(,funSphere,c(-2,-3),c(1,2),control=list(funEvals=11,model=buildKrigingDACE))	
	expect_equal(c(nrow(res$x),length(res$y)), c(11,11))
  res <- spot(,funSphere,c(-2,-3),c(1,2),control=list(optimizer=optimLBFGSB))	
	expect_equal(c(nrow(res$x),length(res$y)), c(20,20))
  res <- spot(,funSphere,c(-2,-3),c(1,2),control=list(design=designUniformRandom))	
	expect_equal(c(nrow(res$x),length(res$y)), c(20,20))
	res <- spot(,funSphere,c(-2,-3),c(1,2),control=list(modelControl=list(target="ei")))
	expect_equal(c(nrow(res$x),length(res$y)), c(20,20))
})