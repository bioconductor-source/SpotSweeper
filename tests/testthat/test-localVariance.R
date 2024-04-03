# run examples from localVariance() function documentation
example(localVariance, echo = FALSE)

test_that("example objects have correct class", {
  expect_s4_class(spe, "SpatialExperiment")
})

test_that("example objects have correct dimensions", {
  expect_equal(dim(spe), c(X, XXXXX))
})

test_that("examples give correct number of colData", {
  expect_equal(length(colnames(colData(spe))), XY)
})

test_that("examples gives correct variance for first 10 spots", {

})
