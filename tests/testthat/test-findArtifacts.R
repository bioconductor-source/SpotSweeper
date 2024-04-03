# run examples from findArtifacts() function documentation
example(findArtifacts, echo = FALSE)

test_that("example objects have correct class", {
  expect_s4_class(spe, "SpatialExperiment")
})

test_that("example object contains artifacts colDta", {
  expect_equal(dim(spe), c(X, XXXXX))
})

test_that("examples give correct number of artifact spots", {
  expect_equal(sum(isTRUE(spe$artifact)), X)
})
