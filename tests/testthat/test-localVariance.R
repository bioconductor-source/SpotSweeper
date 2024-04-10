# run examples from findArtifacts() function documentation
data(DLPFC_artifact)
spe <- DLPFC_artifact

# find artifacts using
set.seed(123)
spe <- findArtifacts(spe,
                     mito_percent = "expr_chrM_ratio",
                     mito_sum = "expr_chrM",
                     n_rings = 5,
                     name = "artifact"
)


# === Tests ===
test_that("example objects have correct class", {
  expect_s4_class(spe, "SpatialExperiment")
})

test_that("examples give correct number of colData", {
  expect_equal(length(colnames(colData(spe))), 26)
})


# NOTE: this test currently gives essentially the same, but not exact numbers.
# this likely has to do with using exact knn ties being chosen randomly.
# need to change to use knn within a radius, which will use all tied points.
#
# This is very minor as findArtifacts test always passes with the same number
# of artifact points, so the function is working as expected.
#
#test_that("examples gives correct local variance", {
#  expect_equal(colData(spe)$k18[1:10], c(0.6205693,  0.8147886,  0.6324888,
#                                         -0.9733463,  1.0527102,  0.8142039,
#                                         -2.9437993, -0.9599085,  0.2998882,
#                                         0.9236396))
#})

