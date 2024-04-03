# run examples from modifiedZ() function documentation
data <- c(1, 2, 3, 4, 5, 100)
z_data <- modifiedZ(data)

# == Tests ==
test_that("example data gives correct Z", {
  expect_type(z_data, "double")
})

test_that("example data gives correct Z", {
  expect_equal(z_data, c(-0.758240, -0.454944, -0.151648,
                         0.151648, 0.454944, 29.268065))
})
