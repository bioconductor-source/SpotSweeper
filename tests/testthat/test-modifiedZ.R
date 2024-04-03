# run examples from modifiedZ() function documentation
example(modifiedZ)

test_that("example data gives correct Z", {
  expect_type(z_data, "numeric")
})
