library(testthat)


test_that("Test about initialization", {
  expect_equal(length(rnorm(50)),50)
})
