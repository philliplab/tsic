context("test-linear_assay_dynamics")

test_that("Test that basics work: diagnostic delay = 20, spread = 0.1", {
  x5 <- linear_assay_dynamics(5, 20, 0.1)
  expect_equal(x5, 0)
  x50 <- linear_assay_dynamics(50, 20, 0.1)
  expect_equal(x50, 1)
  x20 <- linear_assay_dynamics(20, 20, 0.1)
  expect_equal(x20, 0.5)

  x <- linear_assay_dynamics(19, 20, 0.1)
  expect_equal(x, 0)
  x <- linear_assay_dynamics(21, 20, 0.1)
  expect_equal(x, 1)

  x <- linear_assay_dynamics(19.5, 20, 0.1)
  expect_equal(x, 0.25)
  x <- linear_assay_dynamics(20.5, 20, 0.1)
  expect_equal(x, 0.75)
})
