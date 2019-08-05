context("test-step_assay_dynamics")

test_that("Test that the diagnostic_delany works: diagnostic delay = 20", {
  x5 <- step_assay_dynamics(5, 20)
  expect_equal(x5, 0)
  x50 <- step_assay_dynamics(50, 20)
  expect_equal(x50, 1)

  x20 <- step_assay_dynamics(20, 20)
  expect_equal(x20, 0)

  x19.999 <- step_assay_dynamics(19.999, 20)
  expect_equal(x19.999, 0)

  x20.001 <- step_assay_dynamics(20.001, 20)
  expect_equal(x20.001, 1)
})
