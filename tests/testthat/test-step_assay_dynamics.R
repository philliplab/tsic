context("test-step_assay_dynamics")

test_that("Test that the diagnostic_delany works: diagnostic delay = 20", {
  x5 <- step_assay_dynamics(5, 20)
  expect_equal(x5, 0)
  x50 <- step_assay_dynamics(50, 20)
  expect_equal(x50, 1)
})
