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

test_that('step assay dynamics function is vectorized', {
  x_vec <- step_assay_dynamics(1:10, 50)
  expect_equal(unique(x_vec), 0)

  x_vec <- step_assay_dynamics(10:20, 5)
  expect_equal(unique(x_vec), 1)

  x_vec <- step_assay_dynamics(c(1:10, 21:30), 15)
  expect_equal(x_vec, c(rep(0, 10), rep(1, 10)))
})
