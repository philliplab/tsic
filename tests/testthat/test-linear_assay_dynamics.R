context("test-linear_assay_dynamics")

test_that("Test that the spread parameter works: diagnostic delay = 20, spread = 0.1", {
  x5 <- linear_assay_dynamics(5, 20, 0.1)
  expect_equal(x5, 0)
  x50 <- linear_assay_dynamics(50, 20, 0.1)
  expect_equal(x50, 1)
  x20 <- linear_assay_dynamics(20, 20, 0.1)
  expect_equal(x20, 0.5)

  x19 <- linear_assay_dynamics(19, 20, 0.1)
  expect_equal(x19, 0)
  x21 <- linear_assay_dynamics(21, 20, 0.1)
  expect_equal(x21, 1)

  x19.5 <- linear_assay_dynamics(19.5, 20, 0.1)
  expect_equal(x19.5, 0.25)
  x20.5 <- linear_assay_dynamics(20.5, 20, 0.1)
  expect_equal(x20.5, 0.75)
})

test_that("Test that the abs_spread parameter works: diagnostic delay = 10, abs_spread = 5", {
  x5 <- linear_assay_dynamics(5, 10, abs_spread = 5)
  expect_equal(x5, 0)
  x50 <- linear_assay_dynamics(50, 10, abs_spread = 5)
  expect_equal(x50, 1)
  x10 <- linear_assay_dynamics(10, 10, abs_spread = 5)
  expect_equal(x10, 0.5)

  x7.5 <- linear_assay_dynamics(7.5, 10, abs_spread = 5)
  expect_equal(x7.5, 0)
  x12.5 <- linear_assay_dynamics(12.5, 10, abs_spread = 5)
  expect_equal(x12.5, 1)

  x8.75 <- linear_assay_dynamics(8.75, 10, abs_spread = 5)
  expect_equal(x8.75, 0.25)
  x11.25 <- linear_assay_dynamics(11.25, 10, abs_spread = 5)
  expect_equal(x11.25, 0.75)
})

test_that("the linear assay dynamics function is vectorized", {
  x_vec <- linear_assay_dynamics(c(0, 1, 2, 3, 4, 5), 10, abs_spread = 5)
  expect_equal(unique(x_vec), 0)
  
  x_vec <- linear_assay_dynamics(c(20, 21, 22, 23, 24, 25), 10, abs_spread = 5)
  expect_equal(unique(x_vec), 1)

  x_vec <- linear_assay_dynamics(c(0, 1, 2, 3, 4, 5, 20, 21, 22, 23, 24, 25), 10, abs_spread = 5)
  expect_equal(x_vec, c(rep(0, 6), rep(1, 6)))
})

