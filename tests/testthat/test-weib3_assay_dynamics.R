context("test-weib3_assay_dynamics")

test_that("Weib3 assay dynamics works. Aptima from Delany params", {
  # due to rounding of the numbers reported in the paper, there will be small deviations.
  expect_lte(abs(0.5 - weib3_assay_dynamics(11.5, location = 4.8, shape = 1.35, scale = 9)), 0.011)
  expect_lte(abs(0.025 - weib3_assay_dynamics(5.3, location = 4.8, shape = 1.35, scale = 9)), 0.011)
  expect_lte(abs(0.972 - weib3_assay_dynamics(29.1, location = 4.8, shape = 1.35, scale = 9)), 0.011)
  expect_lte(abs(0.99 - weib3_assay_dynamics(33, location = 4.8, shape = 1.35, scale = 9)), 0.011)

  expect_true(all(weib3_assay_dynamics(1:50, location = 4.8, shape = 1.35, scale = 9) >= rep(0, 50)))
  expect_true(all(weib3_assay_dynamics(1:50, location = 4.8, shape = 1.35, scale = 9) <= rep(1, 50)))
})
