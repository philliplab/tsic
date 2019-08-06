context("test-prob_test_result_given_days_to_visit")

test_that("results are correctly converted into probabilities", {
  assay_dynamics <- get_assay_dynamics('step_unit_testing')
  x <- prob_test_result_given_days_to_visit(assay_dynamics = assay_dynamics,
                                            result = '+',
                                            days_to_visit = 20)
  expect_equal(x, 1)
  x <- prob_test_result_given_days_to_visit(assay_dynamics = assay_dynamics,
                                            result = '-',
                                            days_to_visit = 20)
  expect_equal(x, 0)

  x <- prob_test_result_given_days_to_visit(assay_dynamics = assay_dynamics,
                                            result = '+',
                                            days_to_visit = 2)
  expect_equal(x, 0)
  x <- prob_test_result_given_days_to_visit(assay_dynamics = assay_dynamics,
                                            result = '-',
                                            days_to_visit = 2)
  expect_equal(x, 1)


})
