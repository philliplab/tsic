context("test-result_interpretation")

test_that("assay_result_interpreter gets build correctly", {
  assay_dynamics <- get_assay_dynamics(assay = 'step_unit_testing')
  result <- '+'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result)
  expect_type(foo, 'closure')
  expect_equal(foo(20), 1)
  expect_equal(foo(2), 0)
  
  result <- '-'
  foo <- construct_assay_result_interpreter(assay_dynamics = assay_dynamics,
                                            result = result)
  expect_type(foo, 'closure')
  expect_equal(foo(20), 0)
  expect_equal(foo(2), 1)
})
