context('data_processing')

test_that('select_most_informative_results works', {
  # input checking
  ihist <- data.frame(
    ptid = c('p0', 'p1'),
    sample_date = c(as.numeric(as.Date('2016-03-01')) + 0.5, as.numeric(as.Date('2016-09-01')) + 0.5),
    test = c('step_unit_testing', 'step_unit_testing'),
    result = c('-', '+'),
    stringsAsFactors = FALSE
  )
  expect_error(select_most_informative_results(ihist, NULL))

  ihist <- data.frame(
    ptid = c('p0', 'p0'),
    sample_date = c(as.numeric(as.Date('2016-03-01')) + 0.5, as.numeric(as.Date('2016-09-01')) + 0.5),
    test = c('step_unit_testing', 'step_unit_testing'),
    result = c('-', '+'),
    stringsAsFactors = FALSE
  )
  expect_error(select_most_informative_results(ihist, NULL))

  basic_checks_most_infor <- function(ihist, inf_ihist){
    expect_true(class(inf_ihist) == 'list')
    expect_true(all(sort(names(inf_ihist)) == c('kept_ihist', 'rm_ihist')))
    expect_true(all(sort(unique(ihist$sample_date)) == sort(unique(inf_ihist$kept_ihist$sample_date))))
    expect_equal(nrow(ihist), nrow(inf_ihist$kept_ihist) + nrow(inf_ihist$rm_ihist))
    counts <- with(inf_ihist$kept_ihist, tapply(test, list(sample_date, result), length))
    expect_lte(max(counts, na.rm = TRUE), 1)
    expect_gte(min(counts, na.rm = TRUE), 0)
  }

  # simple controlled cases
  ihist <- data.frame(
    ptid = c('p0', 'p0'),
    sample_date = c(as.numeric(as.Date('2016-03-01')) + 0.5, as.numeric(as.Date('2016-03-01')) + 0.5),
    test = c('architect_weib3_delaney', 'taqman_weib3_delaney_and_manufacturer'),
    result = c('+', '+'),
    stringsAsFactors = FALSE
  )
  inf_ihist <- select_most_informative_results(ihist)
  basic_checks_most_infor(ihist, inf_ihist)
  expect_equal(nrow(inf_ihist$kept_ihist), 1)
  expect_equal(inf_ihist$kept_ihist$test, 'architect_weib3_delaney')
  expect_equal(nrow(inf_ihist$rm_ihist), 1)
  expect_equal(inf_ihist$rm_ihist$test, 'taqman_weib3_delaney_and_manufacturer')

  # simple controlled cases
  ihist <- data.frame(
    ptid = c('p0', 'p0'),
    sample_date = c(as.numeric(as.Date('2016-03-01')) + 0.5, as.numeric(as.Date('2016-05-01')) + 0.5),
    test = c('architect_weib3_delaney', 'taqman_weib3_delaney_and_manufacturer'),
    result = c('+', '+'),
    stringsAsFactors = FALSE
  )
  inf_ihist <- select_most_informative_results(ihist)
  basic_checks_most_infor(ihist, inf_ihist)
  expect_equal(nrow(inf_ihist$kept_ihist), 2)
  expect_equal(nrow(inf_ihist$rm_ihist), 0)

  # big realistic case

  ihist <- data.frame(
    ptid = c("p314", "p314", "p314", "p314", "p314", "p314", "p314", "p314", "p314"), 
    sample_date = c(16860.5, 16910.5, 16921.5, 16921.5, 16910.5, 16921.5, 16921.5, 16910.5, 16860.5), 
    test = c("architect_weib3_delaney", "architect_weib3_delaney", "architect_weib3_delaney", 
             "geenius_fr_weib3_delaney", "geenius_indet_weib3_delaney", "geenius_indet_weib3_delaney", 
             "taqman_weib3_delaney_and_manufacturer", "taqman_weib3_delaney_and_manufacturer", 
             "taqman_weib3_delaney_and_manufacturer"), 
    result = c("-", "+", "+", "+", "-", "+", "+", "+", "-"),
    stringsAsFactors = FALSE)
  expected_result <- list(
    kept_ihist = structure(list(
        ptid = c("p314", "p314", "p314", "p314"), 
        sample_date = c(16860.5, 16910.5, 16910.5, 16921.5), 
        test = c("taqman_weib3_delaney_and_manufacturer", "architect_weib3_delaney", "geenius_indet_weib3_delaney", 
                 "geenius_fr_weib3_delaney"), 
        result = c("-", "+", "-", "+")), 
      row.names = c(9L, 2L, 5L, 4L), 
      class = "data.frame"), 
    rm_ihist = structure(list(
        ptid = c("p314", "p314", "p314", "p314", "p314"), 
        sample_date = c(16860.5, 16910.5, 16921.5, 16921.5, 16921.5), 
        test = c("architect_weib3_delaney", "taqman_weib3_delaney_and_manufacturer", "architect_weib3_delaney", 
                 "geenius_indet_weib3_delaney", "taqman_weib3_delaney_and_manufacturer"), 
        result = c("-", "+", "+", "+", "+")), 
      row.names = c(1L, 8L, 3L, 6L, 7L), 
      class = "data.frame"))


  inf_ihist <- select_most_informative_results(ihist)
  expect_equal(inf_ihist, expected_result)
  basic_checks_most_infor(ihist, inf_ihist)
  expect_true(class(inf_ihist) == 'list')
  expect_true(all(sort(names(inf_ihist)) == c('kept_ihist', 'rm_ihist')))
  expect_true(all(sort(unique(ihist$sample_date)) == sort(unique(inf_ihist$kept_ihist$sample_date))))
  expect_equal(nrow(ihist), nrow(inf_ihist$kept_ihist) + nrow(inf_ihist$rm_ihist))
  counts <- with(inf_ihist$kept_ihist, tapply(test, list(sample_date, result), length))
  expect_lte(max(counts, na.rm = TRUE), 1)
  expect_gte(min(counts, na.rm = TRUE), 0)
  expect_equal(inf_ihist, expected_result)
})

test_that('visit_labeller is working', {
  ihist <- data.frame(
    ptid = c("p314", "p314", "p314", "p314", "p314", "p314", "p314", "p314", "p314"), 
    sample_date = c(16860.5, 16910.5, 16921.5, 16921.5, 16910.5, 16921.5, 16921.5, 16910.5, 16860.5), 
    test = c("architect_weib3_delaney", "architect_weib3_delaney", "architect_weib3_delaney", 
             "geenius_fr_weib3_delaney", "geenius_indet_weib3_delaney", "geenius_indet_weib3_delaney", 
             "taqman_weib3_delaney_and_manufacturer", "taqman_weib3_delaney_and_manufacturer", 
             "taqman_weib3_delaney_and_manufacturer"), 
    result = c("-", "+", "+", "+", "-", "+", "+", "+", "-"),
    stringsAsFactors = FALSE)
  ihist <- ihist[order(ihist$sample_date, ihist$test),]
  l_ihist <- visit_labeller(ihist)
  expect_equal(c(-50, 0, 11), sort(unique(l_ihist$rel_sample_date)))
  expect_equal(c(-1, 0, 1), sort(unique(l_ihist$visit_number)))
})
