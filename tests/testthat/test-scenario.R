context('test-scenario')

if (FALSE){
  devtools::load_all()
}

test_that('scenario_recognizer recognizes scenarios', {
  list_of_assays <- c("iscav2_weib3_delaney_and_tosiano", "taqman_weib3_delaney_and_manufacturer",
                      "architect_weib3_delaney", "geenius_indet_weib3_delaney", "geenius_fr_weib3_delaney")
  sc_times <- sim_sc_times(list_of_assays, fix_draw = 0.5)
  
  f_ihist_s1 <- combine_sc_and_visit_times(sc_times,
                                           ((0:2)*28-46),
                                           true_infection_date = as.numeric(as.Date('2019-03-01'))+46.5)
  m_ihist_s1 <- select_most_informative_results(f_ihist_s1)$kept_ihist
  
  f_ihist_s2 <- combine_sc_and_visit_times(sc_times,
                                           ((0:2)*28-42),
                                           true_infection_date = as.numeric(as.Date('2019-03-01'))+42.5)
  m_ihist_s2 <- select_most_informative_results(f_ihist_s2)$kept_ihist
  
  f_ihist_s3 <- combine_sc_and_visit_times(sc_times,
                                           ((0:2)*28-36),
                                           true_infection_date = as.numeric(as.Date('2019-03-01'))+36.5)
  m_ihist_s3 <- select_most_informative_results(f_ihist_s3)$kept_ihist
  
  f_ihist_s4 <- combine_sc_and_visit_times(sc_times,
                                           ((0:2)*28-25),
                                           true_infection_date = as.numeric(as.Date('2019-03-01'))+25.5)
  m_ihist_s4 <- select_most_informative_results(f_ihist_s4)$kept_ihist
  
  f_ihist_s5 <- combine_sc_and_visit_times(sc_times,
                                           ((0:2)*28-20),
                                           true_infection_date = as.numeric(as.Date('2019-03-01'))+20.5)
  m_ihist_s5 <- select_most_informative_results(f_ihist_s5)$kept_ihist

  expect_equal(recognize_scenario(ihist = f_ihist_s1)$matching_scenario, 'Scenario 1')
  expect_equal(recognize_scenario(ihist = f_ihist_s2)$matching_scenario, 'Scenario 2')
  expect_equal(recognize_scenario(ihist = f_ihist_s3)$matching_scenario, 'Scenario 3')
  expect_equal(recognize_scenario(ihist = f_ihist_s4)$matching_scenario, 'Scenario 4')
  expect_equal(recognize_scenario(ihist = f_ihist_s5)$matching_scenario, 'Scenario 5')

  expect_equal(recognize_scenario(ihist = m_ihist_s1)$matching_scenario, 'Scenario 1')
  expect_equal(recognize_scenario(ihist = m_ihist_s2)$matching_scenario, 'Scenario 2')
  expect_equal(recognize_scenario(ihist = m_ihist_s3)$matching_scenario, 'Scenario 3')
  expect_equal(recognize_scenario(ihist = m_ihist_s4)$matching_scenario, 'Scenario 4')
  expect_equal(recognize_scenario(ihist = m_ihist_s5)$matching_scenario, 'Scenario 5')

  expect_equal(recognize_scenario(ihist = f_ihist_s4)$lnfp_gap, 28)
  expect_equal(recognize_scenario(ihist = m_ihist_s5)$lnfp_gap, 28)
})

