testthat::test_that("discovery cohort milestones match the locked analysis", {
  root <- project_root()
  if (!dir.exists(file.path(root, "data", "public", "nhanes"))) {
    download_nhanes(root)
  }
  discovery <- prepare_nhanes_period(root, "2021-2023")
  testthat::expect_equal(nrow(discovery), 11933)
  testthat::expect_equal(sum(discovery$female & discovery$age_eligible, na.rm = TRUE), 1503)
  testthat::expect_equal(sum(discovery$strict_eligible, na.rm = TRUE), 824)
  testthat::expect_equal(sum(primary_complete_indicator(discovery) & positive_weight(discovery$hormone_weight), na.rm = TRUE), 643)
  common <- prepare_common_sample(discovery)
  testthat::expect_equal(nrow(common), 316)
})

testthat::test_that("temporal-comparison strict eligibility matches the locked analysis", {
  root <- project_root()
  replication <- prepare_nhanes_period(root, "2017-2020")
  testthat::expect_equal(sum(replication$strict_eligible, na.rm = TRUE), 1180)
})
