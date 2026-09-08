testthat::test_that("weighted standardization has zero weighted mean and unit weighted spread", {
  x <- c(1, 2, 4, 8)
  w <- c(1, 2, 3, 4)
  z <- weighted_center_scale(x, w)
  testthat::expect_equal(weighted.mean(z, w), 0, tolerance = 1e-12)
  testthat::expect_equal(weighted.mean(z^2, w), 1, tolerance = 1e-12)
})

testthat::test_that("row_mean_min enforces the observed-component minimum", {
  data <- data.frame(a = c(1, 1), b = c(2, NA), c = c(3, NA))
  result <- row_mean_min(data, c("a", "b", "c"), minimum = 2)
  testthat::expect_equal(result[1], 2)
  testthat::expect_true(is.na(result[2]))
})

testthat::test_that("directional independence reference is one eighth", {
  testthat::expect_equal(0.25 * 0.50, 0.125)
})

