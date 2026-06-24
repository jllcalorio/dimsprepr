test_that("run_filtervariance returns correct class and removes low-var features", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100  # clean for variance calc
  # Make one feature constant

  td$x[, 3] <- 42

  res <- run_filtervariance(td$x, percentile = 5, verbose = FALSE)
  expect_s3_class(res, "run_filtervariance")
  expect_false("F3" %in% colnames(res$data))
})

test_that("run_filtervariance preserves class", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  res_mat <- run_filtervariance(td$x, verbose = FALSE)
  expect_true(is.matrix(res_mat$data))

  res_df <- run_filtervariance(as.data.frame(td$x), verbose = FALSE)
  expect_true(is.data.frame(res_df$data))
})

test_that("run_filtervariance rejects bad inputs", {
  expect_error(run_filtervariance("x"), "matrix or data frame")
  td <- make_test_data()
  expect_error(run_filtervariance(td$x, percentile = 200), "between 0 and 100")
})
