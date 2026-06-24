test_that("run_filtermissing returns correct class and dimensions", {
  td <- make_test_data()
  res <- run_filtermissing(td$x, td$metadata, threshold = 0.5, verbose = FALSE)

  expect_s3_class(res, "run_filtermissing")
  expect_equal(nrow(res$data), nrow(td$x))
  expect_true(ncol(res$data) <= ncol(td$x))
  expect_equal(res$n_features_before, ncol(td$x))
  expect_equal(res$n_features_after, ncol(res$data))
  expect_equal(res$n_features_removed, res$n_features_before - res$n_features_after)
})

test_that("run_filtermissing removes fully-missing columns", {
  td <- make_test_data()
  td$x[, 1] <- NA
  res <- run_filtermissing(td$x, td$metadata, threshold = 0.01, verbose = FALSE,
                           zero_as_missing = FALSE, filter_by_group = FALSE)
  expect_false("F1" %in% colnames(res$data))
})

test_that("run_filtermissing treats zeros as missing by default", {
  td <- make_test_data()
  td$x[, 2] <- 0  # all zeros in one feature
  res <- run_filtermissing(td$x, td$metadata, threshold = 0.01,
                           filter_by_group = FALSE, verbose = FALSE)
  expect_false("F2" %in% colnames(res$data))
})

test_that("run_filtermissing preserves matrix class", {
  td <- make_test_data()
  res <- run_filtermissing(td$x, td$metadata, verbose = FALSE)
  expect_true(is.matrix(res$data))

  res_df <- run_filtermissing(as.data.frame(td$x), td$metadata, verbose = FALSE)
  expect_true(is.data.frame(res_df$data))
})

test_that("run_filtermissing rejects bad inputs", {
  td <- make_test_data()
  expect_error(run_filtermissing("not a matrix", td$metadata), "matrix or data frame")
  expect_error(run_filtermissing(td$x, "not a df"), "data frame")
  expect_error(run_filtermissing(td$x, td$metadata[1:5, ]), "nrow")
  expect_error(run_filtermissing(td$x, td$metadata, threshold = 2), "between 0 and 1")
})

test_that("print.run_filtermissing works", {
  td <- make_test_data()
  res <- run_filtermissing(td$x, td$metadata, verbose = FALSE)
  expect_output(print(res), "Missing Value Filtering")
})
