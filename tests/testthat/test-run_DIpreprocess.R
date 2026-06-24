test_that("run_DIpreprocess single combination works end-to-end", {
  td <- make_test_data()

  res <- run_DIpreprocess(
    x                = td$x,
    metadata         = td$metadata,
    normalize_method = "median",
    transform_method = "log10",
    correct_drift    = FALSE,
    rsd_threshold    = 0.5,
    variance_percentile = 5,
    verbose          = FALSE
  )

  expect_s3_class(res, "run_DIpreprocess")
  expect_null(res$error)
  expect_true(nrow(res$data_nonpls) > 0)
  expect_true(ncol(res$data_nonpls) > 0)
  expect_true(length(res$features_final) > 0)
})

test_that("run_DIpreprocess multi-combination mode works", {
  td <- make_test_data()

  res <- run_DIpreprocess(
    x                = td$x,
    metadata         = td$metadata,
    normalize_method = c("sum", "median"),
    transform_method = "log10",
    correct_drift    = FALSE,
    rsd_threshold    = 0.5,
    variance_percentile = 5,
    verbose          = FALSE
  )

  expect_s3_class(res, "run_DIpreprocess_multi")
  expect_equal(length(res), 2L)
  expect_true("sum_log10" %in% names(res))
  expect_true("median_log10" %in% names(res))
})

test_that("run_DIpreprocess rejects bad method names", {
  td <- make_test_data()
  expect_error(
    run_DIpreprocess(td$x, td$metadata, normalize_method = "bogus",
                     correct_drift = FALSE, verbose = FALSE),
    "Unknown normalize_method"
  )
})

test_that("print and summary methods work", {
  td <- make_test_data()
  res <- run_DIpreprocess(
    x = td$x, metadata = td$metadata,
    normalize_method = "median", transform_method = "log10",
    correct_drift = FALSE, rsd_threshold = 0.5,
    variance_percentile = 5, verbose = FALSE
  )
  expect_output(print(res), "run_DIpreprocess")
  expect_output(summary(res), "step-by-step")
})
