test_that("run_filterRSD returns correct class", {
  td <- make_test_data()
  # Remove NAs first so RSD calculation is clean
  td$x[is.na(td$x)] <- 100
  td$x[td$x == 0] <- 100

  res <- run_filterRSD(td$x, td$metadata, max_rsd = 0.5,
                       qc_type = "QC", verbose = FALSE)
  expect_s3_class(res, "run_filterRSD")
  expect_equal(nrow(res$data), nrow(td$x))
  expect_true(ncol(res$data) <= ncol(td$x))
})

test_that("run_filterRSD warns when no QC samples found", {
  td <- make_test_data()
  td$metadata$Group <- "Sample"  # no QC
  expect_warning(
    run_filterRSD(td$x, td$metadata, qc_type = "QC", verbose = FALSE),
    "No QC samples"
  )
})

test_that("run_filterRSD rejects bad inputs", {
  td <- make_test_data()
  expect_error(run_filterRSD("nope", td$metadata), "matrix or data frame")
  expect_error(run_filterRSD(td$x, td$metadata, max_rsd = 5), "between 0 and 1")
})

test_that("print.run_filterRSD works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  res <- run_filterRSD(td$x, td$metadata, verbose = FALSE)
  expect_output(print(res), "RSD Filtering")
})
