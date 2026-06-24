test_that("run_driftBatchCorrect with perform_correction=FALSE returns unchanged", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  res <- run_driftBatchCorrect(td$x, td$metadata,
                                perform_correction = FALSE, verbose = FALSE)
  expect_s3_class(res, "run_driftBatchCorrect")
  expect_false(res$correction_applied)
  expect_equal(res$data, td$x)
})

test_that("run_driftBatchCorrect rejects bad inputs", {
  td <- make_test_data()
  expect_error(run_driftBatchCorrect("x", td$metadata), "matrix or data frame")
  expect_error(run_driftBatchCorrect(td$x, "y"), "data frame")
})

test_that("print.run_driftBatchCorrect works", {
  td <- make_test_data()
  res <- run_driftBatchCorrect(td$x, td$metadata,
                                perform_correction = FALSE, verbose = FALSE)
  expect_output(print(res), "Drift/Batch Correction")
})
