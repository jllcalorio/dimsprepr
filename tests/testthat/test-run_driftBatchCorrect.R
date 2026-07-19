# Helper: minimal data that pmp::QCRSC can actually fit
# 6 QC samples per batch satisfies smooth.spline (needs >= 4 unique x values)
make_qcrsc_data <- function(n_bio = 20, n_qc = 6, p = 5, seed = 4421) {
  set.seed(seed)
  n <- n_bio + n_qc
  x <- matrix(abs(rnorm(n * p, 500, 80)), nrow = n, ncol = p)
  colnames(x) <- paste0("F", seq_len(p))
  rownames(x) <- paste0("S", seq_len(n))
  meta <- data.frame(
    Sample            = rownames(x),
    Group             = c(rep("Sample", n_bio), rep("QC", n_qc)),
    Batch             = rep(1L, n),
    InjectionSequence = seq_len(n),
    stringsAsFactors  = FALSE
  )
  rownames(meta) <- meta$Sample
  list(x = x, metadata = meta)
}

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

test_that("QCRSC tolerates NAs in biological samples", {
  skip_if_not_installed("pmp")
  td <- make_qcrsc_data()
  td$x[1, 3] <- NA   # NA in a biological sample
  td$x[2, 4] <- NA
  res <- run_driftBatchCorrect(
    td$x, td$metadata,
    injection_col = "InjectionSequence",
    spline_smooth_param = 0.75,   # fixed spar avoids CV with small n
    verbose = FALSE
  )
  expect_s3_class(res, "run_driftBatchCorrect")
  expect_true(res$correction_applied)
})

test_that("QCRSC tolerates NAs in QC samples when enough QCs remain", {
  skip_if_not_installed("pmp")
  td <- make_qcrsc_data(n_qc = 8)
  td$x[21, 1] <- NA   # NA in one QC sample (row 21) for feature 1
  res <- run_driftBatchCorrect(
    td$x, td$metadata,
    injection_col = "InjectionSequence",
    spline_smooth_param = 0.75,
    verbose = FALSE
  )
  expect_s3_class(res, "run_driftBatchCorrect")
  expect_true(res$correction_applied)
})

test_that("QCRSC CV crash gives an actionable error message", {
  skip_if_not_installed("pmp")
  # 2 QC samples -> smooth.spline CV needs >= 4 unique x values -> crashes internally
  td <- make_qcrsc_data(n_bio = 10, n_qc = 2)
  expect_error(
    run_driftBatchCorrect(
      td$x, td$metadata,
      injection_col = "InjectionSequence",
      min_QC = 2,
      spline_smooth_param = 0,   # CV mode
      verbose = FALSE
    ),
    regexp = "spline_smooth_param"
  )
})
