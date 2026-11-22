test_that("run_filterRSD returns correct class", {
  td <- make_test_data()
  # Remove NAs first so RSD calculation is clean
  td$x[is.na(td$x)] <- 100
  td$x[td$x == 0] <- 100

  res <- run_filterRSD(td$x, td$metadata, max_rsd = 0.5,
                       verbose = FALSE)
  expect_s3_class(res, "run_filterRSD")
  expect_equal(nrow(res$data), nrow(td$x))
  expect_true(ncol(res$data) <= ncol(td$x))
})

test_that("run_filterRSD warns when no QC samples found", {
  td <- make_test_data()
  td$metadata$Group <- "Sample"  # no QC
  expect_warning(
    run_filterRSD(td$x, td$metadata, verbose = FALSE),
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

test_that("run_filterRSD returns named RSDs and removes undefined RSDs", {
  x <- cbind(stable = rep(100, 6) + 1:6, noisy = c(1, 100, 1, 100, 1, 100),
             empty_qc = c(NA, NA, NA, 5, 6, 7))
  meta <- data.frame(Group = c("QC", "QC", "QC", "A", "A", "A"))
  res <- run_filterRSD(x, meta, max_rsd = 0.3, verbose = FALSE)
  expect_named(res$rsd_values, colnames(x))
  expect_equal(res$features_kept, "stable")
  expect_true(is.na(res$rsd_values[["empty_qc"]]))
})

test_that("run_filterRSD warns on centred (non-positive mean) data", {
  set.seed(11)
  x <- scale(matrix(rnorm(60), 12, 5))
  meta <- data.frame(Group = rep(c("QC", "A"), 6))
  expect_warning(run_filterRSD(x, meta, verbose = FALSE), "<= 0")
})

test_that("run_filterRSD qc_types can select one QC label", {
  x <- cbind(F1 = c(10, 10.1, 9.9, 1, 50, 100, 7, 8))
  meta <- data.frame(Group = c("EQC", "EQC", "EQC", "SQC", "SQC", "SQC", "A", "A"))
  expect_equal(run_filterRSD(x, meta, qc_types = "EQC", verbose = FALSE)$n_features_after, 1L)
  expect_equal(run_filterRSD(x, meta, qc_types = "SQC", verbose = FALSE)$n_features_after, 0L)
})
