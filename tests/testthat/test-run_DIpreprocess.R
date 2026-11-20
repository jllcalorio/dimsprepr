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
  expect_error(run_DIpreprocess(td$x, td$metadata, scale_pls = "zzz", verbose = FALSE),
               "Unknown scale_pls")
  expect_error(run_DIpreprocess(td$x, td$metadata, normalize_method = "pqn_reference",
                                verbose = FALSE), "normalize_ref_sample")
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

test_that("default pipeline runs end-to-end on the bundled data (with QCRSC)", {
  res <- run_DIpreprocess(dims_data, dims_metadata, spline_smooth = 0.5,
                          verbose = FALSE)
  expect_null(res$error)
  expect_equal(nrow(res$data_nonpls), nrow(dims_data))
  expect_false(anyNA(res$data_nonpls))
  # quality filters run on unscaled data, so both branches share one feature set
  expect_identical(colnames(res$data_nonpls), res$features_final)
  expect_identical(colnames(res$data_pls), res$features_final)
  expect_identical(colnames(res$data_transformed), res$features_final)
  # auto-scaled branch has unit variance, pareto branch does not
  expect_equal(unname(matrixStats::colSds(as.matrix(res$data_nonpls))),
               rep(1, length(res$features_final)), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(unname(matrixStats::colSds(as.matrix(res$data_pls))),
                                rep(1, length(res$features_final)))))
  expect_true(res$parameters$correct_drift)
})

test_that("quality filtering matches the step functions on unscaled data", {
  # Regression: filtering used to run after auto-scaling (all variances = 1,
  # QC means ~0), so the kept set was arbitrary.
  res <- run_DIpreprocess(dims_data, dims_metadata, correct_drift = FALSE,
                          variance_percentile = 25, verbose = FALSE)
  norm <- run_normalize(res$data_corrected, res$metadata, method = "median",
                        group_col = "Group_", qc_types = "QC", verbose = FALSE)
  rsd  <- run_filterRSD(norm$data, res$metadata, max_rsd = 0.3, verbose = FALSE)
  trans <- run_transform(norm$data[, rsd$features_kept], "log10", verbose = FALSE)
  vf   <- run_filtervariance(trans$data, percentile = 25, verbose = FALSE)
  expect_identical(res$features_final, vf$features_kept)
})

test_that("rsd_qc_type and track are honoured", {
  td <- make_test_data()
  td$metadata$Group[td$metadata$Group == "QC"] <- c("EQC", "EQC", "EQC", "SQC", "SQC", "SQC")
  res_e <- run_DIpreprocess(td$x, td$metadata, correct_drift = FALSE, rsd_threshold = 0.5,
                            rsd_qc_type = "EQC", variance_percentile = 5, verbose = FALSE)
  expect_null(res_e$error)
  expect_equal(res_e$parameters$rsd_qc_type, "EQC")

  res_t <- run_DIpreprocess(dims_data, dims_metadata, correct_drift = FALSE,
                            track = c("mz_0078", "mz_0001"), verbose = FALSE)
  expect_true("mz_0078" %in% res_t$track$Feature)
  expect_false("mz_0001" %in% res_t$track$Feature)
})

test_that("merge_replicates works with a custom subject column name", {
  meta <- dims_metadata
  meta$Subject <- meta$SubjectID
  meta$SubjectID <- NULL
  meta$Subject[meta$Subject %in% c("BIO02")] <- "BIO01"   # BIO01 measured twice
  res <- run_DIpreprocess(dims_data, meta, subject_id_col = "Subject",
                          correct_drift = FALSE, merge_replicates = TRUE, verbose = FALSE)
  expect_null(res$error)
  expect_equal(nrow(res$data_nonpls_merged), nrow(dims_data) - 1L)
  expect_equal(nrow(res$metadata_merged), nrow(dims_data) - 1L)
})

test_that("errors inside steps are captured in $error, not thrown", {
  res <- suppressMessages(
    run_DIpreprocess(dims_data, dims_metadata, correct_drift = FALSE,
                     rsd_threshold = 0, verbose = FALSE)
  )
  expect_match(res$error, "RSD filter")
  expect_output(print(res), "FAILED")
})

test_that("multi mode names, shares steps 1-5, and prints", {
  multi <- run_DIpreprocess(dims_data, dims_metadata, normalize_method = "median",
                            transform_method = c("log10", "log2"),
                            scale_nonpls = c("auto", "pareto"),
                            correct_drift = FALSE, verbose = FALSE)
  expect_named(multi, c("median_log10_auto", "median_log2_auto",
                        "median_log10_pareto", "median_log2_pareto"))
  expect_identical(multi$median_log10_auto$data_transformed,
                   multi$median_log10_pareto$data_transformed)
  expect_false(is.null(attr(multi, "shared")$data_corrected))
  expect_output(print(multi), "Combinations : 4")
  expect_output(summary(multi), "median_log2_pareto")
})

test_that("failures in steps 1-5 return a well-formed error object", {
  meta <- dims_metadata
  meta$Sample[2] <- meta$Sample[1]          # duplicate sample ID
  res <- suppressMessages(run_DIpreprocess(dims_data, meta, verbose = FALSE))
  expect_s3_class(res, "run_DIpreprocess")
  expect_match(res$error, "Duplicate sample IDs")
  expect_false(anyDuplicated(names(res)) > 0)
  expect_output(print(res), "FAILED")
})
