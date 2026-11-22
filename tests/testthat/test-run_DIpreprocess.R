test_that("run_DIpreprocess single combination works end-to-end", {
  td <- make_test_data()
  res <- run_DIpreprocess(
    x = td$x, metadata = td$metadata,
    normalize_method = "median", transform_method = "log10",
    drift_method = "none", rsd_threshold = 0.5, variance_percentile = 5,
    verbose = FALSE
  )
  expect_s3_class(res, "run_DIpreprocess")
  expect_null(res$error)
  expect_true(nrow(res$data_nonpls) > 0)
  expect_true(length(res$features_final) > 0)
})

test_that("run_DIpreprocess multi-combination mode works", {
  td <- make_test_data()
  res <- run_DIpreprocess(
    x = td$x, metadata = td$metadata,
    normalize_method = c("sum", "median"), transform_method = "log10",
    drift_method = "none", rsd_threshold = 0.5, variance_percentile = 5,
    verbose = FALSE
  )
  expect_s3_class(res, "run_DIpreprocess_multi")
  expect_named(res, c("sum_log10", "median_log10"))
})

test_that("pre-flight checks reject bad arguments before running", {
  td <- make_test_data()
  expect_error(run_DIpreprocess(td$x, td$metadata, normalize_method = "bogus",
                                verbose = FALSE), "Unknown normalize_method")
  expect_error(run_DIpreprocess(td$x, td$metadata, scale_pls = "zzz", verbose = FALSE),
               "Unknown scale_pls")
  expect_error(run_DIpreprocess(td$x, td$metadata, drift_method = "magic", verbose = FALSE),
               "Unknown drift_method")
  expect_error(run_DIpreprocess(td$x, td$metadata, normalize_method = "pqn_reference",
                                verbose = FALSE), "normalize_ref_sample")
  expect_error(run_DIpreprocess(td$x, td$metadata, drift_args = list(spar = 1),
                                verbose = FALSE), "Unknown 'drift_args'")
  expect_error(run_DIpreprocess(td$x, td$metadata, drift_args = list(method = "tGAM"),
                                verbose = FALSE), "not 'drift_args'")
  expect_error(run_DIpreprocess(td$x, td$metadata, rsd_qc_type = "XQC", verbose = FALSE),
               "rsd_qc_type")
})

test_that("print and summary methods work", {
  td <- make_test_data()
  res <- run_DIpreprocess(td$x, td$metadata, drift_method = "none",
                          rsd_threshold = 0.5, variance_percentile = 5, verbose = FALSE)
  expect_output(print(res), "run_DIpreprocess")
  expect_output(summary(res), "step-by-step")
})

test_that("default pipeline runs end-to-end on the bundled data (with QCRSC)", {
  res <- run_DIpreprocess(dims_data, dims_metadata, verbose = FALSE)
  expect_null(res$error)
  expect_equal(res$parameters$drift_method, "QCRSC")
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
})

test_that("drift correction runs on unimputed data and imputation follows it", {
  res <- run_DIpreprocess(dims_data, dims_metadata, verbose = FALSE)
  corr <- as.matrix(res$data_corrected)
  expect_true(anyNA(corr))                              # NAs survive correction
  expect_identical(is.na(corr), is.na(as.matrix(res$data_missing_filtered))[, colnames(corr)])
  expect_false(anyNA(res$data_imputed))
  # observed values equal a direct run_driftBatchCorrect() on the same input
  direct <- run_driftBatchCorrect(res$data_missing_filtered, dims_metadata, verbose = FALSE)
  expect_equal(as.matrix(direct$data), corr)
})

test_that("quality filtering matches the step functions on unscaled data", {
  # Regression: filtering used to run after auto-scaling (all variances = 1,
  # QC means ~0), so the kept set was arbitrary. RSD uses measured values only.
  res <- run_DIpreprocess(dims_data, dims_metadata, drift_method = "none",
                          variance_percentile = 25, verbose = FALSE)
  norm <- run_normalize(res$data_imputed, res$metadata, method = "median",
                        group_col = "Group_", qc_types = "QC", verbose = FALSE)
  measured <- norm$data
  measured[is.na(res$data_corrected)] <- NA
  rsd  <- run_filterRSD(measured, res$metadata, max_rsd = 0.3, verbose = FALSE)
  trans <- run_transform(norm$data[, rsd$features_kept], "log10", verbose = FALSE)
  vf   <- run_filtervariance(trans$data, percentile = 25, verbose = FALSE)
  expect_identical(res$features_final, vf$features_kept)
})

test_that("rsd_qc_type and track are honoured", {
  td <- make_test_data()
  td$metadata$Group[td$metadata$Group == "QC"] <- c("EQC", "EQC", "EQC", "SQC", "SQC", "SQC")
  res_e <- run_DIpreprocess(td$x, td$metadata, drift_method = "none", rsd_threshold = 0.5,
                            rsd_qc_type = "EQC", variance_percentile = 5, verbose = FALSE)
  expect_null(res_e$error)
  expect_equal(res_e$parameters$rsd_qc_type, "EQC")

  res_t <- run_DIpreprocess(dims_data, dims_metadata, drift_method = "none",
                            track = c("mz_0078", "mz_0001"), verbose = FALSE)
  expect_true("mz_0078" %in% res_t$track$Feature)
  expect_false("mz_0001" %in% res_t$track$Feature)
})

test_that("remove_uncorrected drops or keeps features QCRSC cannot correct", {
  x <- dims_data
  b2qc <- which(dims_metadata$Group == "QC" & dims_metadata$Batch == 2)
  x[b2qc[1:3], "mz_0005"] <- 0                         # 3 of 6 QCs left
  drop <- run_DIpreprocess(x, dims_metadata, track = "mz_0005", verbose = FALSE)
  keep <- run_DIpreprocess(x, dims_metadata, remove_uncorrected = FALSE, verbose = FALSE)
  expect_identical(drop$uncorrected_features, "mz_0005")
  expect_false("mz_0005" %in% colnames(drop$data_corrected))
  expect_match(drop$track$Value[drop$track$Feature == "mz_0005"], "batch\\(es\\) 2")
  expect_true("mz_0005" %in% colnames(keep$data_corrected))
  b2 <- dims_metadata$Batch == 2
  expect_equal(keep$data_corrected[b2, "mz_0005"],
               keep$data_missing_filtered[b2, "mz_0005"])  # left uncorrected
})

test_that("robust drift methods run through drift_args", {
  skip_if_not_installed("mgcv")
  res <- run_DIpreprocess(dims_data[, 1:12], dims_metadata, drift_method = "rGAM",
                          drift_args = list(sample_weight = 0.2), verbose = FALSE)
  expect_null(res$error)
  expect_equal(res$parameters$drift_args, list(sample_weight = 0.2))
})

test_that("blank filtering, truly absent zeros, and failed samples are handled", {
  blanks <- matrix(0, 3, ncol(dims_data), dimnames = list(c("B1", "B2", "B3"),
                                                          colnames(dims_data)))
  blanks[, "mz_0011"] <- c(4000, 2500, 6000) *
    stats::median(dims_data[dims_data[, "mz_0011"] > 0, "mz_0011"]) / 1000
  x <- rbind(dims_data, blanks)
  x["S02", ] <- replace(x["S02", ], -(1:5), 0)          # failed injection
  meta <- rbind(dims_metadata, data.frame(Sample = rownames(blanks), Group = "Blank",
                                          Batch = 2L, InjectionSequence = 61:63,
                                          SubjectID = NA))
  res <- run_DIpreprocess(x, meta, blank_types = "Blank", missing_absent_threshold = 0.7,
                          verbose = FALSE)
  expect_null(res$error)
  expect_false(any(c("B1", "B2", "B3") %in% rownames(res$data_nonpls)))
  expect_false("mz_0011" %in% res$features_final)
  expect_identical(res$samples_removed, "S02")
  expect_identical(rownames(res$data_raw), rownames(res$metadata))
  expect_match(capture.output(print(res))[4], "1 removed")
})

test_that("merge_replicates works with a custom subject column name", {
  meta <- dims_metadata
  meta$Subject <- meta$SubjectID
  meta$SubjectID <- NULL
  meta$Subject[meta$Subject %in% c("BIO02")] <- "BIO01"   # BIO01 measured twice
  res <- run_DIpreprocess(dims_data, meta, subject_id_col = "Subject",
                          drift_method = "none", merge_replicates = TRUE, verbose = FALSE)
  expect_null(res$error)
  expect_equal(nrow(res$data_nonpls_merged), nrow(dims_data) - 1L)
  expect_equal(nrow(res$metadata_merged), nrow(dims_data) - 1L)
})

test_that("errors inside steps are captured in $error, not thrown", {
  res <- suppressMessages(
    run_DIpreprocess(dims_data, dims_metadata, drift_method = "none",
                     rsd_threshold = 0, verbose = FALSE)
  )
  expect_match(res$error, "RSD filter")
  expect_output(print(res), "FAILED")
})

test_that("multi mode names, shares steps 1-5, and prints", {
  multi <- run_DIpreprocess(dims_data, dims_metadata, normalize_method = "median",
                            transform_method = c("log10", "log2"),
                            scale_nonpls = c("auto", "pareto"),
                            drift_method = "none", verbose = FALSE)
  expect_named(multi, c("median_log10_auto", "median_log2_auto",
                        "median_log10_pareto", "median_log2_pareto"))
  expect_identical(multi$median_log10_auto$data_transformed,
                   multi$median_log10_pareto$data_transformed)
  expect_false(is.null(attr(multi, "shared")$data_imputed))
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
