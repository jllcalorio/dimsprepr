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

test_that("global filtering also excludes QC samples when include_QC = FALSE", {
  x <- matrix(1, nrow = 10, ncol = 2, dimnames = list(NULL, c("bio_ok", "qc_gap")))
  x[9:10, "qc_gap"] <- NA        # missing only in the two QC rows
  meta <- data.frame(Group = c(rep("A", 8), "QC", "QC"))
  ex <- run_filtermissing(x, meta, threshold = 0.1, filter_by_group = FALSE,
                          include_QC = FALSE, verbose = FALSE)
  inc <- run_filtermissing(x, meta, threshold = 0.1, filter_by_group = FALSE,
                           include_QC = TRUE, verbose = FALSE)
  expect_equal(ex$features_kept, c("bio_ok", "qc_gap"))
  expect_equal(inc$features_kept, "bio_ok")
})

test_that("group-wise filtering works with a single feature", {
  x <- matrix(c(1, 1, NA, NA), ncol = 1, dimnames = list(NULL, "F1"))
  meta <- data.frame(Group = c("A", "A", "B", "B"))
  any_g <- run_filtermissing(x, meta, threshold = 0.5, verbose = FALSE)
  all_g <- run_filtermissing(x, meta, threshold = 0.5, any_group = FALSE, verbose = FALSE)
  expect_equal(any_g$n_features_after, 1L)
  expect_equal(all_g$n_features_after, 0L)
})

test_that("bundled data: almost-empty features are removed", {
  res <- run_filtermissing(dims_data, dims_metadata, verbose = FALSE)
  expect_true(all(c("mz_0078", "mz_0079", "mz_0080") %in% res$features_removed))
})

test_that("a feature observed in exactly 80% of a group is kept (80% rule)", {
  x <- matrix(c(1, 1, 1, 1, NA), ncol = 1, dimnames = list(NULL, "F1"))
  meta <- data.frame(Group = rep("A", 5))
  expect_equal(run_filtermissing(x, meta, threshold = 0.2, verbose = FALSE)$n_features_after, 1L)
})

test_that("absent_threshold keeps mostly-missing features as zeros in that group", {
  x <- cbind(drug = c(5, NA, NA, NA, NA, 4, 6, 5, 7, 8),
             noise = c(1, NA, NA, 2, NA, NA, 3, NA, 4, NA))   # 60% missing in A and B
  meta <- data.frame(Group = rep(c("A", "B"), each = 5))
  res <- run_filtermissing(x, meta, absent_threshold = 0.7, verbose = FALSE)
  expect_identical(res$features_absent, "drug")          # 80% missing in A
  expect_identical(res$features_kept, "drug")            # noise: 60-80% missing
  expect_equal(unname(res$data[1:5, "drug"]), c(5, 0, 0, 0, 0))
  expect_equal(res$n_absent_zeros, 4L)
  expect_error(run_filtermissing(x, meta, absent_threshold = 0.1), "between")
})

test_that("sample_iqr removes samples with abnormally high missingness", {
  x <- dims_data
  x["S02", -(1:5)] <- 0                                   # failed injection
  res <- run_filtermissing(x, dims_metadata, sample_iqr = 3, verbose = FALSE)
  expect_identical(res$samples_removed, "S02")
  expect_equal(nrow(res$data), nrow(x) - 1L)
  expect_identical(rownames(res$data), res$metadata$Sample)
  expect_gt(res$sample_missing_rate[["S02"]], 0.9)
  # off by default, and nothing flagged in clean data
  expect_length(run_filtermissing(x, dims_metadata, verbose = FALSE)$samples_removed, 0L)
  expect_length(run_filtermissing(dims_data, dims_metadata, sample_iqr = 3,
                                  verbose = FALSE)$samples_removed, 0L)
})
