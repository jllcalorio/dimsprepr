# One feature with a known multiplicative drift, QC every 3rd injection,
# two batches with an offset. Returns raw intensities and the true drift.
make_drift_data <- function(seed = 1, n = 90, outlier = FALSE) {
  set.seed(seed)
  inj   <- seq_len(n)
  batch <- rep(1:2, each = n / 2)
  qc    <- inj %% 3 == 1
  drift <- exp(0.3 * sin(inj / 12)) * ifelse(batch == 2, 1.5, 1)
  bio   <- ifelse(qc, 1, exp(rnorm(n, 0, 0.3)))
  y     <- 1000 * drift * bio * exp(rnorm(n, 0, 0.02))
  if (outlier) y[which(qc)[8]] <- y[which(qc)[8]] * 4      # one gross QC outlier
  list(x = matrix(y, ncol = 1, dimnames = list(paste0("S", inj), "F1")),
       metadata = data.frame(Group = ifelse(qc, "QC", "Sample"), Batch = batch,
                             InjectionSequence = inj),
       qc = qc)
}
qc_cv <- function(res, d, drop = integer(0)) {
  v <- as.matrix(res$data)[d$qc, 1]
  if (length(drop)) v <- v[-drop]
  stats::sd(v) / mean(v)
}

test_that("QCRSC matches pmp::QCRSC on complete data", {
  skip_if_not_installed("pmp")
  full <- dims_clean(1:15)
  qc <- dims_metadata$Group == "QC"
  ours <- run_driftBatchCorrect(full, dims_metadata, verbose = FALSE)$data
  ref  <- t(suppressMessages(pmp::QCRSC(
    t(full), order = dims_metadata$InjectionSequence, batch = dims_metadata$Batch,
    classes = ifelse(qc, "QC", "Sample"), spar = 0, minQC = 5, qc_label = "QC")))
  # pmp fits on log((x + sqrt(x^2 + 1)) / 2) instead of log(x): tiny differences
  expect_equal(as.vector(ours), as.vector(ref), tolerance = 1e-3)
})

test_that("every trend method removes a known drift and batch offset", {
  d <- make_drift_data()
  raw_cv <- stats::sd(d$x[d$qc, 1]) / mean(d$x[d$qc, 1])
  methods <- c("QCRSC", "rLOESS", if (requireNamespace("mgcv", quietly = TRUE)) c("rGAM", "tGAM"))
  for (m in methods) {
    res <- run_driftBatchCorrect(d$x, d$metadata, method = m, verbose = FALSE)
    expect_equal(res$method_used, m)
    expect_lt(qc_cv(res, d), raw_cv / 4)
  }
})

test_that("robust methods resist a gross QC outlier better than QCRSC", {
  skip_if_not_installed("mgcv")
  d <- make_drift_data(outlier = TRUE)
  cv <- vapply(c("QCRSC", "rLOESS", "rGAM", "tGAM"), function(m)
    qc_cv(run_driftBatchCorrect(d$x, d$metadata, method = m, verbose = FALSE),
          d, drop = 8), numeric(1))
  expect_true(all(cv[c("rLOESS", "rGAM", "tGAM")] < cv[["QCRSC"]]))
})

test_that("missing values and zeros are ignored by the fit and preserved", {
  x <- dims_data[, 1:8]
  raw <- replace(x, x == 0, NA)
  for (m in c("QCRSC", "rLOESS")) {
    res <- run_driftBatchCorrect(raw, dims_metadata, method = m, verbose = FALSE)
    expect_identical(is.na(as.matrix(res$data)), is.na(raw[, colnames(res$data)]))
  }
  res0 <- run_driftBatchCorrect(x, dims_metadata, verbose = FALSE)
  expect_true(all(res0$data[x[, colnames(res0$data)] == 0] == 0))
})

test_that("features lacking QC values in a batch are flagged, kept as-is, or removed", {
  x <- replace(dims_data[, 1:5], dims_data[, 1:5] == 0, NA)
  b2qc <- which(dims_metadata$Group == "QC" & dims_metadata$Batch == 2)
  x[b2qc[1:3], "mz_0001"] <- NA
  kept <- run_driftBatchCorrect(x, dims_metadata, remove_uncorrected = FALSE, verbose = FALSE)
  expect_equal(kept$uncorrected, data.frame(Feature = "mz_0001", Batch = "2"))
  b2 <- dims_metadata$Batch == 2
  expect_equal(kept$data[b2, "mz_0001"], x[b2, "mz_0001"])
  expect_false(isTRUE(all.equal(kept$data[!b2, "mz_0001"], x[!b2, "mz_0001"])))
  dropped <- run_driftBatchCorrect(x, dims_metadata, verbose = FALSE)
  expect_false("mz_0001" %in% colnames(dropped$data))
  expect_output(print(dropped), "Uncorrected features: 1 \\(removed\\)")
})

test_that("character batch labels are accepted", {
  meta <- dims_metadata
  meta$Batch <- paste0("B", meta$Batch)
  res <- run_driftBatchCorrect(dims_clean(1:5), meta, verbose = FALSE)
  expect_equal(nrow(res$uncorrected), 0L)
  expect_false(anyNA(res$data))
})

test_that("the QC check refits on biological samples when QCs disagree", {
  skip_if_not_installed("mgcv")
  d <- make_drift_data(seed = 3)
  # QCs drift the other way than the samples (unrepresentative QCs)
  d$x[d$qc, 1] <- 1000 * exp(-0.3 * sin(which(d$qc) / 12)) * exp(rnorm(sum(d$qc), 0, 0.02))
  res <- run_driftBatchCorrect(d$x, d$metadata, method = "rGAM", sample_weight = 0.2,
                               qc_check = TRUE, verbose = FALSE)
  expect_true(all(res$qc_check$Sample_based))
  expect_true(all(res$qc_check$p_value < 0.1))
  expect_output(print(res), "QC check")
})

test_that("ComBat corrects batch offsets and handles a single batch", {
  skip_if_not_installed("sva")
  raw <- replace(dims_data[, 1:10], dims_data[, 1:10] == 0, NA)
  res <- run_driftBatchCorrect(raw, dims_metadata, method = "ComBat", verbose = FALSE)
  expect_equal(res$method_used, "ComBat")
  expect_identical(is.na(as.matrix(res$data)), is.na(raw))
  one <- dims_metadata
  one$Batch <- 1L
  expect_warning(
    res1 <- run_driftBatchCorrect(raw, one, method = "ComBat", verbose = FALSE),
    "Only one batch"
  )
  expect_equal(res1$data, raw)
})

test_that("run_driftBatchCorrect validates its inputs", {
  expect_error(run_driftBatchCorrect("x", dims_metadata), "matrix or data frame")
  expect_error(run_driftBatchCorrect(dims_data, "y"), "data frame")
  expect_error(run_driftBatchCorrect(dims_data, dims_metadata, method = "loess"), "one of")
  expect_error(run_driftBatchCorrect(dims_data, dims_metadata, min_QC = 3), ">= 4")
  expect_error(run_driftBatchCorrect(dims_data, dims_metadata, min_QC = 7),
               "fewer than min_QC")
  expect_error(run_driftBatchCorrect(dims_data, dims_metadata, qc_check = TRUE),
               "sample_weight > 0")
  expect_error(run_driftBatchCorrect(dims_data, dims_metadata, sample_weight = 2),
               "between 0 and 1")
  meta <- dims_metadata
  meta$InjectionSequence[2] <- 1L
  expect_error(run_driftBatchCorrect(dims_data, meta), "unique")
})
