test_that("run_normalize sum method works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  td$x[td$x == 0] <- 1
  res <- run_normalize(td$x, td$metadata, method = "sum", verbose = FALSE)

  expect_s3_class(res, "run_normalize")
  expect_equal(dim(res$data), dim(td$x))
})

test_that("run_normalize median method works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  td$x[td$x == 0] <- 1
  res <- run_normalize(td$x, td$metadata, method = "median", verbose = FALSE)
  expect_s3_class(res, "run_normalize")
})

test_that("run_normalize pqn_global method works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  td$x[td$x == 0] <- 1
  res <- run_normalize(td$x, td$metadata, method = "pqn_global", verbose = FALSE)
  expect_s3_class(res, "run_normalize")
})

test_that("run_normalize pqn_group method works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  td$x[td$x == 0] <- 1
  res <- run_normalize(td$x, td$metadata, method = "pqn_group", verbose = FALSE)
  expect_s3_class(res, "run_normalize")
})

test_that("removed duplicate methods are rejected with a clear error", {
  td <- make_test_data()
  for (m in c("quantile", "pqn", "row_rel_abundance", "col_rel_abundance"))
    expect_error(run_normalize(td$x, td$metadata, method = m, verbose = FALSE), "Unknown")
})

test_that("PQN equals the median quotient to the reference spectrum", {
  x <- matrix(c(2, 4, 6, 1, 2, 3, 4, 8, 12), nrow = 3,
              dimnames = list(NULL, c("a", "b", "c")))   # row r = r * (2, 1, 4)
  meta <- data.frame(Group = c("QC", "A", "A"))
  res <- run_normalize(x, meta, method = "pqn_group", verbose = FALSE)
  expect_equal(res$normalization_factors, c(1, 2, 3))
  expect_equal(unname(res$data), matrix(rep(c(2, 1, 4), each = 3), 3))
  mean_ref <- run_normalize(x, meta, method = "pqn_global", reference_method = "mean",
                            verbose = FALSE)
  expect_equal(mean_ref$normalization_factors, c(0.5, 1, 1.5))
})

test_that("run_normalize rejects bad inputs", {
  td <- make_test_data()
  expect_error(run_normalize("x", td$metadata), "matrix or data frame")
  expect_error(run_normalize(td$x, td$metadata, method = "bogus"), "Unknown")
})

test_that("print.run_normalize works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  res <- run_normalize(td$x, td$metadata, method = "sum", verbose = FALSE)
  expect_output(print(res), "Normalization Results")
})

test_that("pqn_group honours group_sample", {
  x <- matrix(c(1, 2, 4, 8, 1, 2, 4, 8) * 10, nrow = 4,
              dimnames = list(NULL, c("a", "b")))
  meta <- data.frame(Group = c("SQC", "EQC", "A", "A"))
  sqc <- run_normalize(x, meta, method = "pqn_group", group_sample = "SQC", verbose = FALSE)
  eqc <- run_normalize(x, meta, method = "pqn_group", group_sample = "EQC", verbose = FALSE)
  expect_equal(sqc$normalization_factors[1], 1)   # SQC row is its own reference
  expect_equal(eqc$normalization_factors[2], 1)
  expect_error(run_normalize(x, meta, method = "pqn_group", group_sample = "zzz",
                             verbose = FALSE), "No 'zzz' samples")
})

test_that("run_normalize validates qc_normalize and reference_method", {
  x <- dims_clean(1:10)
  expect_error(run_normalize(x, dims_metadata, qc_normalize = "bogus"), "should be one of")
  expect_error(run_normalize(x, dims_metadata, reference_method = "mode"), "should be one of")
  res <- run_normalize(x, dims_metadata, method = "pqn_group", verbose = FALSE)
  expect_false(anyNA(res$normalization_factors))
})

test_that("specific_factor divides by the metadata column", {
  x <- matrix(c(2, 4, 6, 8), nrow = 2, dimnames = list(NULL, c("a", "b")))
  meta <- data.frame(Group = c("A", "A"), Osmo = c(2, 4))
  res <- run_normalize(x, meta, method = "Osmo", verbose = FALSE)
  expect_equal(unname(res$data), matrix(c(1, 1, 3, 2), nrow = 2))
})
