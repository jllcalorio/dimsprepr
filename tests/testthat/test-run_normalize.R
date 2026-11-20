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

test_that("run_normalize quantile maps to pqn", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  td$x[td$x == 0] <- 1
  # quantile should not error (it aliases to pqn); may warn if pmp not installed
  skip_if_not_installed("pmp")
  res <- run_normalize(td$x, td$metadata, method = "quantile", verbose = FALSE)
  expect_s3_class(res, "run_normalize")
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

test_that("run_normalize validates qc_normalize and reports pqn factors", {
  x <- dims_clean(1:10)
  expect_error(run_normalize(x, dims_metadata, qc_normalize = "bogus"), "should be one of")
  skip_if_not_installed("pmp")
  res <- run_normalize(x, dims_metadata, method = "pqn", verbose = FALSE)
  expect_false(anyNA(res$normalization_factors))
})

test_that("specific_factor divides by the metadata column", {
  x <- matrix(c(2, 4, 6, 8), nrow = 2, dimnames = list(NULL, c("a", "b")))
  meta <- data.frame(Group = c("A", "A"), Osmo = c(2, 4))
  res <- run_normalize(x, meta, method = "Osmo", verbose = FALSE)
  expect_equal(unname(res$data), matrix(c(1, 1, 3, 2), nrow = 2))
})
