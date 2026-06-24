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
