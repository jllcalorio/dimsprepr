test_that("run_mvimpute deterministic fills NAs", {
  td <- make_test_data()
  # Ensure some NAs exist
  td$x[1, 1] <- NA
  res <- run_mvimpute(td$x, method = 0.2, verbose = FALSE)

  expect_s3_class(res, "run_mvimpute")
  expect_equal(res$n_missing_after, 0L)
  expect_true(!any(is.na(res$data)))
})

test_that("run_mvimpute with method='none' preserves NAs", {
  td <- make_test_data()
  n_na <- sum(is.na(td$x))
  res <- run_mvimpute(td$x, method = "none", verbose = FALSE)
  expect_equal(res$n_missing_after, n_na)
})

test_that("run_mvimpute returns unchanged data when no NAs", {
  mat <- matrix(1:20, nrow = 4)
  colnames(mat) <- paste0("F", 1:5)
  res <- run_mvimpute(mat, method = 0.5, verbose = FALSE)
  expect_equal(res$n_missing_before, 0L)
  expect_equal(res$method_used, "None (no missing values)")
})

test_that("run_mvimpute warns when method > 1", {
  td <- make_test_data()
  td$x[1, 1] <- NA
  expect_warning(run_mvimpute(td$x, method = 2, verbose = FALSE), "exceed")
})

test_that("run_mvimpute rejects bad inputs", {
  expect_error(run_mvimpute("x"), "matrix or data frame")
  td <- make_test_data()
  expect_error(run_mvimpute(td$x, method = -1), "positive")
})

test_that("print.run_mvimpute works", {
  td <- make_test_data()
  res <- run_mvimpute(td$x, method = 0.1, verbose = FALSE)
  expect_output(print(res), "Missing Value Imputation")
})
