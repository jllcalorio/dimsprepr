test_that("run_transform log2 works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  td$x[td$x <= 0] <- 1
  res <- run_transform(td$x, method = "log2", verbose = FALSE)

  expect_s3_class(res, "run_transform")
  expect_equal(dim(res$data), dim(td$x))
  expect_equal(res$shift_applied, 0)
})

test_that("run_transform log10 adds shift for non-positive values", {
  mat <- matrix(c(-5, 0, 10, 20), nrow = 2)
  colnames(mat) <- c("A", "B")
  res <- run_transform(mat, method = "log10", verbose = FALSE)
  expect_true(res$shift_applied > 0)
  expect_true(all(is.finite(res$data)))
})

test_that("run_transform sqrt works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  td$x[td$x < 0] <- 0
  res <- run_transform(td$x, method = "sqrt", verbose = FALSE)
  expect_s3_class(res, "run_transform")
})

test_that("run_transform clr works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  td$x[td$x <= 0] <- 1
  res <- run_transform(td$x, method = "clr", verbose = FALSE)
  expect_s3_class(res, "run_transform")
  # CLR rows should sum to approximately 0

  row_sums <- rowSums(res$data)
  expect_true(all(abs(row_sums) < 1e-10))
})

test_that("run_transform rejects unknown methods", {
  td <- make_test_data()
  expect_error(run_transform(td$x, method = "bogus"), "Unknown")
})

test_that("print.run_transform works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  td$x[td$x <= 0] <- 1
  res <- run_transform(td$x, method = "log2", verbose = FALSE)
  expect_output(print(res), "Transformation Results")
})
