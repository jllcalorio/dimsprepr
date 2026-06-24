test_that("run_scale auto-scaling works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  res <- run_scale(td$x, method = "auto", verbose = FALSE)

  expect_s3_class(res, "run_scale")
  expect_equal(dim(res$data), dim(td$x))
  # Auto-scaled columns should have mean ~0
  col_means <- colMeans(res$data)
  expect_true(all(abs(col_means) < 1e-10))
})

test_that("run_scale pareto works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  res <- run_scale(td$x, method = "pareto", verbose = FALSE)
  expect_s3_class(res, "run_scale")
})

test_that("run_scale mean-centering works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  res <- run_scale(td$x, method = "mean", verbose = FALSE)
  col_means <- colMeans(res$data)
  expect_true(all(abs(col_means) < 1e-10))
})

test_that("run_scale none returns unchanged data", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  res <- run_scale(td$x, method = "none", verbose = FALSE)
  expect_equal(as.numeric(res$data), as.numeric(td$x))
})

test_that("run_scale group bug is fixed (levels(group_vec) not levels(group))", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  res_svast <- run_scale(td$x, method = "svast",
                          metadata = td$metadata, group_col = "Group",
                          verbose = FALSE)
  # Should correctly report group levels, not NULL
  expect_true(!is.null(res_svast$parameters$group_col))
  expect_true(length(res_svast$parameters$group_col) > 0)
})

test_that("run_scale vast works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  res <- run_scale(td$x, method = "vast", verbose = FALSE)
  expect_s3_class(res, "run_scale")
})

test_that("run_scale rejects bad inputs", {
  expect_error(run_scale("x"), "matrix or data frame")
  td <- make_test_data()
  expect_error(run_scale(td$x, method = "bogus"), "Unknown")
  expect_error(run_scale(td$x, method = "svast"), "metadata")
})

test_that("print.run_scale works", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  res <- run_scale(td$x, method = "auto", verbose = FALSE)
  expect_output(print(res), "Scaling Results")
})
