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

test_that("run_scale reports the groups used by supervised methods", {
  td <- make_test_data()
  td$x[is.na(td$x)] <- 100
  res_svast <- run_scale(td$x, method = "svast",
                          metadata = td$metadata, group_col = "Group",
                          verbose = FALSE)
  expect_identical(res_svast$parameters$group_col, "Group")
  expect_identical(res_svast$parameters$groups, sort(unique(td$metadata$Group)))
  expect_error(run_scale(td$x, method = "ln", verbose = FALSE), "Unknown scaling")
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

test_that("svast/xvast need metadata and return finite values", {
  x <- log10(dims_clean(1:10))
  expect_error(run_scale(x, "svast", verbose = FALSE), "metadata")
  for (m in c("svast", "xvast")) {
    res <- run_scale(x, m, metadata = dims_metadata, verbose = FALSE)
    expect_true(all(is.finite(res$data)))
    expect_equal(unname(colMeans(res$data)), rep(0, 10), tolerance = 1e-10)
  }
})
