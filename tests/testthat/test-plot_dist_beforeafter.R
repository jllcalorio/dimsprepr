test_that("plot_dist_beforeafter returns four ggplots and a gtable", {
  before <- log10(dims_clean())
  after  <- run_scale(before, "auto", verbose = FALSE)$data
  pdf(NULL); on.exit(dev.off())
  res <- suppressMessages(plot_dist_beforeafter(before, after, n_random = 10))
  expect_named(res, c("plot_density_before", "plot_density_after",
                      "plot_box_before", "plot_box_after", "plot_combined"))
  expect_s3_class(res$plot_box_after, "ggplot")
  expect_s3_class(res$plot_combined, "gtable")
  # n_random limits the boxes shown
  expect_equal(length(unique(res$plot_box_before$data$id)), 10L)
})

test_that("plot_dist_beforeafter feature mode and run_DIpreprocess input", {
  res <- run_DIpreprocess(dims_data, dims_metadata, correct_drift = FALSE, verbose = FALSE)
  pdf(NULL); on.exit(dev.off())
  out <- suppressMessages(plot_dist_beforeafter(res, group_by = "feature",
                                                plot_what = res$features_final[1:5]))
  expect_equal(length(unique(out$plot_box_after$data$key)), 5L)
})

test_that("plot_dist_beforeafter validates inputs", {
  x <- matrix(1:20, 5, dimnames = list(NULL, letters[1:4]))
  pdf(NULL); on.exit(dev.off())
  expect_error(plot_dist_beforeafter(x, x * 2, group_by = "batch"), "group_by")
  expect_error(plot_dist_beforeafter(x, x * 2, plot_what = "zz"), "not found")
  expect_error(plot_dist_beforeafter(x, x * 2, n_random = 0), "n_random")
  expect_warning(suppressMessages(plot_dist_beforeafter(x, x)), "identical")
})
