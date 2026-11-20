plot_fixture <- function() {
  before <- dims_clean(1:6)
  after  <- run_driftBatchCorrect(before, dims_metadata, spline_smooth_param = 0.5,
                                  verbose = FALSE)$data
  list(before = before, after = after)
}

test_that("plot_beforeafter returns a gtable for chosen features", {
  pf <- plot_fixture()
  pdf(NULL); on.exit(dev.off())
  g <- suppressMessages(plot_beforeafter(pf$before, pf$after, dims_metadata,
                                         plot_what = c("mz_0001", "mz_0002")))
  expect_s3_class(g, "gtable")
})

test_that("plot_beforeafter random selection is reproducible and leaves RNG alone", {
  pf <- plot_fixture()
  pdf(NULL); on.exit(dev.off())
  set.seed(5); before_rng <- runif(1)
  set.seed(5)
  expect_message(plot_beforeafter(pf$before, pf$after, dims_metadata, plot_what = NULL,
                                  seed = 1), "Automatically selected 6")
  expect_equal(runif(1), before_rng)
})

test_that("plot_beforeafter accepts a run_DIpreprocess result", {
  res <- run_DIpreprocess(dims_data[, 1:30], dims_metadata, spline_smooth = 0.5,
                          variance_percentile = 0, verbose = FALSE)
  pdf(NULL); on.exit(dev.off())
  expect_s3_class(suppressMessages(plot_beforeafter(res, plot_what = "mz_0003")), "gtable")
})

test_that("plot_beforeafter validates inputs", {
  pf <- plot_fixture()
  pdf(NULL); on.exit(dev.off())
  expect_error(plot_beforeafter(pf$before, pf$after[, 1:3], dims_metadata), "same dimensions")
  expect_error(plot_beforeafter(pf$before, pf$after, dims_metadata[1:5, ]), "same number of rows")
  expect_error(plot_beforeafter(pf$before, pf$after, dims_metadata, col_batch = "nope"),
               "not found")
  expect_error(suppressMessages(plot_beforeafter(pf$before, pf$after, dims_metadata,
                                                 plot_what = "mz_0001", theme = "neon")),
               "'theme' must be one of")
  expect_warning(suppressMessages(plot_beforeafter(pf$before, pf$after, dims_metadata,
                                                   plot_what = c("mz_0001", "ghost"))),
                 "not found")
})
