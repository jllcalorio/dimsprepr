blank_fixture <- function() {
  set.seed(9)
  x <- cbind(clean   = c(rlnorm(8, 7), 0, 0, 0),            # not in blanks
             stable  = c(rlnorm(8, 4, 0.1), 50, 52, 51),    # contaminant, stable blank
             noisy   = c(rlnorm(8, 4, 0.1), 20, 90, 45),    # contaminant, unstable blank
             high    = c(rlnorm(8, 9), 30, 40, 35))         # well above blank level
  meta <- data.frame(Group = c(rep("Case", 6), "QC", "QC", rep("Blank", 3)))
  list(x = x, metadata = meta)
}

test_that("run_filterblank flags contaminants and removes unstable ones", {
  bf <- blank_fixture()
  res <- run_filterblank(bf$x, bf$metadata, verbose = FALSE)
  expect_s3_class(res, "run_filterblank")
  expect_setequal(res$features_flagged, c("stable", "noisy"))
  expect_identical(res$features_removed, "noisy")
  expect_identical(colnames(res$data), c("clean", "stable", "high"))
  expect_equal(nrow(res$data), 8L)                        # blanks dropped
  expect_false("Blank" %in% res$metadata$Group)
  expect_equal(res$n_blanks, 3L)
  expect_output(print(res), "Blank Filtering Results")
})

test_that("run_filterblank thresholds and inputs are validated", {
  bf <- blank_fixture()
  strict <- run_filterblank(bf$x, bf$metadata, max_blank_rsd = 0, verbose = FALSE)
  expect_setequal(strict$features_removed, c("stable", "noisy"))
  lax <- run_filterblank(bf$x, bf$metadata, fold_change = 0.1, verbose = FALSE)
  expect_length(lax$features_flagged, 0L)
  expect_error(run_filterblank(bf$x, bf$metadata, blank_types = "Solvent"),
               "At least 2 blank")
  expect_error(run_filterblank(bf$x, bf$metadata, fold_change = -1), "positive")
})
