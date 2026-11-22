test_that(".validate_data_matrix rejects non-matrix/df", {
  expect_error(dimsprepr:::.validate_data_matrix("x"), "matrix or data frame")
  expect_error(dimsprepr:::.validate_data_matrix(1:5), "matrix or data frame")
  # Should not error on valid inputs
  expect_silent(dimsprepr:::.validate_data_matrix(matrix(1:4, 2)))
  expect_silent(dimsprepr:::.validate_data_matrix(data.frame(a = 1:3)))
})

test_that(".validate_metadata rejects mismatched rows", {
  x <- matrix(1:6, nrow = 2)
  meta <- data.frame(a = 1:5)
  expect_error(dimsprepr:::.validate_metadata(x, meta), "nrow")
})

test_that(".as_matrix_preserve and .restore_class round-trip", {
  df <- data.frame(a = 1:3, b = 4:6)
  mp <- dimsprepr:::.as_matrix_preserve(df)
  expect_false(mp$was_matrix)
  expect_true(is.matrix(mp$mat))

  restored <- dimsprepr:::.restore_class(mp$mat, mp$was_matrix)
  expect_true(is.data.frame(restored))

  mat <- as.matrix(df)
  mp2 <- dimsprepr:::.as_matrix_preserve(mat)
  expect_true(mp2$was_matrix)
  restored2 <- dimsprepr:::.restore_class(mp2$mat, mp2$was_matrix)
  expect_true(is.matrix(restored2))
})

test_that(".msg returns a function that respects verbose", {
  msg_on  <- dimsprepr:::.msg(TRUE)
  msg_off <- dimsprepr:::.msg(FALSE)
  expect_message(msg_on("hello"), "hello")
  expect_silent(msg_off("hello"))
})

test_that(".require_col errors on missing column", {
  df <- data.frame(a = 1)
  expect_error(dimsprepr:::.require_col(df, "b"), "not found")
  expect_silent(dimsprepr:::.require_col(df, "a"))
})

test_that(".with_seed restores RNG state", {
  # Set a known state
  set.seed(1)
  before <- runif(1)

  set.seed(1)
  # .with_seed should not permanently alter the RNG
  dimsprepr:::.with_seed(999, sample(100, 5))
  after <- runif(1)

  expect_equal(before, after)
})

test_that("method vocabularies are lower-case and unique", {
  for (v in list(dimsprepr:::.normalize_methods, dimsprepr:::.transform_methods,
                 dimsprepr:::.scale_methods, dimsprepr:::.plot_themes)) {
    expect_identical(v, unique(tolower(v)))
  }
})

test_that(".resolve_plot_what keeps found features and falls back when none match", {
  cols <- paste0("F", 1:10)
  expect_equal(dimsprepr:::.resolve_plot_what("F2", cols), "F2")
  expect_warning(out <- suppressMessages(dimsprepr:::.resolve_plot_what("zz", cols)),
                 "Falling back")
  expect_length(out, 4L)
  expect_equal(suppressMessages(dimsprepr:::.resolve_plot_what(NULL, cols[1:3])), cols[1:3])
})

test_that(".build_base_theme accepts every theme and rejects others", {
  for (th in dimsprepr:::.plot_themes)
    expect_s3_class(dimsprepr:::.build_base_theme(th, 11, "sans"), "theme")
  expect_error(dimsprepr:::.build_base_theme("neon", 11, "sans"), "must be one of")
})
