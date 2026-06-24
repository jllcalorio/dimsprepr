#' Impute Missing Values in Metabolomics Data
#'
#' @description
#' Replaces missing values (NAs) using either a simple deterministic method
#' (fraction of minimum), a statistical imputation method (quantile regression
#' imputation of left-censored data, QRILC), or any method available in
#' \code{\link[mice]{mice}} (Multiple Imputation by Chained Equations).
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features
#'   in columns. Missing values should be represented as `NA`.
#' @param method Numeric or character. Controls the imputation strategy:
#'   \itemize{
#'     \item **Numeric** -- fraction of the smallest observed value per feature
#'       used for deterministic imputation (e.g., `0.2` = 1/5th of the
#'       smallest value). By default (`positive_only = TRUE`), the minimum is
#'       computed from strictly positive values only.
#'     \item **`"quantileregression"`** -- uses the QRILC algorithm from the
#'       \pkg{imputeLCMD} package. Requires \pkg{imputeLCMD}.
#'     \item **Any `mice` method string** -- e.g., `"pmm"`, `"norm"`,
#'       `"rf"`, `"cart"`, etc. Requires \pkg{mice}.
#'     \item **`"none"`** -- Skips imputation entirely.
#'   }
#'   Default: `0.2`.
#' @param positive_only Logical. Only used when `method` is numeric. If `TRUE`
#'   (default), per-feature minimum is computed from strictly positive values.
#' @param tune_sigma Numeric. Only used when `method = "quantileregression"`.
#'   Controls the standard deviation of the left-censored distribution.
#'   Default: `1`.
#' @param m Integer. Only used with \pkg{mice} methods. Number of imputations.
#'   Default: `5`.
#' @param maxit Integer. Only used with \pkg{mice} methods. Number of
#'   iterations. Default: `5`.
#' @param seed Integer or `NA`. Only used with \pkg{mice} methods. Random seed.
#'   Default: `NA`.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#' @param ... Additional arguments passed to \code{\link[mice]{mice}}.
#'
#' @return A list of class `"run_mvimpute"` containing:
#'   \item{data}{Matrix or data frame with imputed values (same class as input `x`).}
#'   \item{n_missing_before}{Integer. Number of missing values before imputation.}
#'   \item{n_missing_after}{Integer. Number of missing values after imputation.}
#'   \item{imputed_summary}{Data frame with per-feature imputation statistics.}
#'   \item{parameters}{List of parameters used.}
#'   \item{method_used}{Character string describing the imputation method applied.}
#'
#' @details
#' **Missing Value Imputation Strategies:**
#'
#' 1. **Deterministic (fraction of minimum):** Each missing value is replaced by
#'    a fraction of the per-feature minimum. Standard in metabolomics for MNAR data.
#'
#' 2. **Quantile Regression (QRILC):** Models missing values from a
#'    left-censored distribution. Requires \pkg{imputeLCMD}.
#'
#' 3. **mice:** Multiple Imputation by Chained Equations. Requires \pkg{mice}.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate
#' Imputation by Chained Equations in R. *Journal of Statistical Software*,
#' 45(3), 1-67. \doi{10.18637/jss.v045.i03}
#'
#' Wei, R., et al. (2018). Missing Value Imputation Approach for Mass
#' Spectrometry-based Metabolomics Data. *Scientific Reports*, 8(1), 663.
#' \doi{10.1038/s41598-017-19120-0}
#'
#' @seealso \code{\link{run_DIpreprocess}}, \code{\link[mice]{mice}},
#'   \code{\link[imputeLCMD]{impute.QRILC}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(723)
#' x <- matrix(abs(rnorm(100 * 50, mean = 100)), nrow = 100, ncol = 50)
#' x[sample(length(x), 200)] <- NA
#' colnames(x) <- paste0("Feature", 1:50)
#'
#' result <- run_mvimpute(x, method = 0.2)
#' result
#' }
run_mvimpute <- function(
    x,
    method        = 0.2,
    positive_only = TRUE,
    tune_sigma    = 1,
    m             = 5,
    maxit         = 5,
    seed          = NA,
    verbose       = TRUE,
    ...
) {

  msg <- .msg(verbose)

  .validate_data_matrix(x)
  if (!is.logical(positive_only) || length(positive_only) != 1L)
    stop("'positive_only' must be a single logical value.", call. = FALSE)

  mp <- .as_matrix_preserve(x)
  x_matrix <- mp$mat
  n_missing_before <- sum(is.na(x_matrix))

  if (n_missing_before == 0L) {
    msg("No missing values detected. Returning original data.")
    return(structure(
      list(
        data             = .restore_class(x_matrix, mp$was_matrix),
        n_missing_before = 0L,
        n_missing_after  = 0L,
        imputed_summary  = data.frame(Feature = colnames(x_matrix),
                                       N_Missing = 0L, Imputation_Value = NA_real_),
        parameters  = list(method = method, positive_only = positive_only,
                           tune_sigma = tune_sigma, m = m, maxit = maxit, seed = seed),
        method_used = "None (no missing values)"
      ),
      class = "run_mvimpute"
    ))
  }

  msg(sprintf("Starting missing value imputation (%d missing values)...", n_missing_before))

  # ---------------------------------------------------------------------------
  # Dispatch
  # ---------------------------------------------------------------------------

  if (is.numeric(method)) {

    # -- Deterministic --------------------------------------------------------
    if (length(method) != 1L) stop("'method' when numeric must be a scalar.", call. = FALSE)
    if (method <= 0) stop("'method' when numeric must be positive.", call. = FALSE)
    if (method > 1) warning("'method' > 1: imputed values will exceed the observed minimum.", call. = FALSE)

    msg(sprintf("Using deterministic imputation: %.4f x minimum %s value per feature.",
                method, if (positive_only) "positive" else "observed"))

    na_indices <- is.na(x_matrix)

    # Vectorized min computation via matrixStats
    if (positive_only) {
      x_pos <- x_matrix
      x_pos[x_pos <= 0 | is.na(x_pos)] <- NA
      min_vals <- matrixStats::colMins(x_pos, na.rm = TRUE)
      min_vals[!is.finite(min_vals)] <- 1e-9
    } else {
      min_vals <- matrixStats::colMins(x_matrix, na.rm = TRUE)
      min_vals[!is.finite(min_vals)] <- 1e-9
    }

    imputation_vals <- min_vals * method
    x_matrix[na_indices] <- imputation_vals[col(x_matrix)][na_indices]

    method_description <- sprintf("Deterministic: %.4f x min %s",
                                  method, if (positive_only) "positive" else "observed (all)")
    imputed_summary <- data.frame(
      Feature = colnames(x_matrix), N_Missing = colSums(na_indices),
      Min_Value = min_vals, Imputation_Value = imputation_vals
    )
    params_out <- list(method = method, positive_only = positive_only,
                       tune_sigma = NA, m = NA, maxit = NA, seed = NA)

  } else if (is.character(method)) {

    method_lc <- tolower(trimws(method))

    if (method_lc == "none") {
      msg("Method set to 'none'. No imputation performed.")
      return(structure(
        list(
          data             = .restore_class(x_matrix, mp$was_matrix),
          n_missing_before = n_missing_before,
          n_missing_after  = n_missing_before,
          imputed_summary  = data.frame(Feature = colnames(x_matrix),
                                         N_Missing = colSums(is.na(x_matrix)), Method = "None"),
          parameters  = list(method = method, positive_only = positive_only,
                             tune_sigma = tune_sigma, m = m, maxit = maxit, seed = seed),
          method_used = "None"
        ),
        class = "run_mvimpute"
      ))

    } else if (method_lc == "quantileregression") {

      # -- QRILC --------------------------------------------------------------
      if (!requireNamespace("imputeLCMD", quietly = TRUE))
        stop("Package 'imputeLCMD' required. Install with: BiocManager::install('imputeLCMD')", call. = FALSE)
      if (!is.numeric(tune_sigma) || length(tune_sigma) != 1L || tune_sigma <= 0)
        stop("'tune_sigma' must be a single positive numeric value.", call. = FALSE)

      msg(sprintf("Using QRILC imputation (tune_sigma = %.2f).", tune_sigma))
      na_indices <- is.na(x_matrix)

      tryCatch({
        x_imputed <- imputeLCMD::impute.QRILC(t(x_matrix), tune.sigma = tune_sigma)
        x_matrix  <- t(x_imputed[[1]])
      }, error = function(e) {
        stop("QRILC imputation failed: ", e$message, call. = FALSE)
      })

      method_description <- sprintf("QRILC (tune_sigma = %.2f)", tune_sigma)
      imputed_summary <- data.frame(Feature = colnames(x_matrix),
                                     N_Missing = colSums(na_indices), Method = "QRILC")
      params_out <- list(method = method, positive_only = NA,
                         tune_sigma = tune_sigma, m = NA, maxit = NA, seed = NA)

    } else {

      # -- mice ---------------------------------------------------------------
      if (!requireNamespace("mice", quietly = TRUE))
        stop("Package 'mice' required. Install with: install.packages('mice')", call. = FALSE)
      if (!is.numeric(m) || length(m) != 1L || m < 1L)
        stop("'m' must be a single positive integer.", call. = FALSE)
      if (!is.numeric(maxit) || length(maxit) != 1L || maxit < 1L)
        stop("'maxit' must be a single positive integer.", call. = FALSE)

      msg(sprintf("Using mice imputation (method = '%s', m = %d, maxit = %d).",
                  method, as.integer(m), as.integer(maxit)))
      na_indices <- is.na(x_matrix)

      dots <- list(...)
      if (is.null(dots$predictorMatrix)) {
        dots$predictorMatrix <- tryCatch(
          mice::quickpred(as.data.frame(x_matrix)),
          error = function(e) stop("mice::quickpred() failed: ", e$message, call. = FALSE)
        )
        msg("Predictor matrix built via mice::quickpred().")
      }

      tryCatch({
        mice_obj <- do.call(mice::mice, c(
          list(data = as.data.frame(x_matrix), method = method,
               m = as.integer(m), maxit = as.integer(maxit),
               seed = if (is.na(seed)) NA_integer_ else as.integer(seed),
               printFlag = FALSE),
          dots
        ))
        x_matrix <- as.matrix(mice::complete(mice_obj, action = 1L))
      }, error = function(e) {
        stop("mice imputation failed: ", e$message, call. = FALSE)
      })

      method_description <- sprintf("mice: %s (m = %d, maxit = %d)",
                                    method, as.integer(m), as.integer(maxit))
      imputed_summary <- data.frame(Feature = colnames(x_matrix),
                                     N_Missing = colSums(na_indices), Method = method)
      params_out <- list(method = method, positive_only = NA, tune_sigma = NA,
                         m = as.integer(m), maxit = as.integer(maxit),
                         seed = if (is.na(seed)) NA_integer_ else as.integer(seed))
    }

  } else {
    stop("'method' must be numeric, 'quantileregression', 'none', or a mice method string.", call. = FALSE)
  }

  n_missing_after <- sum(is.na(x_matrix))
  msg(sprintf("Imputation complete. Missing values: %d -> %d.", n_missing_before, n_missing_after))

  structure(
    list(
      data             = .restore_class(x_matrix, mp$was_matrix),
      n_missing_before = n_missing_before,
      n_missing_after  = n_missing_after,
      imputed_summary  = imputed_summary,
      parameters       = params_out,
      method_used      = method_description
    ),
    class = "run_mvimpute"
  )
}

#' @param x Object to print.
#' @param ... Ignored.
#' @rdname run_mvimpute
#' @export
print.run_mvimpute <- function(x, ...) {
  cat("=== Missing Value Imputation Results ===\n")
  cat(sprintf("Method: %s\n", x$method_used))
  cat(sprintf("Missing values before: %d\n", x$n_missing_before))
  cat(sprintf("Missing values after:  %d\n", x$n_missing_after))
  invisible(x)
}
