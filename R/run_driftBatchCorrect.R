#' Correct Signal Drift and Batch Effects Using QC Samples
#'
#' @description
#' Applies Quality Control-based Robust Spline Correction (QCRSC) or batch correction
#' using ComBat to remove systematic signal drift and batch effects in metabolomics data.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#' @param metadata Data frame. Sample metadata with number of rows equal to nrow(x).
#' @param perform_correction Logical. If TRUE, perform correction; if FALSE, return
#'   original data unchanged. Default: TRUE.
#' @param batch_corr_only Logical. If TRUE, perform only ComBat batch correction
#'   (no drift correction). If FALSE, perform QCRSC. Default: FALSE.
#' @param injection_col Character. Column in `metadata` containing injection order.
#'   Default: "InjectionSequence".
#' @param batch_col Character. Column in `metadata` containing batch numbers.
#'   Default: "Batch".
#' @param group_col Character. Column in `metadata` containing sample group labels.
#'   Default: "Group".
#' @param qc_label Character. Group label identifying QC samples for QCRSC. Default: "QC".
#' @param qc_types Character vector. Group labels converted to `qc_label` internally.
#'   Default: c("QC", "SQC", "EQC").
#' @param spline_smooth_param Numeric. Smoothing parameter (0-1). Default: 0.
#' @param min_QC Integer. Minimum QC samples per batch. Default: 5.
#' @param spar_limit Numeric vector of length 2. Spline limits. Default: c(-1.5, 1.5).
#' @param log_scale Logical. Fit on log-transformed data. Default: TRUE.
#' @param use_parametric Logical. Parametric ComBat adjustment. Default: TRUE.
#' @param display_plots Logical. Diagnostic plots for ComBat. Default: FALSE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A list of class `"run_driftBatchCorrect"` containing:
#'   \item{data}{Corrected data (same class as input `x`)}
#'   \item{data_before_correction}{Original data before correction}
#'   \item{correction_applied}{Logical}
#'   \item{method_used}{Character}
#'   \item{uncorrected_features}{Character vector of uncorrected features}
#'   \item{n_uncorrected}{Integer}
#'   \item{parameters}{List of parameters}
#'
#' @details
#' **QCRSC** fits smoothing splines through QC intensities across injection order.
#' Requires \pkg{pmp}.
#'
#' **ComBat** uses empirical Bayes for batch correction. Requires \pkg{sva}.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Kirwan, J.A., et al. (2013). Analytical and Bioanalytical Chemistry, 405, 5147-5157.
#' \doi{10.1007/s00216-013-6856-7}
#'
#' Johnson, W.E., et al. (2007). Biostatistics, 8(1), 118-127.
#' \doi{10.1093/biostatistics/kxj037}
#'
#' @seealso \code{\link{run_DIpreprocess}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(417)
#' n <- 60; p <- 30
#' x <- matrix(abs(rnorm(n * p, 100, 20)), nrow = n, ncol = p)
#' colnames(x) <- paste0("F", seq_len(p))
#'
#' meta <- data.frame(
#'   InjectionSequence = seq_len(n),
#'   Batch = rep(1:2, each = 30),
#'   Group = rep(c(rep("Sample", 9), "QC"), 6)
#' )
#'
#' # Skip correction (returns unchanged)
#' result <- run_driftBatchCorrect(x, meta, perform_correction = FALSE)
#' result
#' }
run_driftBatchCorrect <- function(
    x,
    metadata,
    perform_correction  = TRUE,
    batch_corr_only     = FALSE,
    injection_col       = "InjectionSequence",
    batch_col           = "Batch",
    group_col           = "Group",
    qc_label            = "QC",
    qc_types            = c("QC", "SQC", "EQC"),
    spline_smooth_param = 0,
    min_QC              = 5,
    spar_limit          = c(-1.5, 1.5),
    log_scale           = TRUE,
    use_parametric      = TRUE,
    display_plots       = FALSE,
    verbose             = TRUE
) {

  msg <- .msg(verbose)

  .validate_data_matrix(x)
  .validate_metadata(x, metadata)

  if (!perform_correction) {
    msg("Correction disabled. Returning original data.")
    return(structure(
      list(data = x, data_before_correction = x, correction_applied = FALSE,
           method_used = "None", uncorrected_features = character(0),
           n_uncorrected = 0L, parameters = list(perform_correction = FALSE)),
      class = "run_driftBatchCorrect"
    ))
  }

  mp <- .as_matrix_preserve(x)
  x_matrix <- mp$mat
  x_before <- x_matrix

  .require_col(metadata, group_col)
  metadata <- .add_group_collapsed(metadata, group_col, qc_types, qc_label)

  uncorrected_features <- character(0)
  method_used <- NULL

  if (!batch_corr_only) {

    # QCRSC correction
    msg("Applying QCRSC drift and batch correction...")

    for (col_name in c(injection_col, batch_col, group_col))
      .require_col(metadata, col_name)

    if (!requireNamespace("pmp", quietly = TRUE))
      stop("Package 'pmp' required for QCRSC. Install with: BiocManager::install('pmp')", call. = FALSE)

    if (!is.numeric(spline_smooth_param) || length(spline_smooth_param) != 1 ||
        spline_smooth_param < 0 || spline_smooth_param > 1)
      stop("'spline_smooth_param' must be numeric between 0 and 1.", call. = FALSE)
    if (!is.numeric(spar_limit) || length(spar_limit) != 2)
      stop("'spar_limit' must be a numeric vector of length 2.", call. = FALSE)

    tryCatch({
      # Inner tryCatch to handle sparse-data fallback on the fly
      x_corrected <- tryCatch({
        pmp::QCRSC(
          df       = t(x_matrix),
          order    = as.numeric(metadata[[injection_col]]),
          batch    = as.numeric(metadata[[batch_col]]),
          classes  = as.vector(metadata$Group_),
          spar     = spline_smooth_param,
          log      = log_scale,
          minQC    = min_QC,
          qc_label = qc_label,
          spar_lim = spar_limit
        )
      }, error = function(e) {
        # ponytail: smooth.spline CV crashes when QCs are too sparse (e.g. impute_method="none").
        # Catch the known bug and fallback to a fixed spar instead of dying.
        if (grepl("outp", e$message, fixed = TRUE) && spline_smooth_param == 0) {
          warning("QCRSC CV failed (likely sparse QC data). Falling back to fixed spar = 0.75.", call. = FALSE)
          pmp::QCRSC(
            df       = t(x_matrix),
            order    = as.numeric(metadata[[injection_col]]),
            batch    = as.numeric(metadata[[batch_col]]),
            classes  = as.vector(metadata$Group_),
            spar     = 0.75,
            log      = log_scale,
            minQC    = min_QC,
            qc_label = qc_label,
            spar_lim = spar_limit
          )
        } else {
          stop(e)
        }
      })
      
      x_matrix <- t(x_corrected)
      method_used <- "QCRSC"

      # Identify uncorrected features
      tryCatch({
        if (identical(dim(x_before), dim(x_matrix)) &&
            identical(colnames(x_before), colnames(x_matrix))) {
          tolerance <- .Machine$double.eps^0.5
          diff_matrix <- abs(x_before - x_matrix)
          col_max_diffs <- matrixStats::colMaxs(diff_matrix, na.rm = TRUE)
          uncorrected_indices <- (col_max_diffs <= tolerance) | is.na(col_max_diffs)
          uncorrected_features <- colnames(x_before)[uncorrected_indices]
          if (length(uncorrected_features) > 0)
            msg(sprintf("Warning: %d features could not be corrected", length(uncorrected_features)))
        }
      }, error = function(e) {
        warning("Could not identify uncorrected features: ", e$message, call. = FALSE)
      })

      msg("QCRSC correction completed successfully.")
    }, error = function(e) {
      stop("QCRSC correction failed: ", e$message, call. = FALSE)
    })

  } else {

    # ComBat batch correction
    msg("Applying ComBat batch correction...")
    .require_col(metadata, batch_col)

    if (!requireNamespace("sva", quietly = TRUE))
      stop("Package 'sva' required for ComBat. Install with: BiocManager::install('sva')", call. = FALSE)

    batch_vector <- as.numeric(metadata[[batch_col]])
    if (length(unique(batch_vector)) < 2) {
      warning("Only one batch detected. Returning original data.", call. = FALSE)
      method_used <- "None (single batch)"
    } else {
      tryCatch({
        x_corrected <- sva::ComBat(
          dat = t(x_matrix), batch = batch_vector, mod = NULL,
          par.prior = use_parametric, prior.plots = display_plots,
          mean.only = FALSE, ref.batch = NULL
        )
        x_matrix <- t(x_corrected)
        method_used <- "ComBat"
        msg("ComBat correction completed successfully.")
      }, error = function(e) {
        stop("ComBat correction failed: ", e$message, call. = FALSE)
      })
    }
  }

  new_na_count <- sum(is.na(x_matrix)) - sum(is.na(x_before))
  if (new_na_count > 0)
    msg(sprintf("Warning: Correction introduced %d new missing values", new_na_count))

  structure(
    list(
      data                   = .restore_class(x_matrix, mp$was_matrix),
      data_before_correction = .restore_class(x_before, mp$was_matrix),
      correction_applied     = TRUE,
      method_used            = method_used,
      uncorrected_features   = uncorrected_features,
      n_uncorrected          = length(uncorrected_features),
      parameters = list(
        perform_correction  = perform_correction, batch_corr_only = batch_corr_only,
        injection_col = injection_col, batch_col = batch_col, group_col = group_col,
        qc_label = qc_label, qc_types = qc_types,
        spline_smooth_param = spline_smooth_param, min_QC = min_QC,
        spar_limit = spar_limit, log_scale = log_scale,
        use_parametric = use_parametric, display_plots = display_plots
      )
    ),
    class = "run_driftBatchCorrect"
  )
}

#' @param x Object to print.
#' @param ... Ignored.
#' @rdname run_driftBatchCorrect
#' @export
print.run_driftBatchCorrect <- function(x, ...) {
  cat("=== Drift/Batch Correction Results ===\n")
  cat(sprintf("Correction applied: %s\n", x$correction_applied))
  cat(sprintf("Method used: %s\n", x$method_used))
  if (x$n_uncorrected > 0)
    cat(sprintf("Uncorrected features: %d\n", x$n_uncorrected))
  invisible(x)
}
