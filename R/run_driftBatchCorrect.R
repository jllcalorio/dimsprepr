#' Correct Signal Drift and Batch Effects Using QC Samples
#'
#' @description
#' Applies Quality Control-based Robust Spline Correction (QCRSC) or batch
#' correction using ComBat to remove systematic signal drift and batch effects
#' in metabolomics data.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#'   For `print()`, a `run_driftBatchCorrect` object.
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
#' @param spline_smooth_param Numeric (0-1). Smoothing parameter (`spar`)
#'   of the QC spline. `0` (default) chooses it per feature by
#'   cross-validation; larger values give smoother fits.
#' @param min_QC Integer. Minimum number of QC samples with observed values a
#'   feature needs in a batch to be corrected; features below it are left
#'   uncorrected. Default: 5.
#' @param spar_limit Numeric vector of length 2. Lower and upper limits of
#'   `spar` searched during cross-validation. Default: c(-1.5, 1.5).
#' @param log_scale Logical. Fit the spline on log-transformed data. Default: TRUE.
#' @param use_parametric Logical. Parametric (`TRUE`) or non-parametric
#'   empirical Bayes priors in ComBat. Default: TRUE.
#' @param display_plots Logical. Draw ComBat's prior diagnostic plots. Default: FALSE.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A list of class `"run_driftBatchCorrect"` containing:
#'   \item{data}{Corrected data (same class as input `x`)}
#'   \item{data_before_correction}{Original data before correction}
#'   \item{correction_applied}{Logical}
#'   \item{method_used}{Character: `"QCRSC"`, `"ComBat"`, or `"None"`}
#'   \item{uncorrected_features}{Character vector of features left unchanged
#'     (e.g. too few QC values)}
#'   \item{n_uncorrected}{Integer}
#'   \item{parameters}{List of parameters}
#'
#' @details
#' **QCRSC** (Quality Control-based Robust Spline Correction; Kirwan et al.,
#' 2013) fits, per feature and batch, a smoothing spline through the QC
#' intensities against injection order and removes that trend from every
#' sample, correcting within-batch drift and between-batch offsets in one
#' step. It follows the QC-based correction approach of Dunn et al. (2011). Implemented by [pmp::QCRSC()]. QC samples
#' must be spread through every batch, starting and ending with a QC, with at
#' least 4 (preferably 5 or more) per batch.
#'
#' **ComBat** (`batch_corr_only = TRUE`; Johnson et al., 2007) uses
#' empirical Bayes location/scale adjustment to remove batch effects without
#' using QCs or injection order. Implemented by [sva::ComBat()]. Requires
#' \pkg{sva}. ComBat does not accept missing values, so impute first.
#'
#' Batch labels may be numeric or character.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Dunn, W.B., Broadhurst, D., Begley, P., Zelena, E., Francis-McIntyre, S.,
#' Anderson, N., et al. (2011). Procedures for large-scale metabolic profiling
#' of serum and plasma using gas chromatography and liquid chromatography
#' coupled to mass spectrometry. *Nature Protocols*, 6(7), 1060-1083.
#' \doi{10.1038/nprot.2011.335}
#'
#' Johnson, W.E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in
#' microarray expression data using empirical Bayes methods.
#' *Biostatistics*, 8(1), 118-127. \doi{10.1093/biostatistics/kxj037}
#'
#' Kirwan, J.A., Broadhurst, D.I., Davidson, R.L., & Viant, M.R. (2013).
#' Characterising and correcting batch variation in an automated direct
#' infusion mass spectrometry (DIMS) metabolomics workflow.
#' *Analytical and Bioanalytical Chemistry*, 405(15), 5147-5157.
#' \doi{10.1007/s00216-013-6856-7}
#'
#' Jankevics, A., Lloyd, G.R., & Weber, R.J.M. (2025). *pmp: Peak Matrix
#' Processing and signal batch correction for metabolomics datasets*.
#' R package. \doi{10.18129/B9.bioc.pmp}
#'
#' @family preprocessing steps
#' @seealso [run_DIpreprocess()] for the full pipeline, [plot_beforeafter()]
#'   to inspect the correction along the injection order.
#'
#' @export
#'
#' @examples
#' data(dims_data)
#' data(dims_metadata)
#' raw <- replace(dims_data, dims_data == 0, NA)  # zeros = not detected
#' x <- run_mvimpute(raw[, 1:10], verbose = FALSE)$data
#'
#' # QC-based drift + batch correction (QCRSC)
#' res <- run_driftBatchCorrect(x, dims_metadata, spline_smooth_param = 0.5)
#' res
#'
#' # ComBat batch correction only (requires sva)
#' \donttest{
#' if (requireNamespace("sva", quietly = TRUE)) {
#'   run_driftBatchCorrect(x, dims_metadata, batch_corr_only = TRUE,
#'                         verbose = FALSE)
#' }
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
  .require_col(metadata, batch_col)
  metadata <- .add_group_collapsed(metadata, group_col, qc_types, qc_label)
  # Batch labels may be numbers or strings ("B1"); both become integer codes.
  batch_vector <- as.numeric(factor(metadata[[batch_col]]))

  uncorrected_features <- character(0)
  method_used <- NULL

  if (!batch_corr_only) {

    # QCRSC correction
    msg("Applying QCRSC drift and batch correction...")

    .require_col(metadata, injection_col)

    if (!is.numeric(spline_smooth_param) || length(spline_smooth_param) != 1 ||
        spline_smooth_param < 0 || spline_smooth_param > 1)
      stop("'spline_smooth_param' must be numeric between 0 and 1.", call. = FALSE)
    if (!is.numeric(spar_limit) || length(spar_limit) != 2)
      stop("'spar_limit' must be a numeric vector of length 2.", call. = FALSE)

    .quiet <- if (verbose) identity else suppressMessages   # pmp messages per call
    .qcrsc <- function(spar) .quiet(pmp::QCRSC(
      df       = t(x_matrix),
      order    = as.numeric(metadata[[injection_col]]),
      batch    = batch_vector,
      classes  = as.vector(metadata$Group_),
      spar     = spar,
      log      = log_scale,
      minQC    = min_QC,
      qc_label = qc_label,
      spar_lim = spar_limit
    ))

    tryCatch({
      x_corrected <- tryCatch(.qcrsc(spline_smooth_param), error = function(e) {
        # ponytail: smooth.spline() CV inside pmp crashes ("object 'outp' not
        # found") when a batch has too few QC values; retry once with a fixed spar.
        if (!grepl("outp", conditionMessage(e), fixed = TRUE) || spline_smooth_param != 0)
          stop(e)
        warning("QCRSC CV failed (likely sparse QC data). Falling back to fixed spar = 0.75.",
                call. = FALSE)
        tryCatch(.qcrsc(0.75), error = function(e2) stop(
          "spline fitting failed even with a fixed spar. A batch probably has fewer ",
          "than 4 QC samples with observed values (smooth.spline() needs >= 4). ",
          "Add QC injections, raise 'min_QC' so such features are skipped, ",
          "or set 'spline_smooth_param' > 0. Original error: ", conditionMessage(e2),
          call. = FALSE))
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

    if (!requireNamespace("sva", quietly = TRUE))
      stop("Package 'sva' required for ComBat. Install with: BiocManager::install('sva')", call. = FALSE)

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
