#' Filter Features by Relative Standard Deviation in QC Samples
#'
#' @description
#' Removes features with high relative standard deviation (RSD, also known as
#' coefficient of variation) in quality control samples. High RSD indicates
#' poor analytical reproducibility.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#'   For `print()`, a `run_filterRSD` object.
#' @param metadata Data frame. Sample metadata with number of rows equal to nrow(x).
#'   Must contain column specified by `group_col` to identify QC samples.
#' @param max_rsd Numeric. Maximum allowed RSD (0-1). Features with RSD >= this
#'   threshold in QC samples are removed. Default: 0.3 (30%).
#' @param group_col Character. Name of column in `metadata` containing sample group labels.
#'   Default: "Group".
#' @param qc_type Character. QC sample type to use for RSD calculation. Options:
#'   \itemize{
#'     \item \code{"EQC"}: Extract QC samples (diluted sample pool)
#'     \item \code{"SQC"}: Sample QC samples (undiluted sample pool)
#'     \item \code{"QC"}: Any QC sample, i.e. every label in \code{qc_types}
#'       (use when both EQC and SQC are present, or there is no distinction)
#'   }
#'   Default: "QC".
#' @param qc_types Character vector. Group labels that should be treated as QC samples.
#'   Default: c("QC", "SQC", "EQC").
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list of class `"run_filterRSD"` containing:
#'   \item{data}{Matrix or data frame of filtered data (same class as input `x`)}
#'   \item{features_removed}{Character vector of removed feature names}
#'   \item{features_kept}{Character vector of retained feature names}
#'   \item{rsd_values}{Numeric vector of RSD values for all original features}
#'   \item{n_features_before}{Integer. Number of features before filtering}
#'   \item{n_features_after}{Integer. Number of features after filtering}
#'   \item{n_features_removed}{Integer. Number of features removed}
#'   \item{parameters}{List of parameters used}
#'
#' @details
#' **Relative Standard Deviation (RSD):**
#'
#' RSD is calculated per feature over the selected QC samples as
#' RSD = SD / Mean (missing values ignored). Features whose RSD cannot be
#' computed (e.g. no observed QC values) are removed, matching
#' [pmp::filter_peaks_by_rsd()].
#'
#' **Apply to non-negative intensities**, i.e. before mean-centering or
#' scaling; RSD of centred data is meaningless. A warning is issued if any QC
#' mean is <= 0.
#'
#' Also known as coefficient of variation (CV), RSD expresses variability
#' as a proportion of the mean, making it scale-independent and suitable
#' for comparing features with different magnitudes.
#'
#' **Why Filter by RSD:**
#'
#' Quality control samples should have low variability since they represent
#' the same biological matrix. High RSD in QC samples indicates:
#' - Poor analytical reproducibility
#' - Instrument instability
#' - Ion suppression/enhancement effects
#' - Features near detection limit
#'
#' **Typical RSD Thresholds:**
#'
#' - Stringent: 20% (0.2) - high-quality data only
#' - Standard: 30% (0.3) - recommended for most metabolomics studies
#' - Lenient: 40% (0.4) - when sample size is limited
#'
#' Broadhurst et al. (2018) recommend removing features with QC RSD above
#' 20% (0.2) for untargeted studies, or 30% (0.3) as a commonly used upper
#' limit.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Broadhurst, D.I. (2025). QC:MXP Repeat Injection based Quality Control, Batch Correction,
#' Exploration & Data Cleaning (Version 2.1). Zenodo. \doi{10.5281/zenodo.16824822}.
#' Retrieved from \url{https://github.com/broadhurstdavid/QC-MXP}.
#'
#' Broadhurst, D., Goodacre, R., Reinke, S.N., Kuligowski, J., Wilson, I.D.,
#' Lewis, M.R., & Dunn, W.B. (2018). Guidelines and considerations for the use
#' of system suitability and quality control samples in mass spectrometry
#' assays applied in untargeted clinical metabolomic studies.
#' *Metabolomics*, 14(6), 72. \doi{10.1007/s11306-018-1367-3}
#'
#' Jankevics, A., Lloyd, G.R., & Weber, R.J.M. (2025). *pmp: Peak Matrix
#' Processing and signal batch correction for metabolomics datasets*.
#' R package. \doi{10.18129/B9.bioc.pmp}
#'
#' @family preprocessing steps
#' @seealso [run_DIpreprocess()] for the full pipeline.
#'
#' @export
#'
#' @examples
#' data(dims_data)
#' data(dims_metadata)
#' raw <- replace(dims_data, dims_data == 0, NA)  # zeros = not detected
#' x <- run_mvimpute(raw, verbose = FALSE)$data
#'
#' res <- run_filterRSD(x, dims_metadata, max_rsd = 0.2)
#' res
#' round(head(res$rsd_values), 3)
run_filterRSD <- function(
    x,
    metadata,
    max_rsd = 0.3,
    group_col = "Group",
    qc_type = "QC",
    qc_types = c("QC", "SQC", "EQC"),
    verbose = TRUE
) {

  msg <- .msg(verbose)

  .validate_data_matrix(x)
  .validate_metadata(x, metadata)

  if (!is.numeric(max_rsd) || length(max_rsd) != 1 || max_rsd < 0 || max_rsd > 1)
    stop("'max_rsd' must be a numeric value between 0 and 1.", call. = FALSE)
  .require_col(metadata, group_col)
  if (!qc_type %in% qc_types)
    stop("'qc_type' must be one of: ", paste(qc_types, collapse = ", "), call. = FALSE)

  mp <- .as_matrix_preserve(x)
  x_matrix <- mp$mat
  feature_names <- colnames(x_matrix)
  n_features_before <- ncol(x_matrix)

  msg(sprintf("Filtering features by RSD (threshold: %.1f%%, QC type: %s)...",
              max_rsd * 100, qc_type))

  # Determine QC label for matching

  if (qc_type == "QC") {
    classes <- ifelse(metadata[[group_col]] %in% qc_types, "QC", metadata[[group_col]])
    qc_label <- "QC"
  } else {
    classes <- metadata[[group_col]]
    qc_label <- qc_type
  }

  # Check if QC samples exist
  if (!qc_label %in% classes) {
    warning(sprintf("No %s samples found. Returning original data without filtering.", qc_label),
            call. = FALSE)
    return(structure(
      list(
        data               = .restore_class(x_matrix, mp$was_matrix),
        features_removed   = character(0),
        features_kept      = feature_names,
        rsd_values         = stats::setNames(rep(NA_real_, n_features_before), feature_names),
        n_features_before  = n_features_before,
        n_features_after   = n_features_before,
        n_features_removed = 0L,
        parameters = list(max_rsd = max_rsd, group_col = group_col,
                          qc_type = qc_type, qc_types = qc_types, qc_label = qc_label)
      ),
      class = "run_filterRSD"
    ))
  }

  # Calculate RSD in QC samples
  qc_indices <- which(classes == qc_label)
  x_qc <- x_matrix[qc_indices, , drop = FALSE]

  qc_means <- matrixStats::colMeans2(x_qc, na.rm = TRUE)
  qc_sds   <- matrixStats::colSds(x_qc, na.rm = TRUE)
  if (any(qc_means <= 0, na.rm = TRUE))
    warning("Some QC feature means are <= 0; RSD is only meaningful for ",
            "non-negative intensities (apply RSD filtering before centering/scaling).",
            call. = FALSE)
  rsd_values <- qc_sds / qc_means
  rsd_values[!is.finite(rsd_values)] <- NA
  names(rsd_values) <- feature_names

  # Same rule as pmp::filter_peaks_by_rsd(): features with undefined RSD are removed
  features_kept <- feature_names[!is.na(rsd_values) & rsd_values < max_rsd]

  x_filtered <- x_matrix[, features_kept, drop = FALSE]
  n_features_after <- length(features_kept)
  n_removed <- n_features_before - n_features_after

  msg(sprintf("Removed %d features (%.1f%%) with RSD >= %.1f%%",
              n_removed, (n_removed / n_features_before) * 100, max_rsd * 100))

  structure(
    list(
      data               = .restore_class(x_filtered, mp$was_matrix),
      features_removed   = setdiff(feature_names, features_kept),
      features_kept      = features_kept,
      rsd_values         = rsd_values,
      n_features_before  = n_features_before,
      n_features_after   = n_features_after,
      n_features_removed = n_removed,
      parameters = list(
        max_rsd   = max_rsd,
        group_col = group_col,
        qc_type   = qc_type,
        qc_types  = qc_types,
        qc_label  = qc_label
      )
    ),
    class = "run_filterRSD"
  )
}

#' @param ... Ignored.
#' @rdname run_filterRSD
#' @export
print.run_filterRSD <- function(x, ...) {
  cat("=== RSD Filtering Results ===\n")
  cat(sprintf("Features before: %d\n", x$n_features_before))
  cat(sprintf("Features after:  %d\n", x$n_features_after))
  cat(sprintf("Features removed: %d (%.1f%%)\n",
              x$n_features_removed,
              (x$n_features_removed / x$n_features_before) * 100))
  cat(sprintf("RSD threshold: %.1f%%\n", x$parameters$max_rsd * 100))
  cat(sprintf("QC type used: %s\n", x$parameters$qc_label))
  invisible(x)
}
