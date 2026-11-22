#' Filter Features by Relative Standard Deviation in QC Samples
#'
#' @description
#' Removes features with a high relative standard deviation (RSD, also known
#' as coefficient of variation) across quality control samples. A high QC
#' RSD indicates poor analytical reproducibility.
#'
#' @param x Matrix or data frame. Non-negative intensities with samples in
#'   rows and features in columns.
#'   For `print()`, a `run_filterRSD` object.
#' @param metadata Data frame. Sample metadata with `nrow(metadata) == nrow(x)`.
#' @param max_rsd Numeric in \[0, 1\]. Maximum allowed RSD. Features with a QC
#'   RSD of `max_rsd` or more are removed. Default: `0.3` (30%).
#' @param group_col Character. Column in `metadata` with group labels.
#'   Default: `"Group"`.
#' @param qc_types Character vector. Labels in `group_col` of the QC samples
#'   used to compute the RSD. Give a single label (e.g. `"EQC"`) to use only
#'   that QC type. Default: `c("QC", "SQC", "EQC")`.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return A list of class `"run_filterRSD"` containing:
#'   \item{data}{Filtered data (same class as `x`).}
#'   \item{features_removed, features_kept}{Character vectors of feature names.}
#'   \item{rsd_values}{Named numeric vector of QC RSDs for all input features.}
#'   \item{n_features_before, n_features_after, n_features_removed}{Integers.}
#'   \item{parameters}{List of parameters used.}
#'
#' @details
#' RSD is calculated per feature over the selected QC samples as SD / mean,
#' ignoring missing values. Being scale-independent, it allows features of
#' different magnitudes to be compared. Features whose RSD cannot be computed
#' (e.g. no observed QC values) are removed, as in
#' [pmp::filter_peaks_by_rsd()].
#'
#' Pooled QC samples share one biological matrix, so their variability is
#' technical: poor reproducibility, instrument instability, ion suppression,
#' or proximity to the detection limit.
#'
#' **Apply to measured, non-negative intensities**: before mean-centering or
#' scaling, since the RSD of centred data is meaningless (a warning is issued
#' if any QC mean is <= 0), and without imputed values (set them to `NA`),
#' since near-minimum imputed values among the QCs inflate the RSD.
#' [run_DIpreprocess()] does both.
#'
#' **Thresholds.** Broadhurst et al. (2018) recommend removing features with
#' QC RSD above 20% for untargeted studies, with 30% as a common upper limit;
#' QComics uses RSD < 30% as the acceptance criterion for peak intensities
#' (González-Domínguez et al., 2024). Vynck et al. (2026) caution that a low
#' QC RSD is a weak proxy for good drift correction, as it rewards
#' overfitting.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Broadhurst, D.I. (2025). QC:MXP Repeat Injection based Quality Control,
#' Batch Correction, Exploration & Data Cleaning (Version 2.1). Zenodo.
#' \doi{10.5281/zenodo.16824822}. Retrieved from
#' \url{https://github.com/broadhurstdavid/QC-MXP}.
#'
#' Broadhurst, D., Goodacre, R., Reinke, S.N., Kuligowski, J., Wilson, I.D.,
#' Lewis, M.R., & Dunn, W.B. (2018). Guidelines and considerations for the use
#' of system suitability and quality control samples in mass spectrometry
#' assays applied in untargeted clinical metabolomic studies.
#' *Metabolomics*, 14(6), 72. \doi{10.1007/s11306-018-1367-3}
#'
#' González-Domínguez, A., Estanyol-Torres, N., Brunius, C., Landberg, R., &
#' González-Domínguez, R. (2024). QComics: Recommendations and guidelines for
#' robust, easily implementable and reportable quality control of
#' metabolomics data. *Analytical Chemistry*, 96(3), 1064-1072.
#' \doi{10.1021/acs.analchem.3c03660}
#'
#' Jankevics, A., Lloyd, G.R., & Weber, R.J.M. (2025). *pmp: Peak Matrix
#' Processing and signal batch correction for metabolomics datasets*.
#' R package. \doi{10.18129/B9.bioc.pmp}
#'
#' Vynck, M., Vangeenderhuysen, P., De Paepe, E., Nawrot, T., Plekhova, V., &
#' Vanhaecke, L. (2026). Robust metabolomics data normalization across scales
#' and experimental designs. *Analytical Chemistry*, 98(24), 17627-17637.
#' \doi{10.1021/acs.analchem.5c06841}
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
    max_rsd   = 0.3,
    group_col = "Group",
    qc_types  = c("QC", "SQC", "EQC"),
    verbose   = TRUE
) {

  msg <- .msg(verbose)

  .validate_data_matrix(x)
  .validate_metadata(x, metadata)
  if (!is.numeric(max_rsd) || length(max_rsd) != 1 || is.na(max_rsd) ||
      max_rsd < 0 || max_rsd > 1)
    stop("'max_rsd' must be a numeric value between 0 and 1.", call. = FALSE)
  .require_col(metadata, group_col)
  if (!is.character(qc_types) || length(qc_types) < 1L)
    stop("'qc_types' must be a character vector of QC labels.", call. = FALSE)

  mp <- .as_matrix_preserve(x)
  x_matrix <- mp$mat
  feature_names <- colnames(x_matrix)
  n_features_before <- ncol(x_matrix)
  qc_label <- paste(qc_types, collapse = "/")

  msg(sprintf("Filtering features by RSD (threshold: %.1f%%, QC: %s)...",
              max_rsd * 100, qc_label))

  params <- list(max_rsd = max_rsd, group_col = group_col, qc_types = qc_types)
  qc_rows <- which(metadata[[group_col]] %in% qc_types)

  if (length(qc_rows) == 0L) {
    warning(sprintf("No QC samples (%s) found. Returning original data without filtering.",
                    qc_label), call. = FALSE)
    return(structure(
      list(
        data               = .restore_class(x_matrix, mp$was_matrix),
        features_removed   = character(0),
        features_kept      = feature_names,
        rsd_values         = stats::setNames(rep(NA_real_, n_features_before), feature_names),
        n_features_before  = n_features_before,
        n_features_after   = n_features_before,
        n_features_removed = 0L,
        parameters         = params
      ),
      class = "run_filterRSD"
    ))
  }

  x_qc <- x_matrix[qc_rows, , drop = FALSE]
  qc_means <- matrixStats::colMeans2(x_qc, na.rm = TRUE)
  if (any(qc_means <= 0, na.rm = TRUE))
    warning("Some QC feature means are <= 0; RSD is only meaningful for ",
            "non-negative intensities (apply RSD filtering before centering/scaling).",
            call. = FALSE)
  rsd_values <- matrixStats::colSds(x_qc, na.rm = TRUE) / qc_means
  rsd_values[!is.finite(rsd_values)] <- NA
  names(rsd_values) <- feature_names

  # Same rule as pmp::filter_peaks_by_rsd(): features with undefined RSD are removed
  keep <- !is.na(rsd_values) & rsd_values < max_rsd
  n_removed <- sum(!keep)

  msg(sprintf("Removed %d features (%.1f%%) with RSD >= %.1f%%",
              n_removed, (n_removed / n_features_before) * 100, max_rsd * 100))

  structure(
    list(
      data               = .restore_class(x_matrix[, keep, drop = FALSE], mp$was_matrix),
      features_removed   = feature_names[!keep],
      features_kept      = feature_names[keep],
      rsd_values         = rsd_values,
      n_features_before  = n_features_before,
      n_features_after   = sum(keep),
      n_features_removed = n_removed,
      parameters         = params
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
  cat(sprintf("QC samples used: %s\n", paste(x$parameters$qc_types, collapse = ", ")))
  invisible(x)
}
