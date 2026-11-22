#' Filter Features by Missing Value Threshold
#'
#' @description
#' Removes features (columns) whose proportion of missing values exceeds a
#' threshold, assessed globally or within each sample group, optionally
#' including QC samples. Zeros can be treated as missing, as is usual for
#' below-detection values in metabolomics. Optionally, features that are
#' mostly missing in a group are treated as truly absent there and kept with
#' zeros, following the QComics guidelines.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and
#'   features in columns. Missing values should be `NA`; zeros can be treated
#'   as missing via `zero_as_missing`.
#'   For `print()`, a `run_filtermissing` object.
#' @param metadata Data frame. Sample metadata with `nrow(metadata) == nrow(x)`.
#'   Must contain `group_col` when `filter_by_group = TRUE` or when QC samples
#'   are to be excluded.
#' @param threshold Numeric in \[0, 1\]. Maximum proportion of missing values
#'   allowed. Features with a missing proportion above `threshold` are
#'   removed. Default: `0.2` (the 80% rule).
#' @param filter_by_group Logical. If `TRUE` (default), assess missingness
#'   within each group separately.
#' @param any_group Logical. Only used when `filter_by_group = TRUE`. If `TRUE`
#'   (default), a feature is kept if it passes `threshold` in at least one
#'   group (the "modified 80% rule"; Yang et al., 2015). If `FALSE`, it must
#'   pass in all groups.
#' @param include_QC Logical. If `TRUE`, include QC samples in the
#'   missingness assessment. Default: `FALSE`.
#' @param group_col Character. Column in `metadata` with group labels.
#'   Default: `"Group"`.
#' @param qc_types Character vector. Group labels that identify QC samples.
#'   Default: `c("QC", "SQC", "EQC")`.
#' @param zero_as_missing Logical. If `TRUE` (default), treat zeros as missing
#'   (`NA`) before filtering.
#' @param absent_threshold Numeric in (`threshold`, 1), or `NULL` (default,
#'   off). Within each assessed group, a feature missing in more than this
#'   proportion of samples is considered \emph{truly absent} there: its
#'   missing values in that group are set to 0 instead of being left for
#'   imputation, and the feature is kept even if it fails `threshold`.
#'   QComics uses `0.7` (González-Domínguez et al., 2024).
#' @param sample_iqr Numeric or `NULL` (default, off). QComics sample check:
#'   within each group (QC samples included), a sample whose proportion of
#'   missing values over the commonly detected features (those passing the
#'   feature filter) exceeds the group's third quartile by more than
#'   `sample_iqr` times the interquartile range (taken as at least 0.05) is
#'   removed as a likely failed or low-signal injection, and the feature
#'   filter is recomputed without it. `3` is Tukey's far-out rule.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return A list of class `"run_filtermissing"` containing:
#'   \item{data}{Filtered data (same class as `x`), without the samples in
#'     `samples_removed`.}
#'   \item{metadata}{`metadata` without the samples in `samples_removed`.}
#'   \item{samples_removed}{Samples removed by the `sample_iqr` check.}
#'   \item{features_removed, features_kept}{Character vectors of feature names.}
#'   \item{features_absent}{Features kept as truly absent in at least one group
#'     (empty unless `absent_threshold` is set).}
#'   \item{n_features_before, n_features_after, n_features_removed}{Integers.}
#'   \item{n_zeros_converted}{Number of zeros converted to `NA`.}
#'   \item{n_absent_zeros}{Number of missing values set to 0 as truly absent.}
#'   \item{sample_missing_rate}{Proportion of missing values per input sample
#'     over the commonly detected features, the quantity checked by
#'     `sample_iqr`.}
#'   \item{missingness_summary}{Data frame with the per-feature missing
#'     proportion in each assessed group (one `Missingness` column in global
#'     mode), and logical `Absent` and `Kept` columns.}
#'   \item{parameters}{List of parameters used.}
#'
#' @details
#' **Zeros.** In metabolomics data, zeros usually mean that a peak was below
#' the detection limit, failed integration, or was not detected. With
#' `zero_as_missing = TRUE` (default) they are converted to `NA` first. Set it
#' to `FALSE` only if zeros were already converted, carry quantitative
#' meaning, or the data contain negative values.
#'
#' **Filtering modes.**
#' * Global (`filter_by_group = FALSE`): missingness over all assessed
#'   samples. With `threshold = 0.2` this is the 80% rule of Smilde et al.
#'   (2005): a feature is kept if it is observed in at least 80% of samples.
#' * Group-wise (`filter_by_group = TRUE`): missingness within each group.
#'   With `any_group = TRUE` this is the modified 80% rule (Yang et al.,
#'   2015), which keeps group-specific features; QComics recommends the same
#'   rule (remove features with more than 20% missing in all groups).
#'
#' **QC samples.** With `include_QC = FALSE`, samples labelled with
#' `qc_types` are excluded from the assessment, as their missingness pattern
#' differs from that of biological samples.
#'
#' **Truly absent data (QComics).** Plain filtering removes features that are
#' genuinely present in only some individuals, such as drugs or dietary
#' compounds. Following González-Domínguez et al. (2024), when
#' `absent_threshold` is set, a feature missing in more than that proportion
#' of a group is taken to be absent in that group (real zeros) rather than
#' missing for analytical reasons. Downstream, prefer a transformation that
#' handles zeros (e.g. `"glog"`), since a log transform then needs a shift.
#'
#' **Samples with low detection rates (QComics).** QComics also recommends
#' checking the detection rate of every sample, by group, and discarding
#' samples with abnormally low detection rates, such as failed injections.
#' Their few remaining values are unreliable, and normalization would
#' inflate them. `sample_iqr` applies this check, using the far-out rule that
#' QComics applies to quality markers.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Broadhurst, D.I. (2025). QC:MXP Repeat Injection based Quality Control,
#' Batch Correction, Exploration & Data Cleaning (Version 2.1). Zenodo.
#' \doi{10.5281/zenodo.16824822}. Retrieved from
#' \url{https://github.com/broadhurstdavid/QC-MXP}.
#'
#' González-Domínguez, A., Estanyol-Torres, N., Brunius, C., Landberg, R., &
#' González-Domínguez, R. (2024). QComics: Recommendations and guidelines for
#' robust, easily implementable and reportable quality control of
#' metabolomics data. *Analytical Chemistry*, 96(3), 1064-1072.
#' \doi{10.1021/acs.analchem.3c03660}
#'
#' Smilde, A.K., van der Werf, M.J., Bijlsma, S., van der Werff-van der Vat,
#' B.J., & Jellema, R.H. (2005). Fusion of mass spectrometry-based metabolomics
#' data. *Analytical Chemistry*, 77(20), 6729-6736. \doi{10.1021/ac051080y}
#'
#' Wei, R., Wang, J., Su, M., Jia, E., Chen, S., Chen, T., & Ni, Y. (2018).
#' Missing value imputation approach for mass spectrometry-based metabolomics
#' data. *Scientific Reports*, 8(1), 663. \doi{10.1038/s41598-017-19120-0}
#'
#' Yang, J., Zhao, X., Lu, X., Lin, X., & Xu, G. (2015). A data preprocessing
#' strategy for metabolomics to reduce the mask effect in data analysis.
#' *Frontiers in Molecular Biosciences*, 2, 4. \doi{10.3389/fmolb.2015.00004}
#'
#' @family preprocessing steps
#' @seealso [run_DIpreprocess()] for the full pipeline, [run_mvimpute()] for
#'   the step that usually follows.
#'
#' @export
#'
#' @examples
#' data(dims_data)
#' data(dims_metadata)
#'
#' # Modified 80% rule: keep features observed in >= 80% of any one group
#' res <- run_filtermissing(dims_data, dims_metadata, threshold = 0.2)
#' res
#' res$features_removed
#'
#' # Classic 80% rule across all non-QC samples
#' run_filtermissing(dims_data, dims_metadata, filter_by_group = FALSE,
#'                   verbose = FALSE)
#'
#' # QComics: keep mostly-absent features as zeros instead of removing them
#' run_filtermissing(dims_data, dims_metadata, absent_threshold = 0.7,
#'                   verbose = FALSE)$features_absent
run_filtermissing <- function(
    x,
    metadata,
    threshold        = 0.2,
    filter_by_group  = TRUE,
    any_group        = TRUE,
    include_QC       = FALSE,
    group_col        = "Group",
    qc_types         = c("QC", "SQC", "EQC"),
    zero_as_missing  = TRUE,
    absent_threshold = NULL,
    sample_iqr       = NULL,
    verbose          = TRUE
) {

  msg <- .msg(verbose)

  .validate_data_matrix(x)
  .validate_metadata(x, metadata)

  if (!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) ||
      threshold < 0 || threshold > 1)
    stop("'threshold' must be a numeric value between 0 and 1.", call. = FALSE)
  if (!is.null(absent_threshold) &&
      (!is.numeric(absent_threshold) || length(absent_threshold) != 1 ||
       is.na(absent_threshold) || absent_threshold <= threshold || absent_threshold >= 1))
    stop("'absent_threshold' must be NULL or a number between 'threshold' and 1.",
         call. = FALSE)
  if (!is.null(sample_iqr) &&
      (!is.numeric(sample_iqr) || length(sample_iqr) != 1 || is.na(sample_iqr) ||
       sample_iqr <= 0))
    stop("'sample_iqr' must be NULL or a single positive number.", call. = FALSE)
  if (filter_by_group) .require_col(metadata, group_col)
  if (!is.logical(zero_as_missing) || length(zero_as_missing) != 1)
    stop("'zero_as_missing' must be a single logical value.", call. = FALSE)

  mp <- .as_matrix_preserve(x)
  x_matrix <- mp$mat
  n_features_before <- ncol(x_matrix)
  feature_names <- colnames(x_matrix)

  msg(sprintf("Starting missing value filtering (threshold: %.1f%%)...", threshold * 100))

  n_zeros_converted <- 0L
  if (zero_as_missing) {
    n_zeros_converted <- sum(x_matrix == 0, na.rm = TRUE)
    if (n_zeros_converted > 0L) {
      x_matrix[x_matrix == 0] <- NA
      msg(sprintf("Converted %d zero values to NA (zero_as_missing = TRUE)", n_zeros_converted))
    }
  }
  sample_names <- if (is.null(rownames(x_matrix))) as.character(seq_len(nrow(x_matrix))) else
    rownames(x_matrix)

  # Missing proportions per assessed group, and the feature rule, for the
  # current x_matrix / metadata (global mode = one group).
  .assess <- function() {
    # QC exclusion needs group labels; without group_col all rows count
    a <- if (include_QC || !group_col %in% colnames(metadata)) seq_len(nrow(x_matrix)) else
      which(!metadata[[group_col]] %in% qc_types)
    if (length(a) == 0L)
      stop("No samples available for missingness assessment after QC exclusion.", call. = FALSE)
    g  <- if (filter_by_group) as.character(metadata[[group_col]][a]) else rep("Missingness", length(a))
    gs <- unique(g)
    miss <- matrix(
      vapply(gs, function(k) colMeans(is.na(x_matrix[a[g %in% k], , drop = FALSE])),
             numeric(ncol(x_matrix))),
      nrow = ncol(x_matrix), dimnames = list(feature_names, gs)
    )
    passes <- miss <= threshold
    list(assess = a, grp = g, groups = gs, miss = miss,
         keep = if (any_group) rowSums(passes) > 0L else rowSums(passes) == ncol(passes))
  }
  ev <- .assess()

  # QComics sample check: abnormally high missingness over the commonly
  # detected features (those passing the feature rule), within each group.
  # simplified: Tukey far-out rule with the IQR floored at 0.05, so an almost
  # complete group does not flag a sample for a single missing value.
  sample_missing_rate <- stats::setNames(
    rowMeans(is.na(x_matrix[, ev$keep, drop = FALSE])), sample_names)
  samples_removed <- character(0)
  if (!is.null(sample_iqr)) {
    sg <- if (group_col %in% colnames(metadata)) as.character(metadata[[group_col]]) else
      rep("all", nrow(x_matrix))
    flagged <- logical(nrow(x_matrix))
    for (g in unique(sg)) {
      i <- which(sg %in% g)
      q <- stats::quantile(sample_missing_rate[i], c(0.25, 0.75), names = FALSE, na.rm = TRUE)
      flagged[i] <- sample_missing_rate[i] > q[2] + sample_iqr * max(q[2] - q[1], 0.05)
    }
    if (any(flagged)) {
      samples_removed <- sample_names[flagged]
      x_matrix <- x_matrix[!flagged, , drop = FALSE]
      metadata <- metadata[!flagged, , drop = FALSE]
      msg(sprintf("Removed %d sample(s) with abnormally high missingness: %s",
                  sum(flagged), paste(samples_removed, collapse = ", ")))
      ev <- .assess()
    }
  }
  assess <- ev$assess
  grp    <- ev$grp
  groups <- ev$groups
  missing_by_group <- ev$miss
  keep   <- ev$keep

  # QComics: mostly-missing (feature, group) cells are truly absent -> zeros
  absent <- missing_by_group > (if (is.null(absent_threshold)) 1 else absent_threshold)
  n_absent_zeros <- 0L
  for (g in groups[colSums(absent) > 0]) {
    rows <- assess[grp %in% g]
    cols <- which(absent[, g])
    na_cells <- is.na(x_matrix[rows, cols, drop = FALSE])
    n_absent_zeros <- n_absent_zeros + sum(na_cells)
    x_matrix[rows, cols][na_cells] <- 0
  }
  is_absent <- rowSums(absent) > 0L
  keep <- keep | is_absent

  missingness_summary <- data.frame(Feature = feature_names, missing_by_group,
                                    Absent = is_absent, Kept = keep,
                                    check.names = FALSE, row.names = NULL)

  x_filtered <- x_matrix[, keep, drop = FALSE]
  n_removed <- n_features_before - ncol(x_filtered)

  msg(sprintf("Removed %d features (%.1f%%) with missingness > %.1f%%",
              n_removed, (n_removed / n_features_before) * 100, threshold * 100))
  if (any(is_absent))
    msg(sprintf("Kept %d truly absent feature(s); %d missing value(s) set to 0.",
                sum(is_absent), n_absent_zeros))

  structure(
    list(
      data                = .restore_class(x_filtered, mp$was_matrix),
      features_removed    = feature_names[!keep],
      features_kept       = feature_names[keep],
      features_absent     = feature_names[is_absent],
      n_features_before   = n_features_before,
      n_features_after    = ncol(x_filtered),
      n_features_removed  = n_removed,
      n_zeros_converted   = n_zeros_converted,
      n_absent_zeros      = n_absent_zeros,
      sample_missing_rate = sample_missing_rate,
      samples_removed     = samples_removed,
      metadata            = metadata,
      missingness_summary = missingness_summary,
      parameters = list(
        threshold        = threshold,
        filter_by_group  = filter_by_group,
        any_group        = any_group,
        include_QC       = include_QC,
        group_col        = group_col,
        qc_types         = qc_types,
        zero_as_missing  = zero_as_missing,
        absent_threshold = absent_threshold,
        sample_iqr       = sample_iqr
      )
    ),
    class = "run_filtermissing"
  )
}


#' @param ... Ignored.
#' @rdname run_filtermissing
#' @export
print.run_filtermissing <- function(x, ...) {
  cat("=== Missing Value Filtering Results ===\n")
  cat(sprintf("Features before: %d\n", x$n_features_before))
  cat(sprintf("Features after:  %d\n", x$n_features_after))
  cat(sprintf("Features removed: %d (%.1f%%)\n",
              x$n_features_removed,
              (x$n_features_removed / x$n_features_before) * 100))
  cat(sprintf("Threshold: %.1f%%\n", x$parameters$threshold * 100))
  if (x$parameters$filter_by_group)
    cat(sprintf("Mode: %s group(s) must pass\n",
                if (x$parameters$any_group) "Any" else "All"))
  if (x$parameters$zero_as_missing && x$n_zeros_converted > 0)
    cat(sprintf("Zeros converted to NA: %d\n", x$n_zeros_converted))
  if (length(x$samples_removed) > 0L)
    cat(sprintf("Samples removed (abnormal missingness): %d\n", length(x$samples_removed)))
  if (length(x$features_absent) > 0L)
    cat(sprintf("Truly absent features kept: %d (%d values set to 0)\n",
                length(x$features_absent), x$n_absent_zeros))
  invisible(x)
}
