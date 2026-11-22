#' Filter Background and Contaminant Features Using Procedural Blanks
#'
#' @description
#' Flags features whose signal in study samples is not clearly above the
#' signal in procedural blank samples, removes the flagged features whose
#' blank signal is also unstable, and drops the blank samples. This is the
#' blank inspection step of the QComics guidelines (González-Domínguez et
#' al., 2024).
#'
#' @param x Matrix or data frame of intensities with samples in rows and
#'   features in columns, including the blank samples. Zeros and `NA` are
#'   treated as "not detected".
#'   For `print()`, a `run_filterblank` object.
#' @param metadata Data frame with one row per row of `x`.
#' @param group_col Character. Column in `metadata` with group labels.
#'   Default: `"Group"`.
#' @param blank_types Character vector. Labels in `group_col` that identify
#'   procedural blank samples. Default: `"Blank"`.
#' @param qc_types Character vector. Labels that identify QC samples; they are
#'   not counted as study samples. Default: `c("QC", "SQC", "EQC")`.
#' @param fold_change Numeric > 0. A feature is flagged as a potential
#'   contaminant when its mean in study samples is below `fold_change` times
#'   its mean in blanks. Default: `3` (QComics; mQACC).
#' @param max_blank_rsd Numeric >= 0. Flagged features are removed only if
#'   their relative standard deviation across blanks exceeds this value.
#'   Default: `0.15` (QComics).
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return A list of class `"run_filterblank"` containing:
#'   \item{data}{Data without the blank samples and the removed features
#'     (same class as `x`).}
#'   \item{metadata}{`metadata` without the blank samples.}
#'   \item{features_flagged}{Potential contaminants (mean in samples below
#'     `fold_change` times the mean in blanks).}
#'   \item{features_removed}{Flagged features with blank RSD above
#'     `max_blank_rsd`.}
#'   \item{blank_summary}{Data frame per feature: `Sample_mean`,
#'     `Blank_mean`, `Ratio` (sample mean / blank mean), `Blank_RSD`,
#'     `Flagged`, `Removed`.}
#'   \item{n_blanks}{Number of blank samples removed.}
#'   \item{parameters}{List of parameters used.}
#'
#' @details
#' Means and RSDs are computed over detected (positive, non-missing) values.
#' A feature never detected in the blanks cannot be a blank contaminant.
#' QComics recommends not removing every potential contaminant, because some,
#' such as free fatty acids used as slip agents in plastics, can be
#' biologically relevant. Instead, only those whose blank contribution is
#' variable (blank RSD above 15%) are discarded, since a stable blank
#' contribution does not bias group comparisons. Inject several blanks at the
#' start and at the end of the sequence; at least two are needed here.
#'
#' @author John Lennon L. Calorio
#'
#' @references
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
#' @family preprocessing steps
#' @seealso [run_DIpreprocess()], which runs this step when `blank_types` is
#'   given.
#'
#' @export
#'
#' @examples
#' data(dims_data)
#' data(dims_metadata)
#'
#' # Add three simulated blanks: low background, plus one contaminant
#' blanks <- matrix(0, 3, ncol(dims_data), dimnames = list(
#'   c("B1", "B2", "B3"), colnames(dims_data)))
#' blanks[, "mz_0011"] <- c(4000, 2500, 6000)   # unstable contaminant
#' x <- rbind(dims_data, blanks)
#' meta <- rbind(dims_metadata[, c("Sample", "Group")],
#'               data.frame(Sample = rownames(blanks), Group = "Blank"))
#'
#' res <- run_filterblank(x, meta)
#' res
#' res$features_removed
#' dim(res$data)
run_filterblank <- function(
    x,
    metadata,
    group_col     = "Group",
    blank_types   = "Blank",
    qc_types      = c("QC", "SQC", "EQC"),
    fold_change   = 3,
    max_blank_rsd = 0.15,
    verbose       = TRUE
) {

  msg <- .msg(verbose)

  .validate_data_matrix(x)
  .validate_metadata(x, metadata)
  .require_col(metadata, group_col)
  if (!is.character(blank_types) || length(blank_types) < 1L)
    stop("'blank_types' must be a character vector.", call. = FALSE)
  if (!is.numeric(fold_change) || length(fold_change) != 1L || is.na(fold_change) ||
      fold_change <= 0)
    stop("'fold_change' must be a single positive number.", call. = FALSE)
  if (!is.numeric(max_blank_rsd) || length(max_blank_rsd) != 1L ||
      is.na(max_blank_rsd) || max_blank_rsd < 0)
    stop("'max_blank_rsd' must be a single non-negative number.", call. = FALSE)

  mp <- .as_matrix_preserve(x)
  x_matrix <- mp$mat
  grp <- metadata[[group_col]]
  is_blank  <- grp %in% blank_types
  is_sample <- !is_blank & !grp %in% qc_types
  if (sum(is_blank) < 2L)
    stop(sprintf("At least 2 blank samples (%s) are needed; found %d.",
                 paste(blank_types, collapse = ", "), sum(is_blank)), call. = FALSE)
  if (!any(is_sample))
    stop("No study samples found (all samples are blanks or QCs).", call. = FALSE)

  detected <- replace(x_matrix, is.na(x_matrix) | x_matrix <= 0, NA)
  sample_mean <- colMeans(detected[is_sample, , drop = FALSE], na.rm = TRUE)
  blank_mean  <- colMeans(detected[is_blank, , drop = FALSE], na.rm = TRUE)
  blank_rsd   <- matrixStats::colSds(detected[is_blank, , drop = FALSE], na.rm = TRUE) /
    blank_mean
  sample_mean[is.nan(sample_mean)] <- 0
  blank_mean[is.nan(blank_mean)]   <- 0

  flagged <- blank_mean > 0 & sample_mean < fold_change * blank_mean
  removed <- flagged & !is.na(blank_rsd) & blank_rsd > max_blank_rsd

  msg(sprintf("Blank filter (%d blanks): %d potential contaminant(s) flagged, %d removed (blank RSD > %.0f%%).",
              sum(is_blank), sum(flagged), sum(removed), max_blank_rsd * 100))

  keep_rows <- !is_blank
  structure(
    list(
      data             = .restore_class(x_matrix[keep_rows, !removed, drop = FALSE],
                                        mp$was_matrix),
      metadata         = metadata[keep_rows, , drop = FALSE],
      features_flagged = colnames(x_matrix)[flagged],
      features_removed = colnames(x_matrix)[removed],
      blank_summary    = data.frame(
        Feature = colnames(x_matrix), Sample_mean = sample_mean,
        Blank_mean = blank_mean, Ratio = sample_mean / blank_mean,
        Blank_RSD = blank_rsd, Flagged = flagged, Removed = removed,
        row.names = NULL, stringsAsFactors = FALSE),
      n_blanks         = sum(is_blank),
      parameters       = list(group_col = group_col, blank_types = blank_types,
                              qc_types = qc_types, fold_change = fold_change,
                              max_blank_rsd = max_blank_rsd)
    ),
    class = "run_filterblank"
  )
}

#' @param ... Ignored.
#' @rdname run_filterblank
#' @export
print.run_filterblank <- function(x, ...) {
  cat("=== Blank Filtering Results ===\n")
  cat(sprintf("Blank samples removed: %d\n", x$n_blanks))
  cat(sprintf("Potential contaminants flagged: %d (sample/blank < %g)\n",
              length(x$features_flagged), x$parameters$fold_change))
  cat(sprintf("Features removed (blank RSD > %.0f%%): %d\n",
              x$parameters$max_blank_rsd * 100, length(x$features_removed)))
  invisible(x)
}
