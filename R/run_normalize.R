#' Normalize Metabolomics Data
#'
#' @description
#' Divides each sample by a normalization factor to remove dilution effects
#' and other unwanted sample-to-sample variation in overall concentration.
#'
#' @param x Matrix or data frame. Non-negative intensities with samples in
#'   rows and features in columns. Missing values are ignored when computing
#'   the factors.
#'   For `print()`, a `run_normalize` object.
#' @param metadata Data frame. Sample metadata with `nrow(metadata) == nrow(x)`.
#' @param method Character or numeric vector. Normalization method:
#'   \itemize{
#'     \item `"sum"`: total sum (divide each sample by its row sum)
#'     \item `"median"`: divide each sample by its median
#'     \item `"specific_factor"`: divide by the values in the `factor_col`
#'       metadata column (e.g. osmolality or creatinine)
#'     \item `"pqn_global"`: probabilistic quotient normalization (PQN) with
#'       the reference spectrum computed from all samples
#'     \item `"pqn_reference"`: PQN with one sample (`ref_sample`) as reference
#'     \item `"pqn_group"`: PQN with the reference spectrum computed from the
#'       `group_sample` samples (pooled QCs by default)
#'     \item `"none"`: no normalization
#'     \item a numeric vector: custom factors, one per non-QC sample
#'   }
#'   When `method` is a column name in `metadata`, that column supplies the
#'   factors (equivalent to `"specific_factor"`).
#' @param factor_col Character. Metadata column with normalization factors for
#'   `"specific_factor"`. Default: `"Normalization"`.
#' @param ref_sample Character. Sample identifier (in `sample_id_col`) used as
#'   reference for `"pqn_reference"`. Default: `NULL`.
#' @param group_sample Character. Group label(s) whose samples form the
#'   reference for `"pqn_group"`. Matches either the collapsed QC label `"QC"`
#'   (any of `qc_types`) or a raw label in `group_col` such as `"SQC"`.
#'   Default: `"QC"`.
#' @param qc_normalize Character. How QC samples are scaled when sample factors
#'   are supplied (numeric `method` or `"specific_factor"`), since QCs usually
#'   have no factor of their own: `"mean"` or `"median"` of the sample factors,
#'   or `"none"` (divide by 1). Default: `"none"`.
#' @param group_col Character. Column with group labels. Default: `"Group"`.
#' @param qc_types Character vector. QC group labels.
#'   Default: `c("QC", "SQC", "EQC")`.
#' @param reference_method Character. How the PQN reference spectrum is formed
#'   from its samples: `"median"` (default; Dieterle et al., 2006) or
#'   `"mean"`.
#' @param sample_id_col Character. Column with sample identifiers, used by
#'   `"pqn_reference"`. Default: `"Sample"`.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return A list of class `"run_normalize"` containing:
#'   \item{data}{Normalized data (same class as `x`).}
#'   \item{normalization_factors}{Numeric vector with the factor each sample
#'     was divided by.}
#'   \item{method_used}{Character describing the method.}
#'   \item{parameters}{List of parameters used.}
#'
#' @details
#' Normalization removes unwanted variation in overall sample concentration,
#' such as dilution, which is especially large in urine. Every method divides
#' each sample (row) by one factor.
#'
#' PQN (Dieterle et al., 2006) computes, for each sample, the median of the
#' feature-wise quotients between that sample and a reference spectrum. It is
#' more robust than sum normalization when a few features change strongly.
#' Quotients that are not finite (zero reference) are ignored.
#'
#' Missing values are ignored when computing factors. Any `Inf` created by a
#' zero factor is replaced with `NA`.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Dieterle, F., Ross, A., Schlotterbeck, G., & Senn, H. (2006). Probabilistic
#' quotient normalization as robust method to account for dilution of complex
#' biological mixtures. Application in 1H NMR metabonomics.
#' *Analytical Chemistry*, 78(13), 4281-4290. \doi{10.1021/ac051632c}
#'
#' @family preprocessing steps
#' @seealso [run_DIpreprocess()] for the full pipeline, [run_transform()] for
#'   the step that usually follows.
#'
#' @export
#'
#' @examples
#' data(dims_data)
#' data(dims_metadata)
#' raw <- replace(dims_data, dims_data == 0, NA)  # zeros = not detected
#' x <- run_mvimpute(raw[, 1:20], verbose = FALSE)$data
#'
#' res <- run_normalize(x, dims_metadata, method = "pqn_group")
#' res
#' head(res$normalization_factors)
#'
#' # Sample-specific factors (e.g. osmolality) from a metadata column
#' meta <- dims_metadata
#' meta$Osmolality <- ifelse(meta$Group == "QC", NA, runif(nrow(meta), 0.8, 1.2))
#' run_normalize(x, meta, method = "Osmolality", qc_normalize = "median",
#'               verbose = FALSE)
run_normalize <- function(
    x,
    metadata,
    method           = "median",
    factor_col       = "Normalization",
    ref_sample       = NULL,
    group_sample     = "QC",
    qc_normalize     = "none",
    group_col        = "Group",
    qc_types         = c("QC", "SQC", "EQC"),
    reference_method = "median",
    sample_id_col    = "Sample",
    verbose          = TRUE
) {

  msg <- .msg(verbose)

  .validate_data_matrix(x)
  .validate_metadata(x, metadata)
  qc_normalize     <- match.arg(qc_normalize, c("none", "mean", "median"))
  reference_method <- match.arg(reference_method, c("median", "mean"))

  mp <- .as_matrix_preserve(x)
  x_matrix <- mp$mat

  # A metadata column name as method means "use that column as factors"
  if (is.character(method) && length(method) == 1L && method %in% colnames(metadata)) {
    factor_col <- method
    method <- "specific_factor"
  }
  if (!is.character(method) && !is.numeric(method))
    stop("'method' must be a character string or a numeric vector.", call. = FALSE)
  if (is.character(method) &&
      (length(method) != 1L || !tolower(method) %in% .normalize_methods))
    stop("Unknown normalization method '", paste(method, collapse = ", "),
         "'. Supported: ", paste(.normalize_methods, collapse = ", "),
         ", a numeric vector of factors, or a metadata column name.", call. = FALSE)
  method_lc <- if (is.character(method)) tolower(method) else "custom"

  if (group_col %in% colnames(metadata)) {
    metadata <- .add_group_collapsed(metadata, group_col, qc_types, "QC")
  } else if (method_lc == "pqn_group") {
    stop(sprintf("'%s' column required in metadata for pqn_group.", group_col), call. = FALSE)
  } else {
    metadata$Group_ <- "Sample"
  }
  qc <- metadata$Group_ == "QC"

  msg(sprintf("Applying '%s' normalization...", if (method_lc == "custom") "custom factors" else method_lc))

  # PQN: per-sample median quotient to a reference spectrum
  .pqn <- function(ref_rows) {
    ref_x <- x_matrix[ref_rows, , drop = FALSE]
    ref <- if (reference_method == "median") matrixStats::colMedians(ref_x, na.rm = TRUE) else
      colMeans(ref_x, na.rm = TRUE)
    q <- sweep(x_matrix, 2, ref, "/")
    q[!is.finite(q)] <- NA
    matrixStats::rowMedians(q, na.rm = TRUE)
  }
  # Sample factors for non-QC rows; QC rows get mean/median/1 per qc_normalize
  .with_qc_factor <- function(f_non_qc) {
    f <- numeric(nrow(x_matrix))
    f[!qc] <- f_non_qc
    f[qc] <- switch(qc_normalize,
                    mean   = mean(f_non_qc, na.rm = TRUE),
                    median = stats::median(f_non_qc, na.rm = TRUE),
                    none   = 1)
    f
  }

  factors <- switch(
    method_lc,
    custom = {
      if (length(method) != sum(!qc))
        stop(sprintf("Custom factors must have length %d (non-QC samples), got %d.",
                     sum(!qc), length(method)), call. = FALSE)
      .with_qc_factor(method)
    },
    none   = rep(1, nrow(x_matrix)),
    sum    = matrixStats::rowSums2(x_matrix, na.rm = TRUE),
    median = matrixStats::rowMedians(x_matrix, na.rm = TRUE),
    specific_factor = {
      .require_col(metadata, factor_col)
      f <- as.numeric(metadata[[factor_col]][!qc])
      if (all(is.na(f) | f == 0)) {
        warning("All normalization factors are NA or 0. Using 'sum' instead.", call. = FALSE)
        method_lc <- "sum (fallback from specific_factor)"
        matrixStats::rowSums2(x_matrix, na.rm = TRUE)
      } else {
        .with_qc_factor(f)
      }
    },
    pqn_global = .pqn(seq_len(nrow(x_matrix))),
    pqn_reference = {
      if (is.null(ref_sample))
        stop("'ref_sample' must be specified for pqn_reference.", call. = FALSE)
      .require_col(metadata, sample_id_col)
      ref_idx <- which(metadata[[sample_id_col]] == ref_sample)
      if (length(ref_idx) == 0L)
        stop(sprintf("Reference sample '%s' not found.", ref_sample), call. = FALSE)
      .pqn(ref_idx[1])
    },
    pqn_group = {
      ref_rows <- metadata$Group_ %in% group_sample | metadata[[group_col]] %in% group_sample
      if (!any(ref_rows))
        stop(sprintf("No '%s' samples found in '%s' for group PQN.",
                     paste(group_sample, collapse = ", "), group_col), call. = FALSE)
      .pqn(ref_rows)
    }
  )

  x_matrix <- x_matrix / factors
  if (any(is.infinite(x_matrix))) {
    msg("Warning: Normalization produced Inf values. Replacing with NA.")
    x_matrix[is.infinite(x_matrix)] <- NA
  }

  msg("Normalization complete.")

  structure(
    list(
      data                  = .restore_class(x_matrix, mp$was_matrix),
      normalization_factors = factors,
      method_used           = method_lc,
      parameters = list(
        method           = if (is.numeric(method)) "custom" else method,
        factor_col       = factor_col,
        ref_sample       = ref_sample,
        group_sample     = group_sample,
        qc_normalize     = qc_normalize,
        group_col        = group_col,
        qc_types         = qc_types,
        reference_method = reference_method
      )
    ),
    class = "run_normalize"
  )
}

#' @param ... Ignored.
#' @rdname run_normalize
#' @export
print.run_normalize <- function(x, ...) {
  cat("=== Normalization Results ===\n")
  cat(sprintf("Method: %s\n", x$method_used))
  cat(sprintf("Data dimensions: %d samples x %d features\n",
              nrow(x$data), ncol(x$data)))
  invisible(x)
}
