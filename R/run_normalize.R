#' Normalize Metabolomics Data
#'
#' @description
#' Applies various normalization methods to account for dilution effects and
#' sample-to-sample variation in metabolomics data.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#'   For `print()`, a `run_normalize` object.
#' @param metadata Data frame. Sample metadata with number of rows equal to nrow(x).
#' @param method Character or numeric vector. Normalization method:
#'   \itemize{
#'     \item \code{"sum"}: Total sum normalization (divide each sample by its row sum)
#'     \item \code{"median"}: Median normalization (divide each sample by its row median)
#'     \item \code{"specific_factor"}: Divide by values from the \code{factor_col}
#'       metadata column (e.g. osmolality or creatinine)
#'     \item \code{"pqn_global"}: Probabilistic quotient normalization (PQN) with the
#'       median spectrum of all samples as reference
#'     \item \code{"pqn_reference"}: PQN with one sample (\code{ref_sample}) as reference
#'     \item \code{"pqn_group"}: PQN with the median spectrum of the
#'       \code{group_sample} samples (pooled QCs by default) as reference
#'     \item \code{"pqn"}: PQN via [pmp::pqn_normalisation()] with the QC
#'       samples as reference
#'     \item \code{"quantile"}: Alias for \code{"pqn"}, kept for backward
#'       compatibility. Note this is \emph{not} quantile normalization.
#'     \item \code{"col_rel_abundance"}: Relative abundance per column (feature)
#'     \item \code{"row_rel_abundance"}: Relative abundance per row (sample)
#'     \item \code{"none"}: No normalization
#'     \item Numeric vector: custom factors, one per non-QC sample
#'   }
#'   When \code{method} is a column name in \code{metadata}, that column is used
#'   as normalization factors (equivalent to \code{"specific_factor"}).
#' @param factor_col Character. Metadata column with normalization factors.
#'   Default: "Normalization". Ignored if \code{method} is a metadata column name.
#' @param ref_sample Character. Sample identifier (in \code{sample_id_col}) used
#'   as reference for \code{"pqn_reference"}. Default: NULL.
#' @param group_sample Character. Group label(s) whose median spectrum is the
#'   reference for \code{"pqn_group"}. Matches either the collapsed QC label
#'   \code{"QC"} (any of \code{qc_types}) or a raw label in \code{group_col}
#'   such as \code{"SQC"}. Default: "QC".
#' @param qc_normalize Character. How QC samples are scaled when sample factors
#'   are supplied (numeric \code{method} or \code{"specific_factor"}), since QCs
#'   usually have no factor of their own: \code{"mean"} or \code{"median"} of the
#'   sample factors, or \code{"none"} (divide by 1). Default: "none".
#' @param group_col Character. Column with group labels. Default: "Group".
#' @param qc_types Character vector. QC group labels. Default: c("QC", "SQC", "EQC").
#' @param reference_method Character. Reference method for PQN via pmp: "mean" or "median".
#'   Default: "mean".
#' @param sample_id_col Character. Column with sample identifiers. Default: "Sample".
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A list of class `"run_normalize"` containing:
#'   \item{data}{Normalized data (same class as input `x`)}
#'   \item{normalization_factors}{Numeric vector with the factor each sample
#'     was divided by}
#'   \item{method_used}{Character describing method}
#'   \item{parameters}{List of parameters}
#'
#' @details
#' Normalization removes unwanted sample-to-sample variation such as
#' dilution, which is especially large in urine. Every method divides each
#' sample (row) by one factor.
#'
#' PQN (Dieterle et al., 2006) computes, for each sample, the median of the
#' feature-wise quotients between that sample and a reference spectrum. It is
#' more robust than sum normalization when a few features change strongly.
#' The \code{"pqn"} / \code{"quantile"} options delegate to
#' [pmp::pqn_normalisation()] (Jankevics et al., 2025), falling back to
#' \code{"sum"} with a warning if it fails.
#'
#' Missing values are ignored when computing factors. Any \code{Inf} created
#' by a zero factor is replaced with \code{NA}.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Dieterle, F., Ross, A., Schlotterbeck, G., & Senn, H. (2006). Probabilistic
#' quotient normalization as robust method to account for dilution of complex
#' biological mixtures. Application in 1H NMR metabonomics.
#' *Analytical Chemistry*, 78(13), 4281-4290. \doi{10.1021/ac051632c}
#'
#' Jankevics, A., Lloyd, G.R., & Weber, R.J.M. (2025). *pmp: Peak Matrix
#' Processing and signal batch correction for metabolomics datasets*.
#' R package. \doi{10.18129/B9.bioc.pmp}
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
    method         = "median",
    factor_col     = "Normalization",
    ref_sample     = NULL,
    group_sample   = "QC",
    qc_normalize   = "none",
    group_col      = "Group",
    qc_types       = c("QC", "SQC", "EQC"),
    reference_method = "mean",
    sample_id_col  = "Sample",
    verbose        = TRUE
) {

  msg <- .msg(verbose)

  .validate_data_matrix(x)
  .validate_metadata(x, metadata)

  mp <- .as_matrix_preserve(x)
  x_matrix <- mp$mat
  qc_normalize <- match.arg(qc_normalize, c("none", "mean", "median"))

  # Check if method is a column name in metadata (auto-detect specific_factor)
  if (is.character(method) && length(method) == 1L && method %in% colnames(metadata)) {
    factor_col <- method
    method <- "specific_factor"
  }
  if (is.character(method) &&
      (length(method) != 1L || !tolower(method) %in% .normalize_methods))
    stop("Unknown normalization method '", paste(method, collapse = ", "),
         "'. Supported: ", paste(.normalize_methods, collapse = ", "),
         ", a numeric vector of factors, or a metadata column name.", call. = FALSE)
  if (!is.character(method) && !is.numeric(method))
    stop("'method' must be a character string or a numeric vector.", call. = FALSE)

  # Collapse QC types
  if (group_col %in% colnames(metadata)) {
    metadata <- .add_group_collapsed(metadata, group_col, qc_types, "QC")
  } else if (is.character(method) && tolower(method) %in% c("pqn_group", "pqn", "quantile")) {
    stop(sprintf("'%s' column required in metadata for PQN methods.", group_col), call. = FALSE)
  } else {
    metadata$Group_ <- "Sample"
  }

  qc_indices     <- metadata$Group_ == "QC"
  non_qc_indices <- !qc_indices

  # Initialize all_factors before switch to avoid scope issues
  all_factors <- rep(NA_real_, nrow(x_matrix))

  msg(sprintf("Applying '%s' normalization...",
              if (is.numeric(method)) "custom factors" else method))

  if (is.numeric(method)) {

    n_non_qc <- sum(non_qc_indices)
    if (length(method) != n_non_qc)
      stop(sprintf("Custom factors must have length %d (non-QC samples), got %d.",
                   n_non_qc, length(method)), call. = FALSE)

    norm_factors <- method
    x_matrix[non_qc_indices, ] <- x_matrix[non_qc_indices, ] / norm_factors

    qc_factor <- switch(qc_normalize,
                        "mean"   = mean(norm_factors, na.rm = TRUE),
                        "median" = stats::median(norm_factors, na.rm = TRUE),
                        1)
    x_matrix[qc_indices, ] <- x_matrix[qc_indices, ] / qc_factor

    method_description <- "Custom normalization factors"
    all_factors[non_qc_indices] <- norm_factors
    all_factors[qc_indices] <- qc_factor

  } else {

    method_lc <- tolower(method)
    # Alias: "quantile" -> "pqn" (both use pmp::pqn_normalisation)
    if (method_lc == "quantile") method_lc <- "pqn"

    # ponytail: avoid switch()-as-expr scope issues with <<- inside branches.
    # Set x_matrix and all_factors directly in each branch.
    if (method_lc == "none") {
      msg("No normalization applied.")
      all_factors <- rep(1, nrow(x_matrix))

    } else if (method_lc == "sum") {
      msg("Normalizing by total sum...")
      row_sums <- matrixStats::rowSums2(x_matrix, na.rm = TRUE)
      all_factors <- row_sums
      x_matrix <- x_matrix / row_sums

    } else if (method_lc == "median") {
      msg("Normalizing by median...")
      row_medians <- matrixStats::rowMedians(x_matrix, na.rm = TRUE)
      all_factors <- row_medians
      x_matrix <- x_matrix / row_medians

    } else if (method_lc == "specific_factor") {
      .require_col(metadata, factor_col)
      factor_values <- as.numeric(metadata[[factor_col]][non_qc_indices])

      if (all(is.na(factor_values) | factor_values == 0)) {
        warning("All normalization factors are NA or 0. Using 'sum' instead.", call. = FALSE)
        row_sums <- matrixStats::rowSums2(x_matrix, na.rm = TRUE)
        all_factors <- row_sums
        x_matrix <- x_matrix / row_sums
      } else {
        msg(sprintf("Normalizing using factors from '%s'...", factor_col))
        x_matrix[non_qc_indices, ] <- x_matrix[non_qc_indices, ] / factor_values

        qc_factor <- switch(qc_normalize,
                            "mean"   = mean(factor_values, na.rm = TRUE),
                            "median" = stats::median(factor_values, na.rm = TRUE),
                            1)
        x_matrix[qc_indices, ] <- x_matrix[qc_indices, ] / qc_factor

        all_factors[non_qc_indices] <- factor_values
        all_factors[qc_indices] <- qc_factor
      }

    } else if (method_lc == "pqn_global") {
      msg("Applying PQN (global median reference)...")
      reference_spectrum <- matrixStats::colMedians(x_matrix, na.rm = TRUE)
      quotients <- sweep(x_matrix, 2, reference_spectrum, "/")
      median_quotients <- matrixStats::rowMedians(quotients, na.rm = TRUE)
      all_factors <- median_quotients
      x_matrix <- x_matrix / median_quotients

    } else if (method_lc == "pqn_reference") {
      if (is.null(ref_sample))
        stop("'ref_sample' must be specified for pqn_reference.", call. = FALSE)
      .require_col(metadata, sample_id_col)
      ref_idx <- which(metadata[[sample_id_col]] == ref_sample)
      if (length(ref_idx) == 0)
        stop(sprintf("Reference sample '%s' not found.", ref_sample), call. = FALSE)

      msg(sprintf("Applying PQN (reference: %s)...", ref_sample))
      reference_spectrum <- as.numeric(x_matrix[ref_idx[1], ])
      quotients <- sweep(x_matrix, 2, reference_spectrum, "/")
      median_quotients <- matrixStats::rowMedians(quotients, na.rm = TRUE)
      all_factors <- median_quotients
      x_matrix <- x_matrix / median_quotients

    } else if (method_lc == "pqn_group") {
      # group_sample may be the collapsed "QC" label or one raw label (e.g. "SQC")
      pooled_indices <- metadata$Group_ %in% group_sample |
        metadata[[group_col]] %in% group_sample
      if (sum(pooled_indices) == 0)
        stop(sprintf("No '%s' samples found in '%s' for group PQN.",
                     paste(group_sample, collapse = ", "), group_col), call. = FALSE)

      msg(sprintf("Applying PQN (group: %s)...", group_sample))
      reference_spectrum <- matrixStats::colMedians(
        x_matrix[pooled_indices, , drop = FALSE], na.rm = TRUE
      )
      quotients <- sweep(x_matrix, 2, reference_spectrum, "/")
      median_quotients <- matrixStats::rowMedians(quotients, na.rm = TRUE)
      all_factors <- median_quotients
      x_matrix <- x_matrix / median_quotients

    } else if (method_lc == "pqn") {
      if (!requireNamespace("pmp", quietly = TRUE))
        stop("Package 'pmp' required. Install with: BiocManager::install('pmp')", call. = FALSE)

      tryCatch({
        msg("Applying PQN normalization via pmp...")
        normalized <- pmp::pqn_normalisation(
          t(x_matrix),
          classes    = metadata$Group_,
          qc_label   = "QC",
          ref_mean   = NULL,
          qc_frac    = 0,
          sample_frac = 0,
          ref_method = reference_method
        )
        all_factors <- matrixStats::rowMedians(x_matrix / t(normalized), na.rm = TRUE)
        x_matrix <- t(normalized)
      }, error = function(e) {
        warning("PQN normalization failed: ", e$message, ". Using 'sum' instead.", call. = FALSE)
        row_sums <- matrixStats::rowSums2(x_matrix, na.rm = TRUE)
        all_factors <<- row_sums
        x_matrix <<- x_matrix / row_sums
      })

    } else if (method_lc == "col_rel_abundance") {
      msg("Normalizing by column relative abundance...")
      col_sums <- matrixStats::colSums2(x_matrix, na.rm = TRUE)
      if (any(col_sums == 0))
        warning("Columns sum to 0; relative abundance will produce NaN.", call. = FALSE)
      all_factors <- col_sums
      x_matrix <- sweep(x_matrix, 2, col_sums, "/")

    } else if (method_lc == "row_rel_abundance") {
      msg("Normalizing by row relative abundance...")
      row_sums <- matrixStats::rowSums2(x_matrix, na.rm = TRUE)
      if (any(row_sums == 0))
        warning("Rows sum to 0; relative abundance will produce NaN.", call. = FALSE)
      all_factors <- row_sums
      x_matrix <- x_matrix / row_sums

    }

    method_description <- if (method_lc == "pqn" && tolower(method) == "quantile") {
      "quantile (pqn via pmp)"
    } else {
      method_lc
    }
  }

  # Handle Inf values
  if (any(is.infinite(x_matrix))) {
    msg("Warning: Normalization produced Inf values. Replacing with NA.")
    x_matrix[is.infinite(x_matrix)] <- NA
  }

  msg("Normalization complete.")

  structure(
    list(
      data                  = .restore_class(x_matrix, mp$was_matrix),
      normalization_factors = all_factors,
      method_used           = method_description,
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
