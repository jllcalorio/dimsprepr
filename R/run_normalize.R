#' Normalize Metabolomics Data
#'
#' @description
#' Applies various normalization methods to account for dilution effects and
#' sample-to-sample variation in metabolomics data.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#' @param metadata Data frame. Sample metadata with number of rows equal to nrow(x).
#' @param method Character or numeric vector. Normalization method:
#'   \itemize{
#'     \item \code{"sum"}: Total sum normalization
#'     \item \code{"median"}: Median normalization
#'     \item \code{"specific_factor"}: Use values from a metadata column
#'     \item \code{"pqn_global"}: PQN using global median
#'     \item \code{"pqn_reference"}: PQN using a specific reference sample
#'     \item \code{"pqn_group"}: PQN using pooled QC samples as reference
#'     \item \code{"pqn"}: PQN via the pmp package (requires \pkg{pmp})
#'     \item \code{"quantile"}: Alias for \code{"pqn"} (PQN via pmp)
#'     \item \code{"col_rel_abundance"}: Relative abundance per column (feature)
#'     \item \code{"row_rel_abundance"}: Relative abundance per row (sample)
#'     \item \code{"none"}: No normalization
#'     \item Numeric vector: Custom normalization factors
#'   }
#'   When \code{method} is a column name in \code{metadata}, that column is used
#'   as normalization factors (equivalent to \code{"specific_factor"}).
#' @param factor_col Character. Metadata column with normalization factors.
#'   Default: "Normalization". Ignored if \code{method} is a metadata column name.
#' @param ref_sample Character. Reference sample for PQN. Default: NULL.
#' @param group_sample Character. Group label for PQN group reference.
#'   Default: "QC".
#' @param qc_normalize Character. QC normalization strategy: "mean", "median", "none".
#'   Default: "none".
#' @param group_col Character. Column with group labels. Default: "Group".
#' @param qc_types Character vector. QC group labels. Default: c("QC", "SQC", "EQC").
#' @param reference_method Character. Reference method for PQN via pmp: "mean" or "median".
#'   Default: "mean".
#' @param sample_id_col Character. Column with sample identifiers. Default: "Sample".
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A list of class `"run_normalize"` containing:
#'   \item{data}{Normalized data (same class as input `x`)}
#'   \item{normalization_factors}{Numeric vector of factors used}
#'   \item{method_used}{Character describing method}
#'   \item{parameters}{List of parameters}
#'
#' @details
#' See the package documentation for full descriptions of each method.
#' \code{"pqn"} and \code{"quantile"} both call \code{pmp::pqn_normalisation()}.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Dieterle, F., et al. (2006). Analytical Chemistry, 78(13), 4281-4290.
#' https://doi.org/10.1021/ac051632c
#'
#' @seealso \code{\link{run_DIpreprocess}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(215)
#' x <- matrix(abs(rnorm(80 * 40, mean = 100, sd = 20)), nrow = 80, ncol = 40)
#' colnames(x) <- paste0("Feature", 1:40)
#' meta <- data.frame(
#'   Sample = paste0("S", 1:80),
#'   Group = rep(c("Control", "Treatment", "QC"), c(30, 30, 20))
#' )
#'
#' result <- run_normalize(x, meta, method = "median")
#' result
#' }
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

  # Check if method is a column name in metadata (auto-detect specific_factor)
  if (is.character(method) && length(method) == 1L && method %in% colnames(metadata)) {
    factor_col <- method
    method <- "specific_factor"
  }

  # Collapse QC types
  if (group_col %in% colnames(metadata)) {
    metadata <- .add_group_collapsed(metadata, group_col, qc_types, "QC")
  } else if (is.character(method) && grepl("pqn", tolower(method))) {
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
      pooled_indices <- metadata$Group_ == "QC"
      if (sum(pooled_indices) == 0)
        stop("No QC samples found for group PQN.", call. = FALSE)

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

    } else {
      stop("Unknown normalization method '", method, "'.", call. = FALSE)
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

#' @param x Object to print.
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
