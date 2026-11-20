#' Transform Metabolomics Data
#'
#' @description
#' Applies various transformation methods to stabilize variance and reduce
#' heteroscedasticity in metabolomics data.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#'   For `print()`, a `run_transform` object.
#' @param method Character. Transformation method to apply. Options:
#'   \itemize{
#'     \item \code{"log2"}: Log base 2 transformation
#'     \item \code{"log10"}: Log base 10 transformation
#'     \item \code{"sqrt"}: Square root transformation
#'     \item \code{"cbrt"}: Cube root transformation
#'     \item \code{"clr"}: Centered log-ratio transformation (computed per sample)
#'     \item \code{"arcsin_sqrt"}: Arcsine square root transformation (data must be in \[0, 1\])
#'     \item \code{"vsn"}: Variance Stabilizing Normalization (requires \pkg{vsn})
#'     \item \code{"glog"}: Generalized logarithm transformation (requires \pkg{pmp})
#'     \item \code{"none"}: No transformation
#'   }
#' @param metadata Data frame. Sample metadata, required for glog. Default: NULL.
#' @param group_col Character. Column in `metadata` with group labels.
#'   Default: "Group".
#' @param qc_types Character vector. QC group labels. Default: c("QC", "SQC", "EQC").
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A list of class `"run_transform"` containing:
#'   \item{data}{Transformed data (same class as input `x`)}
#'   \item{method_used}{Character describing method}
#'   \item{shift_applied}{Numeric shift added before transformation}
#'   \item{parameters}{List of parameters}
#'
#' @details
#' Transformations reduce heteroscedasticity and right skew so that fold
#' changes become additive (van den Berg et al., 2006).
#'
#' **log2/log10**: Logarithmic compression. If the data contain non-positive
#' values, a global shift of `abs(min(x)) + 1` is added first (reported in
#' `shift_applied`).
#'
#' **sqrt/cbrt**: Power transformations. A shift of `abs(min(x))` is added for
#' negative values.
#'
#' **clr**: Centered log-ratio (Aitchison, 1982): `log(x)` minus the per-sample
#' mean of `log(x)`. Maps compositional data to real space.
#'
#' **arcsin_sqrt**: For proportion data in \[0, 1\]; errors if the (shifted)
#' maximum exceeds 1.
#'
#' **vsn**: Model-based variance stabilization (Huber et al., 2002) via
#' [vsn::vsn2()]. Requires \pkg{vsn}.
#' Falls back to `log10` with a warning if the fit fails.
#'
#' **glog**: Generalized logarithm (Parsons et al., 2007) via
#' [pmp::glog_transformation()], with the transformation parameter optimised
#' on QC samples. Requires `metadata` with QC labels in `group_col`. Falls back
#' to `log10` with a warning if the fit fails.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Aitchison, J. (1982). The statistical analysis of compositional data.
#' *Journal of the Royal Statistical Society: Series B*, 44(2), 139-177.
#' \doi{10.1111/j.2517-6161.1982.tb01195.x}
#'
#' Huber, W., von Heydebreck, A., Sultmann, H., Poustka, A., & Vingron, M.
#' (2002). Variance stabilization applied to microarray data calibration and
#' to the quantification of differential expression. *Bioinformatics*,
#' 18(Suppl 1), S96-S104. \doi{10.1093/bioinformatics/18.suppl_1.s96}
#'
#' Parsons, H.M., Ludwig, C., Gunther, U.L., & Viant, M.R. (2007). Improved
#' classification accuracy in 1- and 2-dimensional NMR metabolomics data using
#' the variance stabilising generalised logarithm transformation.
#' *BMC Bioinformatics*, 8, 234. \doi{10.1186/1471-2105-8-234}
#'
#' van den Berg, R.A., Hoefsloot, H.C.J., Westerhuis, J.A., Smilde, A.K., &
#' van der Werf, M.J. (2006). Centering, scaling, and transformations:
#' improving the biological information content of metabolomics data.
#' *BMC Genomics*, 7, 142. \doi{10.1186/1471-2164-7-142}
#'
#' Jankevics, A., Lloyd, G.R., & Weber, R.J.M. (2025). *pmp: Peak Matrix
#' Processing and signal batch correction for metabolomics datasets*.
#' R package. \doi{10.18129/B9.bioc.pmp}
#'
#' @family preprocessing steps
#' @seealso [run_DIpreprocess()] for the full pipeline, [run_scale()] for the
#'   step that usually follows.
#'
#' @export
#'
#' @examples
#' data(dims_data)
#' raw <- replace(dims_data, dims_data == 0, NA)  # zeros = not detected
#' x <- run_mvimpute(raw[, 1:20], verbose = FALSE)$data
#'
#' res <- run_transform(x, method = "log10")
#' res
#' range(res$data)
#'
#' # glog uses QC samples, so it needs metadata
#' \donttest{
#' data(dims_metadata)
#' run_transform(x, method = "glog", metadata = dims_metadata, verbose = FALSE)
#' }
run_transform <- function(
    x,
    method    = "log10",
    metadata  = NULL,
    group_col = "Group",
    qc_types  = c("QC", "SQC", "EQC"),
    verbose   = TRUE
) {

  msg <- .msg(verbose)

  .validate_data_matrix(x)

  method <- tolower(method)
  if (length(method) != 1L || !method %in% .transform_methods)
    stop("Unknown transformation method '", paste(method, collapse = ", "),
         "'. Supported: ", paste(.transform_methods, collapse = ", "), ".", call. = FALSE)

  mp <- .as_matrix_preserve(x)
  x_matrix    <- mp$mat
  shift_value <- 0
  method_used <- method

  msg(sprintf("Applying '%s' transformation...", method))

  # Shift data so its minimum is >= 0 (strict = FALSE) or >= 1 (strict = TRUE).
  # ponytail: one global shift for the whole matrix, not per feature.
  .shift <- function(strict) {
    min_val <- min(x_matrix, na.rm = TRUE)
    if (if (strict) min_val <= 0 else min_val < 0) {
      shift_value <<- abs(min_val) + strict
      x_matrix    <<- x_matrix + shift_value
      msg(sprintf("Added shift of %.4g to handle %s values.", shift_value,
                  if (strict) "non-positive" else "negative"))
    }
    x_matrix
  }
  .log10_fallback <- function(e) {
    warning(method, " failed: ", conditionMessage(e), ". Falling back to log10.",
            call. = FALSE)
    method_used <<- sprintf("log10 (fallback from %s)", method)
    log10(.shift(TRUE))
  }

  x_transformed <- switch(
    method,
    "none"  = x_matrix,
    "log2"  = log2(.shift(TRUE)),
    "log10" = log10(.shift(TRUE)),
    "sqrt"  = sqrt(.shift(FALSE)),
    "cbrt"  = .shift(FALSE)^(1 / 3),
    "clr"   = {
      log_x <- log(.shift(TRUE))
      log_x - rowMeans(log_x, na.rm = TRUE)
    },
    "arcsin_sqrt" = {
      xs <- .shift(FALSE)
      if (max(xs, na.rm = TRUE) > 1)
        stop("Arcsine sqrt requires data in [0, 1]. Max after shift is ",
             round(max(xs, na.rm = TRUE), 4), ".", call. = FALSE)
      asin(sqrt(xs))
    },
    "vsn" = {
      if (!requireNamespace("vsn", quietly = TRUE))
        stop("Package 'vsn' required. Install with: BiocManager::install('vsn')", call. = FALSE)
      tryCatch({
        t_matrix <- t(x_matrix)
        t(vsn::predict(vsn::vsn2(t_matrix, verbose = FALSE), newdata = t_matrix))
      }, error = .log10_fallback)
    },
    "glog" = {
      if (is.null(metadata))
        stop("'metadata' is required for glog transformation.", call. = FALSE)
      .validate_metadata(x, metadata)
      .require_col(metadata, group_col)
      metadata <- .add_group_collapsed(metadata, group_col, qc_types, "QC")
      tryCatch(
        t(pmp::glog_transformation(t(x_matrix), classes = metadata$Group_,
                                   qc_label = "QC")),
        error = .log10_fallback
      )
    }
  )

  msg("Transformation complete.")

  structure(
    list(
      data          = .restore_class(x_transformed, mp$was_matrix),
      method_used   = method_used,
      shift_applied = shift_value,
      parameters    = list(
        method    = method,
        group_col = group_col,
        qc_types  = qc_types
      )
    ),
    class = "run_transform"
  )
}

#' @param ... Ignored.
#' @rdname run_transform
#' @export
print.run_transform <- function(x, ...) {
  cat("=== Transformation Results ===\n")
  cat(sprintf("Method: %s\n", x$method_used))
  if (x$shift_applied > 0)
    cat(sprintf("Shift applied: %.4f\n", x$shift_applied))
  cat(sprintf("Data dimensions: %d samples x %d features\n",
              nrow(x$data), ncol(x$data)))
  invisible(x)
}
