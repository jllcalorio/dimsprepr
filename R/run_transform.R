#' Transform Metabolomics Data
#'
#' @description
#' Applies various transformation methods to stabilize variance and reduce
#' heteroscedasticity in metabolomics data.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
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
#'   }
#' @param metadata Data frame. Sample metadata, required for glog. Default: NULL.
#' @param group_col Character. Column in `metadata` with group labels.
#'   Default: "Group".
#' @param qc_types Character vector. QC group labels. Default: c("QC", "SQC", "EQC").
#' @param num_cores Integer or "max". Cores for VSN. Default: 1.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A list of class `"run_transform"` containing:
#'   \item{data}{Transformed data (same class as input `x`)}
#'   \item{method_used}{Character describing method}
#'   \item{shift_applied}{Numeric shift added before transformation}
#'   \item{parameters}{List of parameters}
#'
#' @details
#' **log2/log10**: Logarithmic compression. A shift is added for non-positive values.
#'
#' **sqrt/cbrt**: Power transformations. A shift is added for negative values.
#'
#' **clr**: Centered log-ratio. Maps compositional data to real space.
#'
#' **arcsin_sqrt**: For proportion data in \[0, 1\].
#'
#' **vsn**: Model-based variance stabilization. Requires \pkg{vsn} and \pkg{BiocParallel}.
#'
#' **glog**: Generalized logarithm via \pkg{pmp}. Requires metadata with QC labels.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Huber, W., et al. (2002). Bioinformatics, 18(Suppl 1), S96-S104.
#' \doi{10.1093/bioinformatics/18.suppl_1.s96}
#'
#' Parsons, H.M., et al. (2007). BMC Bioinformatics, 8, 234.
#' \doi{10.1186/1471-2105-8-234}
#'
#' @seealso \code{\link{run_DIpreprocess}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(937)
#' x <- matrix(abs(rnorm(80 * 40, mean = 100, sd = 50)), nrow = 80, ncol = 40)
#' colnames(x) <- paste0("Feature", 1:40)
#'
#' result <- run_transform(x, method = "log2")
#' result
#' }
run_transform <- function(
    x,
    method    = "log2",
    metadata  = NULL,
    group_col = "Group",
    qc_types  = c("QC", "SQC", "EQC"),
    num_cores = 1,
    verbose   = TRUE
) {

  msg <- .msg(verbose)

  .validate_data_matrix(x)

  method <- tolower(method)
  valid_methods <- c("log2", "log10", "sqrt", "cbrt", "clr", "arcsin_sqrt", "vsn", "glog")
  if (!method %in% valid_methods)
    stop("Unknown transformation method '", method, "'. Supported: ",
         paste(valid_methods, collapse = ", "), ".", call. = FALSE)

  mp <- .as_matrix_preserve(x)
  x_matrix    <- mp$mat
  shift_value <- 0

  msg(sprintf("Applying '%s' transformation...", method))

  # VSN core validation
  if (method == "vsn") {
    if (!requireNamespace("parallelly", quietly = TRUE))
      stop("Package 'parallelly' required for VSN core detection.", call. = FALSE)
    max_available_cores <- parallelly::availableCores(omit = 2)
    if (identical(num_cores, "max")) {
      num_cores <- max_available_cores
    } else if (is.numeric(num_cores)) {
      num_cores <- as.integer(num_cores)
      if (num_cores < 1L || num_cores > max_available_cores)
        stop(sprintf("'num_cores' must be between 1 and %d.", max_available_cores), call. = FALSE)
    } else {
      stop("'num_cores' must be \"max\" or a positive integer.", call. = FALSE)
    }
  }

  x_transformed <- switch(
    method,

    "log2" = {
      min_val <- min(x_matrix, na.rm = TRUE)
      if (min_val <= 0) {
        shift_value <- abs(min_val) + 1
        x_matrix <- x_matrix + shift_value
        msg(sprintf("Added shift of %.2f to handle non-positive values.", shift_value))
      }
      log2(x_matrix)
    },

    "log10" = {
      min_val <- min(x_matrix, na.rm = TRUE)
      if (min_val <= 0) {
        shift_value <- abs(min_val) + 1
        x_matrix <- x_matrix + shift_value
        msg(sprintf("Added shift of %.2f to handle non-positive values.", shift_value))
      }
      log10(x_matrix)
    },

    "sqrt" = {
      min_val <- min(x_matrix, na.rm = TRUE)
      if (min_val < 0) {
        shift_value <- abs(min_val)
        x_matrix <- x_matrix + shift_value
        msg(sprintf("Added shift of %.2f to handle negative values.", shift_value))
      }
      sqrt(x_matrix)
    },

    "cbrt" = {
      min_val <- min(x_matrix, na.rm = TRUE)
      if (min_val < 0) {
        shift_value <- abs(min_val)
        x_matrix <- x_matrix + shift_value
        msg(sprintf("Added shift of %.2f to handle negative values.", shift_value))
      }
      x_matrix ^ (1 / 3)
    },

    "clr" = {
      min_val <- min(x_matrix, na.rm = TRUE)
      if (min_val <= 0) {
        shift_value <- abs(min_val) + 1
        x_matrix <- x_matrix + shift_value
        msg(sprintf("Added shift of %.2f to handle non-positive values.", shift_value))
      }
      log_x <- log(x_matrix)
      log_x - rowMeans(log_x, na.rm = TRUE)
    },

    "arcsin_sqrt" = {
      min_val <- min(x_matrix, na.rm = TRUE)
      if (min_val < 0) {
        shift_value <- abs(min_val)
        x_matrix <- x_matrix + shift_value
        msg(sprintf("Added shift of %.2f to handle negative values.", shift_value))
      }
      max_val <- max(x_matrix, na.rm = TRUE)
      if (max_val > 1)
        stop("Arcsine sqrt requires data in [0, 1]. Max after shift is ",
             round(max_val, 4), ".", call. = FALSE)
      asin(sqrt(x_matrix))
    },

    "vsn" = {
      if (!requireNamespace("vsn", quietly = TRUE))
        stop("Package 'vsn' required. Install with: BiocManager::install('vsn')", call. = FALSE)
      if (!requireNamespace("BiocParallel", quietly = TRUE))
        stop("Package 'BiocParallel' required. Install with: BiocManager::install('BiocParallel')", call. = FALSE)

      tryCatch({
        BPPARAM_to_use <- BiocParallel::SnowParam(workers = num_cores)
        BiocParallel::register(BPPARAM_to_use)
        msg(sprintf("Applying VSN transformation (%d core(s))...", num_cores))

        t_matrix <- t(x_matrix)
        vsn_fit <- vsn::vsn2(t_matrix)
        result  <- vsn::predict(vsn_fit, newdata = t_matrix)
        BiocParallel::register(BiocParallel::SerialParam())
        t(result)
      }, error = function(e) {
        tryCatch(BiocParallel::register(BiocParallel::SerialParam()), error = function(e2) NULL)
        warning("VSN failed: ", e$message, ". Falling back to log10.", call. = FALSE)
        min_val <- min(x_matrix, na.rm = TRUE)
        if (min_val <= 0) {
          shift_value <<- abs(min_val) + 1
          x_matrix <<- x_matrix + shift_value
        }
        log10(x_matrix)
      })
    },

    "glog" = {
      if (!requireNamespace("pmp", quietly = TRUE))
        stop("Package 'pmp' required. Install with: BiocManager::install('pmp')", call. = FALSE)
      if (is.null(metadata))
        stop("'metadata' is required for glog transformation.", call. = FALSE)
      .validate_metadata(x, metadata)
      .require_col(metadata, group_col)

      metadata <- .add_group_collapsed(metadata, group_col, qc_types, "QC")

      tryCatch({
        msg("Applying glog transformation...")
        glog_result <- pmp::glog_transformation(
          t(x_matrix), classes = metadata$Group_, qc_label = "QC"
        )
        t(glog_result)
      }, error = function(e) {
        warning("glog failed: ", e$message, ". Falling back to log10.", call. = FALSE)
        min_val <- min(x_matrix, na.rm = TRUE)
        if (min_val <= 0) {
          shift_value <<- abs(min_val) + 1
          x_matrix <<- x_matrix + shift_value
        }
        log10(x_matrix)
      })
    }
  )

  msg("Transformation complete.")

  structure(
    list(
      data          = .restore_class(x_transformed, mp$was_matrix),
      method_used   = method,
      shift_applied = shift_value,
      parameters    = list(
        method    = method,
        num_cores = if (method == "vsn") num_cores else NA,
        group_col = group_col,
        qc_types  = qc_types
      )
    ),
    class = "run_transform"
  )
}

#' @param x Object to print.
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
