#' Scale Metabolomics Data
#'
#' @description
#' Centres and scales features to prepare data for multivariate analysis,
#' with classical scaling methods and the variable stability (VAST) family.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and
#'   features in columns.
#'   For `print()`, a `run_scale` object.
#' @param method Character. Scaling method:
#'   \itemize{
#'     \item `"mean"`: mean-centering only
#'     \item `"auto"`: auto-scaling (mean-centering + unit variance)
#'     \item `"pareto"`: Pareto scaling (mean-centering, divide by sqrt(SD))
#'     \item `"vast"`: VAST, auto-scaling weighted by feature stability
#'     \item `"svast"`: supervised VAST, weighted by the mean within-group
#'       stability
#'     \item `"xvast"`: x-VAST, weighted by the maximum within-group stability
#'     \item `"none"`: no scaling
#'   }
#' @param metadata Data frame. Required for `"svast"` and `"xvast"`.
#'   Default: `NULL`.
#' @param group_col Character. Column in `metadata` with group labels, used by
#'   the supervised methods. Default: `"Group"`.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return A list of class `"run_scale"` containing:
#'   \item{data}{Scaled data (same class as `x`).}
#'   \item{scaling_factors}{Numeric vector of per-feature divisors.}
#'   \item{center_values}{Numeric vector of per-feature means.}
#'   \item{method_used}{Character.}
#'   \item{parameters}{List with `method`, `group_col`, and the group levels
#'     used by supervised methods (`groups`, otherwise `NULL`).}
#'
#' @details
#' With \eqn{\bar{x}_j} and \eqn{s_j} the mean and standard deviation of
#' feature \eqn{j} (van den Berg et al., 2006):
#'
#' * **mean**: \eqn{x - \bar{x}_j}.
#' * **auto** (unit variance): \eqn{(x - \bar{x}_j) / s_j}. Every feature gets
#'   equal weight; recommended for PCA, clustering, t-tests, and ANOVA.
#' * **pareto**: \eqn{(x - \bar{x}_j) / \sqrt{s_j}}. Shrinks large fold changes
#'   less than auto-scaling; recommended for PLS-DA and OPLS-DA.
#' * **vast** (Keun et al., 2003): auto-scaled data multiplied by the
#'   stability \eqn{\bar{x}_j / s_j} (the inverse coefficient of variation),
#'   which down-weights unstable features.
#' * **svast** (supervised VAST; Keun et al., 2003): as VAST, with the
#'   stability averaged over the within-group values
#'   \eqn{\bar{x}_{gj} / s_{gj}}.
#' * **xvast** (x-VAST; Yang et al., 2015): as VAST, with the maximum
#'   within-group stability, so a feature stable in one group keeps a large
#'   weight even when it differs strongly between groups.
#'
#' The VAST weights need positive feature means, i.e. data on a ratio scale.
#' Log-transformed values close to or below zero give meaningless weights, so
#' a warning is issued if any mean used for a weight is <= 0.
#'
#' Zero-variance features get a scaling factor of 1 to avoid division by zero.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Keun, H.C., Ebbels, T.M.D., Antti, H., Bollard, M.E., Beckonert, O.,
#' Holmes, E., Lindon, J.C., & Nicholson, J.K. (2003). Improved analysis of
#' multivariate data by variable stability scaling: application to NMR-based
#' metabolic profiling. *Analytica Chimica Acta*, 490(1-2), 265-276.
#' \doi{10.1016/S0003-2670(03)00094-1}
#'
#' van den Berg, R.A., Hoefsloot, H.C.J., Westerhuis, J.A., Smilde, A.K., &
#' van der Werf, M.J. (2006). Centering, scaling, and transformations:
#' improving the biological information content of metabolomics data.
#' *BMC Genomics*, 7, 142. \doi{10.1186/1471-2164-7-142}
#'
#' Yang, J., Zhao, X., Lu, X., Lin, X., & Xu, G. (2015). A data preprocessing
#' strategy for metabolomics to reduce the mask effect in data analysis.
#' *Frontiers in Molecular Biosciences*, 2, 4. \doi{10.3389/fmolb.2015.00004}
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
#' x <- log10(run_mvimpute(raw[, 1:20], verbose = FALSE)$data)
#'
#' res <- run_scale(x, method = "pareto")
#' res
#' round(head(res$scaling_factors), 3)
#'
#' # Supervised VAST needs group labels
#' run_scale(x, method = "svast", metadata = dims_metadata, verbose = FALSE)
run_scale <- function(
    x,
    method    = "auto",
    metadata  = NULL,
    group_col = "Group",
    verbose   = TRUE
) {

  msg <- .msg(verbose)

  .validate_data_matrix(x)

  method <- tolower(method)
  if (length(method) != 1L || !method %in% .scale_methods)
    stop("Unknown scaling method '", paste(method, collapse = ", "),
         "'. Supported: ", paste(.scale_methods, collapse = ", "), ".", call. = FALSE)

  group_vec <- NULL
  if (method %in% c("svast", "xvast")) {
    if (is.null(metadata))
      stop("'metadata' must be provided for method '", method, "'.", call. = FALSE)
    .validate_metadata(x, metadata)
    .require_col(metadata, group_col)
    group_vec <- as.factor(metadata[[group_col]])
  }

  mp <- .as_matrix_preserve(x)
  x_mat <- mp$mat

  msg(sprintf("Applying '%s' scaling...", method))

  col_means <- colMeans(x_mat, na.rm = TRUE)
  col_sds   <- matrixStats::colSds(x_mat, na.rm = TRUE)
  col_sds[col_sds == 0] <- 1 # guard zero-variance

  # Stability (mean / SD) weights for the VAST family
  .stability <- function() {
    if (method == "vast") {
      st <- col_means / col_sds
    } else {
      st <- vapply(levels(group_vec), function(lv) {
        sub <- x_mat[which(group_vec == lv), , drop = FALSE]
        s <- matrixStats::colSds(sub, na.rm = TRUE)
        s[s == 0] <- 1
        colMeans(sub, na.rm = TRUE) / s
      }, numeric(ncol(x_mat)))
      st <- matrix(st, nrow = ncol(x_mat))
      st <- if (method == "svast") rowMeans(st) else apply(st, 1, max)
    }
    if (any(st <= 0, na.rm = TRUE))
      warning("Some VAST stability weights are <= 0 (feature means <= 0); ",
              "apply VAST scaling to positive, non-log data.", call. = FALSE)
    st[st == 0] <- 1
    st
  }

  sf <- switch(method,
    none   = NULL,
    mean   = rep(1, ncol(x_mat)),
    auto   = col_sds,
    pareto = sqrt(col_sds),
    col_sds / .stability()         # vast, svast, xvast
  )

  if (is.null(sf)) {
    msg("No scaling applied.")
    x_scaled <- x_mat
    sf <- rep(1, ncol(x_mat))
  } else {
    x_scaled <- sweep(sweep(x_mat, 2, col_means, "-"), 2, sf, "/")
  }

  msg("Scaling complete.")

  structure(
    list(
      data            = .restore_class(x_scaled, mp$was_matrix),
      scaling_factors = sf,
      center_values   = col_means,
      method_used     = method,
      parameters      = list(
        method    = method,
        group_col = group_col,
        groups    = if (!is.null(group_vec)) levels(group_vec) else NULL
      )
    ),
    class = "run_scale"
  )
}

#' @param ... Ignored.
#' @rdname run_scale
#' @export
print.run_scale <- function(x, ...) {
  cat("=== Scaling Results ===\n")
  cat(sprintf("Method  : %s\n", x$method_used))
  cat(sprintf("Data    : %d samples x %d features\n",
              nrow(x$data), ncol(x$data)))
  if (!is.null(x$parameters$groups))
    cat(sprintf("Groups  : %s\n", paste(x$parameters$groups, collapse = ", ")))
  invisible(x)
}
