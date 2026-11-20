#' Scale Metabolomics Data
#'
#' @description
#' Applies various scaling methods to standardize features and prepare data
#' for multivariate analysis. Includes classical scaling methods as well as
#' variance-stability (VAST) scaling variants.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#'   For `print()`, a `run_scale` object.
#' @param method Character. Scaling method:
#'   \itemize{
#'     \item \code{"mean"}: Mean-centering only
#'     \item \code{"auto"}: Auto-scaling (mean-centering + unit variance)
#'     \item \code{"pareto"}: Pareto-scaling (mean-centering + sqrt(SD))
#'     \item \code{"ln"}: Natural log transformation
#'     \item \code{"vast"}: VAST scaling (auto-scaling weighted by CV)
#'     \item \code{"svast"}: Supervised VAST (weighted by mean class CV)
#'     \item \code{"xvast"}: Extended supervised VAST (weighted by max class CV)
#'     \item \code{"none"}: No scaling
#'   }
#' @param metadata Data frame. Required for \code{"svast"} and \code{"xvast"}.
#'   Default: NULL.
#' @param group_col Character. Column in \code{metadata} with group labels.
#'   Required for supervised methods. Default: "Group".
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A list of class \code{"run_scale"} containing:
#'   \item{data}{Scaled data (same class as input \code{x})}
#'   \item{scaling_factors}{Numeric vector of per-feature scaling factors}
#'   \item{center_values}{Numeric vector of per-feature means}
#'   \item{method_used}{Character}
#'   \item{parameters}{List of parameters}
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
#' * **vast** (variable stability; Keun et al., 2003): auto-scaled data
#'   multiplied by \eqn{\bar{x}_j / s_j}, i.e. divided by the coefficient of
#'   variation, which down-weights unstable features.
#' * **svast** / **xvast**: supervised VAST variants (Yang et al., 2015) that
#'   use the mean (svast) or maximum (xvast) of the within-group
#'   \eqn{\bar{x} / s} ratios, so they require `metadata` and `group_col`.
#' * **ln**: natural log of the data (a transformation, kept for backward
#'   compatibility; prefer [run_transform()]).
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

  # group vector extraction for supervised VAST methods
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

  # Within-class mean/SD ratio helper for supervised methods
  .class_cv_ratios <- function() {
    lvls <- levels(group_vec)
    ratios <- vapply(lvls, function(lv) {
      sub <- x_mat[which(group_vec == lv), , drop = FALSE]
      m <- colMeans(sub, na.rm = TRUE)
      s <- matrixStats::colSds(sub, na.rm = TRUE)
      s[s == 0] <- 1
      m / s
    }, numeric(ncol(x_mat)))
    t(ratios) # nlevels x nfeatures
  }

  result <- switch(method,

    "none" = {
      msg("No scaling applied.")
      list(x_scaled = x_mat, scaling_factors = rep(1, ncol(x_mat)))
    },

    "mean" = {
      list(x_scaled = scale(x_mat, center = TRUE, scale = FALSE),
           scaling_factors = rep(1, ncol(x_mat)))
    },

    "auto" = {
      list(x_scaled = scale(x_mat, center = TRUE, scale = col_sds),
           scaling_factors = col_sds)
    },

    "pareto" = {
      sf <- sqrt(col_sds)
      list(x_scaled = scale(x_mat, center = TRUE, scale = sf),
           scaling_factors = sf)
    },

    "ln" = {
      list(x_scaled = log(x_mat), scaling_factors = rep(1, ncol(x_mat)))
    },

    "vast" = {
      overall_cv <- col_means / col_sds
      overall_cv[overall_cv == 0] <- 1
      sf <- col_sds / overall_cv
      list(x_scaled = scale(x_mat, center = TRUE, scale = sf),
           scaling_factors = sf)
    },

    "svast" = {
      ratios     <- .class_cv_ratios()
      mean_ratio <- colMeans(ratios)
      mean_ratio[mean_ratio == 0] <- 1
      sf <- col_sds / mean_ratio
      list(x_scaled = scale(x_mat, center = TRUE, scale = sf),
           scaling_factors = sf)
    },

    "xvast" = {
      ratios    <- .class_cv_ratios()
      max_ratio <- apply(ratios, 2, max)
      max_ratio[max_ratio == 0] <- 1
      sf <- col_sds / max_ratio
      list(x_scaled = scale(x_mat, center = TRUE, scale = sf),
           scaling_factors = sf)
    }
  )

  msg("Scaling complete.")

  # Clean scale attributes
  x_scaled <- result$x_scaled
  attributes(x_scaled) <- list(dim = dim(x_scaled), dimnames = dimnames(x_scaled))

  structure(
    list(
      data            = .restore_class(x_scaled, mp$was_matrix),
      scaling_factors = result$scaling_factors,
      center_values   = col_means,
      method_used     = method,
      parameters      = list(
        method    = method,
        # BUG FIX: was `levels(group)` which always returned NULL;
        # now correctly uses `levels(group_vec)`
        group_col = if (!is.null(group_vec)) levels(group_vec) else NULL
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
  if (!is.null(x$parameters$group_col))
    cat(sprintf("Groups  : %s\n", paste(x$parameters$group_col, collapse = ", ")))
  invisible(x)
}
