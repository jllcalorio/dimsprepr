#' Scale Metabolomics Data
#'
#' @description
#' Applies various scaling methods to standardize features and prepare data
#' for multivariate analysis. Includes classical scaling methods as well as
#' variance-stability (VAST) scaling variants.
#'
#' @param x Matrix or data frame. Numeric data with samples in rows and features in columns.
#' @param method Character. Scaling method to apply. Options:
#'   \itemize{
#'     \item \code{"mean"}: Mean-centering only (centers distribution at zero)
#'     \item \code{"auto"}: Auto-scaling (mean-centering + unit variance scaling)
#'     \item \code{"pareto"}: Pareto-scaling (mean-centering + sqrt(SD) scaling)
#'     \item \code{"ln"}: Logarithm (ln) transformation
#'     \item \code{"vast"}: VAST scaling (auto-scaling weighted by CV)
#'     \item \code{"svast"}: Supervised VAST scaling (auto-scaling weighted by
#'       mean class CV)
#'     \item \code{"xvast"}: Extended supervised VAST scaling (auto-scaling
#'       weighted by maximum class CV)
#'     \item \code{"none"}: No scaling method applied
#'   }
#' @param metadata Data frame. Sample metadata with \code{nrow(metadata) == nrow(x)}.
#'   Required when \code{method} is \code{"svast"} or \code{"xvast"}; ignored
#'   otherwise. Default: \code{NULL}.
#' @param group Character. Name of the column in \code{metadata} that contains
#'   group/class labels. Required when \code{method} is \code{"svast"} or
#'   \code{"xvast"}. Default: \code{"Group"}.
#' @param verbose Logical. Print progress messages. Default: \code{TRUE}.
#'
#' @return A list of class \code{"run_scale"} containing:
#'   \item{data}{Matrix or data frame of scaled data (same class as input \code{x})}
#'   \item{scaling_factors}{Numeric vector of scaling factors applied to each feature}
#'   \item{center_values}{Numeric vector of centering values (means) for each feature}
#'   \item{method_used}{Character describing scaling method applied}
#'   \item{parameters}{List of parameters used}
#'
#' @details
#' **Scaling Methods:**
#'
#' - **mean**: Mean-centering only. Each feature is centered to have mean = 0
#'   but retains original variance. Useful when features are already on similar
#'   scales and you want to preserve relative variance information.
#'
#' - **auto**: Auto-scaling (also called standardization or unit variance scaling).
#'   Each feature is mean-centered and divided by its standard deviation,
#'   resulting in mean = 0 and SD = 1 for all features. This gives equal weight
#'   to all features regardless of their original variance.
#'
#'   **Recommended for:** PCA, hierarchical clustering, correlation analysis,
#'   t-tests and ANOVA, and most univariate and multivariate methods.
#'   Note that auto-scaling can amplify noise in low-variance features.
#'
#' - **pareto**: Pareto-scaling. Each feature is mean-centered and divided by
#'   the square root of its standard deviation. This provides a compromise between
#'   no scaling and auto-scaling, reducing the influence of large values while
#'   preserving some of the original variance structure.
#'
#'   **Recommended for:** PLS-DA, OPLS-DA, sPLS-DA, and other PLS-type methods.
#'
#' - **ln**: Logarithm (ln) transformation. Each feature is transformed using
#'   the natural logarithm. This is a non-linear transformation included for
#'   cases where log-transformation is desired as part of the scaling step 
#'   (e.g., to reduce the influence of large values).
#'
#' - **vast**: VAST (Variable Stability) scaling. Extends auto-scaling by
#'   additionally weighting each feature by its overall coefficient of variation
#'   (CV = mean / SD). The scaling factor is \eqn{s_k / ({\bar{x}_k / s_k})} = \eqn{s_k^2 / \bar{x}_k}.
#'   Features with stable measurements (low CV) receive higher weights than
#'   noisy features. Equivalent to dividing auto-scaled data by the overall CV.
#'
#'   **Recommended for:** NMR-based metabolic profiling; data where measurement
#'   stability varies across features.
#'
#' - **svast**: Supervised VAST (s-VAST) scaling. A group-aware extension of
#'   VAST that weights each feature by the average of its within-class CVs
#'   (\eqn{(1/n)\sum_{j=1}^{n} \bar{x}_{jk}/s_{jk}}), rewarding features that
#'   are stable within classes. Requires \code{group}.
#'
#' - **xvast**: Extended supervised VAST (x-VAST) scaling, introduced by
#'   Yang et al. (2015). Weights each feature by the \emph{maximum} within-class
#'   ratio \eqn{\max(\bar{x}_{1k}/s_{1k}, \ldots, \bar{x}_{nk}/s_{nk})}, placing
#'   emphasis on the class in which a feature is most stable. Requires \code{group}.
#'
#' **Choosing a Scaling Method:**
#'
#' The choice depends on your downstream analysis:
#' - For exploratory data analysis (PCA, clustering): **auto-scaling**
#' - For classification/discrimination (PLS-DA): **Pareto-scaling**
#' - When feature magnitudes carry meaning: **mean-centering**
#' - For NMR/metabolomics with variable measurement stability: **VAST**
#' - When class-specific stability is important: **s-VAST** or **x-VAST**
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' van den Berg, R.A., Hoefsloot, H.C., Westerhuis, J.A., Smilde, A.K., & van der Werf, M.J. (2006).
#' Centering, scaling, and transformations: improving the biological information content of
#' metabolomics data. \emph{BMC Genomics}, 7, 142. \doi{10.1186/1471-2164-7-142}
#'
#' Keun, H.C., Ebbels, T.M.D., Antti, H., Bollard, M.E., Beckonert, O., Holmes, E., et al. (2003).
#' Improved analysis of multivariate data by variable stability scaling: application to
#' NMR-based metabolic profiling. \emph{Analytica Chimica Acta}, 490, 265–276.
#' \doi{10.1016/S0003-2670(03)00094-1}
#'
#' Yang, J., Zhao, X., Lu, X., Lin, X., & Xu, G. (2015). A data preprocessing strategy for
#' metabolomics to reduce the mask effect in data analysis. \emph{Frontiers in Molecular
#' Biosciences}, 2, 4. \doi{10.3389/fmolb.2015.00004}
#'
#' Becker, R.A., Chambers, J.M., & Wilks, A.R. (1988). \emph{The New S Language}.
#' Wadsworth & Brooks/Cole. \doi{10.1201/9781351074988}
#'
#' @importFrom matrixStats colSds
#'
#' @seealso \code{\link{run_DIpreprocess}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' x <- matrix(rnorm(100 * 50, mean = 100, sd = 50),
#'             nrow = 100, ncol = 50)
#' colnames(x) <- paste0("Feature", 1:50)
#' group <- rep(c("A", "B"), each = 50)
#'
#' # Auto-scaling for PCA
#' result1 <- run_scale(x, method = "auto")
#'
#' # Pareto-scaling for PLS-DA
#' result2 <- run_scale(x, method = "pareto")
#'
#' # Log (ln) transformation
#' result_ln <- run_scale(x, method = "ln")
#'
#' # Mean-centering only
#' result3 <- run_scale(x, method = "mean")
#'
#' # VAST scaling
#' result4 <- run_scale(x, method = "vast")
#'
#' # Supervised VAST scaling
#' result5 <- run_scale(x, method = "svast", group = group)
#'
#' # Extended supervised VAST scaling (Yang et al., 2015)
#' result6 <- run_scale(x, method = "xvast", group = group)
#' }
run_scale <- function(
    x,
    method   = "auto",
    metadata = NULL,
    group    = "Group",
    verbose  = TRUE
) {

  msg <- function(...) if (verbose) message(...)

  # ---- Input validation -------------------------------------------------------
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }

  method <- tolower(method)
  valid_methods <- c("mean", "auto", "pareto", "ln", "vast", "svast", "xvast", "none")
  if (!method %in% valid_methods) {
    stop("Unknown scaling method '", method, "'. ",
         "Supported: ", paste(valid_methods, collapse = ", "))
  }

  # group vector extraction — required for supervised VAST methods
  group_vec <- NULL
  if (method %in% c("svast", "xvast")) {
    if (is.null(metadata)) {
      stop("'metadata' must be provided for method '", method, "'.")
    }
    if (!is.data.frame(metadata)) {
      stop("'metadata' must be a data frame.")
    }
    if (!group %in% colnames(metadata)) {
      stop("Column '", group, "' not found in 'metadata'.")
    }
    if (nrow(metadata) != nrow(x)) {
      stop("nrow(metadata) (", nrow(metadata), ") must equal nrow(x) (", nrow(x), ").")
    }
    group_vec <- as.factor(metadata[[group]])
  }

  was_matrix <- is.matrix(x)
  x_mat <- as.matrix(x)

  msg(sprintf("Applying '%s' scaling...", method))

  # ---- Always-needed quantities -----------------------------------------------
  col_means <- colMeans(x_mat, na.rm = TRUE)
  col_sds   <- matrixStats::colSds(x_mat, na.rm = TRUE)
  col_sds[col_sds == 0] <- 1  # guard against zero-variance features

  # Helper: compute within-class mean/SD ratio (x̄_jk / s_jk) for each class j
  # Returns a matrix: nlevels(group) × ncol(x_mat)
  .class_cv_ratios <- function() {
    lvls <- levels(group)
    ratios <- vapply(lvls, function(lv) {
      idx  <- which(group == lv)
      sub  <- x_mat[idx, , drop = FALSE]
      m    <- colMeans(sub, na.rm = TRUE)
      s    <- matrixStats::colSds(sub, na.rm = TRUE)
      s[s == 0] <- 1
      m / s
    }, numeric(ncol(x_mat)))
    # vapply returns ncol x nlevels; transpose to nlevels x ncol
    t(ratios)
  }

  # ---- Scaling switch ---------------------------------------------------------
  result <- switch(method,

    "none" = {
      msg("No scaling applied.")
      list(
        x_scaled        = x_mat,
        scaling_factors = rep(1, ncol(x_mat))
      )
    },

    "mean" = {
      msg("Applying mean-centering...")
      list(
        x_scaled        = scale(x_mat, center = TRUE, scale = FALSE),
        scaling_factors = rep(1, ncol(x_mat))
      )
    },

    "auto" = {
      msg("Applying auto-scaling (mean-centering + unit variance)...")
      list(
        x_scaled        = scale(x_mat, center = TRUE, scale = col_sds),
        scaling_factors = col_sds
      )
    },

    "pareto" = {
      msg("Applying Pareto-scaling (mean-centering + sqrt(SD))...")
      sf <- sqrt(col_sds)
      list(
        x_scaled        = scale(x_mat, center = TRUE, scale = sf),
        scaling_factors = sf
      )
    },

    "ln" = {
      msg("Applying natural log (ln) transformation...")
      list(
        x_scaled        = log(x_mat),
        scaling_factors = rep(1, ncol(x_mat))
      )
    },

    # VAST: scale by s_k^2 / x̄_k  (i.e., auto-scale then divide by overall CV)
    # Formula (6) in Yang et al. (2015); originally Keun et al. (2003):
    #   x'_ik = (x̄_k / s_k) * x_ik  →  equivalent divisor: s_k / (x̄_k / s_k) = s_k² / x̄_k
    "vast" = {
      msg("Applying VAST scaling (Keun et al., 2003)...")
      overall_cv <- col_means / col_sds          # x̄_k / s_k  (inverse CV)
      overall_cv[overall_cv == 0] <- 1
      sf <- col_sds / overall_cv                 # = s_k² / x̄_k
      list(
        x_scaled        = scale(x_mat, center = TRUE, scale = sf),
        scaling_factors = sf
      )
    },

    # s-VAST: weight by mean within-class CV ratio  (Formula 7, Yang et al., 2015)
    #   x'_ik = [ (1/n) Σ_j (x̄_jk / s_jk) ] * x_ik
    #   divisor: s_k / mean_j(x̄_jk / s_jk)
    "svast" = {
      msg("Applying supervised VAST scaling / s-VAST (Yang et al., 2015)...")
      ratios      <- .class_cv_ratios()          # nlevels × nfeatures
      mean_ratio  <- colMeans(ratios)            # mean across classes, length = nfeatures
      mean_ratio[mean_ratio == 0] <- 1
      sf <- col_sds / mean_ratio
      list(
        x_scaled        = scale(x_mat, center = TRUE, scale = sf),
        scaling_factors = sf
      )
    },

    # x-VAST: weight by maximum within-class CV ratio  (Formula 5, Yang et al., 2015)
    #   x'_ik = max_j(x̄_jk / s_jk) * x_ik
    #   divisor: s_k / max_j(x̄_jk / s_jk)
    "xvast" = {
      msg("Applying extended supervised VAST scaling / x-VAST (Yang et al., 2015)...")
      ratios    <- .class_cv_ratios()            # nlevels × nfeatures
      max_ratio <- apply(ratios, 2, max)         # max across classes, length = nfeatures
      max_ratio[max_ratio == 0] <- 1
      sf <- col_sds / max_ratio
      list(
        x_scaled        = scale(x_mat, center = TRUE, scale = sf),
        scaling_factors = sf
      )
    }
  )

  msg("Scaling complete.")

  # ---- Clean scale attributes and restore class --------------------------------
  x_scaled <- result$x_scaled
  attributes(x_scaled) <- list(dim = dim(x_scaled), dimnames = dimnames(x_scaled))

  if (!was_matrix) x_scaled <- as.data.frame(x_scaled)

  # ---- Assemble output --------------------------------------------------------
  out <- list(
    data            = x_scaled,
    scaling_factors = result$scaling_factors,
    center_values   = col_means,
    method_used     = method,
    parameters      = list(
      method = method,
      group  = if (!is.null(group)) levels(group) else NULL
    )
  )

  class(out) <- "run_scale"
  out
}


# ---- S3 print method ---------------------------------------------------------
#' @export
print.run_scale <- function(x, ...) {
  cat("=== Scaling Results ===\n")
  cat(sprintf("Method  : %s\n", x$method_used))
  cat(sprintf("Data    : %d samples x %d features\n",
              nrow(x$data), ncol(x$data)))
  if (!is.null(x$parameters$group)) {
    cat(sprintf("Groups  : %s\n", paste(x$parameters$group, collapse = ", ")))
  }
  invisible(x)
}