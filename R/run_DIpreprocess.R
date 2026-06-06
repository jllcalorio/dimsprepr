#' Direct-Injection Metabolomics Preprocessing Pipeline
#'
#' @description
#' Runs a complete, sequential preprocessing pipeline for metabolomics data by
#' orchestrating the individual \code{run_*} functions of \pkg{pondeR}. Steps
#' execute in a fixed, reproducible order regardless of argument arrangement.
#'
#' When \code{normalize_method} and/or \code{transform_method} are supplied as
#' character vectors with more than one element, \strong{all combinations} are
#' computed.  The steps prior to normalization (validation, outlier removal,
#' missing-value filtering, imputation, and drift/batch correction) are executed
#' \strong{only once}; the normalization → transformation → scaling → quality
#' filtering → replicate-merging branch is then repeated for every
#' normalization × transformation pair.  The return value in this case is a
#' plain \code{list} (class \code{"run_DIpreprocess_multi"}) whose elements are
#' named \code{<normalize>_<transform>} and each hold a standard
#' \code{run_DIpreprocess} object.
#'
#' @param x Matrix or data frame. Numeric feature data with \strong{samples in rows}
#'   and \strong{features in columns}. All values should be non-negative raw intensities;
#'   zeros are treated as missing internally.
#' @param metadata Data frame. Sample metadata with \code{nrow(metadata) == nrow(x)}.
#'   At minimum, the group column (default \code{"Group"}) must be present.
#' @param sample_id_col Character. Column in \code{metadata} containing unique
#'   sample identifiers that match \code{rownames(x)}. Default: \code{"Sample"}.
#' @param group_col Character. Column in \code{metadata} containing group labels.
#'   Default: \code{"Group"}.
#' @param qc_types Character vector. Group labels that identify QC samples.
#'   Default: \code{c("QC", "SQC", "EQC")}.
#' @param batch_col Character. Column in \code{metadata} containing batch numbers.
#'   Required for drift correction. Default: \code{"Batch"}.
#' @param injection_col Character. Column in \code{metadata} containing injection
#'   order (integer). Required for drift correction. Default: \code{"InjectionSequence"}.
#' @param norm_factor_col Character. Column in \code{metadata} containing
#'   external normalization factors. Default: \code{"Normalization"}.
#' @param subject_id_col Character. Column in \code{metadata} containing
#'   subject identifiers for technical-replicate merging. Default: \code{"SubjectID"}.
#' @param outliers Character vector or \code{NULL}. Row names of \code{x}
#'   to remove before processing. Default: \code{NULL}.
#' @param missing_threshold Numeric \[0, 1\]. Missing value threshold. Default: \code{0.2}.
#' @param missing_by_group Logical. If \code{TRUE}, assess per group. Default: \code{TRUE}.
#' @param missing_include_qc Logical. If \code{TRUE}, include QC samples. Default: \code{FALSE}.
#' @param impute_fraction Numeric > 0. Fraction of the smallest positive value
#'   per feature used for imputation. Default: \code{0.2}.
#' @param positive_only Logical. If \code{TRUE}, imputation only considers positive values. Default: \code{TRUE}.
#' @param correct_drift Logical. If \code{TRUE}, apply QCRSC correction. Default: \code{TRUE}.
#' @param remove_uncorrected Logical. If \code{TRUE}, remove features QCRSC
#'   could not correct. Default: \code{FALSE}.
#' @param spline_smooth Numeric \[0, 1\]. Smoothing parameter. Default: \code{0}.
#' @param spline_spar_limit Numeric vector. Lower/upper bounds for spline. Default: \code{c(-1.5, 1.5)}.
#' @param correct_on_log Logical. If \code{TRUE}, fit on log-transformed data. Default: \code{TRUE}.
#' @param min_qc_per_batch Integer. Minimum QC samples per batch. Default: \code{5}.
#' @param normalize_method Character or character vector. One or more normalization
#'   methods to apply. When multiple methods are supplied every method is paired
#'   with every element of \code{transform_method} (Cartesian product).
#' \itemize{
#'   \item \code{"sum"}: Normalizes by sum.
#'   \item \code{"median"}: Normalizes by median.
#'   \item \code{"specific_factor"}: Uses external factors in \code{norm_factor_col}.
#'   \item \code{"pqn_global"}: PQN using a global reference.
#'   \item \code{"pqn_reference"}: PQN using \code{normalize_ref_sample}.
#'   \item \code{"pqn_group"}: PQN using pooled QC samples as reference.
#'   \item \code{"quantile"}: Quantile normalization.
#'   \item \code{"none"}: No normalization.
#' }
#' Default: \code{"median"}.
#' @param normalize_ref_sample Character. Reference sample name for PQN. Default: \code{NULL}.
#' @param normalize_qc_method Character. QC normalization when using factors.
#'   One of \code{"mean"}, \code{"median"}, \code{"none"}. Default: \code{"none"}.
#' @param transform_method Character or character vector. One or more transformation
#'   methods to apply. When multiple methods are supplied every method is paired
#'   with every element of \code{normalize_method} (Cartesian product).
#' \itemize{
#'   \item \code{"log2"}, \code{"log10"}, \code{"sqrt"}, \code{"cbrt"},
#'   \code{"clr"}, \code{"vsn"}, \code{"glog"}, or \code{"none"}.
#' }
#' Default: \code{"log10"}.
#' @param vsn_cores Integer or \code{"max"}. Cores for VSN. Default: \code{1}.
#' @param scale_nonpls Character. Scaling for NONPLS branch (\code{"auto"}, \code{"pareto"}, \code{"mean"}).
#' @param scale_pls Character. Scaling for PLS branch (\code{"pareto"}, \code{"auto"}, \code{"mean"}).
#' @param rsd_threshold Numeric \[0, 1\]. QC-RSD threshold. Default: \code{0.3}.
#' @param rsd_qc_type Character. QC type for RSD (\code{"EQC"}, \code{"SQC"}, \code{"QC"}).
#' @param variance_percentile Numeric \[0, 100\]. Percentile to filter. Default: \code{10}.
#' @param scale_filter_ref Character. Harmonisation strategy.
#' \itemize{
#'   \item \code{"auto"}: Keep intersection of both branches.
#'   \item \code{"NONPLS"}: Apply NONPLS filters to both.
#'   \item \code{"PLS"}: Apply PLS filters to both.
#' }
#' @param merge_replicates Logical. Average technical replicates. Default: \code{FALSE}.
#' @param verbose Logical. Print progress. Default: \code{TRUE}.
#'
#' @return
#' \describe{
#'   \item{Single combination (scalar \code{normalize_method} and scalar
#'     \code{transform_method})}{A named list of class \code{run_DIpreprocess}
#'     — identical to the original behaviour.}
#'   \item{Multiple combinations (either argument is a vector of length > 1)}{A
#'     named list of class \code{run_DIpreprocess_multi}.  Each element is named
#'     \code{<normalize>_<transform>} (e.g. \code{sum_log10}) and is itself a
#'     \code{run_DIpreprocess} object.  The shared pre-normalization data
#'     (post-drift-correction) is stored in the \code{shared} attribute.}
#' }
#'
#' @details
#' \strong{Preprocessing Workflow (Fixed Order):}
#' \enumerate{
#'   \item Input validation
#'   \item Outlier removal
#'   \item Missing-value filtering (\code{\link{run_filtermissing}})
#'   \item Missing-value imputation (\code{\link{run_mvimpute}})
#'   \item Signal drift correction (\code{\link{run_driftBatchCorrect}})
#'   \item \emph{\[Per combination\]} Normalization (\code{\link{run_normalize}})
#'   \item \emph{\[Per combination\]} Transformation (\code{\link{run_transform}})
#'   \item \emph{\[Per combination\]} Scaling (\code{\link{run_scale}})
#'   \item \emph{\[Per combination\]} Quality filtering (\code{\link{run_filterRSD}},
#'     \code{\link{run_filtervariance}})
#'   \item \emph{\[Per combination\]} Common-feature harmonisation & replicate merging
#' }
#'
#' @note
#' \strong{Zero Handling:} The pipeline automatically converts all \code{0} values in
#' \code{x} to \code{NA} prior to analysis.
#'
#' \strong{Multi-method pre-flight checks:} When vectors are supplied for
#' \code{normalize_method} or \code{transform_method}, the function validates
#' all method-specific requirements (e.g. \code{norm_factor_col} for
#' \code{"specific_factor"}, \code{normalize_ref_sample} for
#' \code{"pqn_reference"}) \emph{before} executing any preprocessing step, so
#' errors are surfaced early rather than mid-run.
#'
#' @author John Lennon L. Calorio
#'
#' @seealso
#' \code{\link{run_filtermissing}}, \code{\link{run_mvimpute}},
#' \code{\link{run_driftBatchCorrect}}, \code{\link{run_normalize}}
#'
#' @references
#' Kirwan, J.A., et al. (2013). \emph{Analytical and Bioanalytical Chemistry}, 405, 5147-5157.
#'
#' @importFrom matrixStats colMaxs
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n_s <- 80; n_f <- 100
#' x <- matrix(abs(rnorm(n_s * n_f, 500, 150)), nrow = n_s, ncol = n_f)
#' colnames(x) <- paste0("Feature", seq_len(n_f))
#' x[sample(length(x), 400)] <- 0
#'
#' meta <- data.frame(
#'   Sample            = paste0("S", seq_len(n_s)),
#'   Group             = rep(c("Control", "Treatment", "QC"), c(30, 30, 20)),
#'   Batch             = rep(1:2, each = 40),
#'   InjectionSequence = seq_len(n_s),
#'   SubjectID         = c(paste0("BIO", seq_len(60)), rep(NA, 20)),
#'   stringsAsFactors  = FALSE
#' )
#' rownames(x) <- meta$Sample
#'
#' # --- Single combination (original behaviour) ---
#' result_single <- run_DIpreprocess(
#'   x                = x,
#'   metadata         = meta,
#'   normalize_method = "pqn_group",
#'   transform_method = "log2",
#'   correct_drift    = FALSE
#' )
#' class(result_single)          # "run_DIpreprocess"
#' dim(result_single$data_nonpls)
#'
#' # --- Multiple transform methods ---
#' result_multi <- run_DIpreprocess(
#'   x                = x,
#'   metadata         = meta,
#'   normalize_method = "median",
#'   transform_method = c("log10", "vsn"),
#'   correct_drift    = FALSE
#' )
#' class(result_multi)           # "run_DIpreprocess_multi"
#' names(result_multi)           # "median_log10"  "median_vsn"
#'
#' # --- Full Cartesian product ---
#' result_cart <- run_DIpreprocess(
#'   x                = x,
#'   metadata         = meta,
#'   normalize_method = c("sum", "median"),
#'   transform_method = c("log10", "vsn"),
#'   correct_drift    = FALSE
#' )
#' names(result_cart)
#' # "sum_log10"  "sum_vsn"  "median_log10"  "median_vsn"
#' }
run_DIpreprocess <- function(
    x,
    metadata,

    # — Column names ——————————————————————————————————————————————————————————
    sample_id_col    = "Sample",
    group_col        = "Group",
    qc_types         = c("QC", "SQC", "EQC"),
    batch_col        = "Batch",
    injection_col    = "InjectionSequence",
    norm_factor_col  = "Normalization",
    subject_id_col   = "SubjectID",

    # — Outlier removal ———————————————————————————————————————————————————————
    outliers         = NULL,

    # — Missing-value filter ——————————————————————————————————————————————————
    missing_threshold  = 0.2,
    missing_by_group   = TRUE,
    missing_include_qc = FALSE,

    # — Missing-value imputation ——————————————————————————————————————————————
    impute_fraction  = 0.2,
    positive_only    = TRUE,

    # — Drift / batch correction ——————————————————————————————————————————————
    correct_drift        = TRUE,
    remove_uncorrected   = FALSE,
    spline_smooth        = 0,
    spline_spar_limit    = c(-1.5, 1.5),
    correct_on_log       = TRUE,
    min_qc_per_batch     = 5L,

    # — Normalization —————————————————————————————————————————————————————————
    normalize_method     = "median",
    normalize_ref_sample = NULL,
    normalize_qc_method  = "none",

    # — Transformation ————————————————————————————————————————————————————————
    transform_method = "log10",
    vsn_cores        = 1L,

    # — Scaling ———————————————————————————————————————————————————————————————
    scale_nonpls     = "auto",
    scale_pls        = "pareto",

    # — Quality filtering —————————————————————————————————————————————————————
    rsd_threshold       = 0.3,
    rsd_qc_type         = "EQC",
    variance_percentile = 10,
    scale_filter_ref    = "auto",

    # — Replicate merging —————————————————————————————————————————————————————
    merge_replicates = FALSE,

    verbose          = TRUE
) {

  t_start <- proc.time()[["elapsed"]]
  msg <- function(...) if (verbose) message(...)

  # ===========================================================================
  # INTERNAL HELPERS
  # ===========================================================================

  .record_dim <- function(dims_df, label, data) {
    rbind(dims_df, data.frame(
      Step     = label,
      Samples  = nrow(data),
      Features = ncol(data),
      stringsAsFactors = FALSE
    ))
  }

  .empty_dims <- function() {
    data.frame(Step = character(), Samples = integer(), Features = integer(),
               stringsAsFactors = FALSE)
  }

  # ---------------------------------------------------------------------------
  # PRE-FLIGHT: validate method vectors and method-specific requirements
  # ---------------------------------------------------------------------------

  valid_norm_methods <- c("sum", "median", "specific_factor", "pqn_global",
                          "pqn_reference", "pqn_group", "quantile", "none")
  valid_trans_methods <- c("log2", "log10", "sqrt", "cbrt", "clr",
                           "arcsin_sqrt", "vsn", "glog", "none")

  # Coerce to character vector and deduplicate while preserving order
  normalize_method <- unique(as.character(normalize_method))
  transform_method <- unique(as.character(transform_method))

  # Validate every supplied method name up front
  bad_norm <- setdiff(tolower(normalize_method), valid_norm_methods)
  if (length(bad_norm) > 0L)
    stop(sprintf(
      "Unknown normalize_method value(s): %s.\n  Supported: %s.",
      paste(bad_norm, collapse = ", "),
      paste(valid_norm_methods, collapse = ", ")
    ))

  bad_trans <- setdiff(tolower(transform_method), valid_trans_methods)
  if (length(bad_trans) > 0L)
    stop(sprintf(
      "Unknown transform_method value(s): %s.\n  Supported: %s.",
      paste(bad_trans, collapse = ", "),
      paste(valid_trans_methods, collapse = ", ")
    ))

  # Method-specific pre-flight checks -----------------------------------------

  if ("specific_factor" %in% tolower(normalize_method)) {
    if (!norm_factor_col %in% colnames(metadata))
      stop(sprintf(
        paste0("normalize_method = 'specific_factor' requires a metadata column ",
               "named '%s' (set via norm_factor_col), but it was not found.\n",
               "  Available columns: %s"),
        norm_factor_col, paste(colnames(metadata), collapse = ", ")
      ))
    if (!is.numeric(metadata[[norm_factor_col]]))
      stop(sprintf(
        "normalize_method = 'specific_factor': column '%s' must be numeric.",
        norm_factor_col
      ))
  }

  if ("pqn_reference" %in% tolower(normalize_method)) {
    if (is.null(normalize_ref_sample) || !nzchar(normalize_ref_sample))
      stop(paste0(
        "normalize_method = 'pqn_reference' requires a reference sample name ",
        "supplied via normalize_ref_sample, but none was provided."
      ))
    if (sample_id_col %in% colnames(metadata)) {
      if (!normalize_ref_sample %in% metadata[[sample_id_col]])
        stop(sprintf(
          paste0("normalize_method = 'pqn_reference': reference sample '%s' was ",
                 "not found in the '%s' column of metadata."),
          normalize_ref_sample, sample_id_col
        ))
    }
  }

  if ("pqn_group" %in% tolower(normalize_method)) {
    if (group_col %in% colnames(metadata)) {
      qc_present <- any(metadata[[group_col]] %in% qc_types)
      if (!qc_present)
        stop(sprintf(
          paste0("normalize_method = 'pqn_group' requires QC samples (one of: %s) ",
                 "in the '%s' column, but none were found."),
          paste(qc_types, collapse = ", "), group_col
        ))
    }
  }

  if ("glog" %in% tolower(transform_method)) {
    if (!requireNamespace("pmp", quietly = TRUE))
      stop(paste0(
        "transform_method = 'glog' requires the 'pmp' package.\n",
        "  Install it with: BiocManager::install('pmp')"
      ))
  }

  if ("vsn" %in% tolower(transform_method)) {
    if (!requireNamespace("vsn", quietly = TRUE))
      stop(paste0(
        "transform_method = 'vsn' requires the 'vsn' package.\n",
        "  Install it with: BiocManager::install('vsn')"
      ))
    if (!requireNamespace("BiocParallel", quietly = TRUE))
      stop(paste0(
        "transform_method = 'vsn' also requires the 'BiocParallel' package.\n",
        "  Install it with: BiocManager::install('BiocParallel')"
      ))
  }

  # Determine multi-method mode -----------------------------------------------
  is_multi <- length(normalize_method) > 1L || length(transform_method) > 1L

  # Build Cartesian product of combinations
  combos <- expand.grid(
    norm  = normalize_method,
    trans = transform_method,
    stringsAsFactors = FALSE
  )
  combo_names <- paste(combos$norm, combos$trans, sep = "_")

  if (is_multi) {
    msg(sprintf(
      "Multi-method mode: %d combination(s) detected — %s.",
      nrow(combos), paste(combo_names, collapse = ", ")
    ))
    msg("Steps 1-5 (through drift correction) will run once; ",
        "Steps 6-10 will run per combination.")
  }

  # ===========================================================================
  # STEPS 1-5: RUN ONCE
  # ===========================================================================

  shared <- list(
    metadata             = NULL,
    data_raw             = NULL,
    data_missing_filtered = NULL,
    data_imputed         = NULL,
    data_corrected       = NULL,
    uncorrected_features = character(0L),
    dimensions           = .empty_dims(),
    error                = NULL
  )

  tryCatch({

    # =========================================================================
    # 1. INPUT VALIDATION
    # =========================================================================

    msg("Step 1/10: Validating inputs...")

    if (!is.matrix(x) && !is.data.frame(x))
      stop("'x' must be a matrix or data frame.")
    if (!is.data.frame(metadata))
      stop("'metadata' must be a data frame.")
    if (nrow(x) != nrow(metadata))
      stop("nrow(x) != nrow(metadata).")
    if (!is.numeric(missing_threshold) || length(missing_threshold) != 1L ||
        missing_threshold < 0 || missing_threshold > 1)
      stop("'missing_threshold' must be a single numeric value in [0, 1].")
    if (!is.numeric(impute_fraction) || length(impute_fraction) != 1L ||
        impute_fraction <= 0)
      stop("'impute_fraction' must be a single positive numeric value.")
    if (!scale_filter_ref %in% c("auto", "NONPLS", "PLS"))
      stop("'scale_filter_ref' must be one of: 'auto', 'NONPLS', 'PLS'.")
    if (!rsd_qc_type %in% c("QC", "SQC", "EQC"))
      stop("'rsd_qc_type' must be one of: 'QC', 'SQC', 'EQC'.")
    if (!is.numeric(variance_percentile) || length(variance_percentile) != 1L ||
        variance_percentile < 0 || variance_percentile > 100)
      stop("'variance_percentile' must be a numeric value in [0, 100].")
    if (!group_col %in% colnames(metadata))
      stop(sprintf("group_col '%s' not found in metadata.", group_col))

    for (col_name in c(batch_col, injection_col, norm_factor_col)) {
      if (col_name %in% colnames(metadata)) {
        if (!is.numeric(metadata[[col_name]]))
          stop(sprintf(
            "Metadata column '%s' must be numeric.", col_name
          ))
      }
    }

    if (injection_col %in% colnames(metadata)) {
      if (any(duplicated(metadata[[injection_col]], incomparables = NA)))
        stop(sprintf(
          "Duplicate values found in '%s'. Injection sequence orders must be unique.",
          injection_col
        ))
    }

    if (sample_id_col %in% colnames(metadata)) {
      if (any(duplicated(metadata[[sample_id_col]]))) {
        dups <- metadata[[sample_id_col]][duplicated(metadata[[sample_id_col]])]
        stop(sprintf(
          "Duplicate sample IDs in '%s'. Duplicates: %s",
          sample_id_col, paste(head(unique(dups), 5), collapse = ", ")
        ))
      }
    }

    if (norm_factor_col %in% colnames(metadata) && group_col %in% colnames(metadata)) {
      is_non_qc   <- !metadata[[group_col]] %in% qc_types
      is_na_norm  <- is.na(metadata[[norm_factor_col]])
      if (any(is_non_qc & is_na_norm)) {
        bad_samps <- rownames(metadata)[is_non_qc & is_na_norm]
        if (is.null(bad_samps) || length(bad_samps) == 0L)
          bad_samps <- which(is_non_qc & is_na_norm)
        stop(sprintf(
          "Missing values (NA) in '%s' for non-QC samples. Check: %s",
          norm_factor_col, paste(head(bad_samps, 5), collapse = ", ")
        ))
      }
    }

    df <- as.data.frame(as.matrix(x))

    if (sample_id_col %in% colnames(metadata)) {
      rownames(df)       <- as.character(metadata[[sample_id_col]])
      rownames(metadata) <- as.character(metadata[[sample_id_col]])
    } else if (!is.null(rownames(x))) {
      rownames(metadata) <- rownames(x)
    }

    is_numeric_col <- vapply(df, function(col) is.numeric(col) || is.integer(col), logical(1L))
    if (!all(is_numeric_col)) {
      bad_cols <- names(df)[!is_numeric_col]
      stop(sprintf(
        "Input 'x' must contain ONLY numeric columns. Non-numeric: %s.",
        paste(head(bad_cols, 5), collapse = ", ")
      ))
    }

    shared$dimensions <- .record_dim(shared$dimensions, "Original", df)

    # =========================================================================
    # 2. OUTLIER REMOVAL
    # =========================================================================

    msg("Step 2/10: Checking for outliers...")

    if (!is.null(outliers) && length(outliers) > 0L) {
      present   <- intersect(outliers, rownames(df))
      missing_o <- setdiff(outliers,   rownames(df))

      if (length(present) > 0L) {
        df       <- df[!rownames(df)       %in% present, , drop = FALSE]
        metadata <- metadata[!rownames(metadata) %in% present, , drop = FALSE]
        msg(sprintf("  Removed %d outlier(s): %s",
                    length(present), paste(present, collapse = ", ")))
      }
      if (length(missing_o) > 0L)
        warning("Outliers not found in data: ",
                paste(missing_o, collapse = ", "), call. = FALSE)

      shared$dimensions <- .record_dim(shared$dimensions, "After outlier removal", df)
    }

    # =========================================================================
    # METADATA SNAPSHOT (after outlier removal)
    # =========================================================================

    .meta_col <- function(col) {
      if (col %in% colnames(metadata)) metadata[[col]]
      else rep(NA, nrow(metadata))
    }

    meta_snap <- data.frame(
      Sample            = rownames(metadata),
      Group             = .meta_col(group_col),
      Group_            = ifelse(.meta_col(group_col) %in% qc_types,
                                 "QC", .meta_col(group_col)),
      Batch             = suppressWarnings(as.integer(.meta_col(batch_col))),
      InjectionSequence = suppressWarnings(as.numeric(.meta_col(injection_col))),
      SubjectID         = as.character(.meta_col(subject_id_col)),
      stringsAsFactors  = FALSE,
      row.names         = rownames(metadata)
    )
    # Preserve the user-specified normalization factor column under its original
    # name so that run_normalize(factor_col = norm_factor_col) can find it.
    meta_snap[[norm_factor_col]] <- suppressWarnings(
      as.numeric(.meta_col(norm_factor_col))
    )
    shared$metadata <- meta_snap

    qc_rows     <- meta_snap$Group_ == "QC"
    non_qc_rows <- !qc_rows

    # =========================================================================
    # ZEROS → NA
    # =========================================================================

    n_zeros <- sum(df == 0L, na.rm = TRUE)
    if (n_zeros > 0L) {
      df[df == 0L] <- NA
      msg(sprintf("  Converted %d zero(s) to NA (structural missingness).", n_zeros))
    }

    all_na <- colSums(is.na(df)) == nrow(df)
    if (any(all_na)) {
      df <- df[, !all_na, drop = FALSE]
      msg(sprintf("  Removed %d all-NA feature(s). %d remaining.",
                  sum(all_na), ncol(df)))
    }

    shared$data_raw   <- df
    shared$dimensions <- .record_dim(shared$dimensions, "After zero -> NA conversion", df)

    # =========================================================================
    # 3. MISSING-VALUE FILTERING
    # =========================================================================

    msg(sprintf("Step 3/10: Missing-value filtering (threshold = %.0f%%)...",
                missing_threshold * 100))

    filt_miss <- run_filtermissing(
      x               = df,
      metadata        = metadata,
      threshold       = missing_threshold,
      filter_by_group = missing_by_group,
      include_QC      = missing_include_qc,
      group_col       = group_col,
      qc_types        = qc_types,
      zero_as_missing = FALSE,
      verbose         = FALSE
    )

    df                        <- filt_miss$data
    shared$data_missing_filtered <- df
    shared$dimensions         <- .record_dim(shared$dimensions, "After missing filter", df)

    msg(sprintf("  Removed %d feature(s) (%.1f%%). %d remaining.",
                filt_miss$n_features_removed,
                filt_miss$n_features_removed / filt_miss$n_features_before * 100,
                ncol(df)))

    # =========================================================================
    # 4. MISSING-VALUE IMPUTATION
    # =========================================================================

    msg(sprintf("Step 4/10: Imputing missing values (fraction = %.4f)...",
                impute_fraction))

    imp <- run_mvimpute(
      x             = df,
      method        = impute_fraction,
      positive_only = positive_only,
      verbose       = FALSE
    )
    df                   <- imp$data
    shared$data_imputed  <- df
    msg(sprintf("  Imputed %d missing value(s).", imp$n_missing_before))

    # =========================================================================
    # 5. SIGNAL DRIFT AND BATCH CORRECTION
    # =========================================================================

    msg("Step 5/10: Drift/batch correction...")

    uncorrected_features <- character(0L)

    if (correct_drift) {
      required_meta <- c(batch_col, injection_col, group_col)
      missing_meta  <- setdiff(required_meta, colnames(metadata))
      if (length(missing_meta) > 0L) {
        warning("Skipping drift correction: metadata column(s) not found: ",
                paste(missing_meta, collapse = ", "), call. = FALSE)
        correct_drift <- FALSE
      }
    }

    if (correct_drift) {
      corr_result <- run_driftBatchCorrect(
        x                   = df,
        metadata            = metadata,
        perform_correction  = TRUE,
        batch_corr_only     = FALSE,
        injection_sequence  = injection_col,
        batch_numbers       = batch_col,
        groups              = group_col,
        qc_label            = "QC",
        qc_types            = qc_types,
        spline_smooth_param = spline_smooth,
        min_QC              = min_qc_per_batch,
        spar_limit          = spline_spar_limit,
        log_scale           = correct_on_log,
        verbose             = FALSE
      )

      uncorrected_features <- corr_result$uncorrected_features

      if (remove_uncorrected && length(uncorrected_features) > 0L) {
        corr_df <- corr_result$data
        keep    <- setdiff(colnames(corr_df), uncorrected_features)
        df      <- corr_df[, keep, drop = FALSE]
        msg(sprintf("  Removed %d uncorrected feature(s). %d remaining.",
                    length(uncorrected_features), ncol(df)))
      } else {
        df <- corr_result$data
        if (length(uncorrected_features) > 0L)
          msg(sprintf("  %d feature(s) could not be corrected (retained).",
                      length(uncorrected_features)))
      }

      na_before_corr <- sum(is.na(as.matrix(shared$data_imputed)))
      na_after_corr  <- sum(is.na(as.matrix(df)))
      new_na         <- na_after_corr - na_before_corr

      if (na_after_corr > 0L) {
        imp_r <- run_mvimpute(x = df, method = impute_fraction, verbose = FALSE)
        df    <- imp_r$data
        msg(sprintf(
          "  Re-imputed %d NA(s) after QCRSC (%d newly introduced, %d pre-existing).",
          na_after_corr, max(0L, new_na), max(0L, -new_na)
        ))
      }

      msg(sprintf("  QCRSC applied (%s).", corr_result$method_used))
    } else {
      msg("  Skipped (correct_drift = FALSE).")
    }

    shared$data_corrected       <- df
    shared$uncorrected_features <- uncorrected_features
    shared$dimensions           <- .record_dim(shared$dimensions, "After drift correction", df)

  }, error = function(e) {
    message("run_DIpreprocess ERROR (Steps 1-5): ", e$message)
    shared$error <<- e$message
  })

  # If early steps failed, return a minimal error object
  if (!is.null(shared$error)) {
    if (is_multi) {
      out <- structure(list(error = shared$error), class = "run_DIpreprocess_multi")
    } else {
      out <- structure(
        c(shared, list(
          data_normalized  = NULL, data_transformed = NULL,
          data_nonpls      = NULL, data_pls         = NULL,
          data_nonpls_merged = NULL, data_pls_merged  = NULL,
          metadata_merged  = NULL, features_final   = character(0L),
          elapsed_seconds  = proc.time()[["elapsed"]] - t_start,
          parameters       = as.list(match.call()[-1L])
        )),
        class = c("run_DIpreprocess", "list")
      )
    }
    return(out)
  }

  # ===========================================================================
  # STEPS 6-10: PER-COMBINATION HELPER
  # ===========================================================================

  .run_branch <- function(norm_m, trans_m, combo_label) {

    branch_out <- list(
      metadata              = shared$metadata,
      data_raw              = shared$data_raw,
      data_missing_filtered = shared$data_missing_filtered,
      data_imputed          = shared$data_imputed,
      data_corrected        = shared$data_corrected,
      data_normalized       = NULL,
      data_transformed      = NULL,
      data_nonpls           = NULL,
      data_pls              = NULL,
      data_nonpls_merged    = NULL,
      data_pls_merged       = NULL,
      metadata_merged       = NULL,
      features_final        = character(0L),
      uncorrected_features  = shared$uncorrected_features,
      dimensions            = shared$dimensions,
      parameters            = list(
        normalize_method = norm_m,
        transform_method = trans_m,
        scale_nonpls     = scale_nonpls,
        scale_pls        = scale_pls,
        rsd_threshold    = rsd_threshold,
        rsd_qc_type      = rsd_qc_type,
        variance_percentile = variance_percentile,
        scale_filter_ref = scale_filter_ref,
        merge_replicates = merge_replicates
      ),
      elapsed_seconds = NA_real_,
      error           = NULL
    )

    t_branch <- proc.time()[["elapsed"]]
    df       <- shared$data_corrected
    metadata <- as.data.frame(shared$metadata)   # use snapshotted metadata
    meta_snap <- shared$metadata
    qc_rows   <- meta_snap$Group_ == "QC"

    tryCatch({

      # =======================================================================
      # 6. NORMALIZATION
      # =======================================================================

      msg(sprintf("  [%s] Step 6: Normalizing (%s)...", combo_label, norm_m))

      if (norm_m == "none") {
        df_norm <- df
      } else {
        norm_result <- run_normalize(
          x              = df,
          metadata       = metadata,
          method         = norm_m,
          factor_col     = norm_factor_col,
          ref_sample     = normalize_ref_sample,
          qc_normalize   = normalize_qc_method,
          groups         = "Group_",
          qc_types       = "QC",   # meta_snap$Group_ is already recoded to "QC"
          sample_id_col  = "Sample",
          verbose        = FALSE
        )
        df_norm <- norm_result$data
      }

      # =======================================================================
      # 7. TRANSFORMATION
      # =======================================================================

      msg(sprintf("  [%s] Step 7: Transforming (%s)...", combo_label, trans_m))

      if (trans_m == "none") {
        df_trans <- df_norm
      } else {
        trans_result <- run_transform(
          x         = df_norm,
          method    = trans_m,
          metadata  = metadata,
          groups    = "Group_",   # already recoded: all QC types → "QC"
          qc_types  = "QC",
          num_cores = vsn_cores,
          verbose   = FALSE
        )
        df_trans <- trans_result$data
      }

      # =======================================================================
      # 8. SCALING (two parallel branches)
      # =======================================================================

      msg(sprintf("  [%s] Step 8: Scaling (NONPLS='%s', PLS='%s')...",
                  combo_label, scale_nonpls, scale_pls))

      .scale_branch <- function(data, method) {
        if (method == "none") return(as.data.frame(data))
        run_scale(x = data, method = method, verbose = FALSE)$data
      }

      df_nonpls <- .scale_branch(df_trans, scale_nonpls)
      df_pls    <- .scale_branch(df_trans, scale_pls)

      # =======================================================================
      # 9. QUALITY FILTERING
      # =======================================================================

      msg(sprintf(
        "  [%s] Step 9: Quality filtering (RSD<=%.0f%%, var percentile=%g, ref='%s')...",
        combo_label, rsd_threshold * 100, variance_percentile, scale_filter_ref
      ))

      .filter_branch <- function(data, branch_name) {
        # Inside .run_branch, metadata$Group_ has already collapsed all
        # user-supplied qc_types (e.g. "EQC", "SQC") to the single token "QC".
        # Both qc_type and qc_types must therefore use "QC" here so that
        # run_filterRSD's internal validation does not reject the label.
        rsd_res  <- run_filterRSD(
          x        = data,
          metadata = metadata,
          max_rsd  = rsd_threshold,
          groups   = "Group_",
          qc_type  = "QC",
          qc_types = "QC",
          verbose  = FALSE
        )
        data_rsd <- rsd_res$data
        if (is.logical(data_rsd) || is.null(data_rsd) || is.null(ncol(data_rsd)))
          data_rsd <- data[, 0L, drop = FALSE]

        msg(sprintf("    [%s][%s] RSD: removed %d. %d remaining.",
                    combo_label, branch_name, rsd_res$n_features_removed, ncol(data_rsd)))

        if (ncol(data_rsd) == 0L) {
          msg(sprintf("    [%s][%s] Variance: skipped (0 features).",
                      combo_label, branch_name))
          return(list(rsd = data_rsd, final = data_rsd))
        }

        var_res  <- run_filtervariance(
          x          = data_rsd,
          percentile = variance_percentile,
          verbose    = FALSE
        )
        data_var <- var_res$data
        if (is.logical(data_var) || is.null(data_var) || is.null(ncol(data_var)))
          data_var <- data_rsd[, 0L, drop = FALSE]

        msg(sprintf("    [%s][%s] Variance: removed %d. %d remaining.",
                    combo_label, branch_name, var_res$n_features_removed, ncol(data_var)))
        list(rsd = data_rsd, final = data_var)
      }

      if (scale_filter_ref == "auto") {
        filt_nonpls <- .filter_branch(df_nonpls, "NONPLS")
        filt_pls    <- .filter_branch(df_pls,    "PLS")

        branch_out$dimensions <- .record_dim(branch_out$dimensions,
          sprintf("[%s] After NONPLS RSD filter", combo_label), filt_nonpls$rsd)
        branch_out$dimensions <- .record_dim(branch_out$dimensions,
          sprintf("[%s] After NONPLS variance filter", combo_label), filt_nonpls$final)
        branch_out$dimensions <- .record_dim(branch_out$dimensions,
          sprintf("[%s] After PLS RSD filter", combo_label), filt_pls$rsd)
        branch_out$dimensions <- .record_dim(branch_out$dimensions,
          sprintf("[%s] After PLS variance filter", combo_label), filt_pls$final)

        features_final <- intersect(colnames(filt_nonpls$final), colnames(filt_pls$final))
        if (length(features_final) == 0L)
          stop(sprintf(
            "[%s] No features passed quality filtering in both branches. ",
            "Consider relaxing 'rsd_threshold' or 'variance_percentile'.", combo_label
          ))

      } else if (scale_filter_ref == "NONPLS") {
        filt_nonpls <- .filter_branch(df_nonpls, "NONPLS")

        branch_out$dimensions <- .record_dim(branch_out$dimensions,
          sprintf("[%s] After NONPLS RSD filter", combo_label), filt_nonpls$rsd)
        branch_out$dimensions <- .record_dim(branch_out$dimensions,
          sprintf("[%s] After NONPLS variance filter (ref)", combo_label), filt_nonpls$final)

        features_final <- colnames(filt_nonpls$final)
        if (length(features_final) == 0L)
          stop(sprintf(
            "[%s] No features passed quality filtering (NONPLS branch).", combo_label
          ))

      } else {  # "PLS"
        filt_pls <- .filter_branch(df_pls, "PLS")

        branch_out$dimensions <- .record_dim(branch_out$dimensions,
          sprintf("[%s] After PLS RSD filter", combo_label), filt_pls$rsd)
        branch_out$dimensions <- .record_dim(branch_out$dimensions,
          sprintf("[%s] After PLS variance filter (ref)", combo_label), filt_pls$final)

        features_final <- colnames(filt_pls$final)
        if (length(features_final) == 0L)
          stop(sprintf(
            "[%s] No features passed quality filtering (PLS branch).", combo_label
          ))
      }

      msg(sprintf("  [%s] Final feature set: %d feature(s).", combo_label, length(features_final)))

      branch_out$data_nonpls     <- df_nonpls[, features_final, drop = FALSE]
      branch_out$data_pls        <- df_pls[,    features_final, drop = FALSE]
      branch_out$data_normalized <- as.data.frame(df_norm)[,  features_final, drop = FALSE]
      branch_out$data_transformed <- as.data.frame(df_trans)[, features_final, drop = FALSE]
      branch_out$features_final  <- features_final

      branch_out$dimensions <- .record_dim(branch_out$dimensions,
        sprintf("[%s] Final (NONPLS)", combo_label), branch_out$data_nonpls)
      branch_out$dimensions <- .record_dim(branch_out$dimensions,
        sprintf("[%s] Final (PLS)", combo_label), branch_out$data_pls)

      # =======================================================================
      # 10. TECHNICAL-REPLICATE MERGING
      # =======================================================================

      msg(sprintf("  [%s] Step 10: Replicate merging...", combo_label))

      sid <- meta_snap$SubjectID
      do_merge <- merge_replicates &&
        subject_id_col %in% colnames(metadata) &&
        !all(is.na(sid) | sid == "" | sid == "NA")

      if (do_merge) {
        is_qc     <- qc_rows
        has_sid   <- !is.na(sid) & sid != "" & sid != "NA"
        mergeable <- !is_qc & has_sid

        unique_sids <- unique(sid[mergeable])
        has_dups    <- any(tabulate(match(sid[mergeable], unique_sids)) > 1L)

        if (!has_dups) {
          msg(sprintf("  [%s] No technical replicates detected; skipping merge.",
                      combo_label))
        } else {
          .merge_br <- function(data) {
            bio_data  <- data[mergeable,  , drop = FALSE]
            keep_data <- data[!mergeable, , drop = FALSE]
            bio_sid   <- sid[mergeable]
            grp_idx   <- split(seq_len(nrow(bio_data)), bio_sid)
            merged_rows <- do.call(rbind, lapply(grp_idx, function(idx)
              colMeans(bio_data[idx, , drop = FALSE], na.rm = TRUE)))
            merged_df           <- as.data.frame(merged_rows)
            rownames(merged_df) <- names(grp_idx)
            rbind(merged_df, keep_data)
          }

          branch_out$data_nonpls_merged <- .merge_br(branch_out$data_nonpls)
          branch_out$data_pls_merged    <- .merge_br(branch_out$data_pls)

          mergeable_meta <- meta_snap[mergeable, , drop = FALSE]
          keep_meta      <- meta_snap[!mergeable, , drop = FALSE]
          grp_idx_meta   <- split(seq_len(nrow(mergeable_meta)), sid[mergeable])
          merged_meta_rows <- do.call(rbind, lapply(grp_idx_meta, function(idx) {
            r <- mergeable_meta[idx[1L], , drop = FALSE]
            r$InjectionSequence <- min(mergeable_meta$InjectionSequence[idx], na.rm = TRUE)
            r
          }))
          rownames(merged_meta_rows) <- names(grp_idx_meta)

          branch_out$metadata_merged <- rbind(merged_meta_rows, keep_meta)

          n_before <- sum(mergeable)
          n_after  <- nrow(merged_meta_rows)
          msg(sprintf("  [%s] Merged %d rows into %d unique subject(s).",
                      combo_label, n_before, n_after))

          branch_out$dimensions <- .record_dim(branch_out$dimensions,
            sprintf("[%s] After replicate merge (NONPLS)", combo_label),
            branch_out$data_nonpls_merged)
          branch_out$dimensions <- .record_dim(branch_out$dimensions,
            sprintf("[%s] After replicate merge (PLS)", combo_label),
            branch_out$data_pls_merged)
        }

      } else {
        if (merge_replicates && !do_merge)
          msg(sprintf("  [%s] merge_replicates=TRUE but no valid SubjectID; skipping.",
                      combo_label))
        else
          msg(sprintf("  [%s] Skipped (merge_replicates = FALSE).", combo_label))
      }

    }, error = function(e) {
      message(sprintf("run_DIpreprocess ERROR in combination '%s': %s",
                      combo_label, e$message))
      branch_out$error <<- e$message
    })

    branch_out$elapsed_seconds <- proc.time()[["elapsed"]] - t_branch
    structure(branch_out, class = c("run_DIpreprocess", "list"))
  }

  # ===========================================================================
  # DISPATCH: single vs. multi
  # ===========================================================================

  if (!is_multi) {

    # ------------------------------------------------------------------
    # Single combination — original return shape
    # ------------------------------------------------------------------
    msg("Step 6-10: Normalization → Transformation → Scaling → Filtering → Merge...")

    out <- .run_branch(normalize_method, transform_method, combo_names)

    # Back-fill parameters with the full original call
    out$parameters <- as.list(match.call()[-1L])

    out$elapsed_seconds <- proc.time()[["elapsed"]] - t_start
    msg(sprintf("Pipeline completed successfully. Elapsed: %.1f second(s).",
                out$elapsed_seconds))

    return(out)

  } else {

    # ------------------------------------------------------------------
    # Multiple combinations — named list of run_DIpreprocess objects
    # ------------------------------------------------------------------
    results <- vector("list", nrow(combos))
    names(results) <- combo_names

    for (i in seq_len(nrow(combos))) {
      nm  <- combo_names[i]
      msg(sprintf("\n--- Combination %d/%d: %s ---", i, nrow(combos), nm))
      results[[nm]] <- .run_branch(combos$norm[i], combos$trans[i], nm)
    }

    elapsed_total <- proc.time()[["elapsed"]] - t_start

    # Report any failures without aborting the whole return
    failed <- vapply(results, function(r) !is.null(r$error), logical(1L))
    if (any(failed))
      warning(sprintf(
        "%d combination(s) encountered errors: %s\n  Inspect $error on each element for details.",
        sum(failed), paste(names(results)[failed], collapse = ", ")
      ), call. = FALSE)

    msg(sprintf(
      "\nAll combinations completed. Total elapsed: %.1f second(s).", elapsed_total
    ))

    structure(
      results,
      class          = c("run_DIpreprocess_multi", "list"),
      elapsed_seconds = elapsed_total,
      combos         = combos,
      shared         = shared
    )
  }
}


# ==============================================================================
# S3 METHODS — run_DIpreprocess  (single-combination; unchanged)
# ==============================================================================

#' @export
print.run_DIpreprocess <- function(x, ...) {
  cat("=== run_DIpreprocess pipeline result ===\n")
  if (!is.null(x$error)) {
    cat(sprintf("STATUS : FAILED — %s\n", x$error))
    return(invisible(x))
  }
  cat("STATUS : OK\n")
  n_samp <- if (!is.null(x$data_nonpls)) nrow(x$data_nonpls) else NA_integer_
  n_feat <- length(x$features_final)
  cat(sprintf("Data   : %d samples x %d features (final)\n", n_samp, n_feat))
  if (!is.null(x$data_nonpls_merged))
    cat(sprintf("Merged : %d samples x %d features (replicate-merged)\n",
                nrow(x$data_nonpls_merged), ncol(x$data_nonpls_merged)))
  if (length(x$uncorrected_features) > 0L)
    cat(sprintf("QCRSC  : %d uncorrected feature(s) retained\n",
                length(x$uncorrected_features)))
  cat(sprintf("Time   : %.1f second(s)\n", x$elapsed_seconds))
  invisible(x)
}

#' @export
summary.run_DIpreprocess <- function(object, ...) {
  cat("=== run_DIpreprocess — step-by-step dimensions ===\n")
  print(object$dimensions, row.names = FALSE)
  if (!is.null(object$data_nonpls) && !is.null(object$data_pls)) {
    cat("\nFinal datasets:\n")
    cat(sprintf("  data_nonpls      : %d x %d\n",
                nrow(object$data_nonpls), ncol(object$data_nonpls)))
    cat(sprintf("  data_pls         : %d x %d\n",
                nrow(object$data_pls), ncol(object$data_pls)))
    cat(sprintf("  data_normalized  : %d x %d (pre-scale)\n",
                nrow(object$data_normalized), ncol(object$data_normalized)))
    cat(sprintf("  data_transformed : %d x %d (pre-scale)\n",
                nrow(object$data_transformed), ncol(object$data_transformed)))
    if (!is.null(object$data_nonpls_merged))
      cat(sprintf("  data_nonpls_merged: %d x %d\n",
                  nrow(object$data_nonpls_merged), ncol(object$data_nonpls_merged)))
  }
  invisible(object)
}


# ==============================================================================
# S3 METHODS — run_DIpreprocess_multi  (new)
# ==============================================================================

#' @export
print.run_DIpreprocess_multi <- function(x, ...) {
  cat("=== run_DIpreprocess_multi pipeline result ===\n")

  if (!is.null(x$error)) {
    cat(sprintf("STATUS : FAILED — %s\n", x$error))
    return(invisible(x))
  }

  elapsed <- attr(x, "elapsed_seconds")
  combos  <- attr(x, "combos")

  cat(sprintf("Combinations : %d\n", length(x)))
  cat(sprintf("Names        : %s\n", paste(names(x), collapse = ", ")))
  if (!is.null(elapsed))
    cat(sprintf("Total time   : %.1f second(s)\n", elapsed))

  failed <- vapply(x, function(r) !is.null(r$error), logical(1L))
  if (any(failed))
    cat(sprintf("FAILED       : %s\n", paste(names(x)[failed], collapse = ", ")))

  cat("\nPer-combination summary:\n")
  for (nm in names(x)) {
    r <- x[[nm]]
    if (!is.null(r$error)) {
      cat(sprintf("  %-20s  FAILED: %s\n", nm, r$error))
    } else {
      n_samp <- if (!is.null(r$data_nonpls)) nrow(r$data_nonpls) else NA_integer_
      n_feat <- length(r$features_final)
      cat(sprintf("  %-20s  %d samples x %d features\n", nm, n_samp, n_feat))
    }
  }
  invisible(x)
}

#' @export
summary.run_DIpreprocess_multi <- function(object, ...) {
  cat("=== run_DIpreprocess_multi — per-combination summaries ===\n\n")
  for (nm in names(object)) {
    cat(sprintf("--- %s ---\n", nm))
    summary(object[[nm]])
    cat("\n")
  }
  invisible(object)
}