#' Direct-Injection Metabolomics Preprocessing Pipeline
#'
#' @description
#' Runs a complete, sequential preprocessing pipeline for metabolomics data by
#' orchestrating the individual `run_*` functions of \pkg{dimsprepr}. Steps
#' execute in a fixed, reproducible order regardless of argument arrangement.
#'
#' When `normalize_method`, `transform_method`, `scale_nonpls`, and/or
#' `scale_pls` are supplied as character vectors with more than one element,
#' \strong{all combinations} are computed with three levels of deduplication:
#' (1) the steps prior to normalization (validation, outlier removal,
#' missing-value filtering, imputation, and drift/batch correction) run
#' \strong{only once}; (2) normalization, transformation, and quality filtering
#' run \strong{once per unique normalize x transform pair}; (3) scaling and
#' replicate merging run \strong{once per full combination}. The return value
#' in this case is a list of class `"run_DIpreprocess_multi"` whose elements
#' are standard `run_DIpreprocess` objects (see Value for their names).
#'
#' @param x Matrix or data frame. Numeric feature data with \strong{samples in rows}
#'   and \strong{features in columns}. All values should be non-negative raw intensities;
#'   zeros are treated as missing internally.
#'   For `print()`, a `run_DIpreprocess` or `run_DIpreprocess_multi` object.
#' @param metadata Data frame. Sample metadata with `nrow(metadata) == nrow(x)`,
#'   in the same row order as `x`. At minimum, the group column (default
#'   `"Group"`) must be present.
#' @param sample_id_col Character. Column in `metadata` containing unique
#'   sample identifiers; they become the row names of every output. Default:
#'   `"Sample"`.
#' @param group_col Character. Column in `metadata` containing group labels.
#'   Default: `"Group"`.
#' @param qc_types Character vector. Group labels that identify QC samples.
#'   Default: `c("QC", "SQC", "EQC")`.
#' @param batch_col Character. Column in `metadata` containing batch labels
#'   (numeric or character). Required for drift correction. Default: `"Batch"`.
#' @param injection_col Character. Column in `metadata` containing the numeric
#'   injection order (unique values). Required for drift correction. Default:
#'   `"InjectionSequence"`.
#' @param norm_factor_col Character or `NULL`. Column in `metadata`
#'   containing sample-specific normalization factors (e.g. creatinine
#'   concentrations or osmolality for urine). Default: `NULL`.
#'
#'   \strong{Auto-detection behaviour:} When a column name is supplied the
#'   function validates it \emph{before} running any preprocessing:
#'   \itemize{
#'     \item Column must exist in `metadata` and be numeric.
#'     \item Every \emph{non-QC} sample must have a non-`NA` value
#'       (QC samples are permitted to be `NA`).
#'   }
#'   If all checks pass, `normalize_method` is \strong{automatically
#'   overridden} to `"specific_factor"` and the supplied column values
#'   are used for normalization; any value passed in `normalize_method`
#'   is ignored with a message. If the column is absent or not numeric, a
#'   warning is issued and `normalize_method` is used unchanged.
#' @param subject_id_col Character. Column in `metadata` containing
#'   subject identifiers for technical-replicate merging. Default: `"SubjectID"`.
#' @param outliers Character vector or `NULL`. Sample identifiers (row names
#'   after applying `sample_id_col`) to remove before processing. Default: `NULL`.
#' @param missing_threshold Numeric \[0, 1\]. Passed to [run_filtermissing()]
#'   as `threshold`. Default: `0.2`.
#' @param missing_by_group Logical. If `TRUE`, assess missingness per group.
#'   Default: `TRUE`.
#' @param missing_any_group Logical. Only used when `missing_by_group = TRUE`.
#'   If `TRUE` (default), a feature is retained if its missingness is below
#'   the threshold in at least one group (Modified 80% Rule). If `FALSE`,
#'   it must satisfy the threshold in all groups.
#' @param missing_include_qc Logical. If `TRUE`, include QC samples when
#'   assessing missingness. Default: `FALSE`.
#' @param impute_method Imputation method passed to [run_mvimpute()] as
#'   `method`: a fraction (e.g. `0.2` = one fifth of the smallest positive
#'   value per feature), `"quantileregression"`, a \pkg{mice} method, or
#'   `"none"`. Default: `0.2`.
#' @param impute_method_post_drift Imputation method for any missing values
#'   remaining (or introduced) after drift correction. If `NULL` (default),
#'   takes on the value of `impute_method`.
#' @param positive_only Logical. If `TRUE`, the minimum used by fractional
#'   imputation only considers positive values. Default: `TRUE`.
#' @param correct_drift Logical. If `TRUE`, apply QCRSC drift and batch
#'   correction via [run_driftBatchCorrect()]. Default: `TRUE`.
#' @param remove_uncorrected Logical. If `TRUE`, remove features QCRSC
#'   could not correct. Default: `FALSE`.
#' @param spline_smooth Numeric \[0, 1\]. QCRSC spline smoothing parameter;
#'   `0` selects it by cross-validation. Default: `0`.
#' @param spline_spar_limit Numeric vector of length 2. Lower/upper bounds of
#'   the smoothing parameter searched during cross-validation.
#'   Default: `c(-1.5, 1.5)`.
#' @param correct_on_log Logical. If `TRUE`, fit the QC spline on
#'   log-transformed data. Default: `TRUE`.
#' @param min_qc_per_batch Integer. Minimum QC samples per batch for a feature
#'   to be corrected. Default: `5`.
#' @param normalize_method Character or character vector. One or more
#'   normalization methods (see [run_normalize()] for details). When multiple
#'   methods are supplied every method is paired with every element of
#'   `transform_method` (Cartesian product).
#'   \itemize{
#'     \item `"sum"`: Normalizes by sum.
#'     \item `"median"`: Normalizes by median.
#'     \item `"specific_factor"`: Uses external factors in `norm_factor_col`.
#'     \item `"pqn_global"`: PQN using the median of all samples as reference.
#'     \item `"pqn_reference"`: PQN using `normalize_ref_sample` as reference.
#'     \item `"pqn_group"`: PQN using pooled QC samples as reference.
#'     \item `"pqn"`, `"quantile"`: PQN via \pkg{pmp} (`"quantile"` is an
#'       alias kept for backward compatibility; it is not quantile normalization).
#'     \item `"col_rel_abundance"`, `"row_rel_abundance"`: relative abundance.
#'     \item `"none"`: No normalization.
#'   }
#'   Default: `"median"`.
#' @param normalize_ref_sample Character. Reference sample identifier for
#'   `"pqn_reference"`. Default: `NULL`.
#' @param normalize_qc_method Character. How QC samples are scaled when
#'   sample factors are used: `"mean"`, `"median"`, or `"none"`. Default: `"none"`.
#' @param transform_method Character or character vector. One or more
#'   transformation methods (see [run_transform()]): `"log2"`, `"log10"`,
#'   `"sqrt"`, `"cbrt"`, `"clr"`, `"arcsin_sqrt"`, `"vsn"`, `"glog"`, or
#'   `"none"`. When multiple methods are supplied every method is paired with
#'   every element of `normalize_method` (Cartesian product).
#'   Default: `"log10"`.
#' @param scale_nonpls Character or character vector. One or more scaling
#'   methods (see [run_scale()]) for the NONPLS branch, intended for PCA,
#'   clustering, and univariate tests. When multiple methods are supplied,
#'   each is paired with every element of `scale_pls`. Only the varying scale
#'   components are appended to the combination name (e.g.
#'   `median_log10_auto` when only `scale_nonpls` varies). Default: `"auto"`.
#' @param scale_pls Character or character vector. One or more scaling
#'   methods for the PLS branch, intended for PLS-DA/OPLS-DA. Default:
#'   `"pareto"`.
#' @param rsd_threshold Numeric \[0, 1\]. Features with QC RSD at or above this
#'   value are removed ([run_filterRSD()]). Default: `0.3`.
#' @param rsd_qc_type Character. QC type used for the RSD filter: `"QC"` (all
#'   `qc_types`), or one specific label such as `"EQC"` or `"SQC"`.
#'   Default: `"QC"`.
#' @param variance_percentile Numeric \[0, 100\]. Features in this bottom
#'   percentile of variance are removed ([run_filtervariance()]). Default: `10`.
#' @param merge_replicates Logical. Average technical replicates (non-QC rows
#'   sharing a `subject_id_col` value). Default: `FALSE`.
#' @param track Character vector or `NULL`. Feature names to monitor across
#'   preprocessing steps. If a tracked feature is removed, the step and reason
#'   (e.g., missingness percentage, RSD value) are recorded in `$track`.
#'   Default: `NULL`.
#' @param verbose Logical. Print progress. Default: `TRUE`.
#'
#' @return
#' \describe{
#'   \item{Single combination (all four method arguments are scalar)}{A
#'     list of class `run_DIpreprocess` with elements:
#'     \describe{
#'       \item{`data_nonpls`, `data_pls`}{Final data for the NONPLS and PLS
#'         branches (scaled with `scale_nonpls` / `scale_pls`).}
#'       \item{`data_raw`, `data_missing_filtered`, `data_imputed`,
#'         `data_corrected`, `data_normalized`, `data_transformed`}{Data
#'         after each step (`data_normalized`/`data_transformed` are
#'         restricted to `features_final`).}
#'       \item{`data_nonpls_merged`, `data_pls_merged`, `metadata_merged`}{
#'         Replicate-merged data and metadata, or `NULL`.}
#'       \item{`metadata`}{Standardized metadata with columns `Sample`,
#'         `Group`, `Group_` (QC types collapsed to `"QC"`), `Batch`,
#'         `InjectionSequence`, `SubjectID`.}
#'       \item{`features_final`}{Features that passed every filter.}
#'       \item{`uncorrected_features`}{Features QCRSC could not correct.}
#'       \item{`dimensions`}{Data frame of samples x features after each step.}
#'       \item{`track`}{Data frame recording why tracked features were removed.}
#'       \item{`parameters`}{List of parameter values used.}
#'       \item{`elapsed_seconds`, `error`}{Run time, and the error message
#'         (`NULL` on success).}
#'     }}
#'   \item{Multiple combinations (any method argument has length > 1)}{A
#'     named list of class `run_DIpreprocess_multi`. Each element is a
#'     `run_DIpreprocess` object. Element names follow the pattern
#'     `<normalize>_<transform>` when only normalization/transformation vary
#'     (e.g. `sum_log10`), or include scale suffixes when scaling varies
#'     (e.g. `median_log10_auto`, `median_log10_auto_pareto`).
#'     Scale components are appended only for those arguments that have more
#'     than one unique value: `scale_nonpls` first, then `scale_pls`.
#'     The shared pre-normalization data (post-drift-correction) is stored in
#'     the `shared` attribute of the returned list.}
#' }
#'
#' @details
#' \strong{Preprocessing Workflow (Fixed Order):}
#' \enumerate{
#'   \item Input validation
#'   \item Outlier removal
#'   \item Missing-value filtering ([run_filtermissing()])
#'   \item Missing-value imputation ([run_mvimpute()])
#'   \item Signal drift and batch correction ([run_driftBatchCorrect()])
#'   \item \emph{\[Per normalize x transform pair\]} Normalization
#'     ([run_normalize()])
#'   \item \emph{\[Per normalize x transform pair\]} Transformation
#'     ([run_transform()])
#'   \item \emph{\[Per normalize x transform pair\]} Quality filtering:
#'     QC-RSD on the normalized data ([run_filterRSD()]), then variance on the
#'     transformed data ([run_filtervariance()])
#'   \item \emph{\[Per full combination\]} Scaling ([run_scale()]) of the
#'     NONPLS and PLS branches
#'   \item \emph{\[Per full combination\]} Technical-replicate merging
#' }
#'
#' Quality filtering is done \emph{before} scaling because the RSD of
#' mean-centred data and the variance of auto-scaled data (always 1) carry no
#' information. Both branches therefore share the same feature set,
#' `features_final`.
#'
#' The workflow follows the QC-based guidelines of Broadhurst et al. (2018)
#' and the DIMS batch-correction strategy of Kirwan et al. (2013).
#'
#' @note
#' \strong{Zero Handling:} The pipeline automatically converts all `0` values in
#' `x` to `NA` prior to analysis.
#'
#' \strong{Pre-flight checks:} Method names, and method-specific requirements
#' (e.g. `norm_factor_col` for `"specific_factor"`, `normalize_ref_sample` for
#' `"pqn_reference"`), are validated \emph{before} any preprocessing step, so
#' errors surface early rather than mid-run. Errors raised inside a step are
#' caught: the function returns an object whose `$error` holds the message.
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
#' Kirwan, J.A., Broadhurst, D.I., Davidson, R.L., & Viant, M.R. (2013).
#' Characterising and correcting batch variation in an automated direct
#' infusion mass spectrometry (DIMS) metabolomics workflow.
#' *Analytical and Bioanalytical Chemistry*, 405(15), 5147-5157.
#' \doi{10.1007/s00216-013-6856-7}
#'
#' van den Berg, R.A., Hoefsloot, H.C.J., Westerhuis, J.A., Smilde, A.K., &
#' van der Werf, M.J. (2006). Centering, scaling, and transformations:
#' improving the biological information content of metabolomics data.
#' *BMC Genomics*, 7, 142. \doi{10.1186/1471-2164-7-142}
#'
#' See the individual step functions for method-specific references.
#'
#' @family preprocessing steps
#' @seealso [plot_beforeafter()] and [plot_dist_beforeafter()] accept a
#'   `run_DIpreprocess` result directly; [dims_data] and [dims_metadata] are
#'   example inputs.
#'
#'
#' @export
#'
#' @examples
#' data(dims_data)
#' data(dims_metadata)
#'
#' res <- run_DIpreprocess(dims_data, dims_metadata, correct_drift = FALSE,
#'                         verbose = FALSE)
#' res
#' summary(res)
#' dim(res$data_nonpls)
#'
#' # Full default pipeline, including QCRSC drift correction (slower: loads pmp)
#' \donttest{
#' res_drift <- run_DIpreprocess(dims_data, dims_metadata, verbose = FALSE)
#' res_drift$parameters$correct_drift
#' }
#'
#' # Multi-method mode: steps 1-5 run once, then every combination
#' multi <- run_DIpreprocess(
#'   dims_data, dims_metadata,
#'   normalize_method = c("median", "pqn_group"),
#'   transform_method = "log10",
#'   scale_nonpls     = c("auto", "pareto"),
#'   correct_drift    = FALSE,
#'   verbose          = FALSE
#' )
#' names(multi)
#' multi
run_DIpreprocess <- function(
    x,
    metadata,

    # - Column names ----------------------------------------------------------
    sample_id_col    = "Sample",
    group_col        = "Group",
    qc_types         = c("QC", "SQC", "EQC"),
    batch_col        = "Batch",
    injection_col    = "InjectionSequence",
    norm_factor_col  = NULL,
    subject_id_col   = "SubjectID",

    # - Outlier removal -------------------------------------------------------
    outliers         = NULL,

    # - Missing-value filter --------------------------------------------------
    missing_threshold  = 0.2,
    missing_by_group   = TRUE,
    missing_any_group  = TRUE,
    missing_include_qc = FALSE,

    # - Missing-value imputation ----------------------------------------------
    impute_method  = 0.2,
    impute_method_post_drift = NULL,
    positive_only    = TRUE,

    # - Drift / batch correction ----------------------------------------------
    correct_drift        = TRUE,
    remove_uncorrected   = FALSE,
    spline_smooth        = 0,
    spline_spar_limit    = c(-1.5, 1.5),
    correct_on_log       = TRUE,
    min_qc_per_batch     = 5L,

    # - Normalization ---------------------------------------------------------
    normalize_method     = "median",
    normalize_ref_sample = NULL,
    normalize_qc_method  = "none",

    # - Transformation --------------------------------------------------------
    transform_method = "log10",

    # - Scaling ---------------------------------------------------------------
    scale_nonpls     = "auto",
    scale_pls        = "pareto",

    # - Quality filtering -----------------------------------------------------
    rsd_threshold       = 0.3,
    rsd_qc_type         = "QC",
    variance_percentile = 10,

    # - Replicate merging -----------------------------------------------------
    merge_replicates = FALSE,
    track            = NULL,

    verbose          = TRUE
) {

  t_start <- proc.time()[["elapsed"]]
  msg <- .msg(verbose)

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

  .track_init <- function() {
    data.frame(Feature = character(), Step = character(), Value = character(),
               stringsAsFactors = FALSE)
  }

  # Parameters actually used (evaluated at call time, so e.g. a correct_drift
  # switched off for missing columns is reported as FALSE).
  .params <- function(norm_m, trans_m, scale_np, scale_p) list(
    sample_id_col        = sample_id_col,
    group_col            = group_col,
    qc_types             = qc_types,
    batch_col            = batch_col,
    injection_col        = injection_col,
    norm_factor_col      = norm_factor_col,
    subject_id_col       = subject_id_col,
    outliers             = outliers,
    missing_threshold    = missing_threshold,
    missing_by_group     = missing_by_group,
    missing_any_group    = missing_any_group,
    missing_include_qc   = missing_include_qc,
    impute_method        = impute_method,
    impute_method_post_drift = impute_method_post_drift,
    positive_only        = positive_only,
    correct_drift        = correct_drift,
    remove_uncorrected   = remove_uncorrected,
    spline_smooth        = spline_smooth,
    spline_spar_limit    = spline_spar_limit,
    correct_on_log       = correct_on_log,
    min_qc_per_batch     = min_qc_per_batch,
    normalize_method     = norm_m,
    normalize_ref_sample = normalize_ref_sample,
    normalize_qc_method  = normalize_qc_method,
    transform_method     = trans_m,
    scale_nonpls         = scale_np,
    scale_pls            = scale_p,
    rsd_threshold        = rsd_threshold,
    rsd_qc_type          = rsd_qc_type,
    variance_percentile  = variance_percentile,
    merge_replicates     = merge_replicates,
    track                = track,
    verbose              = verbose
  )

  # ---------------------------------------------------------------------------
  # PRE-FLIGHT: validate method vectors and method-specific requirements
  # ---------------------------------------------------------------------------

  if (!is.null(track) && !is.character(track))
    stop("'track' must be a character vector of feature names.", call. = FALSE)

  # Coerce to lower-case character vectors and deduplicate, preserving order.
  normalize_method <- unique(tolower(as.character(normalize_method)))
  transform_method <- unique(tolower(as.character(transform_method)))
  scale_nonpls     <- unique(tolower(as.character(scale_nonpls)))
  scale_pls        <- unique(tolower(as.character(scale_pls)))

  # Fail fast on unknown method names (vocabularies live in globals.R).
  .check_methods <- function(m, valid, arg) {
    bad <- setdiff(m, valid)
    if (length(bad) > 0L)
      stop(sprintf("Unknown %s: %s. Supported: %s.", arg,
                   paste(bad, collapse = ", "), paste(valid, collapse = ", ")),
           call. = FALSE)
  }
  .check_methods(normalize_method, .normalize_methods, "normalize_method")
  .check_methods(transform_method, .transform_methods, "transform_method")
  .check_methods(scale_nonpls,     .scale_methods,     "scale_nonpls")
  .check_methods(scale_pls,        .scale_methods,     "scale_pls")
  if ("pqn_reference" %in% normalize_method && is.null(normalize_ref_sample))
    stop("normalize_method = 'pqn_reference' requires 'normalize_ref_sample'.", call. = FALSE)

  # Method-specific pre-flight checks -----------------------------------------

  # norm_factor_col auto-detection: validate the column and, if good, override
  # normalize_method to "specific_factor". QC samples may have NA values.
  if (!is.null(norm_factor_col)) {
    if (!norm_factor_col %in% colnames(metadata)) {
      warning(sprintf(
        "norm_factor_col '%s' not found in metadata; ignoring it and using normalize_method = '%s'.",
        norm_factor_col, paste(normalize_method, collapse = ", ")
      ), call. = FALSE)
      norm_factor_col <- NULL
    } else if (!is.numeric(metadata[[norm_factor_col]])) {
      warning(sprintf(
        "norm_factor_col '%s' is not numeric; ignoring it and using normalize_method = '%s'.",
        norm_factor_col, paste(normalize_method, collapse = ", ")
      ), call. = FALSE)
      norm_factor_col <- NULL
    } else {
      # Require all non-QC samples to have a value
      if (group_col %in% colnames(metadata)) {
        is_non_qc  <- !metadata[[group_col]] %in% qc_types
        is_na_norm <- is.na(metadata[[norm_factor_col]])
        if (any(is_non_qc & is_na_norm)) {
          bad_samps <- rownames(metadata)[is_non_qc & is_na_norm]
          if (length(bad_samps) == 0L) bad_samps <- which(is_non_qc & is_na_norm)
          stop(sprintf(
            "norm_factor_col '%s' has missing values (NA) for non-QC sample(s). Check: %s",
            norm_factor_col, paste(utils::head(bad_samps, 5), collapse = ", ")
          ), call. = FALSE)
        }
      }
      # Valid: auto-switch to specific_factor normalization
      if (!all(normalize_method == "specific_factor")) {
        message(sprintf(
          "norm_factor_col '%s' provided and valid: overriding normalize_method ('%s') to 'specific_factor'.",
          norm_factor_col, paste(normalize_method, collapse = ", ")
        ))
        normalize_method <- "specific_factor"
      }
    }
  }

  # Guard: if "specific_factor" was requested explicitly, norm_factor_col must exist
  if ("specific_factor" %in% normalize_method && is.null(norm_factor_col))
    stop("normalize_method = 'specific_factor' requires 'norm_factor_col' to be specified.",
         call. = FALSE)

  # Determine multi-method mode -----------------------------------------------
  is_multi <- length(normalize_method) > 1L || length(transform_method) > 1L ||
              length(scale_nonpls) > 1L     || length(scale_pls) > 1L

  # Build Cartesian product of all four method dimensions
  combos <- expand.grid(
    norm     = normalize_method,
    trans    = transform_method,
    scale_np = scale_nonpls,
    scale_p  = scale_pls,
    stringsAsFactors = FALSE
  )

  # Only append a scale component to the name when it has > 1 unique value,
  # mirroring how norm/trans are always both included.
  name_parts <- list(combos$norm, combos$trans)
  if (length(scale_nonpls) > 1L) name_parts <- c(name_parts, list(combos$scale_np))
  if (length(scale_pls)    > 1L) name_parts <- c(name_parts, list(combos$scale_p))
  combo_names <- do.call(paste, c(name_parts, list(sep = "_")))

  if (is_multi) {
    msg(sprintf(
      "Multi-method mode: %d combination(s) detected - %s.",
      nrow(combos), paste(combo_names, collapse = ", ")
    ))
    msg("Steps 1-5 (through drift correction) will run once; ",
        "Steps 6-8 (norm + transform + quality filter) once per unique norm x transform pair; ",
        "Steps 9-10 (scaling + merging) once per full combination.")
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
    track                = .track_init(),
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
    if (!group_col %in% colnames(metadata))
      stop(sprintf("group_col '%s' not found in metadata.", group_col))

    for (col_name in c(injection_col, norm_factor_col)) {
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
          sample_id_col, paste(utils::head(unique(dups), 5), collapse = ", ")
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
        paste(utils::head(bad_cols, 5), collapse = ", ")
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
      Batch             = .meta_col(batch_col),
      InjectionSequence = suppressWarnings(as.numeric(.meta_col(injection_col))),
      SubjectID         = as.character(.meta_col(subject_id_col)),
      stringsAsFactors  = FALSE,
      row.names         = rownames(metadata)
    )
    # Preserve the normalization factor column (if provided) so that
    # run_normalize(factor_col = norm_factor_col) can find it in meta_snap.
    if (!is.null(norm_factor_col))
      meta_snap[[norm_factor_col]] <- suppressWarnings(
        as.numeric(.meta_col(norm_factor_col))
      )
    shared$metadata <- meta_snap

    # =========================================================================
    # ZEROS -> NA
    # =========================================================================

    n_zeros <- sum(df == 0L, na.rm = TRUE)
    if (n_zeros > 0L) {
      df[df == 0L] <- NA
      msg(sprintf("  Converted %d zero(s) to NA (structural missingness).", n_zeros))
    }

    all_na <- colSums(is.na(df)) == nrow(df)
    if (any(all_na)) {
      if (!is.null(track)) {
        lost_na <- intersect(track, names(all_na)[all_na])
        for (f in lost_na) {
          shared$track <- rbind(shared$track, data.frame(
            Feature = f, Step = "All-NA removal", Value = "100% missing",
            stringsAsFactors = FALSE
          ))
        }
      }
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

    df_before_miss <- df
    filt_miss <- run_filtermissing(
      x               = df,
      metadata        = metadata,
      threshold       = missing_threshold,
      filter_by_group = missing_by_group,
      any_group       = missing_any_group,
      include_QC      = missing_include_qc,
      group_col       = group_col,
      qc_types        = qc_types,
      zero_as_missing = FALSE,
      verbose         = FALSE
    )

    df                        <- filt_miss$data

    if (!is.null(track)) {
      lost_miss <- setdiff(intersect(track, colnames(df_before_miss)), colnames(df))
      for (f in lost_miss) {
        vals <- df_before_miss[[f]]
        gps  <- as.character(metadata[[group_col]])
        unique_gps <- sort(unique(gps))
        miss_info <- vapply(unique_gps, function(g) {
          sub <- vals[gps == g]
          sprintf("%s=%.0f%%", g, (sum(is.na(sub)) / length(sub)) * 100)
        }, character(1L))

        shared$track <- rbind(shared$track, data.frame(
          Feature = f, Step = "Missing value", 
          Value = paste(miss_info, collapse = ", "),
          stringsAsFactors = FALSE
        ))
      }
    }

    shared$data_missing_filtered <- df
    shared$dimensions         <- .record_dim(shared$dimensions, "After missing filter", df)

    msg(sprintf("  Removed %d feature(s) (%.1f%%). %d remaining.",
                filt_miss$n_features_removed,
                filt_miss$n_features_removed / filt_miss$n_features_before * 100,
                ncol(df)))

    # =========================================================================
    # 4. MISSING-VALUE IMPUTATION
    # =========================================================================

    msg(sprintf("Step 4/10: Imputing missing values (method = %s)...",
                as.character(impute_method)))

    imp <- run_mvimpute(
      x             = df,
      method        = impute_method,
      positive_only = positive_only,
      verbose       = FALSE
    )
    df                   <- imp$data
    shared$data_imputed  <- df
    msg(sprintf("  Imputed %d missing value(s).", imp$n_missing_before - imp$n_missing_after))

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
        injection_col       = injection_col,
        batch_col           = batch_col,
        group_col           = group_col,
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
        if (!is.null(track)) {
          lost_drift <- intersect(track, uncorrected_features)
          for (f in lost_drift) {
            shared$track <- rbind(shared$track, data.frame(
              Feature = f, Step = "Drift correction", Value = "Uncorrected",
              stringsAsFactors = FALSE
            ))
          }
        }
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
        # Fall back to first method if second is not explicitly set
        post_method <- if (is.null(impute_method_post_drift)) impute_method else impute_method_post_drift

        if (is.character(post_method) && tolower(trimws(post_method)) == "none") {
          msg(sprintf(
            "  %d NA(s) present after QCRSC (%d newly introduced), but post-drift imputation is set to 'none'. Skipping.",
            na_after_corr, max(0L, new_na)
          ))
        } else {
          imp_r <- run_mvimpute(x = df, method = post_method, positive_only = positive_only, verbose = FALSE)
          df    <- imp_r$data
          msg(sprintf(
            "  Re-imputed %d NA(s) after QCRSC (%d newly introduced, %d pre-existing) using method: %s.",
            na_after_corr, max(0L, new_na), na_after_corr - max(0L, new_na), as.character(post_method)
          ))
        }
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
          parameters       = .params(normalize_method, transform_method,
                                     scale_nonpls, scale_pls)
        )),
        class = c("run_DIpreprocess", "list")
      )
    }
    return(out)
  }

  # ===========================================================================
  # STEPS 6-8: PER NORM x TRANSFORM HELPER (cached across scale combinations)
  # ===========================================================================

  # Normalization -> transformation -> quality filtering. Filtering runs on
  # UNSCALED data: RSD needs non-negative intensities (normalized data) and
  # variance ranking is meaningless after auto-scaling (every variance = 1).
  # Returns df_norm/df_trans restricted to the surviving features.
  .run_norm_trans <- function(norm_m, trans_m) {
    metadata  <- shared$metadata
    nt_label  <- paste(norm_m, trans_m, sep = "_")
    track_add <- .track_init()
    dims_add  <- .empty_dims()
    .add_track <- function(f, step, value)
      track_add <<- rbind(track_add, data.frame(Feature = f, Step = step, Value = value,
                                                stringsAsFactors = FALSE))

    tryCatch({

      # 6. NORMALIZATION
      msg(sprintf("  [%s] Step 6: Normalizing (%s)...", nt_label, norm_m))
      df_norm <- run_normalize(
        x             = shared$data_corrected,
        metadata      = metadata,
        method        = norm_m,
        factor_col    = norm_factor_col,
        ref_sample    = normalize_ref_sample,
        qc_normalize  = normalize_qc_method,
        group_col     = "Group_",
        qc_types      = "QC",   # meta_snap$Group_ is already recoded to "QC"
        sample_id_col = "Sample",
        verbose       = FALSE
      )$data

      # 7. TRANSFORMATION
      msg(sprintf("  [%s] Step 7: Transforming (%s)...", nt_label, trans_m))
      df_trans <- run_transform(
        x         = df_norm,
        method    = trans_m,
        metadata  = metadata,
        group_col = "Group_",
        qc_types  = "QC",
        verbose   = FALSE
      )$data

      # 8. QUALITY FILTERING: QC-RSD on normalized data, variance on transformed
      msg(sprintf(
        "  [%s] Step 8: Quality filtering (QC RSD < %.0f%%, variance percentile = %g)...",
        nt_label, rsd_threshold * 100, variance_percentile
      ))

      rsd_res <- run_filterRSD(
        x         = df_norm,
        metadata  = metadata,
        max_rsd   = rsd_threshold,
        group_col = "Group",      # raw labels, so rsd_qc_type can pick EQC/SQC
        qc_type   = rsd_qc_type,
        qc_types  = qc_types,
        verbose   = FALSE
      )
      for (f in setdiff(intersect(track, colnames(df_norm)), rsd_res$features_kept)) {
        r <- rsd_res$rsd_values[[f]]
        .add_track(f, "RSD filter", if (is.na(r)) "NA" else sprintf("%.1f%%", r * 100))
      }
      dims_add <- .record_dim(dims_add, sprintf("[%s] After RSD filter", nt_label),
                              rsd_res$data)
      msg(sprintf("    RSD: removed %d. %d remaining.",
                  rsd_res$n_features_removed, rsd_res$n_features_after))
      if (rsd_res$n_features_after == 0L)
        stop("No features passed the RSD filter. Consider relaxing 'rsd_threshold'.",
             call. = FALSE)

      var_res <- run_filtervariance(
        x          = df_trans[, rsd_res$features_kept, drop = FALSE],
        percentile = variance_percentile,
        verbose    = FALSE
      )
      for (f in setdiff(intersect(track, rsd_res$features_kept), var_res$features_kept))
        .add_track(f, "Variance filter", sprintf("Var=%.2e", var_res$variance_values[[f]]))
      dims_add <- .record_dim(dims_add, sprintf("[%s] After variance filter", nt_label),
                              var_res$data)
      msg(sprintf("    Variance: removed %d. %d remaining.",
                  var_res$n_features_removed, var_res$n_features_after))

      features_final <- var_res$features_kept
      if (length(features_final) == 0L)
        stop("No features passed quality filtering. ",
             "Consider relaxing 'rsd_threshold' or 'variance_percentile'.", call. = FALSE)

      list(df_norm        = df_norm[,  features_final, drop = FALSE],
           df_trans       = df_trans[, features_final, drop = FALSE],
           features_final = features_final,
           track = track_add, dims = dims_add,
           norm_m = norm_m, trans_m = trans_m, error = NULL)

    }, error = function(e) {
      message(sprintf("run_DIpreprocess ERROR in norm+transform '%s': %s",
                      nt_label, conditionMessage(e)))
      list(track = track_add, dims = dims_add,
           norm_m = norm_m, trans_m = trans_m, error = conditionMessage(e))
    })
  }

  # ===========================================================================
  # STEPS 9-10: PER FULL-COMBINATION HELPER
  # ===========================================================================

  # Scaling (two branches) -> replicate merging for one combination.
  .run_scale_merge <- function(nt, scale_np, scale_p, combo_label) {

    branch_out <- list(
      metadata              = shared$metadata,
      data_raw              = shared$data_raw,
      data_missing_filtered = shared$data_missing_filtered,
      data_imputed          = shared$data_imputed,
      data_corrected        = shared$data_corrected,
      data_normalized       = nt$df_norm,
      data_transformed      = nt$df_trans,
      data_nonpls           = NULL,
      data_pls              = NULL,
      data_nonpls_merged    = NULL,
      data_pls_merged       = NULL,
      metadata_merged       = NULL,
      features_final        = if (is.null(nt$features_final)) character(0L) else nt$features_final,
      uncorrected_features  = shared$uncorrected_features,
      track                 = rbind(shared$track, nt$track),
      dimensions            = rbind(shared$dimensions, nt$dims),
      parameters            = .params(nt$norm_m, nt$trans_m, scale_np, scale_p),
      elapsed_seconds       = NA_real_,
      error                 = nt$error
    )

    # Propagate norm+transform+filter failures immediately
    if (!is.null(nt$error))
      return(structure(branch_out, class = c("run_DIpreprocess", "list")))

    t_branch  <- proc.time()[["elapsed"]]
    meta_snap <- shared$metadata

    tryCatch({

      # 9. SCALING (two parallel branches)
      msg(sprintf("  [%s] Step 9: Scaling (NONPLS='%s', PLS='%s')...",
                  combo_label, scale_np, scale_p))

      .scale_branch <- function(method)
        run_scale(x = nt$df_trans, method = method, metadata = meta_snap,
                  group_col = "Group_", verbose = FALSE)$data

      branch_out$data_nonpls <- .scale_branch(scale_np)
      branch_out$data_pls    <- .scale_branch(scale_p)

      branch_out$dimensions <- .record_dim(branch_out$dimensions,
        sprintf("[%s] Final (NONPLS)", combo_label), branch_out$data_nonpls)
      branch_out$dimensions <- .record_dim(branch_out$dimensions,
        sprintf("[%s] Final (PLS)", combo_label), branch_out$data_pls)

      # 10. TECHNICAL-REPLICATE MERGING
      msg(sprintf("  [%s] Step 10: Replicate merging...", combo_label))

      # meta_snap$SubjectID is filled from subject_id_col (all NA if absent)
      sid      <- meta_snap$SubjectID
      has_sid  <- !is.na(sid) & sid != "" & sid != "NA"
      do_merge <- merge_replicates && any(has_sid)

      if (do_merge) {
        mergeable <- meta_snap$Group_ != "QC" & has_sid

        if (!anyDuplicated(sid[mergeable])) {
          msg(sprintf("  [%s] No technical replicates detected; skipping merge.",
                      combo_label))
        } else {
          .merge_br <- function(data) {
            bio_data  <- data[mergeable,  , drop = FALSE]
            keep_data <- data[!mergeable, , drop = FALSE]
            grp_idx   <- split(seq_len(nrow(bio_data)), sid[mergeable])
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

          msg(sprintf("  [%s] Merged %d rows into %d unique subject(s).",
                      combo_label, sum(mergeable), nrow(merged_meta_rows)))

          branch_out$dimensions <- .record_dim(branch_out$dimensions,
            sprintf("[%s] After replicate merge (NONPLS)", combo_label),
            branch_out$data_nonpls_merged)
          branch_out$dimensions <- .record_dim(branch_out$dimensions,
            sprintf("[%s] After replicate merge (PLS)", combo_label),
            branch_out$data_pls_merged)
        }

      } else if (merge_replicates) {
        msg(sprintf("  [%s] merge_replicates=TRUE but no valid SubjectID; skipping.",
                    combo_label))
      } else {
        msg(sprintf("  [%s] Skipped (merge_replicates = FALSE).", combo_label))
      }

    }, error = function(e) {
      message(sprintf("run_DIpreprocess ERROR in combination '%s': %s",
                      combo_label, conditionMessage(e)))
      branch_out$error <<- conditionMessage(e)
    })

    branch_out$elapsed_seconds <- proc.time()[["elapsed"]] - t_branch
    structure(branch_out, class = c("run_DIpreprocess", "list"))
  }

  # ===========================================================================
  # DISPATCH: single vs. multi
  # ===========================================================================

  if (!is_multi) {

    # Single combination - original return shape
    msg("Steps 6-10: Normalization -> Transformation -> Filtering -> Scaling -> Merge...")

    nt  <- .run_norm_trans(normalize_method, transform_method)
    out <- .run_scale_merge(nt, scale_nonpls, scale_pls, combo_names)

    out$elapsed_seconds <- proc.time()[["elapsed"]] - t_start
    if (is.null(out$error))
      msg(sprintf("Pipeline completed successfully. Elapsed: %.1f second(s).",
                  out$elapsed_seconds))

    return(out)

  } else {

    # Multiple combinations - named list of run_DIpreprocess objects.
    # Norm+transform is cached per unique (norm, trans) pair; scale/filter/merge
    # runs once per full (norm, trans, scale_np, scale_p) combination.

    nt_grid <- unique(combos[, c("norm", "trans"), drop = FALSE])
    nt_keys <- paste(nt_grid$norm, nt_grid$trans, sep = "_")

    msg(sprintf(
      "\nRunning %d norm+transform pair(s), then fanning out to %d full combination(s)...",
      nrow(nt_grid), nrow(combos)
    ))

    # Cache norm+transform to avoid re-running for each scale variant
    nt_cache        <- vector("list", nrow(nt_grid))
    names(nt_cache) <- nt_keys

    for (i in seq_len(nrow(nt_grid))) {
      msg(sprintf("\n--- Norm+Transform %d/%d: %s ---", i, nrow(nt_grid), nt_keys[i]))
      nt_cache[[nt_keys[i]]] <- .run_norm_trans(nt_grid$norm[i], nt_grid$trans[i])
    }

    results        <- vector("list", nrow(combos))
    names(results) <- combo_names

    for (i in seq_len(nrow(combos))) {
      nm     <- combo_names[i]
      nt_key <- paste(combos$norm[i], combos$trans[i], sep = "_")
      msg(sprintf("\n--- Combination %d/%d: %s ---", i, nrow(combos), nm))
      results[[nm]] <- .run_scale_merge(
        nt_cache[[nt_key]], combos$scale_np[i], combos$scale_p[i], nm
      )
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
      class           = c("run_DIpreprocess_multi", "list"),
      elapsed_seconds = elapsed_total,
      combos          = combos,
      shared          = shared
    )
  }
}


# ==============================================================================
# S3 METHODS - run_DIpreprocess
# ==============================================================================

#' @param object A `run_DIpreprocess` or `run_DIpreprocess_multi` object.
#' @param ... Ignored.
#' @rdname run_DIpreprocess
#' @export
print.run_DIpreprocess <- function(x, ...) {
  cat("=== run_DIpreprocess pipeline result ===\n")
  if (!is.null(x$error)) {
    cat(sprintf("STATUS : FAILED - %s\n", x$error))
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
  if (!is.null(x$track) && nrow(x$track) > 0L)
    cat(sprintf("Track  : %d features tracked in removal summary\n", nrow(x$track)))
  cat(sprintf("Time   : %.1f second(s)\n", x$elapsed_seconds))
  invisible(x)
}

#' @rdname run_DIpreprocess
#' @export
summary.run_DIpreprocess <- function(object, ...) {
  cat("=== run_DIpreprocess - step-by-step dimensions ===\n")
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
# S3 METHODS - run_DIpreprocess_multi
# ==============================================================================

#' @rdname run_DIpreprocess
#' @export
print.run_DIpreprocess_multi <- function(x, ...) {
  cat("=== run_DIpreprocess_multi pipeline result ===\n")

  if (!is.null(x$error)) {
    cat(sprintf("STATUS : FAILED - %s\n", x$error))
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

#' @rdname run_DIpreprocess
#' @export
summary.run_DIpreprocess_multi <- function(object, ...) {
  cat("=== run_DIpreprocess_multi - per-combination summaries ===\n\n")
  for (nm in names(object)) {
    cat(sprintf("--- %s ---\n", nm))
    summary(object[[nm]])
    cat("\n")
  }
  invisible(object)
}