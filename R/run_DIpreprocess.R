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
#' (1) the steps up to and including imputation run \strong{only once};
#' (2) normalization, transformation, and quality filtering run
#' \strong{once per unique normalize x transform pair}; (3) scaling and
#' replicate merging run \strong{once per full combination}. The return value
#' is then a list of class `"run_DIpreprocess_multi"` whose elements are
#' standard `run_DIpreprocess` objects (see Value for their names).
#'
#' @param x Matrix or data frame. Numeric feature data with \strong{samples in
#'   rows} and \strong{features in columns}. All values should be non-negative
#'   raw intensities; zeros are treated as missing.
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
#'   injection order (unique values). Required for QC-based drift correction.
#'   Default: `"InjectionSequence"`.
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
#' @param blank_types Character vector or `NULL`. Labels in `group_col` of
#'   procedural blank samples. When given, potential contaminants are flagged
#'   and the unstable ones removed with [run_filterblank()] (QComics rules),
#'   and the blank samples are then dropped. Default: `NULL` (no blank step).
#' @param missing_threshold Numeric \[0, 1\]. Maximum proportion of missing
#'   values allowed, passed to [run_filtermissing()] as `threshold`.
#'   Default: `0.2`.
#' @param missing_by_group Logical. If `TRUE`, assess missingness per group.
#'   Default: `TRUE`.
#' @param missing_any_group Logical. Only used when `missing_by_group = TRUE`.
#'   If `TRUE` (default), a feature is retained if it passes the threshold in
#'   at least one group (modified 80% rule). If `FALSE`, it must pass in all
#'   groups.
#' @param missing_include_qc Logical. If `TRUE`, include QC samples when
#'   assessing missingness. Default: `FALSE`.
#' @param missing_absent_threshold Numeric or `NULL`. Passed to
#'   [run_filtermissing()] as `absent_threshold`: features missing in more
#'   than this proportion of a group are kept as truly absent (zeros) there
#'   (QComics uses `0.7`). Default: `NULL` (off).
#' @param missing_sample_iqr Numeric or `NULL`. Passed to
#'   [run_filtermissing()] as `sample_iqr`: samples whose missingness over the
#'   commonly detected features is abnormally high within their group (more
#'   than `missing_sample_iqr` interquartile ranges above the third quartile)
#'   are removed as likely failed injections, as QComics recommends. `NULL`
#'   disables the check. Default: `3`.
#' @param drift_method Character. Signal drift and batch correction method
#'   passed to [run_driftBatchCorrect()]: `"QCRSC"` (default), `"tGAM"`,
#'   `"rGAM"`, `"rLOESS"`, `"ComBat"`, or `"none"` to skip the step.
#' @param drift_args Named list of further arguments for
#'   [run_driftBatchCorrect()], e.g.
#'   `list(spline_smooth_param = 0.5, min_QC = 4)` or
#'   `list(sample_weight = 0.2, qc_check = TRUE)`. Default: `list()`.
#' @param remove_uncorrected Logical. If `TRUE` (default), remove features
#'   that drift correction could not correct in every batch (too few QC
#'   values). If `FALSE`, keep them with their original values in those
#'   batches. See [run_driftBatchCorrect()].
#' @param impute_method Imputation method passed to [run_mvimpute()] as
#'   `method`: a fraction (e.g. `0.2` = one fifth of the smallest positive
#'   value per feature), `"quantileregression"`, a \pkg{mice} method, or
#'   `"none"`. Imputation runs once, after drift correction. Default: `0.2`.
#' @param positive_only Logical. If `TRUE`, the minimum used by fractional
#'   imputation only considers positive values. Default: `TRUE`.
#' @param normalize_method Character or character vector. One or more
#'   normalization methods (see [run_normalize()] for details). When multiple
#'   methods are supplied every method is paired with every element of
#'   `transform_method` (Cartesian product).
#'   \itemize{
#'     \item `"sum"`: Normalizes by sum.
#'     \item `"median"`: Normalizes by median.
#'     \item `"specific_factor"`: Uses external factors in `norm_factor_col`.
#'     \item `"pqn_global"`: PQN using all samples as reference.
#'     \item `"pqn_reference"`: PQN using `normalize_ref_sample` as reference.
#'     \item `"pqn_group"`: PQN using pooled QC samples as reference.
#'     \item `"none"`: No normalization.
#'   }
#'   Default: `"median"`.
#' @param normalize_ref_sample Character. Reference sample identifier for
#'   `"pqn_reference"`. Default: `NULL`.
#' @param normalize_qc_method Character. How QC samples are scaled when
#'   sample factors are used: `"mean"`, `"median"`, or `"none"`. Default: `"none"`.
#' @param transform_method Character or character vector. One or more
#'   transformation methods (see [run_transform()]): `"ln"`, `"log2"`,
#'   `"log10"`, `"sqrt"`, `"cbrt"`, `"clr"`, `"arcsin_sqrt"`, `"vsn"`,
#'   `"glog"`, or `"none"`. When multiple methods are supplied every method is
#'   paired with every element of `normalize_method` (Cartesian product).
#'   Default: `"log10"`.
#' @param rsd_threshold Numeric \[0, 1\]. Features with QC RSD at or above this
#'   value are removed ([run_filterRSD()]). Default: `0.3`.
#' @param rsd_qc_type Character. QC samples used for the RSD filter: `"QC"`
#'   (all `qc_types`), or one or more specific labels such as `"EQC"`.
#'   Default: `"QC"`.
#' @param variance_percentile Numeric \[0, 100\]. Features in this bottom
#'   percentile of variance are removed ([run_filtervariance()]). Default: `10`.
#' @param scale_nonpls Character or character vector. One or more scaling
#'   methods (see [run_scale()]) for the NONPLS branch, intended for PCA,
#'   clustering, and univariate tests. When multiple methods are supplied,
#'   each is paired with every element of `scale_pls`. Only the varying scale
#'   components are appended to the combination name (e.g.
#'   `median_log10_auto` when only `scale_nonpls` varies). Default: `"auto"`.
#' @param scale_pls Character or character vector. One or more scaling
#'   methods for the PLS branch, intended for PLS-DA/OPLS-DA. Default:
#'   `"pareto"`.
#' @param merge_replicates Logical. Average technical replicates (non-QC rows
#'   sharing a `subject_id_col` value). Default: `FALSE`.
#' @param track Character vector or `NULL`. Feature names to monitor across
#'   preprocessing steps. If a tracked feature is removed, the step and reason
#'   (e.g., missingness percentage, RSD value) are recorded in `$track`.
#'   Default: `NULL`.
#' @param verbose Logical. Print progress. Default: `TRUE`.
#' @param object A `run_DIpreprocess` or `run_DIpreprocess_multi` object.
#' @param ... Ignored.
#'
#' @return
#' \describe{
#'   \item{Single combination (all four method arguments are scalar)}{A
#'     list of class `run_DIpreprocess` with elements:
#'     \describe{
#'       \item{`data_nonpls`, `data_pls`}{Final data for the NONPLS and PLS
#'         branches (scaled with `scale_nonpls` / `scale_pls`).}
#'       \item{`data_raw`, `data_missing_filtered`, `data_corrected`,
#'         `data_imputed`, `data_normalized`, `data_transformed`}{Data after
#'         each step (`data_normalized`/`data_transformed` are restricted to
#'         `features_final`). `data_raw` has zeros as `NA`; missing values
#'         remain up to `data_corrected`.}
#'       \item{`data_nonpls_merged`, `data_pls_merged`, `metadata_merged`}{
#'         Replicate-merged data and metadata, or `NULL`.}
#'       \item{`metadata`}{Standardized metadata with columns `Sample`,
#'         `Group`, `Group_` (QC types collapsed to `"QC"`), `Batch`,
#'         `InjectionSequence`, `SubjectID`.}
#'       \item{`features_final`}{Features that passed every filter.}
#'       \item{`samples_removed`}{Samples removed for abnormally high
#'         missingness (see `missing_sample_iqr`).}
#'       \item{`uncorrected_features`}{Features drift correction could not
#'         correct in every batch.}
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
#'     The shared data (up to imputation) is stored in the `shared` attribute
#'     of the returned list.}
#' }
#'
#' @details
#' \strong{Preprocessing Workflow (Fixed Order):}
#' \enumerate{
#'   \item Input validation
#'   \item Sample removal: `outliers`, then blank filtering and blank removal
#'     ([run_filterblank()]) when `blank_types` is given
#'   \item Zeros to `NA`, then missing-value handling ([run_filtermissing()]):
#'     removal of samples with abnormally low detection rates, feature
#'     filtering, and optional truly-absent zeros
#'   \item Signal drift and batch correction ([run_driftBatchCorrect()]),
#'     on data that still contain missing values
#'   \item Missing-value imputation ([run_mvimpute()])
#'   \item \emph{\[Per normalize x transform pair\]} Normalization
#'     ([run_normalize()])
#'   \item \emph{\[Per normalize x transform pair\]} Transformation
#'     ([run_transform()])
#'   \item \emph{\[Per normalize x transform pair\]} Quality filtering:
#'     QC-RSD on the measured (non-imputed) normalized values
#'     ([run_filterRSD()]), then variance on the transformed data
#'     ([run_filtervariance()])
#'   \item \emph{\[Per full combination\]} Scaling ([run_scale()]) of the
#'     NONPLS and PLS branches
#'   \item \emph{\[Per full combination\]} Technical-replicate merging
#' }
#'
#' \strong{Why this order.} Drift correction comes before imputation because
#' QC-based correction ignores missing values, whereas imputed QC values
#' would bias the fitted trends (Kirwan et al., 2013); imputation therefore
#' runs once, on corrected data. The QC-RSD leaves imputed values out: it
#' measures the precision of measured values, and imputed near-minimum values
#' among the QCs would inflate it. Quality filtering comes \emph{before}
#' scaling because the RSD of mean-centred data and the variance of
#' auto-scaled data (always 1) carry no information. Both branches therefore
#' share the same feature set, `features_final`.
#'
#' The workflow follows the QC-based guidelines of Broadhurst et al. (2018)
#' and QComics (González-Domínguez et al., 2024) and the DIMS
#' batch-correction strategy of Kirwan et al. (2013). The outlier-robust
#' drift correction methods of Vynck et al. (2026) are available through
#' `drift_method`.
#'
#' @note
#' \strong{Zero Handling:} The pipeline converts all `0` values in `x` to
#' `NA` before analysis (unless `missing_absent_threshold` marks them as
#' truly absent).
#'
#' \strong{Pre-flight checks:} Method names, and method-specific requirements
#' (e.g. `norm_factor_col` for `"specific_factor"`, `normalize_ref_sample` for
#' `"pqn_reference"`, `drift_args` names), are validated \emph{before} any
#' preprocessing step, so errors surface early rather than mid-run. Errors
#' raised inside a step are caught: the function returns an object whose
#' `$error` holds the message.
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
#' González-Domínguez, A., Estanyol-Torres, N., Brunius, C., Landberg, R., &
#' González-Domínguez, R. (2024). QComics: Recommendations and guidelines for
#' robust, easily implementable and reportable quality control of
#' metabolomics data. *Analytical Chemistry*, 96(3), 1064-1072.
#' \doi{10.1021/acs.analchem.3c03660}
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
#' Vynck, M., Vangeenderhuysen, P., De Paepe, E., Nawrot, T., Plekhova, V., &
#' Vanhaecke, L. (2026). Robust metabolomics data normalization across scales
#' and experimental designs. *Analytical Chemistry*, 98(24), 17627-17637.
#' \doi{10.1021/acs.analchem.5c06841}
#'
#' See the individual step functions for method-specific references.
#'
#' @family preprocessing steps
#' @seealso [plot_beforeafter()] and [plot_dist_beforeafter()] accept a
#'   `run_DIpreprocess` result directly; [dims_data] and [dims_metadata] are
#'   example inputs.
#'
#' @export
#'
#' @examples
#' data(dims_data)
#' data(dims_metadata)
#'
#' # Default pipeline, including QC-RSC drift correction
#' res <- run_DIpreprocess(dims_data, dims_metadata, verbose = FALSE)
#' res
#' summary(res)
#' dim(res$data_nonpls)
#'
#' # Multi-method mode: steps 1-5 run once, then every combination
#' multi <- run_DIpreprocess(
#'   dims_data, dims_metadata,
#'   normalize_method = c("median", "pqn_group"),
#'   transform_method = "log10",
#'   scale_nonpls     = c("auto", "pareto"),
#'   drift_method     = "none",
#'   verbose          = FALSE
#' )
#' names(multi)
#' multi
#'
#' # Outlier-robust drift correction (Vynck et al., 2026; needs mgcv)
#' \donttest{
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   run_DIpreprocess(dims_data[, 1:20], dims_metadata, drift_method = "tGAM",
#'                    drift_args = list(sample_weight = 0.2), verbose = FALSE)
#' }
#' }
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

    # - Sample removal --------------------------------------------------------
    outliers         = NULL,
    blank_types      = NULL,

    # - Missing-value filter --------------------------------------------------
    missing_threshold        = 0.2,
    missing_by_group         = TRUE,
    missing_any_group        = TRUE,
    missing_include_qc       = FALSE,
    missing_absent_threshold = NULL,
    missing_sample_iqr       = 3,

    # - Drift / batch correction ----------------------------------------------
    drift_method       = "QCRSC",
    drift_args         = list(),
    remove_uncorrected = TRUE,

    # - Missing-value imputation ----------------------------------------------
    impute_method    = 0.2,
    positive_only    = TRUE,

    # - Normalization ---------------------------------------------------------
    normalize_method     = "median",
    normalize_ref_sample = NULL,
    normalize_qc_method  = "none",

    # - Transformation --------------------------------------------------------
    transform_method = "log10",

    # - Quality filtering -----------------------------------------------------
    rsd_threshold       = 0.3,
    rsd_qc_type         = "QC",
    variance_percentile = 10,

    # - Scaling ---------------------------------------------------------------
    scale_nonpls     = "auto",
    scale_pls        = "pareto",

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

  # Parameters actually used (evaluated at call time, so e.g. a drift_method
  # switched to "none" for missing columns is reported as such).
  .params <- function(norm_m, trans_m, scale_np, scale_p) list(
    sample_id_col            = sample_id_col,
    group_col                = group_col,
    qc_types                 = qc_types,
    batch_col                = batch_col,
    injection_col            = injection_col,
    norm_factor_col          = norm_factor_col,
    subject_id_col           = subject_id_col,
    outliers                 = outliers,
    blank_types              = blank_types,
    missing_threshold        = missing_threshold,
    missing_by_group         = missing_by_group,
    missing_any_group        = missing_any_group,
    missing_include_qc       = missing_include_qc,
    missing_absent_threshold = missing_absent_threshold,
    missing_sample_iqr       = missing_sample_iqr,
    drift_method             = drift_method,
    drift_args               = drift_args,
    remove_uncorrected       = remove_uncorrected,
    impute_method            = impute_method,
    positive_only            = positive_only,
    normalize_method         = norm_m,
    normalize_ref_sample     = normalize_ref_sample,
    normalize_qc_method      = normalize_qc_method,
    transform_method         = trans_m,
    rsd_threshold            = rsd_threshold,
    rsd_qc_type              = rsd_qc_type,
    variance_percentile      = variance_percentile,
    scale_nonpls             = scale_np,
    scale_pls                = scale_p,
    merge_replicates         = merge_replicates,
    track                    = track,
    verbose                  = verbose
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
  if (!is.character(drift_method) || length(drift_method) != 1L)
    stop("'drift_method' must be a single character string.", call. = FALSE)
  .check_methods(drift_method, c("none", .drift_methods), "drift_method")
  if ("pqn_reference" %in% normalize_method && is.null(normalize_ref_sample))
    stop("normalize_method = 'pqn_reference' requires 'normalize_ref_sample'.", call. = FALSE)

  # drift_args: named run_driftBatchCorrect() arguments not set elsewhere
  if (!is.list(drift_args) ||
      (length(drift_args) > 0L && (is.null(names(drift_args)) || any(names(drift_args) == ""))))
    stop("'drift_args' must be a named list of run_driftBatchCorrect() arguments.",
         call. = FALSE)
  reserved <- c("x", "metadata", "method", "injection_col", "batch_col", "group_col",
                "qc_types", "remove_uncorrected", "verbose")
  bad_args <- intersect(names(drift_args), reserved)
  if (length(bad_args) > 0L)
    stop("Set these through run_DIpreprocess() arguments, not 'drift_args': ",
         paste(bad_args, collapse = ", "), ".", call. = FALSE)
  unknown_args <- setdiff(names(drift_args), names(formals(run_driftBatchCorrect)))
  if (length(unknown_args) > 0L)
    stop("Unknown 'drift_args' for run_driftBatchCorrect(): ",
         paste(unknown_args, collapse = ", "), ".", call. = FALSE)

  # RSD filter QC samples: "QC" means all qc_types
  if (!is.character(rsd_qc_type) || !all(rsd_qc_type %in% c("QC", qc_types)))
    stop("'rsd_qc_type' must be \"QC\" or labels from 'qc_types'.", call. = FALSE)
  rsd_qc <- if (identical(rsd_qc_type, "QC")) qc_types else rsd_qc_type

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
    msg("Steps 1-5 (through imputation) will run once; ",
        "Steps 6-8 (norm + transform + quality filter) once per unique norm x transform pair; ",
        "Steps 9-10 (scaling + merging) once per full combination.")
  }

  # ===========================================================================
  # STEPS 1-5: RUN ONCE
  # ===========================================================================

  shared <- list(
    metadata              = NULL,
    data_raw              = NULL,
    data_missing_filtered = NULL,
    data_corrected        = NULL,
    data_imputed          = NULL,
    samples_removed       = character(0L),
    uncorrected_features  = character(0L),
    dimensions            = .empty_dims(),
    track                 = .track_init(),
    error                 = NULL
  )
  .add_shared_track <- function(f, step, value)
    shared$track <<- rbind(shared$track, data.frame(Feature = f, Step = step, Value = value,
                                                    stringsAsFactors = FALSE))

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
          stop(sprintf("Metadata column '%s' must be numeric.", col_name))
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

    is_numeric_col <- vapply(df, is.numeric, logical(1L))
    if (!all(is_numeric_col)) {
      bad_cols <- names(df)[!is_numeric_col]
      stop(sprintf(
        "Input 'x' must contain ONLY numeric columns. Non-numeric: %s.",
        paste(utils::head(bad_cols, 5), collapse = ", ")
      ))
    }

    shared$dimensions <- .record_dim(shared$dimensions, "Original", df)

    # =========================================================================
    # 2. SAMPLE REMOVAL: OUTLIERS, THEN BLANKS
    # =========================================================================

    msg("Step 2/10: Removing outliers and blanks...")

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

    if (!is.null(blank_types)) {
      if (!any(metadata[[group_col]] %in% blank_types)) {
        warning("No blank samples (", paste(blank_types, collapse = ", "),
                ") found; skipping the blank filter.", call. = FALSE)
      } else {
        bf <- run_filterblank(df, metadata, group_col = group_col,
                              blank_types = blank_types, qc_types = qc_types,
                              verbose = FALSE)
        for (f in intersect(track, bf$features_removed)) {
          r <- bf$blank_summary[bf$blank_summary$Feature == f, ]
          .add_shared_track(f, "Blank filter",
                            sprintf("sample/blank = %.2f, blank RSD = %.0f%%",
                                    r$Ratio, r$Blank_RSD * 100))
        }
        df       <- bf$data
        metadata <- bf$metadata
        msg(sprintf("  Blank filter: %d potential contaminant(s) flagged, %d removed; %d blank sample(s) dropped.",
                    length(bf$features_flagged), length(bf$features_removed), bf$n_blanks))
        shared$dimensions <- .record_dim(shared$dimensions, "After blank filter", df)
      }
    }

    # =========================================================================
    # METADATA SNAPSHOT (after sample removal)
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
    # 3. ZEROS -> NA AND MISSING-VALUE FILTERING
    # =========================================================================

    msg(sprintf("Step 3/10: Missing values (zeros as NA; threshold = %.0f%%)...",
                missing_threshold * 100))

    n_zeros <- sum(df == 0, na.rm = TRUE)
    if (n_zeros > 0L) {
      df[df == 0] <- NA
      msg(sprintf("  Converted %d zero(s) to NA (structural missingness).", n_zeros))
    }

    all_na <- colSums(is.na(df)) == nrow(df)
    if (any(all_na)) {
      for (f in intersect(track, names(all_na)[all_na]))
        .add_shared_track(f, "All-NA removal", "100% missing")
      df <- df[, !all_na, drop = FALSE]
      msg(sprintf("  Removed %d all-NA feature(s). %d remaining.",
                  sum(all_na), ncol(df)))
    }

    shared$data_raw   <- df
    shared$dimensions <- .record_dim(shared$dimensions, "After zero -> NA conversion", df)

    df_before_miss <- df
    filt_miss <- run_filtermissing(
      x                = df,
      metadata         = metadata,
      threshold        = missing_threshold,
      filter_by_group  = missing_by_group,
      any_group        = missing_any_group,
      include_QC       = missing_include_qc,
      group_col        = group_col,
      qc_types         = qc_types,
      zero_as_missing  = FALSE,
      absent_threshold = missing_absent_threshold,
      sample_iqr       = missing_sample_iqr,
      verbose          = FALSE
    )
    df <- filt_miss$data

    if (length(filt_miss$samples_removed) > 0L) {
      # keep metadata and earlier data aligned with the remaining samples
      kept_rows <- !rownames(metadata) %in% filt_miss$samples_removed
      metadata        <- metadata[kept_rows, , drop = FALSE]
      shared$metadata <- shared$metadata[kept_rows, , drop = FALSE]
      shared$data_raw <- shared$data_raw[kept_rows, , drop = FALSE]
      shared$samples_removed <- filt_miss$samples_removed
      msg(sprintf("  Removed %d sample(s) with abnormally high missingness: %s",
                  length(filt_miss$samples_removed),
                  paste(filt_miss$samples_removed, collapse = ", ")))
    }

    for (f in intersect(track, filt_miss$features_removed)) {
      vals <- df_before_miss[[f]]
      gps  <- as.character(metadata[[group_col]])
      miss_info <- vapply(sort(unique(gps)), function(g) {
        sub <- vals[gps == g]
        sprintf("%s=%.0f%%", g, mean(is.na(sub)) * 100)
      }, character(1L))
      .add_shared_track(f, "Missing value", paste(miss_info, collapse = ", "))
    }

    shared$data_missing_filtered <- df
    shared$dimensions <- .record_dim(shared$dimensions, "After missing filter", df)

    msg(sprintf("  Removed %d feature(s) (%.1f%%). %d remaining.",
                filt_miss$n_features_removed,
                filt_miss$n_features_removed / filt_miss$n_features_before * 100,
                ncol(df)))
    if (length(filt_miss$features_absent) > 0L)
      msg(sprintf("  Kept %d truly absent feature(s) with %d zero(s).",
                  length(filt_miss$features_absent), filt_miss$n_absent_zeros))

    # =========================================================================
    # 4. SIGNAL DRIFT AND BATCH CORRECTION (before imputation)
    # =========================================================================

    msg(sprintf("Step 4/10: Drift/batch correction (%s)...", drift_method))

    if (drift_method != "none") {
      required_meta <- c(batch_col, group_col, if (drift_method != "ComBat") injection_col)
      missing_meta  <- setdiff(required_meta, colnames(metadata))
      if (length(missing_meta) > 0L) {
        warning("Skipping drift correction: metadata column(s) not found: ",
                paste(missing_meta, collapse = ", "), call. = FALSE)
        drift_method <- "none"
      }
    }

    if (drift_method != "none") {
      corr <- do.call(run_driftBatchCorrect, c(list(
        x                  = df,
        metadata           = metadata,
        method             = drift_method,
        injection_col      = injection_col,
        batch_col          = batch_col,
        group_col          = group_col,
        qc_types           = qc_types,
        remove_uncorrected = remove_uncorrected,
        verbose            = FALSE
      ), drift_args))

      shared$uncorrected_features <- corr$uncorrected_features
      if (remove_uncorrected) {
        for (f in intersect(track, corr$uncorrected_features))
          .add_shared_track(f, "Drift correction", sprintf(
            "Uncorrected in batch(es) %s",
            paste(corr$uncorrected$Batch[corr$uncorrected$Feature == f], collapse = ", ")))
      }
      df <- corr$data
      msg(sprintf("  %s applied. %d feature(s) not correctable in every batch (%s). %d remaining.",
                  drift_method, corr$n_uncorrected,
                  if (remove_uncorrected) "removed" else "kept uncorrected", ncol(df)))
    } else {
      msg("  Skipped (drift_method = 'none').")
    }

    shared$data_corrected <- df
    shared$dimensions     <- .record_dim(shared$dimensions, "After drift correction", df)

    # =========================================================================
    # 5. MISSING-VALUE IMPUTATION (once, on corrected data)
    # =========================================================================

    msg(sprintf("Step 5/10: Imputing missing values (method = %s)...",
                as.character(impute_method)))

    imp <- run_mvimpute(
      x             = df,
      method        = impute_method,
      positive_only = positive_only,
      verbose       = FALSE
    )
    df <- imp$data
    shared$data_imputed <- df
    msg(sprintf("  Imputed %d missing value(s).", imp$n_missing_before - imp$n_missing_after))

  }, error = function(e) {
    message("run_DIpreprocess ERROR (Steps 1-5): ", conditionMessage(e))
    shared$error <<- conditionMessage(e)
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
        x             = shared$data_imputed,
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

      # RSD measures the precision of measured QC values, so imputed values
      # (still NA in data_corrected) are left out
      df_obs <- df_norm
      df_obs[is.na(shared$data_corrected)] <- NA
      rsd_res <- run_filterRSD(
        x         = df_obs,
        metadata  = metadata,
        max_rsd   = rsd_threshold,
        group_col = "Group",      # raw labels, so rsd_qc_type can pick EQC/SQC
        qc_types  = rsd_qc,
        verbose   = FALSE
      )
      for (f in intersect(track, rsd_res$features_removed)) {
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
      for (f in intersect(track, var_res$features_removed))
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
      data_corrected        = shared$data_corrected,
      data_imputed          = shared$data_imputed,
      data_normalized       = nt$df_norm,
      samples_removed       = shared$samples_removed,
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
    # Norm+transform+filter is cached per unique (norm, trans) pair;
    # scale/merge runs once per full (norm, trans, scale_np, scale_p) combination.

    nt_grid <- unique(combos[, c("norm", "trans"), drop = FALSE])
    nt_keys <- paste(nt_grid$norm, nt_grid$trans, sep = "_")

    msg(sprintf(
      "\nRunning %d norm+transform pair(s), then fanning out to %d full combination(s)...",
      nrow(nt_grid), nrow(combos)
    ))

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
  if (length(x$samples_removed) > 0L)
    cat(sprintf("Samples: %d removed for abnormally high missingness\n",
                length(x$samples_removed)))
  if (!is.null(x$data_nonpls_merged))
    cat(sprintf("Merged : %d samples x %d features (replicate-merged)\n",
                nrow(x$data_nonpls_merged), ncol(x$data_nonpls_merged)))
  if (!is.null(x$parameters$drift_method))
    cat(sprintf("Drift  : %s", x$parameters$drift_method),
        if (length(x$uncorrected_features) > 0L)
          sprintf("(%d uncorrected feature(s) %s)", length(x$uncorrected_features),
                  if (isTRUE(x$parameters$remove_uncorrected)) "removed" else "kept"),
        "\n")
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

  if (!is.null(x[["error"]])) {
    cat(sprintf("STATUS : FAILED - %s\n", x[["error"]]))
    return(invisible(x))
  }

  elapsed <- attr(x, "elapsed_seconds")

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
