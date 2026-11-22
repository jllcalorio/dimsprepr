#' Correct Signal Drift and Batch Effects Using QC Samples
#'
#' @description
#' Removes within-batch signal drift and between-batch intensity differences.
#' For every feature and batch, a smooth trend of intensity against injection
#' order is fitted through the quality control (QC) samples and divided out.
#' Four trend models are available: the classic QC-based robust spline
#' correction (QC-RSC; Kirwan et al., 2013) and the outlier-robust rLOESS,
#' rGAM and tGAM models of Vynck et al. (2026). ComBat (Johnson et al., 2007)
#' corrects batch offsets without using QCs or injection order.
#'
#' @param x Matrix or data frame of raw, non-negative intensities with samples
#'   in rows and features in columns. Missing values (`NA`) are allowed: they
#'   are ignored when fitting and stay missing, so impute \emph{after}
#'   correction. Zeros are also ignored when fitting and stay zero.
#'   For `print()`, a `run_driftBatchCorrect` object.
#' @param metadata Data frame with one row per row of `x`, in the same order.
#' @param method Character. Correction method: `"QCRSC"` (default), `"tGAM"`,
#'   `"rGAM"`, `"rLOESS"`, or `"ComBat"`. See Details.
#' @param injection_col Character. Column in `metadata` with the numeric,
#'   unique injection order. Not used by ComBat. Default: `"InjectionSequence"`.
#' @param batch_col Character. Column in `metadata` with batch labels (numeric
#'   or character). Default: `"Batch"`.
#' @param group_col Character. Column in `metadata` with sample group labels.
#'   Default: `"Group"`.
#' @param qc_types Character vector. Labels in `group_col` that identify QC
#'   samples. Default: `c("QC", "SQC", "EQC")`.
#' @param min_QC Integer, at least 4. Minimum number of observed (positive,
#'   non-missing) QC values a feature needs within a batch to be corrected in
#'   that batch. Default: `5`.
#' @param log_scale Logical. Fit the trend on log intensities (`TRUE`,
#'   default) or on the raw scale. The correction is multiplicative either
#'   way. ComBat is applied to log intensities when `TRUE`.
#' @param remove_uncorrected Logical. If `TRUE` (default), features that could
#'   not be corrected in every batch (too few QC values, or a failed fit) are
#'   removed. If `FALSE` they are kept, with their original values in the
#'   batches that could not be corrected.
#' @param spline_smooth_param Numeric in \[0, 1\]. QCRSC only. Smoothing
#'   parameter (`spar`) of [stats::smooth.spline()]; `0` (default) selects it
#'   per feature and batch by leave-one-out cross-validation.
#' @param spar_limit Numeric vector of length 2. QCRSC only. Search range for
#'   `spar` when `spline_smooth_param = 0`. Default: `c(-1.5, 1.5)`.
#' @param gam_k Integer, at least 3. tGAM and rGAM only. Maximum basis
#'   dimension of the thin plate regression spline; it is lowered
#'   automatically when a batch has few fitting points. Default: `10`.
#' @param sample_weight Numeric in \[0, 1\]. tGAM and rGAM only. Weight of the
#'   biological samples relative to the QC samples (weight 1) in the trend
#'   fit. `0` (default) fits QCs only; `0.2` is the 5:1 QC-to-sample weighting
#'   used by Vynck et al. (2026).
#' @param qc_check Logical. tGAM and rGAM only, and requires
#'   `sample_weight > 0`. If `TRUE`, test for every feature and batch whether
#'   the QC samples follow the same trend as the biological samples; where
#'   they do not (p < `qc_check_p`), fit the trend on the biological samples
#'   only. Default: `FALSE`.
#' @param qc_check_p Numeric in (0, 1). p-value cutoff of the QC check.
#'   Default: `0.1`.
#' @param use_parametric Logical. ComBat only. Use parametric (`TRUE`,
#'   default) or non-parametric empirical Bayes priors.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return A list of class `"run_driftBatchCorrect"` containing:
#'   \item{data}{Corrected data (same class as `x`), without the uncorrected
#'     features when `remove_uncorrected = TRUE`.}
#'   \item{data_before_correction}{The input data as a matrix-like object of
#'     the same class.}
#'   \item{method_used}{Character: the method applied.}
#'   \item{uncorrected}{Data frame with one row per feature and batch that
#'     could not be corrected (columns `Feature`, `Batch`).}
#'   \item{uncorrected_features}{Character vector of features not corrected
#'     in at least one batch.}
#'   \item{n_uncorrected}{Integer, `length(uncorrected_features)`.}
#'   \item{qc_check}{When `qc_check = TRUE`, a data frame with the QC check
#'     p-value per feature and batch and whether the trend was refitted on
#'     biological samples only (`Sample_based`); otherwise `NULL`.}
#'   \item{parameters}{List of parameters used.}
#'
#' @details
#' **Trend models.** For each feature and batch, a trend \eqn{f(t)} of the
#' (log) intensity against injection order \eqn{t} is fitted and each value
#' is corrected multiplicatively as \eqn{x \cdot m / f(t)}, where \eqn{m} is
#' the median intensity of the feature over all samples (on the log scale,
#' \eqn{f} is back-transformed first). Fitting each batch separately removes
#' both the within-batch drift and the between-batch offsets.
#'
#' * `"QCRSC"`: a cubic smoothing spline through the QC samples
#'   (Kirwan et al., 2013), the spline version of the QC-RLSC approach of
#'   Dunn et al. (2011). The smoothing parameter is chosen by leave-one-out
#'   cross-validation unless `spline_smooth_param` fixes it. This matches
#'   [pmp::QCRSC()] (Jankevics et al., 2025), except that values in
#'   batches that cannot be corrected are kept rather than set to `NA`.
#' * `"rLOESS"`: robust LOESS through the QC samples, with outliers
#'   down-weighted by Tukey's biweight (`loess(family = "symmetric")`) and the
#'   span chosen by generalized cross-validation (Vynck et al., 2026).
#'   Samples injected before the first or after the last QC of a batch get
#'   the trend value of the nearest QC.
#' * `"rGAM"`: a Gaussian additive model with a thin plate regression spline,
#'   fitted by REML, then refitted once with Tukey biweight weights
#'   \eqn{(1 - (r_i/c)^2)^2} computed from the residuals divided by their
#'   median absolute deviation, with \eqn{c = 4.685} (Vynck et al., 2026).
#' * `"tGAM"`: the same additive model with a scaled-t error distribution
#'   ([mgcv::scat()]), whose heavy tails make the fit resistant to outlying
#'   QC values. Vynck et al. (2026) found tGAM best overall and proposed it as
#'   the default for LC-MS data; it is the slowest method.
#' * `"ComBat"`: empirical Bayes location and scale adjustment of batches
#'   ([sva::ComBat()]); no drift correction. Requires \pkg{sva}.
#'
#' rGAM and tGAM require \pkg{mgcv}. By default they are fitted on QC samples
#' only, like QCRSC and rLOESS. With `sample_weight > 0` the biological
#' samples also inform the trend, which helps when QCs are sparse, and
#' `qc_check = TRUE` adds the QC representativeness check of Vynck et al.
#' (2026). That check compares, with a likelihood ratio test on
#' maximum-likelihood fits, one common trend against separate trends for QC
#' and biological samples.
#'
#' **QC layout.** QC samples should be spread through each batch, ideally
#' starting and ending it (González-Domínguez et al., 2024), with at least
#' `min_QC` per batch. The function stops if a batch has fewer QC samples
#' than `min_QC`. A feature with fewer than `min_QC` observed QC values in a
#' batch is not corrected in that batch; see `remove_uncorrected`.
#'
#' **Order of operations.** Correction runs on data with missing values, as
#' QC-RSC was designed to (Kirwan et al., 2013): impute afterwards, otherwise
#' imputed QC values distort the trend. [run_DIpreprocess()] follows this
#' order.
#'
#' @author John Lennon L. Calorio
#'
#' @references
#' Dunn, W.B., Broadhurst, D., Begley, P., Zelena, E., Francis-McIntyre, S.,
#' Anderson, N., et al. (2011). Procedures for large-scale metabolic profiling
#' of serum and plasma using gas chromatography and liquid chromatography
#' coupled to mass spectrometry. *Nature Protocols*, 6(7), 1060-1083.
#' \doi{10.1038/nprot.2011.335}
#'
#' González-Domínguez, A., Estanyol-Torres, N., Brunius, C., Landberg, R., &
#' González-Domínguez, R. (2024). QComics: Recommendations and guidelines for
#' robust, easily implementable and reportable quality control of
#' metabolomics data. *Analytical Chemistry*, 96(3), 1064-1072.
#' \doi{10.1021/acs.analchem.3c03660}
#'
#' Jankevics, A., Lloyd, G.R., & Weber, R.J.M. (2025). *pmp: Peak Matrix
#' Processing and signal batch correction for metabolomics datasets*.
#' R package. \doi{10.18129/B9.bioc.pmp}
#'
#' Johnson, W.E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in
#' microarray expression data using empirical Bayes methods.
#' *Biostatistics*, 8(1), 118-127. \doi{10.1093/biostatistics/kxj037}
#'
#' Kirwan, J.A., Broadhurst, D.I., Davidson, R.L., & Viant, M.R. (2013).
#' Characterising and correcting batch variation in an automated direct
#' infusion mass spectrometry (DIMS) metabolomics workflow.
#' *Analytical and Bioanalytical Chemistry*, 405(15), 5147-5157.
#' \doi{10.1007/s00216-013-6856-7}
#'
#' Vynck, M., Vangeenderhuysen, P., De Paepe, E., Nawrot, T., Plekhova, V., &
#' Vanhaecke, L. (2026). Robust metabolomics data normalization across scales
#' and experimental designs. *Analytical Chemistry*, 98(24), 17627-17637.
#' \doi{10.1021/acs.analchem.5c06841}
#'
#' Wood, S.N. (2011). Fast stable restricted maximum likelihood and marginal
#' likelihood estimation of semiparametric generalized linear models.
#' *Journal of the Royal Statistical Society: Series B*, 73(1), 3-36.
#' \doi{10.1111/j.1467-9868.2010.00749.x}
#'
#' @family preprocessing steps
#' @seealso [run_DIpreprocess()] for the full pipeline, [plot_beforeafter()]
#'   to inspect the correction along the injection order.
#'
#' @export
#'
#' @examples
#' data(dims_data)
#' data(dims_metadata)
#' raw <- replace(dims_data, dims_data == 0, NA)  # zeros = not detected
#'
#' # QC-based robust spline correction (missing values are fine)
#' res <- run_driftBatchCorrect(raw[, 1:10], dims_metadata,
#'                              spline_smooth_param = 0.5)
#' res
#'
#' # Outlier-robust tGAM using QCs and, with a 5:1 weight, biological samples
#' if (requireNamespace("mgcv", quietly = TRUE)) {
#'   run_driftBatchCorrect(raw[, 1:2], dims_metadata, method = "tGAM",
#'                         sample_weight = 0.2, verbose = FALSE)
#' }
#'
#' # ComBat (batch offsets only, needs sva) is called the same way:
#' # run_driftBatchCorrect(raw, dims_metadata, method = "ComBat")
run_driftBatchCorrect <- function(
    x,
    metadata,
    method              = "QCRSC",
    injection_col       = "InjectionSequence",
    batch_col           = "Batch",
    group_col           = "Group",
    qc_types            = c("QC", "SQC", "EQC"),
    min_QC              = 5,
    log_scale           = TRUE,
    remove_uncorrected  = TRUE,
    spline_smooth_param = 0,
    spar_limit          = c(-1.5, 1.5),
    gam_k               = 10,
    sample_weight       = 0,
    qc_check            = FALSE,
    qc_check_p          = 0.1,
    use_parametric      = TRUE,
    verbose             = TRUE
) {

  msg <- .msg(verbose)

  .validate_data_matrix(x)
  .validate_metadata(x, metadata)
  if (!is.character(method) || length(method) != 1L || !method %in% .drift_methods)
    stop("'method' must be one of: ", paste(.drift_methods, collapse = ", "), ".",
         call. = FALSE)
  .require_col(metadata, batch_col)
  .require_col(metadata, group_col)
  for (arg in c("log_scale", "remove_uncorrected", "qc_check", "use_parametric", "verbose"))
    if (!isTRUE(get(arg)) && !isFALSE(get(arg)))
      stop(sprintf("'%s' must be TRUE or FALSE.", arg), call. = FALSE)
  if (anyNA(metadata[[batch_col]]))
    stop(sprintf("Column '%s' contains missing values.", batch_col), call. = FALSE)

  mp       <- .as_matrix_preserve(x)
  x_matrix <- mp$mat
  storage.mode(x_matrix) <- "double"
  batch    <- factor(metadata[[batch_col]])
  is_qc    <- metadata[[group_col]] %in% qc_types

  params <- list(method = method, injection_col = injection_col,
                 batch_col = batch_col, group_col = group_col, qc_types = qc_types,
                 min_QC = min_QC, log_scale = log_scale,
                 remove_uncorrected = remove_uncorrected,
                 spline_smooth_param = spline_smooth_param, spar_limit = spar_limit,
                 gam_k = gam_k, sample_weight = sample_weight, qc_check = qc_check,
                 qc_check_p = qc_check_p, use_parametric = use_parametric)
  .result <- function(corrected, uncorrected, qc_tab = NULL) {
    unc_feat <- unique(uncorrected$Feature)
    if (remove_uncorrected && length(unc_feat) > 0L)
      corrected <- corrected[, !colnames(corrected) %in% unc_feat, drop = FALSE]
    structure(list(
      data                   = .restore_class(corrected, mp$was_matrix),
      data_before_correction = .restore_class(x_matrix, mp$was_matrix),
      method_used            = method,
      uncorrected            = uncorrected,
      uncorrected_features   = unc_feat,
      n_uncorrected          = length(unc_feat),
      qc_check               = qc_tab,
      parameters             = params
    ), class = "run_driftBatchCorrect")
  }
  no_uncorrected <- data.frame(Feature = character(), Batch = character(),
                               stringsAsFactors = FALSE)

  # Values usable for fitting: positive and observed; zeros and NA are skipped.
  usable <- !is.na(x_matrix) & x_matrix > 0
  z <- x_matrix
  z[!usable] <- NA
  if (log_scale) z <- log(z)

  # ---------------------------------------------------------------------------
  # ComBat: batch location/scale adjustment, no QCs needed
  # ---------------------------------------------------------------------------
  if (method == "ComBat") {
    if (!requireNamespace("sva", quietly = TRUE))
      stop("Package 'sva' is required for method = 'ComBat'. ",
           "Install with: BiocManager::install('sva')", call. = FALSE)
    if (nlevels(batch) < 2L) {
      warning("Only one batch detected. Returning original data.", call. = FALSE)
      return(.result(x_matrix, no_uncorrected))
    }
    msg(sprintf("Applying ComBat batch correction to %d batches...", nlevels(batch)))
    quiet <- if (verbose) identity else suppressMessages
    zc <- tryCatch(
      t(quiet(sva::ComBat(dat = t(z), batch = batch, mod = NULL,
                          par.prior = use_parametric, prior.plots = FALSE))),
      error = function(e) stop("ComBat correction failed: ", conditionMessage(e),
                               call. = FALSE))
    corrected <- if (log_scale) exp(zc) else zc
    corrected[!usable] <- x_matrix[!usable]   # zeros stay zero, NA stay NA
    dimnames(corrected) <- dimnames(x_matrix)
    msg("ComBat correction complete.")
    return(.result(corrected, no_uncorrected))
  }

  # ---------------------------------------------------------------------------
  # Trend-based (QC) methods: validation
  # ---------------------------------------------------------------------------
  .require_col(metadata, injection_col)
  inj <- metadata[[injection_col]]
  if (!is.numeric(inj) || anyNA(inj) || anyDuplicated(inj))
    stop(sprintf("'%s' must be numeric, complete, and unique.", injection_col),
         call. = FALSE)
  if (!is.numeric(min_QC) || length(min_QC) != 1L || is.na(min_QC) || min_QC < 4)
    stop("'min_QC' must be a single number >= 4 (smoothing needs at least 4 QC points).",
         call. = FALSE)
  if (!is.numeric(spline_smooth_param) || length(spline_smooth_param) != 1L ||
      is.na(spline_smooth_param) || spline_smooth_param < 0 || spline_smooth_param > 1)
    stop("'spline_smooth_param' must be a number between 0 and 1.", call. = FALSE)
  if (!is.numeric(spar_limit) || length(spar_limit) != 2L || anyNA(spar_limit))
    stop("'spar_limit' must be a numeric vector of length 2.", call. = FALSE)
  if (!is.numeric(gam_k) || length(gam_k) != 1L || is.na(gam_k) || gam_k < 3)
    stop("'gam_k' must be a single number >= 3.", call. = FALSE)
  if (!is.numeric(sample_weight) || length(sample_weight) != 1L ||
      is.na(sample_weight) || sample_weight < 0 || sample_weight > 1)
    stop("'sample_weight' must be a number between 0 and 1.", call. = FALSE)
  if (!is.numeric(qc_check_p) || length(qc_check_p) != 1L || is.na(qc_check_p) ||
      qc_check_p <= 0 || qc_check_p >= 1)
    stop("'qc_check_p' must be a number between 0 and 1.", call. = FALSE)
  is_gam <- method %in% c("tGAM", "rGAM")
  if (is_gam && !requireNamespace("mgcv", quietly = TRUE))
    stop("Package 'mgcv' is required for method = '", method, "'.", call. = FALSE)
  if (qc_check && !(is_gam && sample_weight > 0))
    stop("'qc_check = TRUE' needs method 'tGAM' or 'rGAM' and sample_weight > 0.",
         call. = FALSE)

  n_qc_batch <- tapply(is_qc, batch, sum)
  if (any(n_qc_batch < min_QC))
    stop(sprintf(
      "Batch(es) %s have fewer than min_QC = %d QC samples, so QC-based correction is impossible. Lower 'min_QC' (>= 4), merge batches, or use method = 'ComBat'.",
      paste(names(n_qc_batch)[n_qc_batch < min_QC], collapse = ", "), min_QC
    ), call. = FALSE)

  msg(sprintf("Applying %s drift/batch correction to %d features in %d batch(es)...",
              method, ncol(x_matrix), nlevels(batch)))

  # Robust GAM fit (tGAM: scaled-t errors; rGAM: one Tukey-biweight refit).
  .fit_gam <- function(d, k, fit_method = "REML") {
    if (method == "tGAM")
      return(mgcv::gam(y ~ s(x, k = k), data = d, weights = d$w,
                       family = mgcv::scat(), method = fit_method))
    fit <- mgcv::gam(y ~ s(x, k = k), data = d, weights = d$w, method = fit_method)
    r <- stats::residuals(fit, type = "response")
    s <- stats::mad(r)                      # mad() already includes 1.4826
    if (is.finite(s) && s > 0) {
      u  <- r / (4.685 * s)
      bw <- (1 - u^2)^2 * (abs(u) < 1)
      if (sum(bw > 0) >= k + 3) {
        d$w <- d$w * bw
        fit <- mgcv::gam(y ~ s(x, k = k), data = d, weights = d$w, method = fit_method)
      }
    }
    fit
  }

  # Trend of one feature in one batch at every injection order `t_all`, or NULL.
  # simplified: one loop over features x batches in R, no parallelism; ~1 ms per
  # spline/LOESS fit and ~0.1 s per GAM fit. Upgrade path: parallel::mclapply.
  .trend <- function(t_all, y_all, qc) {
    fit_qc <- qc & !is.na(y_all)
    if (sum(fit_qc) < min_QC) return(NULL)
    tq <- t_all[fit_qc]
    yq <- y_all[fit_qc]
    out <- list(trend = NULL, p = NA_real_, sample_based = FALSE)

    if (method == "QCRSC") {
      # smooth.spline()'s C code prints "spar-finding" notes straight to stderr;
      # capture them, as pmp::QCRSC() does.
      utils::capture.output(type = "message", sp <- if (spline_smooth_param == 0) {
        stats::smooth.spline(tq, yq, cv = TRUE,
                             control.spar = list(low = spar_limit[1], high = spar_limit[2]))
      } else {
        stats::smooth.spline(tq, yq, spar = spline_smooth_param)
      })
      out$trend <- stats::predict(sp, t_all)$y

    } else if (method == "rLOESS") {
      n <- length(tq)
      # neighbourhoods of >= 4 points; at most 20 candidate spans
      spans <- unique(round(seq(4, n, length.out = min(n - 3, 20)))) / n
      .lo <- function(span) stats::loess(yq ~ tq, span = span, degree = 2,
                                         family = "symmetric")
      gcv <- vapply(spans, function(sp) {
        f <- tryCatch(.lo(sp), error = function(e) NULL)
        if (is.null(f) || n - f$trace.hat <= 0) return(Inf)
        n * sum(stats::residuals(f)^2) / (n - f$trace.hat)^2
      }, numeric(1))
      if (all(!is.finite(gcv))) return(NULL)
      f <- .lo(spans[which.min(gcv)])
      # simplified: loess cannot extrapolate, so injections outside the QC range
      # get the nearest QC's trend value (flat ends); fine when batches start
      # and end with QCs, as recommended.
      out$trend <- stats::predict(f, data.frame(tq = pmin(pmax(t_all, min(tq)), max(tq))))

    } else {
      fit_bio <- !qc & !is.na(y_all)
      w_bio   <- if (sample_weight > 0) sample_weight else 0
      use     <- fit_qc | (fit_bio & w_bio > 0)
      d <- data.frame(x = t_all[use], y = y_all[use],
                      w = ifelse(qc[use], 1, w_bio), qc = factor(qc[use]))
      k <- min(gam_k, nrow(d) - 3, length(unique(d$x)) - 1)
      if (k < 3) return(NULL)
      fit <- .fit_gam(d, k)
      k_chk <- min(k, sum(fit_qc) - 3, sum(fit_bio) - 3)
      if (qc_check && k_chk >= 3) {
        # Likelihood ratio test (ML fits): one trend vs separate QC/sample trends.
        # A failed test keeps the combined fit.
        fam <- if (method == "tGAM") mgcv::scat() else stats::gaussian()
        out$p <- tryCatch({
          m0 <- mgcv::gam(y ~ s(x, k = k_chk), data = d, family = fam, method = "ML")
          m1 <- mgcv::gam(y ~ qc + s(x, k = k_chk, by = qc), data = d,
                          family = fam, method = "ML")
          stats::anova(m0, m1, test = "Chisq")[["Pr(>Chi)"]][2]
        }, error = function(e) NA_real_)
        if (!is.na(out$p) && out$p < qc_check_p) {
          ds <- d[d$qc == "FALSE", , drop = FALSE]
          ds$w <- 1
          fit <- .fit_gam(ds, min(k, nrow(ds) - 3))
          out$sample_based <- TRUE
        }
      }
      out$trend <- as.numeric(stats::predict(fit, newdata = data.frame(x = t_all)))
    }
    out
  }

  corrected <- x_matrix
  ref <- matrixStats::colMedians(replace(x_matrix, !usable, NA), na.rm = TRUE)
  failed <- matrix(FALSE, ncol(x_matrix), nlevels(batch),
                   dimnames = list(colnames(x_matrix), levels(batch)))
  qc_rows <- list()

  for (b in levels(batch)) {
    rows <- which(batch == b)
    for (j in seq_len(ncol(x_matrix))) {
      fit <- tryCatch(suppressWarnings(.trend(inj[rows], z[rows, j], is_qc[rows])),
                      error = function(e) NULL)
      f <- if (is.null(fit$trend)) NULL else if (log_scale) exp(fit$trend) else fit$trend
      if (is.null(f) || any(!is.finite(f)) || any(f <= 0) || !is.finite(ref[j])) {
        failed[j, b] <- TRUE
        next
      }
      corrected[rows, j] <- x_matrix[rows, j] * ref[j] / f
      if (qc_check)
        qc_rows[[length(qc_rows) + 1L]] <- data.frame(
          Feature = colnames(x_matrix)[j], Batch = b, p_value = fit$p,
          Sample_based = fit$sample_based, stringsAsFactors = FALSE)
    }
  }

  idx <- which(failed, arr.ind = TRUE)
  uncorrected <- data.frame(Feature = rownames(failed)[idx[, 1]],
                            Batch   = colnames(failed)[idx[, 2]],
                            stringsAsFactors = FALSE)
  uncorrected <- uncorrected[order(match(uncorrected$Feature, colnames(x_matrix))), ,
                             drop = FALSE]
  rownames(uncorrected) <- NULL
  qc_tab <- if (qc_check) do.call(rbind, c(qc_rows, list(data.frame(
    Feature = character(), Batch = character(), p_value = numeric(),
    Sample_based = logical(), stringsAsFactors = FALSE)))) else NULL

  out <- .result(corrected, uncorrected, qc_tab)
  if (out$n_uncorrected > 0L)
    msg(sprintf("%d feature(s) could not be corrected in at least one batch (%s).",
                out$n_uncorrected, if (remove_uncorrected) "removed" else "kept"))
  if (qc_check)
    msg(sprintf("QC check: trend refitted on biological samples in %d feature-batch fit(s).",
                sum(qc_tab$Sample_based)))
  msg(sprintf("%s correction complete.", method))
  out
}

#' @param ... Ignored.
#' @rdname run_driftBatchCorrect
#' @export
print.run_driftBatchCorrect <- function(x, ...) {
  cat("=== Drift/Batch Correction Results ===\n")
  cat(sprintf("Method used: %s\n", x$method_used))
  if (x$n_uncorrected > 0L)
    cat(sprintf("Uncorrected features: %d (%s)\n", x$n_uncorrected,
                if (isTRUE(x$parameters$remove_uncorrected)) "removed" else "kept"))
  if (!is.null(x$qc_check))
    cat(sprintf("QC check: %d of %d fits refitted on biological samples\n",
                sum(x$qc_check$Sample_based), nrow(x$qc_check)))
  cat(sprintf("Data: %d samples x %d features\n", nrow(x$data), ncol(x$data)))
  invisible(x)
}
