# dimsprepr 1.0.0

First CRAN release. Changes relative to the earlier GitHub versions:

## New methods

* `run_driftBatchCorrect()` gains `method`: `"QCRSC"` (default), the
  outlier-robust `"rLOESS"`, `"rGAM"`, and `"tGAM"` of Vynck et al. (2026),
  and `"ComBat"`. rGAM and tGAM can also use the biological samples, with a
  lower weight (`sample_weight`), and test whether the QCs are representative
  (`qc_check`, `qc_check_p`). They need the suggested package mgcv.
* QComics (González-Domínguez et al., 2024):
  * New `run_filterblank()` flags potential contaminants (signal in samples
    below 3 times the blank signal), removes those with an unstable blank
    signal, and drops the blanks (`blank_types` in `run_DIpreprocess()`).
  * `run_filtermissing(sample_iqr)` removes samples with abnormally high
    missingness, such as failed injections, before filtering features
    (`missing_sample_iqr = 3` in `run_DIpreprocess()`). It also reports
    `sample_missing_rate`.
  * `run_filtermissing(absent_threshold)` keeps features that are mostly
    missing in a group as truly absent (zeros) instead of removing them
    (`missing_absent_threshold` in `run_DIpreprocess()`).
* `run_transform()` gains `"ln"` (natural log), moved from `run_scale()`.

## Pipeline changes in `run_DIpreprocess()`

* New step order: drift and batch correction now runs on data that still
  contain missing values, and imputation follows it, once. Imputed QC
  values no longer bias the QC trends, so `impute_method_post_drift` was
  removed.
* Drift options are now `drift_method` (`"none"` replaces
  `correct_drift = FALSE`) and `drift_args`, a list of further
  `run_driftBatchCorrect()` arguments. This replaces `correct_drift`,
  `spline_smooth`, `spline_spar_limit`, `correct_on_log`, and
  `min_qc_per_batch`.
* `remove_uncorrected` now defaults to `TRUE` and works per batch: a feature
  with too few measured QC values in any batch is removed or, with `FALSE`,
  kept with its original values in that batch. Previously the batch was
  silently set to `NA` and re-imputed, creating artificial low values.
* Quality filtering (QC-RSD and variance) runs before scaling: RSD on the
  measured (non-imputed) normalized values, and variance on the transformed
  data. Before, it ran on scaled data, where QC means are near 0 and every
  auto-scaled feature has variance 1. Both branches share one feature set,
  so `scale_filter_ref` was removed.
* `rsd_qc_type` is honoured and may name several QC labels.
* Unknown method names, invalid `drift_args`, and `"pqn_reference"` without
  `normalize_ref_sample` are rejected before any processing starts.
* `merge_replicates = TRUE` works when `subject_id_col` is not named
  `"SubjectID"`; supervised scaling (`"svast"`, `"xvast"`) receives the group
  labels; batch labels may be character.
* The result gains `samples_removed`, and its data elements follow the step
  order: `data_raw`, `data_missing_filtered`, `data_corrected`,
  `data_imputed`, and so on.

## Step functions

* `run_driftBatchCorrect()`: QC-RSC is now implemented directly (it matches
  `pmp::QCRSC()`), so a failed spline fit affects one feature and batch only
  instead of aborting the run. Removed `perform_correction`,
  `batch_corr_only` (use `method = "ComBat"`), `qc_label`, and
  `display_plots`. ComBat runs on log intensities when `log_scale = TRUE`.
* `run_filtermissing()`: a feature observed in exactly 80% of samples is now
  kept, as the 80% rule intends (it was removed); `include_QC = FALSE` also
  applies in global mode; fixed a single-feature shape error.
* `run_filterRSD()`: `qc_type` was removed; `qc_types` selects the QC samples.
  Features with undefined RSD are removed (as in `pmp::filter_peaks_by_rsd()`);
  a warning is given when QC means are <= 0.
* `run_normalize()`: removed the duplicate methods `"pqn"` and `"quantile"`
  (use `"pqn_group"`, with `reference_method = "mean"` for a mean reference)
  and `"row_rel_abundance"` (identical to `"sum"`), and `"col_rel_abundance"`
  (a feature scaling, not a normalization). `reference_method` now applies to
  every PQN variant and defaults to `"median"`. `group_sample` is used by
  `"pqn_group"`.
* `run_mvimpute()`: removed `m` (a single completed data set is returned, so
  more imputations only cost time); mice now copes with non-syntactic feature
  names.
* `run_scale()`: removed `"ln"` (a transformation; see `run_transform()`);
  VAST weights warn when feature means are <= 0; `parameters` reports
  `group_col` and the group levels separately.
* `run_transform()`: removed `num_cores` (vsn never used the parallel
  backend that was registered); glog now needs pmp explicitly.

## Plotting

* `base_size` defaults to 11 (ggplot2's default, which suits the 7 x 5 inch
  figures of R Markdown and Quarto) and all text sizes follow it:
  `title_size`, `xlab_size`, and `ylab_size` default to `base_size`, and
  `subtitle_size` to `1.2 * base_size`. `global_font_size` was removed.
* Like ggplots, the plots are now drawn when printed (automatically when
  not assigned) instead of on every call, so `p <- plot_beforeafter(...)`
  followed by `p` no longer draws twice. Inside loops, use `print()`.
  `plot_beforeafter()` selects 4 random features by default, and
  `plot_dist_beforeafter()` shows 20 box plots by default.
* `plot_beforeafter()` on a `run_DIpreprocess` result compares the data
  before and after drift correction.

## Packaging

* pmp, BiocParallel, and parallelly are no longer required: the default
  pipeline installs without Bioconductor. imputeLCMD, mgcv, mice, pmp, sva,
  and vsn are suggested for optional methods.
* Example data has QC samples in every batch, so all examples run, including
  drift correction.
* References with DOIs for every method, runnable examples, and tests for
  every exported function.
