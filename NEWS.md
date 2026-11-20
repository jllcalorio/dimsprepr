# dimsprepr 1.0.0

First CRAN release.

## Pipeline changes in `run_DIpreprocess()`

* Quality filtering (QC-RSD and variance) now runs **before** scaling: RSD on
  the normalized data and variance on the transformed data. Before, filtering
  ran on scaled data, where QC means are near 0 (so RSD is meaningless) and
  every auto-scaled feature has variance 1 (so the variance filter removed an
  arbitrary set of features). Both branches now share one feature set, so the
  `scale_filter_ref` argument was removed.
* `rsd_qc_type` is now honoured (it was ignored before). Its default is `"QC"`
  (all `qc_types`), which matches the previous effective behaviour.
* Unknown method names and `"pqn_reference"` without `normalize_ref_sample`
  are rejected before any processing starts.
* `merge_replicates = TRUE` now works when `subject_id_col` is not named
  `"SubjectID"`.
* Supervised scaling (`"svast"`, `"xvast"`) now receives the group labels.
* Batch labels may be character (e.g. `"B1"`).
* `vsn_cores` was removed (see `run_transform()`).

## Step functions

* `run_filtermissing()`: `include_QC = FALSE` now also excludes QC samples in
  global mode (`filter_by_group = FALSE`); fixed a shape error with a single
  feature.
* `run_filterRSD()`: the `qc_type` default is now `"QC"`. Features with
  undefined RSD are removed (as in `pmp::filter_peaks_by_rsd()`), `rsd_values`
  is named, and a warning is given when QC means are <= 0 (centred data).
* `run_normalize()`: `group_sample` is now used by `"pqn_group"`;
  `qc_normalize` is validated; `"pqn"` reports its normalization factors;
  documented that `"quantile"` is an alias of PQN, not quantile normalization.
* `run_transform()`: removed `num_cores`. `vsn` never used the BiocParallel
  backend it registered, so a SNOW cluster was started for nothing.
  `method_used` now reports a fallback to log10.
* `run_driftBatchCorrect()`: an actionable error is raised when QCRSC cannot
  fit splines (too few QCs per batch); pmp messages respect `verbose`.

## Plotting

* `plot_beforeafter()` on a `run_DIpreprocess` result now compares imputed
  with drift-corrected data (its injection-order view is a drift diagnostic);
  previously it log-transformed auto-scaled data.
* Both plot functions now share the theme and feature-selection helpers.

## Packaging

* `BiocParallel` and `parallelly` are no longer dependencies; `imputeLCMD`,
  `mice`, `sva`, and `vsn` moved to Suggests.
* Example data now has QC samples in every batch, so all examples run,
  including drift correction.
* Added references (with DOIs) for every method, runnable examples, and tests
  for every exported function.
