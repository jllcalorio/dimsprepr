# Correct Signal Drift and Batch Effects Using QC Samples

Applies Quality Control-based Robust Spline Correction (QCRSC) or batch
correction using ComBat to remove systematic signal drift and batch
effects in metabolomics data.

## Usage

``` r
run_driftBatchCorrect(
  x,
  metadata,
  perform_correction = TRUE,
  batch_corr_only = FALSE,
  injection_col = "InjectionSequence",
  batch_col = "Batch",
  group_col = "Group",
  qc_label = "QC",
  qc_types = c("QC", "SQC", "EQC"),
  spline_smooth_param = 0,
  min_QC = 5,
  spar_limit = c(-1.5, 1.5),
  log_scale = TRUE,
  use_parametric = TRUE,
  display_plots = FALSE,
  verbose = TRUE
)

# S3 method for class 'run_driftBatchCorrect'
print(x, ...)
```

## Arguments

- x:

  Object to print.

- metadata:

  Data frame. Sample metadata with number of rows equal to nrow(x).

- perform_correction:

  Logical. If TRUE, perform correction; if FALSE, return original data
  unchanged. Default: TRUE.

- batch_corr_only:

  Logical. If TRUE, perform only ComBat batch correction (no drift
  correction). If FALSE, perform QCRSC. Default: FALSE.

- injection_col:

  Character. Column in `metadata` containing injection order. Default:
  "InjectionSequence".

- batch_col:

  Character. Column in `metadata` containing batch numbers. Default:
  "Batch".

- group_col:

  Character. Column in `metadata` containing sample group labels.
  Default: "Group".

- qc_label:

  Character. Group label identifying QC samples for QCRSC. Default:
  "QC".

- qc_types:

  Character vector. Group labels converted to `qc_label` internally.
  Default: c("QC", "SQC", "EQC").

- spline_smooth_param:

  Numeric. Smoothing parameter (0-1). Default: 0.

- min_QC:

  Integer. Minimum QC samples per batch. Default: 5.

- spar_limit:

  Numeric vector of length 2. Spline limits. Default: c(-1.5, 1.5).

- log_scale:

  Logical. Fit on log-transformed data. Default: TRUE.

- use_parametric:

  Logical. Parametric ComBat adjustment. Default: TRUE.

- display_plots:

  Logical. Diagnostic plots for ComBat. Default: FALSE.

- verbose:

  Logical. Print progress. Default: TRUE.

- ...:

  Ignored.

## Value

A list of class `"run_driftBatchCorrect"` containing:

- data:

  Corrected data (same class as input `x`)

- data_before_correction:

  Original data before correction

- correction_applied:

  Logical

- method_used:

  Character

- uncorrected_features:

  Character vector of uncorrected features

- n_uncorrected:

  Integer

- parameters:

  List of parameters

## Details

**QCRSC** fits smoothing splines through QC intensities across injection
order. Requires pmp.

**ComBat** uses empirical Bayes for batch correction. Requires sva.

## References

Kirwan, J.A., et al. (2013). Analytical and Bioanalytical Chemistry,
405, 5147-5157.
[doi:10.1007/s00216-013-6856-7](https://doi.org/10.1007/s00216-013-6856-7)

Johnson, W.E., et al. (2007). Biostatistics, 8(1), 118-127.
[doi:10.1093/biostatistics/kxj037](https://doi.org/10.1093/biostatistics/kxj037)

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.md)

## Author

John Lennon L. Calorio

## Examples

``` r
# \donttest{
set.seed(417)
n <- 60; p <- 30
x <- matrix(abs(rnorm(n * p, 100, 20)), nrow = n, ncol = p)
colnames(x) <- paste0("F", seq_len(p))

meta <- data.frame(
  InjectionSequence = seq_len(n),
  Batch = rep(1:2, each = 30),
  Group = rep(c(rep("Sample", 9), "QC"), 6)
)

# Skip correction (returns unchanged)
result <- run_driftBatchCorrect(x, meta, perform_correction = FALSE)
#> Correction disabled. Returning original data.
result
#> === Drift/Batch Correction Results ===
#> Correction applied: FALSE
#> Method used: None
# }
```
