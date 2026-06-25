# Filter Features by Low Variance

Removes features with low variance across samples. Features with minimal
variation provide little information for downstream analysis and can be
safely removed to reduce dimensionality.

## Usage

``` r
run_filtervariance(x, percentile = 10, verbose = TRUE)

# S3 method for class 'run_filtervariance'
print(x, ...)
```

## Arguments

- x:

  Object to print.

- percentile:

  Numeric. Percentile threshold (0-100) for variance filtering. Features
  in the bottom X percentile of variance are removed. For example, 10
  removes the 10% of features with lowest variance. Default: 10.

- verbose:

  Logical. Print progress messages. Default: TRUE.

- ...:

  Ignored.

## Value

A list of class `"run_filtervariance"` containing:

- data:

  Matrix or data frame of filtered data (same class as input `x`)

- features_removed:

  Character vector of removed feature names

- features_kept:

  Character vector of retained feature names

- variance_values:

  Numeric vector of variance values for all original features

- variance_threshold:

  Numeric. The variance value corresponding to the percentile cutoff

- n_features_before:

  Integer. Number of features before filtering

- n_features_after:

  Integer. Number of features after filtering

- n_features_removed:

  Integer. Number of features removed

- parameters:

  List of parameters used

## Details

**Why Filter by Variance:**

Features with very low variance across samples:

- Provide minimal discriminatory power for classification

- Contribute little to multivariate models (PCA, PLS-DA)

- May represent technical noise or detection artifacts

- Increase computational burden without adding information

**Choosing the Percentile:**

- **Conservative** (5-10%): Remove only the least variable features

- **Moderate** (10-20%): Standard approach for most studies

- **Aggressive** (20-30%): When high dimensionality is a concern

## References

Broadhurst, D.I. (2025). QC:MXP Repeat Injection based Quality Control,
Batch Correction, Exploration & Data Cleaning (Version 2.1) Zendono.
[doi:10.5281/zenodo.16824822](https://doi.org/10.5281/zenodo.16824822) .
Retrieved from <https://github.com/broadhurstdavid/QC-MXP>.

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.md)

## Author

John Lennon L. Calorio

## Examples

``` r
# \donttest{
set.seed(519)
x <- matrix(rnorm(100 * 50, mean = 100, sd = seq(0.1, 25, length.out = 50)),
            nrow = 100, ncol = 50, byrow = TRUE)
colnames(x) <- paste0("Feature", 1:50)

result <- run_filtervariance(x, percentile = 10)
#> Filtering features by variance (10th percentile cutoff)...
#> Removed 5 features (10.0%) with variance <= 6.884355
result
#> === Variance Filtering Results ===
#> Features before: 50
#> Features after:  45
#> Features removed: 5 (10.0%)
#> Percentile cutoff: 10%
#> Variance threshold: 6.884355
# }
```
