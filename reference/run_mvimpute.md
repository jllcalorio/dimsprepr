# Impute Missing Values in Metabolomics Data

Replaces missing values (NAs) using either a simple deterministic method
(fraction of minimum), a statistical imputation method (quantile
regression imputation of left-censored data, QRILC), or any method
available in [`mice`](https://amices.org/mice/reference/mice.html)
(Multiple Imputation by Chained Equations).

## Usage

``` r
run_mvimpute(
  x,
  method = 0.2,
  positive_only = TRUE,
  tune_sigma = 1,
  m = 5,
  maxit = 5,
  seed = NA,
  verbose = TRUE,
  ...
)

# S3 method for class 'run_mvimpute'
print(x, ...)
```

## Arguments

- x:

  Object to print.

- method:

  Numeric or character. Controls the imputation strategy:

  - **Numeric** – fraction of the smallest observed value per feature
    used for deterministic imputation (e.g., `0.2` = 1/5th of the
    smallest value). By default (`positive_only = TRUE`), the minimum is
    computed from strictly positive values only.

  - **`"quantileregression"`** – uses the QRILC algorithm from the
    imputeLCMD package. Requires imputeLCMD.

  - **Any `mice` method string** – e.g., `"pmm"`, `"norm"`, `"rf"`,
    `"cart"`, etc. Requires mice.

  - **`"none"`** – Skips imputation entirely.

  Default: `0.2`.

- positive_only:

  Logical. Only used when `method` is numeric. If `TRUE` (default),
  per-feature minimum is computed from strictly positive values.

- tune_sigma:

  Numeric. Only used when `method = "quantileregression"`. Controls the
  standard deviation of the left-censored distribution. Default: `1`.

- m:

  Integer. Only used with mice methods. Number of imputations. Default:
  `5`.

- maxit:

  Integer. Only used with mice methods. Number of iterations. Default:
  `5`.

- seed:

  Integer or `NA`. Only used with mice methods. Random seed. Default:
  `NA`.

- verbose:

  Logical. Print progress messages. Default: `TRUE`.

- ...:

  Ignored.

## Value

A list of class `"run_mvimpute"` containing:

- data:

  Matrix or data frame with imputed values (same class as input `x`).

- n_missing_before:

  Integer. Number of missing values before imputation.

- n_missing_after:

  Integer. Number of missing values after imputation.

- imputed_summary:

  Data frame with per-feature imputation statistics.

- parameters:

  List of parameters used.

- method_used:

  Character string describing the imputation method applied.

## Details

**Missing Value Imputation Strategies:**

1.  **Deterministic (fraction of minimum):** Each missing value is
    replaced by a fraction of the per-feature minimum. Standard in
    metabolomics for MNAR data.

2.  **Quantile Regression (QRILC):** Models missing values from a
    left-censored distribution. Requires imputeLCMD.

3.  **mice:** Multiple Imputation by Chained Equations. Requires mice.

## References

Van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate
Imputation by Chained Equations in R. *Journal of Statistical Software*,
45(3), 1-67.
[doi:10.18637/jss.v045.i03](https://doi.org/10.18637/jss.v045.i03)

Wei, R., et al. (2018). Missing Value Imputation Approach for Mass
Spectrometry-based Metabolomics Data. *Scientific Reports*, 8(1), 663.
[doi:10.1038/s41598-017-19120-0](https://doi.org/10.1038/s41598-017-19120-0)

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.md),
[`mice`](https://amices.org/mice/reference/mice.html),
[`impute.QRILC`](https://rdrr.io/pkg/imputeLCMD/man/impute.QRILC.html)

## Author

John Lennon L. Calorio

## Examples

``` r
# \donttest{
set.seed(723)
x <- matrix(abs(rnorm(100 * 50, mean = 100)), nrow = 100, ncol = 50)
x[sample(length(x), 200)] <- NA
colnames(x) <- paste0("Feature", 1:50)

result <- run_mvimpute(x, method = 0.2)
#> Starting missing value imputation (200 missing values)...
#> Using deterministic imputation: 0.2000 x minimum positive value per feature.
#> Imputation complete. Missing values: 200 -> 0.
result
#> === Missing Value Imputation Results ===
#> Method: Deterministic: 0.2000 x min positive
#> Missing values before: 200
#> Missing values after:  0
# }
```
