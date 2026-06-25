# Scale Metabolomics Data

Applies various scaling methods to standardize features and prepare data
for multivariate analysis. Includes classical scaling methods as well as
variance-stability (VAST) scaling variants.

## Usage

``` r
run_scale(
  x,
  method = "auto",
  metadata = NULL,
  group_col = "Group",
  verbose = TRUE
)

# S3 method for class 'run_scale'
print(x, ...)
```

## Arguments

- x:

  Object to print.

- method:

  Character. Scaling method:

  - `"mean"`: Mean-centering only

  - `"auto"`: Auto-scaling (mean-centering + unit variance)

  - `"pareto"`: Pareto-scaling (mean-centering + sqrt(SD))

  - `"ln"`: Natural log transformation

  - `"vast"`: VAST scaling (auto-scaling weighted by CV)

  - `"svast"`: Supervised VAST (weighted by mean class CV)

  - `"xvast"`: Extended supervised VAST (weighted by max class CV)

  - `"none"`: No scaling

- metadata:

  Data frame. Required for `"svast"` and `"xvast"`. Default: NULL.

- group_col:

  Character. Column in `metadata` with group labels. Required for
  supervised methods. Default: "Group".

- verbose:

  Logical. Print progress. Default: TRUE.

- ...:

  Ignored.

## Value

A list of class `"run_scale"` containing:

- data:

  Scaled data (same class as input `x`)

- scaling_factors:

  Numeric vector of per-feature scaling factors

- center_values:

  Numeric vector of per-feature means

- method_used:

  Character

- parameters:

  List of parameters

## Details

**auto**: Recommended for PCA, clustering, t-tests, ANOVA.

**pareto**: Recommended for PLS-DA, OPLS-DA.

**vast/svast/xvast**: Variance-stability methods from Keun et al. (2003)
and Yang et al. (2015).

## References

van den Berg, R.A., et al. (2006). BMC Genomics, 7, 142.
[doi:10.1186/1471-2164-7-142](https://doi.org/10.1186/1471-2164-7-142)

Yang, J., et al. (2015). Frontiers in Molecular Biosciences, 2, 4.
[doi:10.3389/fmolb.2015.00004](https://doi.org/10.3389/fmolb.2015.00004)

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.md)

## Author

John Lennon L. Calorio

## Examples

``` r
# \donttest{
set.seed(662)
x <- matrix(rnorm(100 * 50, mean = 100, sd = 50), nrow = 100, ncol = 50)
colnames(x) <- paste0("Feature", 1:50)

result <- run_scale(x, method = "auto")
#> Applying 'auto' scaling...
#> Scaling complete.
result
#> === Scaling Results ===
#> Method  : auto
#> Data    : 100 samples x 50 features
# }
```
