# Normalize Metabolomics Data

Applies various normalization methods to account for dilution effects
and sample-to-sample variation in metabolomics data.

## Usage

``` r
run_normalize(
  x,
  metadata,
  method = "sum",
  factor_col = "Normalization",
  ref_sample = NULL,
  group_sample = "QC",
  qc_normalize = "median",
  group_col = "Group",
  qc_types = c("QC", "SQC", "EQC"),
  reference_method = "mean",
  sample_id_col = "Sample",
  verbose = TRUE
)

# S3 method for class 'run_normalize'
print(x, ...)
```

## Arguments

- x:

  Object to print.

- metadata:

  Data frame. Sample metadata with number of rows equal to nrow(x).

- method:

  Character or numeric vector. Normalization method:

  - `"sum"`: Total sum normalization

  - `"median"`: Median normalization

  - `"specific_factor"`: Use values from a metadata column

  - `"pqn_global"`: PQN using global median

  - `"pqn_reference"`: PQN using a specific reference sample

  - `"pqn_group"`: PQN using pooled QC samples as reference

  - `"pqn"`: PQN via the pmp package (requires pmp)

  - `"quantile"`: Alias for `"pqn"` (PQN via pmp)

  - `"col_rel_abundance"`: Relative abundance per column (feature)

  - `"row_rel_abundance"`: Relative abundance per row (sample)

  - Numeric vector: Custom normalization factors

- factor_col:

  Character. Metadata column with normalization factors. Default:
  "Normalization".

- ref_sample:

  Character. Reference sample for PQN. Default: NULL.

- group_sample:

  Character. Group label for PQN group reference. Default: "QC".

- qc_normalize:

  Character. QC normalization strategy: "mean", "median", "none".
  Default: "median".

- group_col:

  Character. Column with group labels. Default: "Group".

- qc_types:

  Character vector. QC group labels. Default: c("QC", "SQC", "EQC").

- reference_method:

  Character. Reference method for PQN via pmp: "mean" or "median".
  Default: "mean".

- sample_id_col:

  Character. Column with sample identifiers. Default: "Sample".

- verbose:

  Logical. Print progress. Default: TRUE.

- ...:

  Ignored.

## Value

A list of class `"run_normalize"` containing:

- data:

  Normalized data (same class as input `x`)

- normalization_factors:

  Numeric vector of factors used

- method_used:

  Character describing method

- parameters:

  List of parameters

## Details

See the package documentation for full descriptions of each method.
`"pqn"` and `"quantile"` both call
[`pmp::pqn_normalisation()`](https://rdrr.io/pkg/pmp/man/pqn_normalisation.html).

## References

Dieterle, F., et al. (2006). Analytical Chemistry, 78(13), 4281-4290.
[doi:10.1021/ac051632c](https://doi.org/10.1021/ac051632c)

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.md)

## Author

John Lennon L. Calorio

## Examples

``` r
# \donttest{
set.seed(215)
x <- matrix(abs(rnorm(80 * 40, mean = 100, sd = 20)), nrow = 80, ncol = 40)
colnames(x) <- paste0("Feature", 1:40)
meta <- data.frame(
  Sample = paste0("S", 1:80),
  Group = rep(c("Control", "Treatment", "QC"), c(30, 30, 20))
)

result <- run_normalize(x, meta, method = "median")
#> Applying 'median' normalization...
#> Normalizing by median...
#> Normalization complete.
result
#> === Normalization Results ===
#> Method: median
#> Data dimensions: 80 samples x 40 features
# }
```
