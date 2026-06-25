# Transform Metabolomics Data

Applies various transformation methods to stabilize variance and reduce
heteroscedasticity in metabolomics data.

## Usage

``` r
run_transform(
  x,
  method = "log2",
  metadata = NULL,
  group_col = "Group",
  qc_types = c("QC", "SQC", "EQC"),
  num_cores = 1,
  verbose = TRUE
)

# S3 method for class 'run_transform'
print(x, ...)
```

## Arguments

- x:

  Object to print.

- method:

  Character. Transformation method to apply. Options:

  - `"log2"`: Log base 2 transformation

  - `"log10"`: Log base 10 transformation

  - `"sqrt"`: Square root transformation

  - `"cbrt"`: Cube root transformation

  - `"clr"`: Centered log-ratio transformation (computed per sample)

  - `"arcsin_sqrt"`: Arcsine square root transformation (data must be in
    \[0, 1\])

  - `"vsn"`: Variance Stabilizing Normalization (requires vsn)

  - `"glog"`: Generalized logarithm transformation (requires pmp)

- metadata:

  Data frame. Sample metadata, required for glog. Default: NULL.

- group_col:

  Character. Column in `metadata` with group labels. Default: "Group".

- qc_types:

  Character vector. QC group labels. Default: c("QC", "SQC", "EQC").

- num_cores:

  Integer or "max". Cores for VSN. Default: 1.

- verbose:

  Logical. Print progress. Default: TRUE.

- ...:

  Ignored.

## Value

A list of class `"run_transform"` containing:

- data:

  Transformed data (same class as input `x`)

- method_used:

  Character describing method

- shift_applied:

  Numeric shift added before transformation

- parameters:

  List of parameters

## Details

**log2/log10**: Logarithmic compression. A shift is added for
non-positive values.

**sqrt/cbrt**: Power transformations. A shift is added for negative
values.

**clr**: Centered log-ratio. Maps compositional data to real space.

**arcsin_sqrt**: For proportion data in \[0, 1\].

**vsn**: Model-based variance stabilization. Requires vsn and
BiocParallel.

**glog**: Generalized logarithm via pmp. Requires metadata with QC
labels.

## References

Huber, W., et al. (2002). Bioinformatics, 18(Suppl 1), S96-S104.
[doi:10.1093/bioinformatics/18.suppl_1.s96](https://doi.org/10.1093/bioinformatics/18.suppl_1.s96)

Parsons, H.M., et al. (2007). BMC Bioinformatics, 8, 234.
[doi:10.1186/1471-2105-8-234](https://doi.org/10.1186/1471-2105-8-234)

## See also

[`run_DIpreprocess`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.md)

## Author

John Lennon L. Calorio

## Examples

``` r
# \donttest{
set.seed(937)
x <- matrix(abs(rnorm(80 * 40, mean = 100, sd = 50)), nrow = 80, ncol = 40)
colnames(x) <- paste0("Feature", 1:40)

result <- run_transform(x, method = "log2")
#> Applying 'log2' transformation...
#> Transformation complete.
result
#> === Transformation Results ===
#> Method: log2
#> Data dimensions: 80 samples x 40 features
# }
```
