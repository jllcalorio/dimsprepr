# Example Direct-Injection Metabolomics Data Matrix

A simulated numeric matrix of 50 samples x 80 features representing raw
peak intensities from a direct-injection mass spectrometry experiment.
Rows are samples; columns are features (m/z bins). Zero values represent
below-detection-limit measurements. Features 78-80 are intentionally
sparse to exercise missing-value filtering.

## Usage

``` r
dims_data
```

## Format

A numeric matrix with 50 rows and 80 columns. Row names are sample
identifiers (`S01`–`S50`); column names are feature labels
(`mz_0001`–`mz_0080`).

## Source

Simulated data. See `data-raw/dims_example.R` for the generation script.

## See also

[dims_metadata](https://jllcalorio.github.io/dimsprepr/reference/dims_metadata.md),
[`run_DIpreprocess()`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.md)

## Examples

``` r
data(dims_data)
dim(dims_data)
#> [1] 50 80
```
