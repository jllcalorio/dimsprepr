# Example Sample Metadata for `dims_data`

A data frame of sample-level metadata matching the rows of
[`dims_data`](https://jllcalorio.github.io/dimsprepr/reference/dims_data.md).
Contains group labels, batch assignments, injection order, and subject
identifiers.

## Usage

``` r
dims_metadata
```

## Format

A data frame with 50 rows and 5 columns:

- Sample:

  Character. Unique sample identifier (`S01`–`S50`).

- Group:

  Character. Sample group: `"Control"`, `"Treatment"`, or `"QC"`.

- Batch:

  Integer. Batch number (1 or 2).

- InjectionSequence:

  Integer. Order of injection (1–50).

- SubjectID:

  Character. Biological subject ID for non-QC samples; `NA` for QC.

## Source

Simulated data. See `data-raw/dims_example.R`.

## See also

[dims_data](https://jllcalorio.github.io/dimsprepr/reference/dims_data.md),
[`run_DIpreprocess()`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.md)

## Examples

``` r
data(dims_metadata)
table(dims_metadata$Group)
#> 
#>   Control        QC Treatment 
#>        20        10        20 
```
