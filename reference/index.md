# Package index

## Direct-Injection Metabolomics Preprocessing Pipeline

Function for running a sequential preprocessing pipeline for
metabolomics data by orchestrating the individual run\_\* functions of
dimsprepr.

- [`run_DIpreprocess()`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.md)
  : Direct-Injection Metabolomics Preprocessing Pipeline

## Data Cleaning & Quality Control

Filtering features based on missingness, variance, or RSD thresholds.

- [`run_filtermissing()`](https://jllcalorio.github.io/dimsprepr/reference/run_filtermissing.md)
  [`print(`*`<run_filtermissing>`*`)`](https://jllcalorio.github.io/dimsprepr/reference/run_filtermissing.md)
  : Filter Features by Missing Value Threshold
- [`run_filterRSD()`](https://jllcalorio.github.io/dimsprepr/reference/run_filterRSD.md)
  [`print(`*`<run_filterRSD>`*`)`](https://jllcalorio.github.io/dimsprepr/reference/run_filterRSD.md)
  : Filter Features by Relative Standard Deviation in QC Samples
- [`run_filtervariance()`](https://jllcalorio.github.io/dimsprepr/reference/run_filtervariance.md)
  [`print(`*`<run_filtervariance>`*`)`](https://jllcalorio.github.io/dimsprepr/reference/run_filtervariance.md)
  : Filter Features by Low Variance

## Missing Value Imputation

Addressing missing data using various imputation algorithms.

- [`run_mvimpute()`](https://jllcalorio.github.io/dimsprepr/reference/run_mvimpute.md)
  [`print(`*`<run_mvimpute>`*`)`](https://jllcalorio.github.io/dimsprepr/reference/run_mvimpute.md)
  : Impute Missing Values in Metabolomics Data

## Signal Drift and Batch Effects Correction

Correcting for technical variation and instrument signal drift.

- [`run_driftBatchCorrect()`](https://jllcalorio.github.io/dimsprepr/reference/run_driftBatchCorrect.md)
  [`print(`*`<run_driftBatchCorrect>`*`)`](https://jllcalorio.github.io/dimsprepr/reference/run_driftBatchCorrect.md)
  : Correct Signal Drift and Batch Effects Using QC Samples

## Data Normalization, Transformation, and Scaling

Transforming and scaling data to ensure comparability across samples.

- [`run_normalize()`](https://jllcalorio.github.io/dimsprepr/reference/run_normalize.md)
  [`print(`*`<run_normalize>`*`)`](https://jllcalorio.github.io/dimsprepr/reference/run_normalize.md)
  : Normalize Metabolomics Data
- [`run_transform()`](https://jllcalorio.github.io/dimsprepr/reference/run_transform.md)
  [`print(`*`<run_transform>`*`)`](https://jllcalorio.github.io/dimsprepr/reference/run_transform.md)
  : Transform Metabolomics Data
- [`run_scale()`](https://jllcalorio.github.io/dimsprepr/reference/run_scale.md)
  [`print(`*`<run_scale>`*`)`](https://jllcalorio.github.io/dimsprepr/reference/run_scale.md)
  : Scale Metabolomics Data

## Plotting Functions of Two Data Frames

Plots two data frames at a time, can be before and after implementing a
method.

- [`plot_beforeafter()`](https://jllcalorio.github.io/dimsprepr/reference/plot_beforeafter.md)
  : Plot Pairwise Before-and-After Comparison for Selected Features
- [`plot_dist_beforeafter()`](https://jllcalorio.github.io/dimsprepr/reference/plot_dist_beforeafter.md)
  : Plot Distribution Comparison Before and After Data Transformation

## Example Datasets

Simulated direct-injection metabolomics data included with dimsprepr for
testing and demonstrating the preprocessing pipeline.

- [`dims_data`](https://jllcalorio.github.io/dimsprepr/reference/dims_data.md)
  : Example Direct-Injection Metabolomics Data Matrix

- [`dims_metadata`](https://jllcalorio.github.io/dimsprepr/reference/dims_metadata.md)
  :

  Example Sample Metadata for `dims_data`

## Internal functions

Internal functions used by dimsprepr.
