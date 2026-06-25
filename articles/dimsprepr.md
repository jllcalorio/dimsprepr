# Getting Started with dimsprepr

## Introduction

**dimsprepr** provides a complete, reproducible preprocessing pipeline
for direct-injection mass spectrometry (DIMS) metabolomics data. It
covers: missing-value filtering, imputation, signal-drift and batch
correction, normalization, transformation, scaling, and quality
filtering — all orchestrated through a single function call.

## Quick Example

``` r
library(dimsprepr)

# Load the bundled example dataset
data(dims_data, package = "dimsprepr")
data(dims_metadata, package = "dimsprepr")

dim(dims_data)
#> [1] 50 80
table(dims_metadata$Group)
#> 
#>   Control        QC Treatment 
#>        20        10        20
```

### Run the full pipeline

The
[`run_DIpreprocess()`](https://jllcalorio.github.io/dimsprepr/reference/run_DIpreprocess.md)
function executes all preprocessing steps in a fixed, validated order.
Here we skip drift correction (since this is a small simulated dataset
without a realistic QC injection pattern):

``` r
result <- run_DIpreprocess(
  x                = dims_data,
  metadata         = dims_metadata,
  normalize_method = "median",
  transform_method = "log10",
  correct_drift    = FALSE,
  rsd_threshold    = 0.5,
  variance_percentile = 5,
  verbose          = FALSE
)

result
#> === run_DIpreprocess pipeline result ===
#> STATUS : OK
#> Data   : 50 samples x 36 features (final)
#> Time   : 4.8 second(s)
```

The result is a named list with intermediate and final data matrices:

``` r
# Final preprocessed matrix (auto-scaled for general analysis)
dim(result$data_nonpls)
#> [1] 50 36

# Step-by-step dimension tracking
result$dimensions
#>                                           Step Samples Features
#> 1                                     Original      50       80
#> 2                  After zero -> NA conversion      50       78
#> 3                         After missing filter      50       77
#> 4                       After drift correction      50       77
#> 5       [median_log10] After NONPLS RSD filter      50       44
#> 6  [median_log10] After NONPLS variance filter      50       38
#> 7          [median_log10] After PLS RSD filter      50       44
#> 8     [median_log10] After PLS variance filter      50       41
#> 9                [median_log10] Final (NONPLS)      50       36
#> 10                  [median_log10] Final (PLS)      50       36
```

## Using Individual Functions

Each preprocessing step is also available as a standalone function for
maximum flexibility.

### 1. Missing-Value Filtering

``` r
filt <- run_filtermissing(
  dims_data,
  dims_metadata,
  threshold = 0.3,
  verbose   = FALSE
)
filt
#> === Missing Value Filtering Results ===
#> Features before: 80
#> Features after:  77
#> Features removed: 3 (3.8%)
#> Threshold: 30.0%
#> Mode: Any group(s) must pass
#> Zeros converted to NA: 342
```

### 2. Imputation

``` r
imp <- run_mvimpute(filt$data, method = 0.2, verbose = FALSE)
imp
#> === Missing Value Imputation Results ===
#> Method: Deterministic: 0.2000 x min positive
#> Missing values before: 195
#> Missing values after:  0
```

### 3. Normalization

``` r
norm <- run_normalize(imp$data, dims_metadata, method = "median", verbose = FALSE)
```

### 4. Transformation

``` r
trans <- run_transform(norm$data, method = "log10", verbose = FALSE)
```

### 5. Scaling

``` r
scaled <- run_scale(trans$data, method = "auto", verbose = FALSE)
scaled
#> === Scaling Results ===
#> Method  : auto
#> Data    : 50 samples x 77 features
```

## Cartesian-Product Pipelines

Supply vectors to `normalize_method` and/or `transform_method` to
compute all combinations automatically. Expensive upstream steps
(filtering, imputation, drift correction) run only once:

``` r
multi <- run_DIpreprocess(
  x                = dims_data,
  metadata         = dims_metadata,
  normalize_method = c("sum", "median"),
  transform_method = c("log10", "log2"),
  correct_drift    = FALSE,
  verbose          = FALSE
)
names(multi)
# [1] "sum_log10"    "sum_log2"     "median_log10" "median_log2"
```

## Session Info

``` r
sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dimsprepr_1.0.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] SummarizedExperiment_1.42.0 gtable_0.3.6               
#>  [3] impute_1.86.0               xfun_0.59                  
#>  [5] bslib_0.11.0                ggplot2_4.0.3              
#>  [7] htmlwidgets_1.6.4           Biobase_2.72.0             
#>  [9] lattice_0.22-9              vctrs_0.7.3                
#> [11] tools_4.6.1                 Rdpack_2.6.6               
#> [13] generics_0.1.4              stats4_4.6.1               
#> [15] parallel_4.6.1              missForest_1.6.1           
#> [17] tibble_3.3.1                pkgconfig_2.0.3            
#> [19] Matrix_1.7-5                RColorBrewer_1.1-3         
#> [21] S7_0.2.2                    desc_1.4.3                 
#> [23] S4Vectors_0.50.1            rngtools_1.5.2             
#> [25] lifecycle_1.0.5             stringr_1.6.0              
#> [27] compiler_4.6.1              farver_2.1.2               
#> [29] textshaping_1.0.5           Seqinfo_1.2.0              
#> [31] codetools_0.2-20            htmltools_0.5.9            
#> [33] sass_0.4.10                 yaml_2.3.12                
#> [35] pillar_1.11.1               pkgdown_2.2.0              
#> [37] jquerylib_0.1.4             cachem_1.1.0               
#> [39] DelayedArray_0.38.2         doRNG_1.8.6.3              
#> [41] iterators_1.0.14            abind_1.4-8                
#> [43] foreach_1.5.2               pcaMethods_2.4.0           
#> [45] tidyselect_1.2.1            digest_0.6.39              
#> [47] stringi_1.8.7               reshape2_1.4.5             
#> [49] dplyr_1.2.1                 fastmap_1.2.0              
#> [51] grid_4.6.1                  cli_3.6.6                  
#> [53] SparseArray_1.12.2          magrittr_2.0.5             
#> [55] S4Arrays_1.12.0             randomForest_4.7-1.2       
#> [57] scales_1.4.0                pmp_1.24.0                 
#> [59] rmarkdown_2.31              XVector_0.52.0             
#> [61] matrixStats_1.5.0           otel_0.2.0                 
#> [63] gridExtra_2.3               ranger_0.18.0              
#> [65] ragg_1.5.2                  evaluate_1.0.5             
#> [67] knitr_1.51                  rbibutils_2.4.1            
#> [69] GenomicRanges_1.64.0        IRanges_2.46.0             
#> [71] rlang_1.2.0                 itertools_0.1-3            
#> [73] Rcpp_1.1.1-1.1              glue_1.8.1                 
#> [75] BiocGenerics_0.58.1         jsonlite_2.0.0             
#> [77] plyr_1.8.9                  R6_2.6.1                   
#> [79] MatrixGenerics_1.24.0       systemfonts_1.3.2          
#> [81] fs_2.1.0
```
