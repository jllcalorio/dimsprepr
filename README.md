# dimsprepr

End-to-end preprocessing of direct-injection mass spectrometry (DIMS)
metabolomics data in a single reproducible pipeline.

`dimsprepr` runs the preprocessing steps of a DIMS metabolomics study in a
fixed, validated order:

1. sample removal (listed outliers, and procedural blanks after blank
   filtering);
2. missing-value handling (failed injections, the modified 80% rule, and
   optionally truly absent data);
3. QC-based signal drift and batch correction;
4. imputation;
5. normalization;
6. transformation;
7. QC-RSD and variance filtering;
8. scaling;
9. replicate merging.

Its quality control steps follow the QComics guidelines
([González-Domínguez et al., 2024](https://doi.org/10.1021/acs.analchem.3c03660)).
Drift correction offers the QC-RSC spline
([Kirwan et al., 2013](https://doi.org/10.1007/s00216-013-6856-7)) and the
outlier-robust rLOESS, rGAM, and tGAM methods
([Vynck et al., 2026](https://doi.org/10.1021/acs.analchem.5c06841)).

Several normalization, transformation, and scaling methods can be supplied
at once: the pipeline computes every combination while running the shared
upstream steps only once, for side-by-side comparison of preprocessing
strategies.

## Installation

Install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("jllcalorio/dimsprepr")
```

The default pipeline needs no Bioconductor packages. Optional methods use
`mgcv` (robust tGAM/rGAM drift correction), `imputeLCMD` and `mice`
(imputation), `sva` (ComBat), `pmp` (glog transformation), and `vsn` (VSN
transformation). Install them if you need those methods; the Bioconductor
ones (`pmp`, `sva`, `vsn`, `imputeLCMD`'s dependencies) via `BiocManager`:

```r
install.packages("BiocManager")
BiocManager::install(c("pmp", "sva", "vsn"))
```

## Quick start

```r
library(dimsprepr)

res <- run_DIpreprocess(dims_data, dims_metadata)
res                 # summary of what was removed and kept
res$data_nonpls     # auto-scaled: PCA, clustering, univariate tests
res$data_pls        # Pareto-scaled: PLS-DA / OPLS-DA

plot_beforeafter(res)          # drift correction along the injection order
plot_dist_beforeafter(res)     # distributions before and after
```

See `vignette("dimsprepr")` for a walk-through.

## Contributing

Contributions are welcome. Please review the Code of Conduct before
contributing.

## Reporting bugs

Please report bugs and feature requests on the
[GitHub issues page](https://github.com/jllcalorio/dimsprepr/issues).

## License

MIT. See `LICENSE`.

## Acknowledgements

Thanks to the R Project for Statistical Computing, Positron, RStudio, and the
developers of ggplot2, mgcv, pmp, sva, mice, vsn, imputeLCMD, and all the
other dependencies.
