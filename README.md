# dimsprepr: An R package for end-to-end preprocessing of direct-injection metabolomics data — missing value filtering, imputation, signal drift and batch correction, normalization, transformation, scaling, and quality filtering — in a single reproducible pipeline.

`dimsprepr` provides a complete, reproducible preprocessing pipeline for direct-injection mass spectrometry (DIMS) metabolomics data. Built around a single orchestrating function, the package executes preprocessing steps — outlier removal, missing value filtering, imputation, signal drift and batch correction, normalization, transformation, scaling, and quality filtering — in a fixed, validated order, ensuring methodological consistency across studies.

A key feature is its support for Cartesian product pipelines: users may supply multiple normalization and transformation methods simultaneously, with computationally expensive upstream steps executed only once and shared across all method combinations. This enables systematic, side-by-side comparison of preprocessing strategies within a single function call.

`dimsprepr` is designed for metabolomics researchers who need a structured, auditable preprocessing workflow ready for downstream statistical analysis.

## Installation

You can install the development version of `dimsprepr` directly from GitHub:

```R
# Install pak if you haven't already
install.packages("pak")

# Install the GitHub package
pak::pak("username/repository")

```
## Contributing

We welcome contributions to `dimsprepr`! Please review our Code of Conduct before contributing.

## Reporting Bugs

If you encounter any bugs or have feature requests, please open an issue on our GitHub Issues page.

## License

This project is licensed under the MIT License.

## Acknowledgements

Huge thanks to The R Project for Statistical Computing, Positron, RStudio, the developers of the dplyr, ggplot2, everything in the tidyverse package, pmp, sva, mice, and all the dependencies I fail to mention. Your work makes `dimsprepr` possible.

## Future Development

More functions are continuously being added to `dimsprepr` to further enhance its capabilities and ease of use. Stay tuned for updates!
