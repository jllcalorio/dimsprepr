#' Example Direct-Injection Metabolomics Data Matrix
#'
#' A simulated numeric matrix of 50 samples x 80 features representing
#' raw peak intensities from a direct-injection mass spectrometry experiment.
#' Rows are samples; columns are features (m/z bins). Zero values represent
#' below-detection-limit measurements. Features 78-80 are intentionally
#' sparse to exercise missing-value filtering.
#'
#' @format A numeric matrix with 50 rows and 80 columns. Row names are
#'   sample identifiers (`S01`--`S50`); column names are feature labels
#'   (`mz_0001`--`mz_0080`).
#'
#' @source Simulated data. See `data-raw/dims_example.R` for the generation
#'   script.
#'
#' @seealso [dims_metadata], [run_DIpreprocess()]
#'
#' @examples
#' data(dims_data)
#' dim(dims_data)
#'
#' @keywords datasets
"dims_data"

#' Example Sample Metadata for \code{dims_data}
#'
#' A data frame of sample-level metadata matching the rows of
#' \code{\link{dims_data}}. Contains group labels, batch assignments,
#' injection order, and subject identifiers.
#'
#' @format A data frame with 50 rows and 5 columns:
#' \describe{
#'   \item{Sample}{Character. Unique sample identifier (`S01`--`S50`).}
#'   \item{Group}{Character. Sample group: `"Control"`, `"Treatment"`, or `"QC"`.}
#'   \item{Batch}{Integer. Batch number (1 or 2).}
#'   \item{InjectionSequence}{Integer. Order of injection (1--50).}
#'   \item{SubjectID}{Character. Biological subject ID for non-QC samples; `NA` for QC.}
#' }
#'
#' @source Simulated data. See `data-raw/dims_example.R`.
#'
#' @seealso [dims_data], [run_DIpreprocess()]
#'
#' @examples
#' data(dims_metadata)
#' table(dims_metadata$Group)
#'
#' @keywords datasets
"dims_metadata"
