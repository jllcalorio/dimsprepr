# Internal validation helpers – shared across run_* functions
# ponytail: single source of truth for repeated validation boilerplate

#' Validate that x is a matrix or data frame
#' @noRd
.validate_data_matrix <- function(x, name = "x") {
  if (!is.matrix(x) && !is.data.frame(x))
    stop(sprintf("'%s' must be a matrix or data frame.", name), call. = FALSE)
}

#' Validate metadata is a data frame with matching rows
#' @noRd
.validate_metadata <- function(x, metadata) {
  if (!is.data.frame(metadata))
    stop("'metadata' must be a data frame.", call. = FALSE)
  if (nrow(x) != nrow(metadata))
    stop(sprintf(
      "nrow(x) (%d) must equal nrow(metadata) (%d).",
      nrow(x), nrow(metadata)
    ), call. = FALSE)
}

#' Convert to matrix preserving original class info
#' @return list(mat, was_matrix)
#' @noRd
.as_matrix_preserve <- function(x) {
  list(mat = as.matrix(x), was_matrix = is.matrix(x))
}

#' Restore original class (matrix or data.frame)
#' @noRd
.restore_class <- function(x_matrix, was_matrix) {
  if (was_matrix) x_matrix else as.data.frame(x_matrix)
}

#' Verbose message factory
#' @noRd
.msg <- function(verbose) {
  function(...) if (verbose) message(...)
}

#' Collapse QC types into a single label in a Group_ column
#' @noRd
.add_group_collapsed <- function(metadata, group_col, qc_types, qc_label = "QC") {
  metadata$Group_ <- ifelse(
    metadata[[group_col]] %in% qc_types, qc_label, metadata[[group_col]]
  )
  metadata
}

#' Require a column to exist in a data frame
#' @noRd
.require_col <- function(df, col, df_name = "metadata") {
  if (!col %in% colnames(df))
    stop(sprintf("Column '%s' not found in %s.", col, df_name), call. = FALSE)
}
